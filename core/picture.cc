
#include <math.h>
#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "fmo.h"
#include "bitstream_nal.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"

#include "sei.h"
#include "output.h"
#include "neighbour.h"
#include "memalloc.h"
#include "macroblock.h"
#include "mb.h"

#include "intra_prediction.h"
#include "deblock.h"

#include "biaridecod.h"

#include "erc_api.h"
#include "dpb_common.h"
#include "dpb_mvc.h"



static void init_picture_decoding(VideoParameters *p_Vid)
{
    Slice *pSlice = p_Vid->ppSliceList[0];
    int j, i, iDeblockMode = 1;

    if (p_Vid->iSliceNumOfCurrPic >= MAX_NUM_SLICES)
        error ("Maximum number of supported slices exceeded. \nPlease recompile with increased value for MAX_NUM_SLICES", 200);

    if (p_Vid->pNextPPS->Valid && p_Vid->pNextPPS->pic_parameter_set_id == pSlice->pic_parameter_set_id) {
        pps_t tmpPPS;
        memcpy(&tmpPPS, &(p_Vid->PicParSet[pSlice->pic_parameter_set_id]), sizeof (pps_t));
        (p_Vid->PicParSet[pSlice->pic_parameter_set_id]).slice_group_id = NULL;
        MakePPSavailable(p_Vid, p_Vid->pNextPPS->pic_parameter_set_id, p_Vid->pNextPPS);
        memcpy(p_Vid->pNextPPS, &tmpPPS, sizeof (pps_t));
        tmpPPS.slice_group_id = NULL;
    }

    UseParameterSet(pSlice);
    if (pSlice->idr_flag)
        p_Vid->number = 0;

    p_Vid->PicHeightInMbs = p_Vid->FrameHeightInMbs / ( 1 + pSlice->field_pic_flag );
    p_Vid->PicSizeInMbs   = p_Vid->PicWidthInMbs * p_Vid->PicHeightInMbs;
    p_Vid->FrameSizeInMbs = p_Vid->PicWidthInMbs * p_Vid->FrameHeightInMbs;
    p_Vid->structure = pSlice->structure;

    fmo_init(p_Vid, pSlice);

#if (MVC_EXTENSION_ENABLE)
    if (pSlice->layer_id > 0 && pSlice->svc_extension_flag == 0 && pSlice->NaluHeaderMVCExt.non_idr_flag == 0)
        idr_memory_management(p_Vid->p_Dpb_layer[pSlice->layer_id], p_Vid->dec_picture);
    update_ref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_ltref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_pic_num(pSlice);
    i = pSlice->view_id;
#endif
    init_Deblock(p_Vid, pSlice->mb_aff_frame_flag);
    //init mb_data;
    for (j = 0; j < p_Vid->iSliceNumOfCurrPic; j++) {
        if (p_Vid->ppSliceList[j]->disable_deblocking_filter_idc != 1)
            iDeblockMode = 0;
#if (MVC_EXTENSION_ENABLE)
        assert(p_Vid->ppSliceList[j]->view_id == i);
#endif
    }
    p_Vid->iDeblockMode = iDeblockMode;
}


//! gives CBP value from codeword number, both for intra and inter
static const byte NCBP[2][48][2] = {
  {  // 0      1        2       3       4       5       6       7       8       9      10      11
    {15, 0},{ 0, 1},{ 7, 2},{11, 4},{13, 8},{14, 3},{ 3, 5},{ 5,10},{10,12},{12,15},{ 1, 7},{ 2,11},
    { 4,13},{ 8,14},{ 6, 6},{ 9, 9},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},
    { 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},
    { 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0}
  },
  {
    {47, 0},{31,16},{15, 1},{ 0, 2},{23, 4},{27, 8},{29,32},{30, 3},{ 7, 5},{11,10},{13,12},{14,15},
    {39,47},{43, 7},{45,11},{46,13},{16,14},{ 3, 6},{ 5, 9},{10,31},{12,35},{19,37},{21,42},{26,44},
    {28,33},{35,34},{37,36},{42,40},{44,39},{ 1,43},{ 2,45},{ 4,46},{ 8,17},{17,18},{18,20},{20,24},
    {24,19},{ 6,21},{ 9,26},{22,28},{25,23},{32,27},{33,29},{34,30},{36,22},{40,25},{38,38},{41,41}
  }
};

static void linfo_cbp_intra_normal(int len,int info,int *cbp, int *dummy)
{
  int cbp_idx;

  linfo_ue(len, info, &cbp_idx, dummy);
  *cbp=NCBP[1][cbp_idx][0];
}

static void linfo_cbp_intra_other(int len,int info,int *cbp, int *dummy)
{
  int cbp_idx;

  linfo_ue(len, info, &cbp_idx, dummy);
  *cbp=NCBP[0][cbp_idx][0];
}

static void linfo_cbp_inter_normal(int len,int info,int *cbp, int *dummy)
{
  int cbp_idx;

  linfo_ue(len, info, &cbp_idx, dummy);
  *cbp=NCBP[1][cbp_idx][1];
}

static void linfo_cbp_inter_other(int len,int info,int *cbp, int *dummy)
{
  int cbp_idx;

  linfo_ue(len, info, &cbp_idx, dummy);
  *cbp=NCBP[0][cbp_idx][1];
}


/*!
 ************************************************************************
 * \brief
 *    write the encoding mode and motion vectors of current
 *    MB to the buffer of the error concealment module.
 ************************************************************************
 */
static void ercWriteMBMODEandMV(Macroblock *currMB)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  int i, ii, jj, currMBNum = currMB->mbAddrX; //p_Vid->currentSlice->current_mb_nr;
  StorablePicture *dec_picture = p_Vid->dec_picture;
  int mbx = xPosMB(currMBNum, dec_picture->size_x), mby = yPosMB(currMBNum, dec_picture->size_x);
  objectBuffer_t *currRegion, *pRegion;

  currRegion = p_Vid->erc_object_list + (currMBNum<<2);

  if(p_Vid->type != B_SLICE) //non-B frame
  {
    for (i=0; i<4; ++i)
    {
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
        currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  :
        currMB->b8mode[i]==0      ? REGMODE_INTER_COPY :
        currMB->b8mode[i]==1      ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);
      if (currMB->b8mode[i]==0 || currMB->b8mode[i]==IBLOCK)  // INTRA OR COPY
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        ii              = 4*mbx + (i & 0x01)*2;// + BLOCK_SIZE;
        jj              = 4*mby + (i >> 1  )*2;
        if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS
        {
          pRegion->mv[0]  = (dec_picture->mv_info[jj][ii].mv[LIST_0].mv_x + dec_picture->mv_info[jj][ii + 1].mv[LIST_0].mv_x + dec_picture->mv_info[jj + 1][ii].mv[LIST_0].mv_x + dec_picture->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_x + 2)/4;
          pRegion->mv[1]  = (dec_picture->mv_info[jj][ii].mv[LIST_0].mv_y + dec_picture->mv_info[jj][ii + 1].mv[LIST_0].mv_y + dec_picture->mv_info[jj + 1][ii].mv[LIST_0].mv_y + dec_picture->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_y + 2)/4;
        }
        else // 16x16, 16x8, 8x16, 8x8
        {
          pRegion->mv[0]  = dec_picture->mv_info[jj][ii].mv[LIST_0].mv_x;
          pRegion->mv[1]  = dec_picture->mv_info[jj][ii].mv[LIST_0].mv_y;
        }
        currMB->p_Slice->erc_mvperMB      += iabs(pRegion->mv[0]) + iabs(pRegion->mv[1]);
        pRegion->mv[2]    = dec_picture->mv_info[jj][ii].ref_idx[LIST_0];
      }
    }
  }
  else  //B-frame
  {
    for (i=0; i<4; ++i)
    {
      ii                  = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
      jj                  = 4*mby + (i/2)*2;
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
        currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
      if (currMB->mb_type==I16MB || currMB->b8mode[i]==IBLOCK)  // INTRA
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        int idx = (dec_picture->mv_info[jj][ii].ref_idx[0] < 0) ? 1 : 0;
        pRegion->mv[0]    = (dec_picture->mv_info[jj][ii].mv[idx].mv_x + 
          dec_picture->mv_info[jj][ii+1].mv[idx].mv_x + 
          dec_picture->mv_info[jj+1][ii].mv[idx].mv_x + 
          dec_picture->mv_info[jj+1][ii+1].mv[idx].mv_x + 2)/4;
        pRegion->mv[1]    = (dec_picture->mv_info[jj][ii].mv[idx].mv_y + 
          dec_picture->mv_info[jj][ii+1].mv[idx].mv_y + 
          dec_picture->mv_info[jj+1][ii].mv[idx].mv_y + 
          dec_picture->mv_info[jj+1][ii+1].mv[idx].mv_y + 2)/4;
        currMB->p_Slice->erc_mvperMB      += iabs(pRegion->mv[0]) + iabs(pRegion->mv[1]);

        pRegion->mv[2]  = (dec_picture->mv_info[jj][ii].ref_idx[idx]);
      }
    }
  }
}

static void reorder_lists(Slice *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;

    if (currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE) {
        if (currSlice->ref_pic_list_reordering_flag[LIST_0])
            reorder_ref_pic_list(currSlice, LIST_0);
        if (p_Vid->no_reference_picture == currSlice->listX[0][currSlice->num_ref_idx_l0_active_minus1]) {
            if (p_Vid->non_conforming_stream)
                printf("RefPicList0[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_l0_active_minus1);
            else
                error("RefPicList0[ num_ref_idx_l0_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
        }
        // that's a definition
        currSlice->listXsize[0] = (char) currSlice->num_ref_idx_l0_active_minus1 + 1;
    }

    if (currSlice->slice_type == B_SLICE) {
        if (currSlice->ref_pic_list_reordering_flag[LIST_1])
            reorder_ref_pic_list(currSlice, LIST_1);
        if (p_Vid->no_reference_picture == currSlice->listX[1][currSlice->num_ref_idx_l1_active_minus1]) {
            if (p_Vid->non_conforming_stream)
                printf("RefPicList1[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_l1_active_minus1);
            else
                error("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
        }
        // that's a definition
        currSlice->listXsize[1] = (char) currSlice->num_ref_idx_l1_active_minus1 + 1;
    }

    free_ref_pic_list_reordering_buffer(currSlice);
}


static void init_slice(VideoParameters *p_Vid, Slice *currSlice)
{
    int i;

    p_Vid->active_sps = currSlice->active_sps;
    p_Vid->active_pps = currSlice->active_pps;

    init_lists(currSlice);

#if (MVC_EXTENSION_ENABLE)
    if (currSlice->svc_extension_flag == 0 || currSlice->svc_extension_flag == 1)
        reorder_lists_mvc(currSlice, currSlice->ThisPOC);
    else
        reorder_lists(currSlice);

    if (currSlice->fs_listinterview0) {
        free(currSlice->fs_listinterview0);
        currSlice->fs_listinterview0 = NULL;
    }
    if (currSlice->fs_listinterview1) {
        free(currSlice->fs_listinterview1);
        currSlice->fs_listinterview1 = NULL;
    }
#endif

    if (currSlice->structure == FRAME)
        init_mbaff_lists(p_Vid, currSlice);

    // update reference flags and set current p_Vid->ref_flag
    if (!(currSlice->redundant_pic_cnt != 0 && p_Vid->previous_frame_num == currSlice->frame_num)) {
        for (i = 16; i > 0; i--)
            currSlice->ref_flag[i] = currSlice->ref_flag[i-1];
    }
    currSlice->ref_flag[0] = currSlice->redundant_pic_cnt == 0 ? p_Vid->Is_primary_correct
                                                               : p_Vid->Is_redundant_correct;

    if (currSlice->active_sps->chroma_format_idc == 0 || currSlice->active_sps->chroma_format_idc == 3) {
        currSlice->linfo_cbp_intra = linfo_cbp_intra_other;
        currSlice->linfo_cbp_inter = linfo_cbp_inter_other;
    } else {
        currSlice->linfo_cbp_intra = linfo_cbp_intra_normal;
        currSlice->linfo_cbp_inter = linfo_cbp_inter_normal;
    }

    if (currSlice->active_pps->entropy_coding_mode_flag) {
        init_contexts(currSlice);
        currSlice->last_dquant = 0;
    }

    if ((currSlice->active_pps->weighted_bipred_idc > 0 && currSlice->slice_type == B_SLICE) ||
        (currSlice->active_pps->weighted_pred_flag && currSlice->slice_type != I_SLICE))
        fill_wp_params(currSlice);
}

// this is intended to make get_block_luma faster by doing this at a more appropriate level
// i.e. per slice rather than per MB
static void init_cur_imgy(Slice *currSlice, VideoParameters *p_Vid)
{
    int i, j;
    if (p_Vid->separate_colour_plane_flag != 0) {
        StorablePicture *vidref = p_Vid->no_reference_picture;
        int noref = (currSlice->framepoc < p_Vid->recovery_poc);
        if (currSlice->colour_plane_id == 0) {
            for (j = 0; j < 6; j++) {
                for (i = 0; i < MAX_LIST_SIZE; i++) {
                    StorablePicture *curr_ref = currSlice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = curr_ref->imgY;
                    }
                }
            }
        }
    } else {
        StorablePicture *vidref = p_Vid->no_reference_picture;
        int noref = (currSlice->framepoc < p_Vid->recovery_poc);
        int total_lists = currSlice->mb_aff_frame_flag ? 6 :
                          currSlice->slice_type == B_SLICE ? 2 : 1;
        for (j = 0; j < total_lists; j++) {
            // note that if we always set this to MAX_LIST_SIZE, we avoid crashes with invalid ref_idx being set
            // since currently this is done at the slice level, it seems safe to do so.
            // Note for some reason I get now a mismatch between version 12 and this one in cabac. I wonder why.
            for (i = 0; i < MAX_LIST_SIZE; i++) {
                StorablePicture *curr_ref = currSlice->listX[j][i];
                if (curr_ref) {
                    curr_ref->no_ref = noref && (curr_ref == vidref);
                    curr_ref->cur_imgY = curr_ref->imgY;
                }
            }
        }
    }
}

/*!
 ************************************************************************
 * \brief
 *    decodes one slice
 ************************************************************************
 */
static void decode_one_slice(Slice *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    int current_header = currSlice->current_header;

    bool end_of_slice = 0;
    Macroblock *currMB = NULL;
    currSlice->cod_counter = -1;

    init_slice(p_Vid, currSlice);

    if ((current_header != SOP && current_header != SOS) || currSlice->ei_flag != 0)
        return;

    if (p_Vid->separate_colour_plane_flag != 0)
        change_plane_JV(p_Vid, currSlice->colour_plane_id, currSlice);
    else {
        currSlice->mb_data     = p_Vid->mb_data;
        currSlice->dec_picture = p_Vid->dec_picture;
        currSlice->siblock     = p_Vid->siblock;
        currSlice->ipredmode   = p_Vid->ipredmode;
        currSlice->intra_block = p_Vid->intra_block;
    }

    if (currSlice->slice_type == B_SLICE)
        compute_colocated(currSlice, currSlice->listX);

    if (currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE)
        init_cur_imgy(currSlice, p_Vid); 

    while (!end_of_slice) { // loop over macroblocks
        // Initializes the current macroblock
        start_macroblock(currSlice, &currMB);
        // Get the syntax elements from the NAL
        currSlice->read_one_macroblock(currMB);
        decode_one_macroblock(currMB, currSlice->dec_picture);

        if (currSlice->mb_aff_frame_flag && currMB->mb_field) {
            currSlice->num_ref_idx_l0_active_minus1 = ((currSlice->num_ref_idx_l0_active_minus1 + 1) >> 1) - 1;
            currSlice->num_ref_idx_l1_active_minus1 = ((currSlice->num_ref_idx_l1_active_minus1 + 1) >> 1) - 1;
        }

#if (DISABLE_ERC == 0)
        ercWriteMBMODEandMV(currMB);
#endif

        end_of_slice = exit_macroblock(currSlice, (!currSlice->mb_aff_frame_flag|| currSlice->current_mb_nr%2));
    }
}


void decode_picture(VideoParameters *p_Vid)
{
    Slice **ppSliceList = p_Vid->ppSliceList;
    int iSliceNo;

    p_Vid->num_dec_mb = 0;

    init_picture_decoding(p_Vid);

    for (iSliceNo = 0; iSliceNo < p_Vid->iSliceNumOfCurrPic; iSliceNo++) {
        Slice *currSlice = ppSliceList[iSliceNo];

        assert(currSlice->current_header != EOS);
        assert(currSlice->current_slice_nr == iSliceNo);

        decode_one_slice(currSlice);

        p_Vid->num_dec_mb  += currSlice->num_dec_mb;
        p_Vid->erc_mvperMB += currSlice->erc_mvperMB;
    }
}
