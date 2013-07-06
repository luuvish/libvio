
/*!
 ***********************************************************************
 * \file image.c
 *
 * \brief
 *    Decode a Slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Gabi Blaettermann
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *    - Alexis Tourapis                 <alexismt@ieee.org>
 *    - Jill Boyce                      <jill.boyce@thomson.net>
 *    - Saurav K Bandyopadhyay          <saurav@ieee.org>
 *    - Zhenyu Wu                       <Zhenyu.Wu@thomson.net
 *    - Purvin Pandit                   <Purvin.Pandit@thomson.net>
 *
 ***********************************************************************
 */

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


static inline float psnr(int max_sample_sq, int samples, float sse_distortion ) 
{
  return (float) (sse_distortion == 0.0 ? 0.0 : (10.0 * log10(max_sample_sq * (double) ((double) samples / sse_distortion))));
}

static inline int iabs2(int x) 
{
  return (x) * (x);
}

static inline void reset_mbs(Macroblock *currMB)
{
  currMB->slice_nr = -1; 
  currMB->ei_flag  =  1;
  currMB->dpl_flag =  0;
}

static void setup_buffers(VideoParameters *p_Vid, int layer_id)
{
  CodingParameters *cps = p_Vid->p_EncodePar[layer_id];
  int i;

  if(p_Vid->last_dec_layer_id != layer_id)
  {
    p_Vid->imgY_ref = cps->imgY_ref;
    p_Vid->imgUV_ref = cps->imgUV_ref;
    if(p_Vid->active_sps->separate_colour_plane_flag)
    {
     for( i=0; i<MAX_PLANE; i++ )
     {
       p_Vid->mb_data_JV[i] = cps->mb_data_JV[i];
       p_Vid->intra_block_JV[i] = cps->intra_block_JV[i];
       p_Vid->ipredmode_JV[i] = cps->ipredmode_JV[i];
       p_Vid->siblock_JV[i] = cps->siblock_JV[i];
     }
     p_Vid->mb_data = NULL;
     p_Vid->intra_block = NULL;
     p_Vid->ipredmode = NULL;
     p_Vid->siblock = NULL;
    }
    else
    {
      p_Vid->mb_data = cps->mb_data;
      p_Vid->intra_block = cps->intra_block;
      p_Vid->ipredmode = cps->ipredmode;
      p_Vid->siblock = cps->siblock;
    }
    p_Vid->PicPos = cps->PicPos;
    p_Vid->nz_coeff = cps->nz_coeff;
    p_Vid->qp_per_matrix = cps->qp_per_matrix;
    p_Vid->qp_rem_matrix = cps->qp_rem_matrix;
    p_Vid->img2buf = cps->img2buf;
    p_Vid->last_dec_layer_id = layer_id;
  }
}

#if MVC_EXTENSION_ENABLE
static void init_mvc_picture(Slice *currSlice)
{
  int i;
  VideoParameters *p_Vid = currSlice->p_Vid;
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb_layer[0];

  StorablePicture *p_pic = NULL;

  // find BL reconstructed picture
  if (!currSlice->field_pic_flag)
  {
    for (i = 0; i < (int)p_Dpb->used_size/*size*/; i++)
    {
      FrameStore *fs = p_Dpb->fs[i];
      if ((fs->frame->view_id == 0) && (fs->frame->frame_poc == currSlice->framepoc))
      {
        p_pic = fs->frame;
        break;
      }
    }
  }
  else if (!currSlice->bottom_field_flag)
  {
    for (i = 0; i < (int)p_Dpb->used_size/*size*/; i++)
    {
      FrameStore *fs = p_Dpb->fs[i];
      if ((fs->top_field->view_id == 0) && (fs->top_field->top_poc == currSlice->toppoc))
      {
        p_pic = fs->top_field;
        break;
      }
    }
  }
  else
  {
    for (i = 0; i < (int)p_Dpb->used_size/*size*/; i++)
    {
      FrameStore *fs = p_Dpb->fs[i];
      if ((fs->bottom_field->view_id == 0) && (fs->bottom_field->bottom_poc == currSlice->bottompoc))
      {
        p_pic = fs->bottom_field;
        break;
      }
    }
  }
  if(!p_pic)
  {
    p_Vid->bFrameInit = 0;
  }
  else
  {
    process_picture_in_dpb_s(p_Vid, p_pic);
    store_proc_picture_in_dpb (currSlice->p_Dpb, clone_storable_picture(p_Vid, p_pic));
  }
}
#endif


static void copy_dec_picture_JV( VideoParameters *p_Vid, StorablePicture *dst, StorablePicture *src )
{
  dst->top_poc              = src->top_poc;
  dst->bottom_poc           = src->bottom_poc;
  dst->frame_poc            = src->frame_poc;
  dst->qp                   = src->qp;
  dst->slice_qp_delta       = src->slice_qp_delta;
  dst->chroma_qp_offset[0]  = src->chroma_qp_offset[0];
  dst->chroma_qp_offset[1]  = src->chroma_qp_offset[1];

  dst->poc                  = src->poc;

  dst->slice_type           = src->slice_type;
  dst->used_for_reference   = src->used_for_reference;
  dst->idr_flag             = src->idr_flag;
  dst->no_output_of_prior_pics_flag = src->no_output_of_prior_pics_flag;
  dst->long_term_reference_flag = src->long_term_reference_flag;
  dst->adaptive_ref_pic_buffering_flag = src->adaptive_ref_pic_buffering_flag;

  dst->dec_ref_pic_marking_buffer = src->dec_ref_pic_marking_buffer;

  dst->mb_aff_frame_flag    = src->mb_aff_frame_flag;
  dst->PicWidthInMbs        = src->PicWidthInMbs;
  dst->pic_num              = src->pic_num;
  dst->frame_num            = src->frame_num;
  dst->recovery_frame       = src->recovery_frame;
  dst->coded_frame          = src->coded_frame;

  dst->chroma_format_idc    = src->chroma_format_idc;

  dst->frame_mbs_only_flag  = src->frame_mbs_only_flag;
  dst->frame_cropping_flag  = src->frame_cropping_flag;

  dst->frame_crop_left_offset   = src->frame_crop_left_offset;
  dst->frame_crop_right_offset  = src->frame_crop_right_offset;
  dst->frame_crop_top_offset    = src->frame_crop_top_offset;
  dst->frame_crop_bottom_offset = src->frame_crop_bottom_offset;

#if (ENABLE_OUTPUT_TONEMAPPING)
  // store the necessary tone mapping sei into StorablePicture structure
  dst->seiHasTone_mapping = src->seiHasTone_mapping;

  dst->seiHasTone_mapping    = src->seiHasTone_mapping;
  dst->tone_mapping_model_id = src->tone_mapping_model_id;
  dst->tonemapped_bit_depth  = src->tonemapped_bit_depth;
  if( src->tone_mapping_lut )
  {
    int coded_data_bit_max = (1 << p_Vid->seiToneMapping->coded_data_bit_depth);
    dst->tone_mapping_lut      = (imgpel *)malloc(sizeof(int) * coded_data_bit_max);
    if (NULL == dst->tone_mapping_lut)
    {
      no_mem_exit("copy_dec_picture_JV: tone_mapping_lut");
    }
    memcpy(dst->tone_mapping_lut, src->tone_mapping_lut, sizeof(imgpel) * coded_data_bit_max);
  }
#endif
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new picture
 ************************************************************************
 */
void init_picture(VideoParameters *p_Vid, Slice *currSlice, InputParameters *p_Inp)
{
  int i;
  int nplane;
  StorablePicture *dec_picture = NULL;
  sps_t *sps = p_Vid->active_sps;
  pps_t *pps = p_Vid->active_pps;
  DecodedPictureBuffer *p_Dpb = currSlice->p_Dpb;

  int PicSizeInMbs = sps->PicWidthInMbs * (sps->FrameHeightInMbs / (1 + currSlice->field_pic_flag));

  p_Vid->bFrameInit = 1;
  if (p_Vid->dec_picture)
  {
    // this may only happen on slice loss
    exit_picture(p_Vid, &p_Vid->dec_picture);
  }
  p_Vid->dpb_layer_id = currSlice->layer_id;
  //set buffers;
  setup_buffers(p_Vid, currSlice->layer_id);

  if (p_Vid->recovery_point)
    p_Vid->recovery_frame_num = (currSlice->frame_num + p_Vid->recovery_frame_cnt) % sps->MaxFrameNum;

  if (currSlice->idr_flag)
    p_Vid->recovery_frame_num = currSlice->frame_num;

  if (p_Vid->recovery_point == 0 &&
    currSlice->frame_num != p_Vid->pre_frame_num &&
    currSlice->frame_num != (p_Vid->pre_frame_num + 1) % sps->MaxFrameNum)
  {
    if (sps->gaps_in_frame_num_value_allowed_flag == 0)
    {
#if (DISABLE_ERC == 0)
      // picture error concealment
      if(p_Inp->conceal_mode !=0)
      {
        if((currSlice->frame_num) < ((p_Vid->pre_frame_num + 1) % sps->MaxFrameNum))
        {
          /* Conceal lost IDR frames and any frames immediately
             following the IDR. Use frame copy for these since
             lists cannot be formed correctly for motion copy*/
          p_Vid->conceal_mode = 1;
          p_Vid->IDR_concealment_flag = 1;
          conceal_lost_frames(p_Dpb, currSlice);
          //reset to original concealment mode for future drops
          p_Vid->conceal_mode = p_Inp->conceal_mode;
        }
        else
        {
          //reset to original concealment mode for future drops
          p_Vid->conceal_mode = p_Inp->conceal_mode;

          p_Vid->IDR_concealment_flag = 0;
          conceal_lost_frames(p_Dpb, currSlice);
        }
      }
      else
#endif
      {   /* Advanced Error Concealment would be called here to combat unintentional loss of pictures. */
        error("An unintentional loss of pictures occurs! Exit\n", 100);
      }
    }
    if(p_Vid->conceal_mode == 0)
      fill_frame_num_gap(p_Vid, currSlice);
  }

  if(currSlice->nal_ref_idc)
  {
    p_Vid->pre_frame_num = currSlice->frame_num;
  }

  //calculate POC
  decode_poc(p_Vid, currSlice);

  if (p_Vid->recovery_frame_num == (int) currSlice->frame_num && p_Vid->recovery_poc == 0x7fffffff)
    p_Vid->recovery_poc = currSlice->framepoc;

  if(currSlice->nal_ref_idc)
    p_Vid->last_ref_pic_poc = currSlice->framepoc;

  if (!currSlice->field_pic_flag || !currSlice->bottom_field_flag)
  {
    gettime (&(p_Vid->start_time));             // start time
  }

  dec_picture = p_Vid->dec_picture = alloc_storable_picture (p_Vid, currSlice->structure, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr, 1);
  dec_picture->top_poc=currSlice->toppoc;
  dec_picture->bottom_poc=currSlice->bottompoc;
  dec_picture->frame_poc=currSlice->framepoc;
  dec_picture->qp = currSlice->SliceQpY;
  dec_picture->slice_qp_delta = currSlice->slice_qp_delta;
  dec_picture->chroma_qp_offset[0] = pps->chroma_qp_index_offset;
  dec_picture->chroma_qp_offset[1] = pps->second_chroma_qp_index_offset;
  dec_picture->iCodingType = !currSlice->field_pic_flag ? (currSlice->MbaffFrameFlag? FRAME_MB_PAIR_CODING:FRAME_CODING): FIELD_CODING;
  dec_picture->layer_id = currSlice->layer_id;
#if (MVC_EXTENSION_ENABLE)
  dec_picture->view_id         = currSlice->view_id;
  dec_picture->inter_view_flag = currSlice->inter_view_flag;
  dec_picture->anchor_pic_flag = currSlice->anchor_pic_flag;
  if (dec_picture->view_id == 1)
  {
    if((p_Vid->profile_idc == MVC_HIGH) || (p_Vid->profile_idc == STEREO_HIGH))
      init_mvc_picture(currSlice);
  }
#endif

  // reset all variables of the error concealment instance before decoding of every frame.
  // here the third parameter should, if perfectly, be equal to the number of slices per frame.
  // using little value is ok, the code will allocate more memory if the slice number is larger
#if (DISABLE_ERC == 0)
  ercReset(p_Vid->erc_errorVar, PicSizeInMbs, PicSizeInMbs, dec_picture->size_x);
#endif
  p_Vid->erc_mvperMB = 0;

    if (!currSlice->field_pic_flag)
        dec_picture->poc = currSlice->framepoc;
    else if (!currSlice->bottom_field_flag) {
        dec_picture->poc = currSlice->toppoc;
        p_Vid->number *= 2;
    } else {
        dec_picture->poc = currSlice->bottompoc;
        p_Vid->number = p_Vid->number * 2 + 1;
    }

  if (p_Vid->type > SI_SLICE)
  {
    p_Vid->type = P_SLICE;  // concealed element
  }

  // CAVLC init
  if (!pps->entropy_coding_mode_flag)
  {
    memset(p_Vid->nz_coeff[0][0][0], -1, PicSizeInMbs * 48 *sizeof(byte)); // 3 * 4 * 4
  }

  // Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  // TO set Macroblock Map (mark all MBs as 'have to be concealed')
  if(sps->separate_colour_plane_flag)
  {
    for( nplane=0; nplane<MAX_PLANE; ++nplane )
    {      
      Macroblock *currMB = p_Vid->mb_data_JV[nplane];
      char *intra_block = p_Vid->intra_block_JV[nplane];
      for(i=0; i < PicSizeInMbs; ++i)
      {
        reset_mbs(currMB++);
      }
      memset(p_Vid->ipredmode_JV[nplane][0], DC_PRED, 16 * sps->FrameHeightInMbs * sps->PicWidthInMbs * sizeof(char));
      if(pps->constrained_intra_pred_flag)
      {
        for (i=0; i<PicSizeInMbs; ++i)
        {
          intra_block[i] = 1;
        }
      }
    }
  }
  else
  {
    Macroblock *currMB = p_Vid->mb_data;
    for(i=0; i<PicSizeInMbs; ++i)
      reset_mbs(currMB++);
    if(pps->constrained_intra_pred_flag)
    {
      for (i=0; i<PicSizeInMbs; ++i)
      {
        p_Vid->intra_block[i] = 1;
      }
    }
    memset(p_Vid->ipredmode[0], DC_PRED, 16 * sps->FrameHeightInMbs * sps->PicWidthInMbs * sizeof(char));
  }  

  dec_picture->slice_type = p_Vid->type;
  dec_picture->used_for_reference = (currSlice->nal_ref_idc != 0);
  dec_picture->idr_flag = currSlice->idr_flag;
  dec_picture->no_output_of_prior_pics_flag = currSlice->no_output_of_prior_pics_flag;
  dec_picture->long_term_reference_flag     = currSlice->long_term_reference_flag;
  dec_picture->adaptive_ref_pic_buffering_flag = currSlice->adaptive_ref_pic_marking_mode_flag;

  dec_picture->dec_ref_pic_marking_buffer = currSlice->dec_ref_pic_marking_buffer;
  currSlice->dec_ref_pic_marking_buffer   = NULL;

  dec_picture->mb_aff_frame_flag = currSlice->MbaffFrameFlag;
  dec_picture->PicWidthInMbs     = sps->PicWidthInMbs;

  p_Vid->get_mb_block_pos = dec_picture->mb_aff_frame_flag ? get_mb_block_pos_mbaff : get_mb_block_pos_normal;
  p_Vid->getNeighbour     = dec_picture->mb_aff_frame_flag ? getAffNeighbour : getNonAffNeighbour;

  dec_picture->pic_num   = currSlice->frame_num;
  dec_picture->frame_num = currSlice->frame_num;

  dec_picture->recovery_frame = (unsigned int) ((int) currSlice->frame_num == p_Vid->recovery_frame_num);

  dec_picture->coded_frame = !currSlice->field_pic_flag;

  dec_picture->chroma_format_idc = sps->chroma_format_idc;

  dec_picture->frame_mbs_only_flag = sps->frame_mbs_only_flag;
  dec_picture->frame_cropping_flag = sps->frame_cropping_flag;

  if (dec_picture->frame_cropping_flag)
  {
    dec_picture->frame_crop_left_offset   = sps->frame_crop_left_offset;
    dec_picture->frame_crop_right_offset  = sps->frame_crop_right_offset;
    dec_picture->frame_crop_top_offset    = sps->frame_crop_top_offset;
    dec_picture->frame_crop_bottom_offset = sps->frame_crop_bottom_offset;
  }

#if (ENABLE_OUTPUT_TONEMAPPING)
  // store the necessary tone mapping sei into StorablePicture structure
  if (p_Vid->seiToneMapping->seiHasTone_mapping)
  {
    int coded_data_bit_max = (1 << p_Vid->seiToneMapping->coded_data_bit_depth);
    dec_picture->seiHasTone_mapping    = 1;
    dec_picture->tone_mapping_model_id = p_Vid->seiToneMapping->model_id;
    dec_picture->tonemapped_bit_depth  = p_Vid->seiToneMapping->sei_bit_depth;
    dec_picture->tone_mapping_lut      = (imgpel *)malloc(coded_data_bit_max * sizeof(int));
    if (NULL == dec_picture->tone_mapping_lut)
    {
      no_mem_exit("init_picture: tone_mapping_lut");
    }
    memcpy(dec_picture->tone_mapping_lut, p_Vid->seiToneMapping->lut, sizeof(imgpel) * coded_data_bit_max);
    update_tone_mapping_sei(p_Vid->seiToneMapping);
  }
  else
    dec_picture->seiHasTone_mapping = 0;
#endif

  if(sps->separate_colour_plane_flag)
  {
    p_Vid->dec_picture_JV[0] = p_Vid->dec_picture;
    p_Vid->dec_picture_JV[1] = alloc_storable_picture (p_Vid, (PictureStructure) currSlice->structure, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr, 1);
    copy_dec_picture_JV( p_Vid, p_Vid->dec_picture_JV[1], p_Vid->dec_picture_JV[0] );
    p_Vid->dec_picture_JV[2] = alloc_storable_picture (p_Vid, (PictureStructure) currSlice->structure, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr, 1);
    copy_dec_picture_JV( p_Vid, p_Vid->dec_picture_JV[2], p_Vid->dec_picture_JV[0] );
  }
}


static void Error_tracking(VideoParameters *p_Vid, Slice *currSlice)
{
    if (currSlice->redundant_pic_cnt == 0)
        p_Vid->Is_primary_correct = p_Vid->Is_redundant_correct = 1;

    if (currSlice->redundant_pic_cnt == 0 && p_Vid->type != I_SLICE) {
        for (int i = 0; i < currSlice->num_ref_idx_l0_active_minus1 + 1 ; i++) {
            if (currSlice->ref_flag[i] == 0)  // any reference of primary slice is incorrect
                p_Vid->Is_primary_correct = 0; // primary slice is incorrect
        }
    } else if (currSlice->redundant_pic_cnt != 0 && p_Vid->type != I_SLICE) {
        int redundant_slice_ref_idx = currSlice->abs_diff_pic_num_minus1[0][0] + 1;
        if (currSlice->ref_flag[redundant_slice_ref_idx] == 0)  // reference of redundant slice is incorrect
            p_Vid->Is_redundant_correct = 0;  // redundant slice is incorrect
    }
}

static Slice *malloc_slice(InputParameters *p_Inp, VideoParameters *p_Vid)
{
    int i, j, memory_size = 0;
    Slice *currSlice;

    currSlice = (Slice *)calloc(1, sizeof(Slice));
    if (!currSlice) {
        snprintf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", p_Inp->FileFormat);
        error(errortext,100);
    }

    // create all context models
    currSlice->mot_ctx = create_contexts_MotionInfo();
    currSlice->tex_ctx = create_contexts_TextureInfo();

    currSlice->max_part_nr = 3;  //! assume data partitioning (worst case) for the following mallocs()
    currSlice->partArr = AllocPartition(currSlice->max_part_nr);

    memory_size += get_mem3Dint(&(currSlice->wp_weight), 2, MAX_REFERENCE_PICTURES, 3);
    memory_size += get_mem3Dint(&(currSlice->wp_offset), 6, MAX_REFERENCE_PICTURES, 3);
    memory_size += get_mem4Dint(&(currSlice->wbp_weight), 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);

    memory_size += get_mem3Dpel(&(currSlice->mb_pred), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dpel(&(currSlice->mb_rec ), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&(currSlice->mb_rres), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&(currSlice->cof    ), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    memory_size += get_mem2Dpel(&currSlice->tmp_block_l0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dpel(&currSlice->tmp_block_l1, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dpel(&currSlice->tmp_block_l2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dpel(&currSlice->tmp_block_l3, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dint(&currSlice->tmp_res, MB_BLOCK_SIZE + 5, MB_BLOCK_SIZE + 5);

#if (MVC_EXTENSION_ENABLE)
    currSlice->view_id = MVC_INIT_VIEW_ID;
    currSlice->inter_view_flag = 0;
    currSlice->anchor_pic_flag = 0;
#endif
    // reference flag initialization
    for (i = 0; i < 17; i++)
        currSlice->ref_flag[i] = 1;
    for (i = 0; i < 6; i++) {
        currSlice->listX[i] = (StorablePicture **)calloc(MAX_LIST_SIZE, sizeof (StorablePicture*)); // +1 for reordering
        if (!currSlice->listX[i])
            no_mem_exit("malloc_slice: currSlice->listX[i]");
    }
    for (j = 0; j < 6; j++) {
        for (i = 0; i < MAX_LIST_SIZE; i++)
            currSlice->listX[j][i] = NULL;
        currSlice->listXsize[j]=0;
    }

    return currSlice;
}

static void copy_slice_info(Slice *currSlice, OldSliceParams *p_old_slice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;

    p_old_slice->pps_id         = currSlice->pic_parameter_set_id;
    p_old_slice->frame_num      = currSlice->frame_num;
    p_old_slice->field_pic_flag = currSlice->field_pic_flag;

    if (currSlice->field_pic_flag)
        p_old_slice->bottom_field_flag = currSlice->bottom_field_flag;

    p_old_slice->nal_ref_idc = currSlice->nal_ref_idc;
    p_old_slice->idr_flag    = (byte) currSlice->idr_flag;

    if (currSlice->idr_flag)
        p_old_slice->idr_pic_id = currSlice->idr_pic_id;

    if (p_Vid->active_sps->pic_order_cnt_type == 0) {
        p_old_slice->pic_oder_cnt_lsb          = currSlice->pic_order_cnt_lsb;
        p_old_slice->delta_pic_oder_cnt_bottom = currSlice->delta_pic_order_cnt_bottom;
    }

    if (p_Vid->active_sps->pic_order_cnt_type == 1) {
        p_old_slice->delta_pic_order_cnt[0] = currSlice->delta_pic_order_cnt[0];
        p_old_slice->delta_pic_order_cnt[1] = currSlice->delta_pic_order_cnt[1];
    }
#if (MVC_EXTENSION_ENABLE)
    p_old_slice->view_id = currSlice->view_id;
    p_old_slice->inter_view_flag = currSlice->inter_view_flag; 
    p_old_slice->anchor_pic_flag = currSlice->anchor_pic_flag;
#endif
    p_old_slice->layer_id = currSlice->layer_id;
}

/*!
 ***********************************************************************
 * \brief
 *    decodes one I- or P-frame
 *
 ***********************************************************************
 */
int decode_one_frame(DecoderParams *pDecoder)
{
    VideoParameters *p_Vid = pDecoder->p_Vid;
    InputParameters *p_Inp = p_Vid->p_Inp;
    int current_header, iRet;
    Slice *currSlice;
    Slice **ppSliceList = p_Vid->ppSliceList;
    //read one picture first;
    current_header = 0; 
    p_Vid->iSliceNumOfCurrPic = 0;
    p_Vid->num_dec_mb = 0;

    if (p_Vid->newframe) {
        if (p_Vid->pNextPPS->Valid) {
            MakePPSavailable(p_Vid, p_Vid->pNextPPS->pic_parameter_set_id, p_Vid->pNextPPS);
            p_Vid->pNextPPS->Valid = 0;
        }

        //get the first slice from currentslice;
        assert(ppSliceList[p_Vid->iSliceNumOfCurrPic]);
        currSlice = p_Vid->pNextSlice;
        p_Vid->pNextSlice = ppSliceList[p_Vid->iSliceNumOfCurrPic];
        ppSliceList[p_Vid->iSliceNumOfCurrPic] = currSlice;
        assert(currSlice->current_slice_nr == 0);

        UseParameterSet(currSlice);

        init_picture(p_Vid, currSlice, p_Inp);

        p_Vid->iSliceNumOfCurrPic++;
        current_header = SOS;
    }

    while (current_header != SOP && current_header != EOS) {
        //no pending slices;
        assert(p_Vid->iSliceNumOfCurrPic < p_Vid->iNumOfSlicesAllocated);
        if (!ppSliceList[p_Vid->iSliceNumOfCurrPic])
            ppSliceList[p_Vid->iSliceNumOfCurrPic] = malloc_slice(p_Inp, p_Vid);
        currSlice = ppSliceList[p_Vid->iSliceNumOfCurrPic];
        currSlice->p_Vid = p_Vid;
        currSlice->p_Inp = p_Inp;
        currSlice->p_Dpb = p_Vid->p_Dpb_layer[0]; //set default value;

        current_header = read_new_slice(currSlice);
        //init;
        currSlice->current_header = current_header;

        // error tracking of primary and redundant slices.
        Error_tracking(p_Vid, currSlice);
        // If primary and redundant are received and primary is correct, discard the redundant
        // else, primary slice will be replaced with redundant slice.
        if (currSlice->frame_num == p_Vid->previous_frame_num &&
            currSlice->redundant_pic_cnt != 0 &&
            p_Vid->Is_primary_correct != 0 && current_header != EOS)
            continue;

        if ((current_header != SOP && current_header != EOS) ||
            (p_Vid->iSliceNumOfCurrPic == 0 && current_header == SOP)) {
            currSlice->current_slice_nr = (short) p_Vid->iSliceNumOfCurrPic;
            p_Vid->dec_picture->max_slice_id = (short) imax(currSlice->current_slice_nr, p_Vid->dec_picture->max_slice_id);
            if (p_Vid->iSliceNumOfCurrPic > 0) {
                currSlice->framepoc  = (*ppSliceList)->framepoc;
                currSlice->toppoc    = (*ppSliceList)->toppoc;
                currSlice->bottompoc = (*ppSliceList)->bottompoc;  
                currSlice->ThisPOC   = (*ppSliceList)->ThisPOC;
            }
            p_Vid->iSliceNumOfCurrPic++;
            if (p_Vid->iSliceNumOfCurrPic >= p_Vid->iNumOfSlicesAllocated) {
                Slice **tmpSliceList = (Slice **)realloc(p_Vid->ppSliceList, (p_Vid->iNumOfSlicesAllocated+MAX_NUM_DECSLICES)*sizeof(Slice*));
                if (!tmpSliceList) {
                    tmpSliceList = (Slice **)calloc((p_Vid->iNumOfSlicesAllocated+MAX_NUM_DECSLICES), sizeof(Slice*));
                    memcpy(tmpSliceList, p_Vid->ppSliceList, p_Vid->iSliceNumOfCurrPic*sizeof(Slice*));
                    //free;
                    free(p_Vid->ppSliceList);
                    ppSliceList = p_Vid->ppSliceList = tmpSliceList;
                } else {
                    ppSliceList = p_Vid->ppSliceList = tmpSliceList;
                    memset(p_Vid->ppSliceList+p_Vid->iSliceNumOfCurrPic, 0, sizeof(Slice*)*MAX_NUM_DECSLICES);
                }
                p_Vid->iNumOfSlicesAllocated += MAX_NUM_DECSLICES;
            }
            current_header = SOS;       
        } else {
            p_Vid->newframe = 1;
            currSlice->current_slice_nr = 0;
            //keep it in currentslice;
            ppSliceList[p_Vid->iSliceNumOfCurrPic] = p_Vid->pNextSlice;
            p_Vid->pNextSlice = currSlice; 
        }

        copy_slice_info(currSlice, p_Vid->old_slice);
    }
    iRet = current_header;

    decode_picture(p_Vid);
    return (iRet);
}

/*!
 ************************************************************************
 * \brief
 *    Convert file read buffer to source picture structure
 * \param imgX
 *    Pointer to image plane
 * \param buf
 *    Buffer for file output
 * \param size_x
 *    horizontal image size in pixel
 * \param size_y
 *    vertical image size in pixel
 * \param symbol_size_in_bytes
 *    number of bytes used per pel
 ************************************************************************
 */
static void buffer2img (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes)
{
  int i,j;

  uint16 tmp16, ui16;
  unsigned long  tmp32, ui32;

  if (symbol_size_in_bytes> sizeof(imgpel))
  {
    error ("Source picture has higher bit depth than imgpel data type. \nPlease recompile with larger data type for imgpel.", 500);
  }

  if (( sizeof(char) == sizeof (imgpel)) && ( sizeof(char) == symbol_size_in_bytes))
  {
    // imgpel == pixel_in_file == 1 byte -> simple copy
    memcpy(&imgX[0][0], buf, size_x * size_y);
  }
  else
  {
    // sizeof (imgpel) > sizeof(char)
    if (testEndian())
    {
      // big endian
      switch (symbol_size_in_bytes)
      {
      case 1:
        {
          for(j = 0; j < size_y; ++j)
            for(i = 0; i < size_x; ++i)
            {
              imgX[j][i]= buf[i+j*size_x];
            }
          break;
        }
      case 2:
        {
          for(j=0;j<size_y;++j)
            for(i=0;i<size_x;++i)
            {
              memcpy(&tmp16, buf+((i+j*size_x)*2), 2);
              ui16  = (uint16) ((tmp16 >> 8) | ((tmp16&0xFF)<<8));
              imgX[j][i] = (imgpel) ui16;
            }
          break;
        }
      case 4:
        {
          for(j=0;j<size_y;++j)
            for(i=0;i<size_x;++i)
            {
              memcpy(&tmp32, buf+((i+j*size_x)*4), 4);
              ui32  = ((tmp32&0xFF00)<<8) | ((tmp32&0xFF)<<24) | ((tmp32&0xFF0000)>>8) | ((tmp32&0xFF000000)>>24);
              imgX[j][i] = (imgpel) ui32;
            }
        }
      default:
        {
           error ("reading only from formats of 8, 16 or 32 bit allowed on big endian architecture", 500);
           break;
        }
      }

    }
    else
    {
      // little endian
      if (symbol_size_in_bytes == 1)
      {
        for (j=0; j < size_y; ++j)
        {
          for (i=0; i < size_x; ++i)
          {
            imgX[j][i]=*(buf++);
          }
        }
      }
      else
      {
        for (j=0; j < size_y; ++j)
        {
          int jpos = j*size_x;
          for (i=0; i < size_x; ++i)
          {
            imgX[j][i]=0;
            memcpy(&(imgX[j][i]), buf +((i+jpos)*symbol_size_in_bytes), symbol_size_in_bytes);
          }
        }
      }

    }
  }
}


/*!
 ***********************************************************************
 * \brief
 *    compute generic SSE
 ***********************************************************************
 */
static int64 compute_SSE(imgpel **imgRef, imgpel **imgSrc, int xRef, int xSrc, int ySize, int xSize)
{
  int i, j;
  imgpel *lineRef, *lineSrc;
  int64 distortion = 0;

  for (j = 0; j < ySize; j++)
  {
    lineRef = &imgRef[j][xRef];    
    lineSrc = &imgSrc[j][xSrc];

    for (i = 0; i < xSize; i++)
      distortion += iabs2( *lineRef++ - *lineSrc++ );
  }
  return distortion;
}

/*!
 ************************************************************************
 * \brief
 *    Calculate the value of frame_no
 ************************************************************************
*/
void calculate_frame_no(VideoParameters *p_Vid, StorablePicture *p)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  // calculate frame number
  int  psnrPOC = p_Vid->active_sps->mb_adaptive_frame_field_flag ? p->poc /(p_Inp->poc_scale) : p->poc/(p_Inp->poc_scale);
  
  if (psnrPOC==0)// && p_Vid->psnr_number)
  {
    p_Vid->idr_psnr_number = p_Vid->g_nFrame * p_Vid->ref_poc_gap/(p_Inp->poc_scale);
  }
  p_Vid->psnr_number = imax(p_Vid->psnr_number, p_Vid->idr_psnr_number+psnrPOC);

  p_Vid->frame_no = p_Vid->idr_psnr_number + psnrPOC;
}


/*!
************************************************************************
* \brief
*    Find PSNR for all three components.Compare decoded frame with
*    the original sequence. Read p_Inp->jumpd frames to reflect frame skipping.
* \param p_Vid
*      video encoding parameters for current picture
* \param p
*      picture to be compared
* \param p_ref
*      file pointer piont to reference YUV reference file
************************************************************************
*/
void find_snr(VideoParameters *p_Vid, 
              StorablePicture *p,
              int *p_ref)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  SNRParameters   *snr   = p_Vid->snr;
  sps_t *sps = p_Vid->active_sps;

  int k;
  int ret;
  int64 diff_comp[3] = {0};
  int64  status;
  int symbol_size_in_bytes = (p_Vid->pic_unit_bitsize_on_disk >> 3);
  int comp_size_x[3], comp_size_y[3];
  int64 framesize_in_bytes;

  Boolean rgb_output = (Boolean) (p_Vid->active_sps->vui_parameters.matrix_coefficients==0);
  unsigned char *buf;
  imgpel **cur_ref [3];
  imgpel **cur_comp[3]; 
  // picture error concealment
  char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};

  cur_ref[0]  = p_Vid->imgY_ref;
  cur_ref[1]  = p->chroma_format_idc != YUV400 ? p_Vid->imgUV_ref[0] : NULL;
  cur_ref[2]  = p->chroma_format_idc != YUV400 ? p_Vid->imgUV_ref[1] : NULL;

  cur_comp[0] = p->imgY;
  cur_comp[1] = p->chroma_format_idc != YUV400 ? p->imgUV[0]  : NULL;
  cur_comp[2] =  p->chroma_format_idc!= YUV400 ? p->imgUV[1]  : NULL; 

  comp_size_x[0] = p_Inp->source.width[0];
  comp_size_y[0] = p_Inp->source.height[0];
  comp_size_x[1] = comp_size_x[2] = p_Inp->source.width[1];
  comp_size_y[1] = comp_size_y[2] = p_Inp->source.height[1];

  framesize_in_bytes = (((int64) comp_size_x[0] * comp_size_y[0]) + ((int64) comp_size_x[1] * comp_size_y[1] ) * 2) * symbol_size_in_bytes;

  // KS: this buffer should actually be allocated only once, but this is still much faster than the previous version
  buf = (unsigned char *)malloc ( comp_size_x[0] * comp_size_y[0] * symbol_size_in_bytes );

  if (NULL == buf)
  {
    no_mem_exit("find_snr: buf");
  }

  status = lseek (*p_ref, framesize_in_bytes * p_Vid->frame_no, SEEK_SET);
  if (status == -1)
  {
    fprintf(stderr, "Warning: Could not seek to frame number %d in reference file. Shown PSNR might be wrong.\n", p_Vid->frame_no);
    free (buf);
    return;
  }

  if(rgb_output)
    lseek (*p_ref, framesize_in_bytes/3, SEEK_CUR);

  for (k = 0; k < ((p->chroma_format_idc != YUV400) ? 3 : 1); ++k)
  {

    if(rgb_output && k == 2)
      lseek (*p_ref, -framesize_in_bytes, SEEK_CUR);

    ret = read(*p_ref, buf, comp_size_x[k] * comp_size_y[k] * symbol_size_in_bytes);
    if (ret != comp_size_x[k] * comp_size_y[k] * symbol_size_in_bytes)
    {
      printf ("Warning: could not read from reconstructed file\n");
      memset (buf, 0, comp_size_x[k] * comp_size_y[k] * symbol_size_in_bytes);
      close(*p_ref);
      *p_ref = -1;
      break;
    }
    buffer2img(cur_ref[k], buf, comp_size_x[k], comp_size_y[k], symbol_size_in_bytes);

    // Compute SSE
    diff_comp[k] = compute_SSE(cur_ref[k], cur_comp[k], 0, 0, comp_size_y[k], comp_size_x[k]);

    int max_pix_value_sqd = iabs2((1 << (k > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1);
    // Collecting SNR statistics
    snr->snr[k] = psnr( max_pix_value_sqd, comp_size_x[k] * comp_size_y[k], (float) diff_comp[k]);   

    if (snr->frame_ctr == 0) // first
    {
      snr->snra[k] = snr->snr[k];                                                        // keep luma snr for first frame
    }
    else
    {
      snr->snra[k] = (float)(snr->snra[k]*(snr->frame_ctr)+snr->snr[k])/(snr->frame_ctr + 1); // average snr chroma for all frames
    }
  }

  if(rgb_output)
    lseek (*p_ref, framesize_in_bytes * 2 / 3, SEEK_CUR);

  free (buf);

  // picture error concealment
  if(p->concealed_pic)
  {
    fprintf(stdout,"%04d(P)  %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      p_Vid->frame_no, p->frame_poc, p->pic_num, p->qp,
      snr->snr[0], snr->snr[1], snr->snr[2], yuv_types[p->chroma_format_idc], 0);
  }
}









#if (MVC_EXTENSION_ENABLE)
int GetVOIdx(VideoParameters *p_Vid, int iViewId)
{
  int iVOIdx = -1;
  int *piViewIdMap;
  if(p_Vid->active_subset_sps)
  {
    piViewIdMap = p_Vid->active_subset_sps->view_id;
    for(iVOIdx = p_Vid->active_subset_sps->num_views_minus1; iVOIdx>=0; iVOIdx--)
      if(piViewIdMap[iVOIdx] == iViewId)
        break;
  }
  else
  {
    sub_sps_t *curr_subset_sps;
    int i;

    curr_subset_sps = p_Vid->SubsetSeqParSet;
    for(i=0; i<MAXSPS; i++)
    {
      if(curr_subset_sps->num_views_minus1>=0 && curr_subset_sps->sps.Valid)
      {
        break;
      }
      curr_subset_sps++;
    }

    if( i < MAXSPS )
    {
      p_Vid->active_subset_sps = curr_subset_sps;

      piViewIdMap = p_Vid->active_subset_sps->view_id;
      for(iVOIdx = p_Vid->active_subset_sps->num_views_minus1; iVOIdx>=0; iVOIdx--)
        if(piViewIdMap[iVOIdx] == iViewId)
          break;

      return iVOIdx;
    }
    else
    {
      iVOIdx = 0;
    }
  }

  return iVOIdx;
}

#endif
