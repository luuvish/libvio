#include "global.h"
#include "input_parameters.h"
#include "report.h"

#include "slice.h"
#include "image.h"
#include "fmo.h"
#include "data_partition.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"

#include "sei.h"
#include "output.h"
#include "memalloc.h"
#include "macroblock.h"
#include "neighbour.h"

using vio::h264::mb_t;

#include "erc_api.h"
#include "dpb.h"


static inline int is_BL_profile(unsigned int profile_idc) 
{
  return ( profile_idc == FREXT_CAVLC444 || profile_idc == BASELINE || profile_idc == MAIN || profile_idc == EXTENDED ||
           profile_idc == FREXT_HP || profile_idc == FREXT_Hi10P || profile_idc == FREXT_Hi422 || profile_idc == FREXT_Hi444);
}

static inline int is_EL_profile(unsigned int profile_idc) 
{
  return ( (profile_idc == MVC_HIGH) || (profile_idc == STEREO_HIGH)
           );
}

static inline int is_MVC_profile(unsigned int profile_idc)
{
  return ( (0)
#if (MVC_EXTENSION_ENABLE)
  || (profile_idc == MVC_HIGH) || (profile_idc == STEREO_HIGH)
#endif
  );
}


static int init_global_buffers(VideoParameters *p_Vid, int layer_id)
{
    int memory_size=0;
    int i;
    CodingParameters *cps = p_Vid->p_EncodePar[layer_id];
    sps_t *sps = p_Vid->active_sps;
    int FrameSizeInMbs = sps->PicWidthInMbs * sps->FrameHeightInMbs;

    if (p_Vid->global_init_done[layer_id])
        free_layer_buffers(p_Vid, layer_id);

    // allocate memory in structure p_Vid
    if (sps->separate_colour_plane_flag) {
        for (i = 0; i < 3; i++)
            cps->mb_data_JV[i] = new mb_t[FrameSizeInMbs];
        cps->mb_data = NULL;
    } else {
        cps->mb_data = new mb_t[FrameSizeInMbs];
    }

    p_Vid->global_init_done[layer_id] = 1;

    return memory_size;
}

static void setup_layer_info(VideoParameters *p_Vid, sps_t *sps, LayerParameters *p_Lps)
{
    int layer_id = p_Lps->layer_id;
    p_Lps->p_Vid = p_Vid;
    p_Lps->p_Cps = p_Vid->p_EncodePar[layer_id];
    p_Lps->p_SPS = sps;
    p_Lps->p_Dpb = p_Vid->p_Dpb_layer[layer_id];
}

void activate_sps(VideoParameters *p_Vid, sps_t *sps)
{
    if (p_Vid->active_sps != sps) {
        int prev_profile_idc = p_Vid->active_sps ? p_Vid->active_sps->profile_idc : 0;

        if (p_Vid->dec_picture)
            // this may only happen on slice loss
            exit_picture(p_Vid, &p_Vid->dec_picture);
        p_Vid->active_sps = sps;

        if (p_Vid->dpb_layer_id == 0 && is_BL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[0]->init_done) {
            setup_layer_info(p_Vid, sps, p_Vid->p_LayerPar[0]);
        } else if (p_Vid->dpb_layer_id == 1 && is_EL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[1]->init_done) {
            setup_layer_info(p_Vid, sps, p_Vid->p_LayerPar[1]);
        }

#if (MVC_EXTENSION_ENABLE)
        if (prev_profile_idc != sps->profile_idc) {
            if (is_BL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[0]->init_done) {
                init_global_buffers(p_Vid, 0);

                if (!p_Vid->no_output_of_prior_pics_flag) {
                    flush_dpb(p_Vid->p_Dpb_layer[0]);
                    flush_dpb(p_Vid->p_Dpb_layer[1]);
                }
                init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);

            } else if ((is_MVC_profile(prev_profile_idc) ||
                        is_MVC_profile(sps->profile_idc)) &&
                       !p_Vid->p_Dpb_layer[1]->init_done) {

                assert(p_Vid->p_Dpb_layer[0]->init_done);
                if (p_Vid->p_Dpb_layer[0]->init_done) {
                    free_dpb(p_Vid->p_Dpb_layer[0]);
                    init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);
                }

                init_global_buffers(p_Vid, 1);
                // for now lets re_init both buffers. Later, we should only re_init appropriate one
                // Note that we seem to be doing this for every frame which seems not good.
                init_dpb(p_Vid, p_Vid->p_Dpb_layer[1], 2);
            }
        }
#endif

#if (DISABLE_ERC == 0)
        ercInit(p_Vid, sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16, 1);
        if (p_Vid->dec_picture) {
            slice_t *currSlice = p_Vid->ppSliceList[0];
            ercReset(p_Vid->erc_errorVar, currSlice->PicSizeInMbs, currSlice->PicSizeInMbs, p_Vid->dec_picture->size_x);
            p_Vid->erc_mvperMB = 0;
        }
#endif
    }
}

void activate_pps(VideoParameters *p_Vid, pps_t *pps)
{  
    if (p_Vid->active_pps != pps) {
        if (p_Vid->dec_picture) {
            // this may only happen on slice loss
            exit_picture(p_Vid, &p_Vid->dec_picture);
        }

        p_Vid->active_pps = pps;
    }
}



void UseParameterSet (slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    int PicParsetId = currSlice->pic_parameter_set_id;  
    pps_t *pps = &p_Vid->PicParSet[PicParsetId];
    sps_t *sps = &p_Vid->SeqParSet[pps->seq_parameter_set_id];

    if (!pps->Valid)
        printf ("Trying to use an invalid (uninitialized) Picture Parameter Set with ID %d, expect the unexpected...\n", PicParsetId);
#if (MVC_EXTENSION_ENABLE)
    if (currSlice->svc_extension_flag == -1) {
        if (!sps->Valid)
            printf ("PicParset %d references an invalid (uninitialized) Sequence Parameter Set with ID %d, expect the unexpected...\n", 
        PicParsetId, (int) pps->seq_parameter_set_id);
    } else {
        // Set SPS to the subset SPS parameters
        p_Vid->active_subset_sps = p_Vid->SubsetSeqParSet + pps->seq_parameter_set_id;
        sps = &(p_Vid->active_subset_sps->sps);
        if (!p_Vid->SubsetSeqParSet[pps->seq_parameter_set_id].Valid)
            printf ("PicParset %d references an invalid (uninitialized) Subset Sequence Parameter Set with ID %d, expect the unexpected...\n", 
                    PicParsetId, (int) pps->seq_parameter_set_id);
    }
#endif

    // In theory, and with a well-designed software, the lines above
    // are everything necessary.  In practice, we need to patch many values
    // in p_Vid-> (but no more in p_Inp-> -- these have been taken care of)

    // Set Sequence Parameter Stuff first
    if ((int) sps->pic_order_cnt_type < 0 || sps->pic_order_cnt_type > 2) {
        printf("invalid sps->pic_order_cnt_type = %d\n", (int) sps->pic_order_cnt_type);
        error("pic_order_cnt_type != 1", -1000);
    }

    if (sps->pic_order_cnt_type == 1) {
        if ((int)sps->num_ref_frames_in_pic_order_cnt_cycle >= MAX_NUM_REF_FRAMES)
            error("num_ref_frames_in_pic_order_cnt_cycle too large",-1011);
    }
    p_Vid->dpb_layer_id = currSlice->layer_id;
    activate_sps(p_Vid, sps);
    activate_pps(p_Vid, pps);

    p_Vid->type = currSlice->slice_type;
}

static void init_picture_decoding(VideoParameters *p_Vid)
{
    slice_t *pSlice = p_Vid->ppSliceList[0];

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

    p_Vid->structure = pSlice->structure;

    fmo_init(p_Vid, pSlice);

#if (MVC_EXTENSION_ENABLE)
    if (pSlice->layer_id > 0 && pSlice->svc_extension_flag == 0 && pSlice->NaluHeaderMVCExt.non_idr_flag == 0)
        idr_memory_management(p_Vid->p_Dpb_layer[pSlice->layer_id], p_Vid->dec_picture);
    update_ref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_ltref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_pic_num(pSlice);
#endif
}


void decode_picture(VideoParameters *p_Vid)
{
    slice_t **ppSliceList = p_Vid->ppSliceList;
    int iSliceNo;

    p_Vid->num_dec_mb = 0;

    init_picture_decoding(p_Vid);

    for (iSliceNo = 0; iSliceNo < p_Vid->iSliceNumOfCurrPic; iSliceNo++) {
        slice_t *currSlice = ppSliceList[iSliceNo];

        assert(currSlice->current_header != EOS);
        assert(currSlice->current_slice_nr == iSliceNo);

        if (!currSlice->init())
            continue;
        currSlice->decode();

        p_Vid->num_dec_mb  += currSlice->num_dec_mb;
        p_Vid->erc_mvperMB += currSlice->erc_mvperMB;
    }

    exit_picture(p_Vid, &p_Vid->dec_picture);

    p_Vid->previous_frame_num = ppSliceList[0]->frame_num;
}


void pad_buf(imgpel *pImgBuf, int iWidth, int iHeight, int iStride, int iPadX, int iPadY)
{
    int j;
    imgpel *pLine0 = pImgBuf - iPadX, *pLine;
    int i;
    for (i = -iPadX; i < 0; i++)
        pImgBuf[i] = *pImgBuf;
    for (i = 0; i < iPadX; i++)
        pImgBuf[i+iWidth] = *(pImgBuf+iWidth-1);

    for (j = -iPadY; j < 0; j++)
        memcpy(pLine0+j*iStride, pLine0, iStride*sizeof(imgpel));
    for (j = 1; j < iHeight; j++) {
        pLine = pLine0 + j*iStride;
        for (i = 0; i < iPadX; i++)
            pLine[i] = pLine[iPadX];
        pLine += iPadX+iWidth-1;
        for (i = 1; i < iPadX + 1; i++)
            pLine[i] = *pLine;
    }
    pLine = pLine0 + (iHeight - 1) * iStride;
    for (j = iHeight; j < iHeight + iPadY; j++)
        memcpy(pLine0+j*iStride,  pLine, iStride*sizeof(imgpel));
}

void pad_dec_picture(VideoParameters *p_Vid, storable_picture *dec_picture)
{
    int iPadX = MCBUF_LUMA_PAD_X;
    int iPadY = MCBUF_LUMA_PAD_Y;
    int iWidth = dec_picture->size_x;
    int iHeight = dec_picture->size_y;
    int iStride = dec_picture->iLumaStride;

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (dec_picture->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (dec_picture->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

    pad_buf(*dec_picture->imgY, iWidth, iHeight, iStride, iPadX, iPadY);

    if (dec_picture->chroma_format_idc != YUV400) {
        iPadX = iChromaPadX;
        iPadY = iChromaPadY;
        iWidth = dec_picture->size_x_cr;
        iHeight = dec_picture->size_y_cr;
        iStride = dec_picture->iChromaStride;
        pad_buf(*dec_picture->imgUV[0], iWidth, iHeight, iStride, iPadX, iPadY);
        pad_buf(*dec_picture->imgUV[1], iWidth, iHeight, iStride, iPadX, iPadY);
    }
}

void exit_picture(VideoParameters *p_Vid, storable_picture **dec_picture)
{
    sps_t *sps = p_Vid->active_sps;
    slice_t *currSlice = p_Vid->ppSliceList[0];
    int PicSizeInMbs = sps->PicWidthInMbs * (sps->FrameHeightInMbs / (1 + currSlice->field_pic_flag));

    // return if the last picture has already been finished
    if (*dec_picture == NULL ||
        (p_Vid->num_dec_mb != PicSizeInMbs &&
         (sps->chroma_format_idc != YUV444 || !sps->separate_colour_plane_flag)))
        return;

#if (DISABLE_ERC == 0)
    erc_picture(p_Vid, dec_picture);
#endif

    currSlice->decoder.deblock_filter(*currSlice);

    if (p_Vid->structure != FRAME)
        p_Vid->number /= 2;
#if (MVC_EXTENSION_ENABLE)
    if ((*dec_picture)->used_for_reference || ((*dec_picture)->inter_view_flag == 1))
        pad_dec_picture(p_Vid, *dec_picture);
#endif
#if MVC_EXTENSION_ENABLE
    store_picture_in_dpb(p_Vid->p_Dpb_layer[(*dec_picture)->view_id], *dec_picture);
#endif

    if (p_Vid->last_has_mmco_5)
        p_Vid->pre_frame_num = 0;

    p_Vid->status(dec_picture);

    *dec_picture = NULL;
}
