#include "global.h"
#include "input_parameters.h"
#include "h264decoder.h"
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
#include "ref_list.h"


static inline float psnr(int max_sample_sq, int samples, float sse_distortion ) 
{
  return (float) (sse_distortion == 0.0 ? 0.0 : (10.0 * log10(max_sample_sq * (double) ((double) samples / sse_distortion))));
}

static inline void reset_mbs(mb_t *currMB)
{
    currMB->slice_nr = -1;
    currMB->ei_flag  =  1;
    currMB->dpl_flag =  0;
}

static void setup_buffers(VideoParameters* p_Vid, int layer_id)
{
    CodingParameters* cps = p_Vid->p_EncodePar[layer_id];

    if (p_Vid->last_dec_layer_id != layer_id) {
        if (p_Vid->active_sps->separate_colour_plane_flag) {
            for (int i = 0; i < 3; ++i)
                p_Vid->mb_data_JV[i] = cps->mb_data_JV[i];
            p_Vid->mb_data = NULL;
        } else
            p_Vid->mb_data = cps->mb_data;
        p_Vid->last_dec_layer_id = layer_id;
    }
}

#if MVC_EXTENSION_ENABLE
static void init_mvc_picture(slice_t* currSlice)
{
    VideoParameters* p_Vid = currSlice->p_Vid;
    dpb_t* p_Dpb = p_Vid->p_Dpb_layer[0];

    storable_picture* p_pic = NULL;

    // find BL reconstructed picture
    if (!currSlice->field_pic_flag) {
        for (int i = 0; i < (int)p_Dpb->used_size/*size*/; ++i) {
            frame_store* fs = p_Dpb->fs[i];
            if (fs->frame->slice.view_id == 0 && fs->frame->frame_poc == currSlice->PicOrderCnt) {
                p_pic = fs->frame;
                break;
            }
        }
    } else if (!currSlice->bottom_field_flag) {
        for (int i = 0; i < (int)p_Dpb->used_size/*size*/; ++i) {
            frame_store* fs = p_Dpb->fs[i];
            if (fs->top_field->slice.view_id == 0 && fs->top_field->top_poc == currSlice->TopFieldOrderCnt) {
                p_pic = fs->top_field;
                break;
            }
        }
    } else {
        for (int i = 0; i < (int)p_Dpb->used_size/*size*/; ++i) {
            frame_store* fs = p_Dpb->fs[i];
            if (fs->bottom_field->slice.view_id == 0 && fs->bottom_field->bottom_poc == currSlice->BottomFieldOrderCnt) {
                p_pic = fs->bottom_field;
                break;
            }
        }
    }

    if (p_pic)
        picture_in_dpb(currSlice, p_Vid, p_pic);
}
#endif


static void copy_dec_picture_JV(VideoParameters *p_Vid, storable_picture *dst, storable_picture *src)
{
    dst->top_poc                         = src->top_poc;
    dst->bottom_poc                      = src->bottom_poc;
    dst->frame_poc                       = src->frame_poc;

    dst->poc                             = src->poc;
    dst->used_for_reference              = src->used_for_reference;

    dst->slice.coded_frame                     = src->slice.coded_frame;
    dst->slice.slice_type                      = src->slice.slice_type;
    dst->slice.idr_flag                        = src->slice.idr_flag;
    dst->slice.mb_aff_frame_flag               = src->slice.mb_aff_frame_flag;
    dst->slice.no_output_of_prior_pics_flag    = src->slice.no_output_of_prior_pics_flag;
    dst->slice.long_term_reference_flag        = src->slice.long_term_reference_flag;
    dst->slice.adaptive_ref_pic_buffering_flag = src->slice.adaptive_ref_pic_buffering_flag;
    dst->slice.dec_ref_pic_marking_buffer      = src->slice.dec_ref_pic_marking_buffer;

    dst->PicWidthInMbs                   = src->PicWidthInMbs;
    dst->PicNum                          = src->PicNum;
    dst->frame_num                       = src->frame_num;
    dst->recovery_frame                  = src->recovery_frame;

    // store the necessary tone mapping sei into storable_picture structure
    dst->seiHasTone_mapping    = src->seiHasTone_mapping;
    dst->tone_mapping_model_id = src->tone_mapping_model_id;
    dst->tonemapped_bit_depth  = src->tonemapped_bit_depth;
    if (src->tone_mapping_lut) {
        int coded_data_bit_max = (1 << p_Vid->seiToneMapping->coded_data_bit_depth);
        dst->tone_mapping_lut = (imgpel *)malloc(sizeof(int) * coded_data_bit_max);
        if (NULL == dst->tone_mapping_lut)
            no_mem_exit("copy_dec_picture_JV: tone_mapping_lut");
        memcpy(dst->tone_mapping_lut, src->tone_mapping_lut, sizeof(imgpel) * coded_data_bit_max);
    }
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new picture
 ************************************************************************
 */
void init_picture(VideoParameters *p_Vid, slice_t *currSlice, InputParameters *p_Inp)
{
    int i;
    int nplane;
    storable_picture* dec_picture = NULL;
    sps_t* sps = p_Vid->active_sps;
    dpb_t* p_Dpb = currSlice->p_Dpb;

    if (p_Vid->dec_picture)
        // this may only happen on slice loss
        exit_picture(p_Vid, &p_Vid->dec_picture);
    p_Vid->dpb_layer_id = currSlice->layer_id;
    //set buffers;
    setup_buffers(p_Vid, currSlice->layer_id);

    if (p_Vid->recovery_point)
        p_Vid->recovery_frame_num = (currSlice->frame_num + p_Vid->recovery_frame_cnt) % sps->MaxFrameNum;

    if (currSlice->idr_flag)
        p_Vid->recovery_frame_num = currSlice->frame_num;

    if (p_Vid->recovery_point == 0 &&
        currSlice->frame_num != p_Vid->pre_frame_num &&
        currSlice->frame_num != (p_Vid->pre_frame_num + 1) % sps->MaxFrameNum) {
        if (sps->gaps_in_frame_num_value_allowed_flag == 0) {
#if (DISABLE_ERC == 0)
            // picture error concealment
            if (p_Inp->conceal_mode != 0) {
                if (currSlice->frame_num < (p_Vid->pre_frame_num + 1) % sps->MaxFrameNum) {
                    /* Conceal lost IDR frames and any frames immediately
                       following the IDR. Use frame copy for these since
                       lists cannot be formed correctly for motion copy*/
                    p_Vid->conceal_mode = 1;
                    p_Vid->IDR_concealment_flag = 1;
                    conceal_lost_frames(p_Dpb, currSlice);
                    //reset to original concealment mode for future drops
                    p_Vid->conceal_mode = p_Inp->conceal_mode;
                } else {
                    //reset to original concealment mode for future drops
                    p_Vid->conceal_mode = p_Inp->conceal_mode;

                    p_Vid->IDR_concealment_flag = 0;
                    conceal_lost_frames(p_Dpb, currSlice);
                }
            } else
#endif
                /* Advanced Error Concealment would be called here to combat unintentional loss of pictures. */
                error("An unintentional loss of pictures occurs! Exit\n", 100);
        }
        if (p_Vid->conceal_mode == 0)
            fill_frame_num_gap(p_Vid, currSlice);
    }

    if (currSlice->nal_ref_idc)
        p_Vid->pre_frame_num = currSlice->frame_num;

    //calculate POC
    decode_poc(p_Vid, currSlice);

    if (p_Vid->recovery_frame_num == (int) currSlice->frame_num && p_Vid->recovery_poc == 0x7fffffff)
        p_Vid->recovery_poc = currSlice->PicOrderCnt;

    if (currSlice->nal_ref_idc)
        p_Vid->last_ref_pic_poc = currSlice->PicOrderCnt;

    if (!currSlice->field_pic_flag || !currSlice->bottom_field_flag)
        p_Vid->snr->start_time = std::chrono::system_clock::now();

    dec_picture = p_Vid->dec_picture = new storable_picture(p_Vid, currSlice->structure,
        sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
        sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);
    dec_picture->top_poc     = currSlice->TopFieldOrderCnt;
    dec_picture->bottom_poc  = currSlice->BottomFieldOrderCnt;
    dec_picture->frame_poc   = currSlice->PicOrderCnt;
    dec_picture->slice.iCodingType = !currSlice->field_pic_flag ?
                                     (currSlice->MbaffFrameFlag ? FRAME_MB_PAIR_CODING : FRAME_CODING) : FIELD_CODING;
    dec_picture->slice.layer_id    = currSlice->layer_id;
#if (MVC_EXTENSION_ENABLE)
    dec_picture->slice.view_id         = currSlice->view_id;
    dec_picture->slice.inter_view_flag = currSlice->inter_view_flag;
    dec_picture->slice.anchor_pic_flag = currSlice->anchor_pic_flag;
    if (dec_picture->slice.view_id == 1) {
        if ((p_Vid->profile_idc == MVC_HIGH) || (p_Vid->profile_idc == STEREO_HIGH))
            init_mvc_picture(currSlice);
    }
#endif

    // reset all variables of the error concealment instance before decoding of every frame.
    // here the third parameter should, if perfectly, be equal to the number of slices per frame.
    // using little value is ok, the code will allocate more memory if the slice number is larger
#if (DISABLE_ERC == 0)
    ercReset(p_Vid->erc_errorVar, currSlice->PicSizeInMbs, currSlice->PicSizeInMbs, dec_picture->size_x);
#endif
    p_Vid->erc_mvperMB = 0;

    if (!currSlice->field_pic_flag)
        dec_picture->poc = currSlice->PicOrderCnt;
    else if (!currSlice->bottom_field_flag) {
        dec_picture->poc = currSlice->TopFieldOrderCnt;
        p_Vid->number *= 2;
    } else {
        dec_picture->poc = currSlice->BottomFieldOrderCnt;
        p_Vid->number = p_Vid->number * 2 + 1;
    }

    if (p_Vid->type > SI_slice)
        p_Vid->type = P_slice;  // concealed element

    // TO set mb_t Map (mark all MBs as 'have to be concealed')
    if (sps->separate_colour_plane_flag) {
        for (nplane = 0; nplane < 3; ++nplane) {
            mb_t* currMB = p_Vid->mb_data_JV[nplane];
            for (i = 0; i < currSlice->PicSizeInMbs; ++i)
                reset_mbs(currMB++);
        }
    } else {
        mb_t* currMB = p_Vid->mb_data;
        for (i = 0; i < currSlice->PicSizeInMbs; ++i)
            reset_mbs(currMB++);
    }

    dec_picture->used_for_reference              = currSlice->nal_ref_idc != 0;

    dec_picture->slice.coded_frame                     = !currSlice->field_pic_flag;
    dec_picture->slice.idr_flag                        = currSlice->idr_flag;
    dec_picture->slice.slice_type                      = p_Vid->type;
    dec_picture->slice.no_output_of_prior_pics_flag    = currSlice->no_output_of_prior_pics_flag;
    dec_picture->slice.long_term_reference_flag        = currSlice->long_term_reference_flag;
    dec_picture->slice.adaptive_ref_pic_buffering_flag = currSlice->adaptive_ref_pic_marking_mode_flag;
    dec_picture->slice.dec_ref_pic_marking_buffer      = currSlice->dec_ref_pic_marking_buffer;
    dec_picture->slice.mb_aff_frame_flag               = currSlice->MbaffFrameFlag;
    currSlice->dec_ref_pic_marking_buffer        = NULL;

    dec_picture->PicWidthInMbs                   = sps->PicWidthInMbs;

    dec_picture->PicNum                          = currSlice->frame_num;
    dec_picture->frame_num                       = currSlice->frame_num;
    dec_picture->recovery_frame                  = (unsigned int) ((int) currSlice->frame_num == p_Vid->recovery_frame_num);

    // store the necessary tone mapping sei into storable_picture structure
    if (p_Vid->seiToneMapping->seiHasTone_mapping) {
        int coded_data_bit_max = (1 << p_Vid->seiToneMapping->coded_data_bit_depth);
        dec_picture->seiHasTone_mapping    = 1;
        dec_picture->tone_mapping_model_id = p_Vid->seiToneMapping->model_id;
        dec_picture->tonemapped_bit_depth  = p_Vid->seiToneMapping->sei_bit_depth;
        dec_picture->tone_mapping_lut      = (imgpel *)malloc(coded_data_bit_max * sizeof(int));
        if (NULL == dec_picture->tone_mapping_lut)
            no_mem_exit("init_picture: tone_mapping_lut");
        memcpy(dec_picture->tone_mapping_lut, p_Vid->seiToneMapping->lut, sizeof(imgpel) * coded_data_bit_max);
        update_tone_mapping_sei(p_Vid->seiToneMapping);
    } else
        dec_picture->seiHasTone_mapping = 0;

    if (sps->separate_colour_plane_flag) {
        p_Vid->dec_picture_JV[0] = p_Vid->dec_picture;
        p_Vid->dec_picture_JV[1] = new storable_picture(
            p_Vid, (PictureStructure) currSlice->structure,
            sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
            sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);
        copy_dec_picture_JV( p_Vid, p_Vid->dec_picture_JV[1], p_Vid->dec_picture_JV[0] );
        p_Vid->dec_picture_JV[2] = new storable_picture(
            p_Vid, (PictureStructure) currSlice->structure,
            sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
            sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);
        copy_dec_picture_JV( p_Vid, p_Vid->dec_picture_JV[2], p_Vid->dec_picture_JV[0] );
    }
}


static void Error_tracking(VideoParameters *p_Vid, slice_t *currSlice)
{
    if (currSlice->redundant_pic_cnt == 0)
        p_Vid->Is_primary_correct = p_Vid->Is_redundant_correct = 1;

    if (currSlice->redundant_pic_cnt == 0 && p_Vid->type != I_slice) {
        for (int i = 0; i < currSlice->num_ref_idx_l0_active_minus1 + 1 ; i++) {
            if (currSlice->ref_flag[i] == 0)  // any reference of primary slice is incorrect
                p_Vid->Is_primary_correct = 0; // primary slice is incorrect
        }
    } else if (currSlice->redundant_pic_cnt != 0 && p_Vid->type != I_slice) {
        int redundant_slice_ref_idx = currSlice->abs_diff_pic_num_minus1[0][0] + 1;
        if (currSlice->ref_flag[redundant_slice_ref_idx] == 0)  // reference of redundant slice is incorrect
            p_Vid->Is_redundant_correct = 0;  // redundant slice is incorrect
    }
}


int DecoderParams::decode_one_frame()
{
    VideoParameters *p_Vid = this->p_Vid;
    InputParameters *p_Inp = p_Vid->p_Inp;
    int current_header, iRet;
    slice_t *currSlice;
    slice_t **ppSliceList = p_Vid->ppSliceList;
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
            ppSliceList[p_Vid->iSliceNumOfCurrPic] = new slice_t;
        currSlice = ppSliceList[p_Vid->iSliceNumOfCurrPic];
        currSlice->p_Vid = p_Vid;
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
            if (p_Vid->iSliceNumOfCurrPic > 0) {
                currSlice->PicOrderCnt  = (*ppSliceList)->PicOrderCnt;
                currSlice->TopFieldOrderCnt    = (*ppSliceList)->TopFieldOrderCnt;
                currSlice->BottomFieldOrderCnt = (*ppSliceList)->BottomFieldOrderCnt;  
                currSlice->PicOrderCnt   = (*ppSliceList)->PicOrderCnt;
            }
            p_Vid->iSliceNumOfCurrPic++;
            if (p_Vid->iSliceNumOfCurrPic >= p_Vid->iNumOfSlicesAllocated) {
                slice_t **tmpSliceList = (slice_t **)realloc(p_Vid->ppSliceList, (p_Vid->iNumOfSlicesAllocated+MAX_NUM_DECSLICES)*sizeof(slice_t*));
                if (!tmpSliceList) {
                    tmpSliceList = (slice_t **)calloc((p_Vid->iNumOfSlicesAllocated+MAX_NUM_DECSLICES), sizeof(slice_t*));
                    memcpy(tmpSliceList, p_Vid->ppSliceList, p_Vid->iSliceNumOfCurrPic*sizeof(slice_t*));
                    //free;
                    free(p_Vid->ppSliceList);
                    ppSliceList = p_Vid->ppSliceList = tmpSliceList;
                } else {
                    ppSliceList = p_Vid->ppSliceList = tmpSliceList;
                    memset(p_Vid->ppSliceList+p_Vid->iSliceNumOfCurrPic, 0, sizeof(slice_t*)*MAX_NUM_DECSLICES);
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

        *(p_Vid->old_slice) = *currSlice;
    }
    iRet = current_header;

    decode_picture(p_Vid);
    return (iRet);
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
