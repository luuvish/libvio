
#include <math.h>
#include <limits.h>

#include "global.h"
#include "image.h"
#include "fmo.h"
#include "data_partition.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"
#include "slice.h"

#include "sei.h"
#include "output.h"
#include "neighbour.h"
#include "memalloc.h"
#include "macroblock.h"

#include "deblock.h"

#include "quantization.h"

#include "erc_api.h"
#include "dpb.h"


using namespace arrow::video::h264;


/*!
 ************************************************************************
 * \brief
 *    detect if current slice is "first VCL NAL unit of a picture"
 ************************************************************************
 */
static int is_new_picture(StorablePicture *dec_picture, slice_t *currSlice, OldSliceParams *p_old_slice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;

    int result = 0;

    result |= (NULL==dec_picture);
    result |= (p_old_slice->pps_id != currSlice->pic_parameter_set_id);
    result |= (p_old_slice->frame_num != currSlice->frame_num);
    result |= (p_old_slice->field_pic_flag != currSlice->field_pic_flag);

    if (currSlice->field_pic_flag && p_old_slice->field_pic_flag)
        result |= (p_old_slice->bottom_field_flag != currSlice->bottom_field_flag);

    result |= (p_old_slice->nal_ref_idc != currSlice->nal_ref_idc) && ((p_old_slice->nal_ref_idc == 0) || (currSlice->nal_ref_idc == 0));
    result |= (p_old_slice->idr_flag    != currSlice->idr_flag);

    if (currSlice->idr_flag && p_old_slice->idr_flag)
        result |= (p_old_slice->idr_pic_id != currSlice->idr_pic_id);

    if (p_Vid->active_sps->pic_order_cnt_type == 0) {
        result |= (p_old_slice->pic_oder_cnt_lsb          != currSlice->pic_order_cnt_lsb);
        if (p_Vid->active_pps->bottom_field_pic_order_in_frame_present_flag  ==  1 &&  !currSlice->field_pic_flag )
            result |= (p_old_slice->delta_pic_oder_cnt_bottom != currSlice->delta_pic_order_cnt_bottom);
    }

    if (p_Vid->active_sps->pic_order_cnt_type == 1) {
        if (!p_Vid->active_sps->delta_pic_order_always_zero_flag) {
            result |= (p_old_slice->delta_pic_order_cnt[0] != currSlice->delta_pic_order_cnt[0]);
            if (p_Vid->active_pps->bottom_field_pic_order_in_frame_present_flag  ==  1 &&  !currSlice->field_pic_flag )
                result |= (p_old_slice->delta_pic_order_cnt[1] != currSlice->delta_pic_order_cnt[1]);
        }
    }

#if (MVC_EXTENSION_ENABLE)
    result |= (currSlice->view_id != p_old_slice->view_id);
    result |= (currSlice->inter_view_flag != p_old_slice->inter_view_flag);
    result |= (currSlice->anchor_pic_flag != p_old_slice->anchor_pic_flag);
#endif
    result |= (currSlice->layer_id != p_old_slice->layer_id);
    return result;
}


static int parse_idr(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    InputParameters *p_Inp = currSlice->p_Inp;
    nalu_t *nalu = p_Vid->nalu; 
    int current_header = 0;

    if (p_Vid->recovery_point || nalu->nal_unit_type == NALU_TYPE_IDR) {
        if (p_Vid->recovery_point_found == 0) {
            if (nalu->nal_unit_type != NALU_TYPE_IDR) {
                printf("Warning: Decoding does not start with an IDR picture.\n");
                p_Vid->non_conforming_stream = 1;
            } else
                p_Vid->non_conforming_stream = 0;
        }
        p_Vid->recovery_point_found = 1;
    }

    if (p_Vid->recovery_point_found == 0)
        return current_header;

    currSlice->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
    currSlice->nal_ref_idc = nalu->nal_ref_idc;
    currSlice->dp_mode = PAR_DP_1;
    currSlice->max_part_nr = 1;
#if (MVC_EXTENSION_ENABLE)
    if (currSlice->svc_extension_flag != 0)
#endif
        currSlice->partArr[0].init(nalu);

#if (MVC_EXTENSION_ENABLE)
    if (currSlice->svc_extension_flag == 0) {
        currSlice->view_id         = currSlice->NaluHeaderMVCExt.view_id;
        currSlice->inter_view_flag = currSlice->NaluHeaderMVCExt.inter_view_flag;
        currSlice->anchor_pic_flag = currSlice->NaluHeaderMVCExt.anchor_pic_flag;
    } else if (currSlice->svc_extension_flag == -1) { //SVC and the normal AVC;
        if (p_Vid->active_subset_sps == NULL) {
            currSlice->view_id = GetBaseViewId(p_Vid, &p_Vid->active_subset_sps);
            if (currSlice->NaluHeaderMVCExt.iPrefixNALU > 0) {
                assert(currSlice->view_id == currSlice->NaluHeaderMVCExt.view_id);
                currSlice->inter_view_flag = currSlice->NaluHeaderMVCExt.inter_view_flag;
                currSlice->anchor_pic_flag = currSlice->NaluHeaderMVCExt.anchor_pic_flag;
            } else {
                currSlice->inter_view_flag = 1;
                currSlice->anchor_pic_flag = currSlice->idr_flag;
            }
        } else {
            assert(p_Vid->active_subset_sps->num_views_minus1 >=0);
            // prefix NALU available
            if (currSlice->NaluHeaderMVCExt.iPrefixNALU > 0) {
                currSlice->view_id = currSlice->NaluHeaderMVCExt.view_id;
                currSlice->inter_view_flag = currSlice->NaluHeaderMVCExt.inter_view_flag;
                currSlice->anchor_pic_flag = currSlice->NaluHeaderMVCExt.anchor_pic_flag;
            } else { //no prefix NALU;
                currSlice->view_id = p_Vid->active_subset_sps->view_id[0];
                currSlice->inter_view_flag = 1;
                currSlice->anchor_pic_flag = currSlice->idr_flag;
            }
        }
    }
    currSlice->layer_id = currSlice->view_id = GetVOIdx( p_Vid, currSlice->view_id );
#endif

    // Some syntax of the slice_t Header depends on the parameter set, which depends on
    // the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
    // of the slice header first, then setup the active parameter sets, and then read
    // the rest of the slice header
    slice_header(currSlice);
#if (MVC_EXTENSION_ENABLE)
    if (currSlice->view_id >= 0)
        currSlice->p_Dpb = p_Vid->p_Dpb_layer[currSlice->view_id];
#endif

    assign_quant_params (currSlice);        

    // if primary slice is replaced with redundant slice, set the correct image type
    if (currSlice->redundant_pic_cnt && p_Vid->Is_primary_correct == 0 && p_Vid->Is_redundant_correct)
        p_Vid->dec_picture->slice_type = p_Vid->type;

    if (is_new_picture(p_Vid->dec_picture, currSlice, p_Vid->old_slice)) {
        if (p_Vid->iSliceNumOfCurrPic==0)
            init_picture(p_Vid, currSlice, p_Inp);

        current_header = SOP;
        //check zero_byte if it is also the first NAL unit in the access unit
        p_Vid->bitstream.CheckZeroByteVCL(nalu);
    } else
        current_header = SOS;

    // From here on, p_Vid->active_sps, p_Vid->active_pps and the slice header are valid
    if (currSlice->MbaffFrameFlag)
        currSlice->current_mb_nr = currSlice->first_mb_in_slice << 1;
    else
        currSlice->current_mb_nr = currSlice->first_mb_in_slice;

    if (p_Vid->active_pps->entropy_coding_mode_flag)
        currSlice->partArr[0].de_cabac.init(&currSlice->partArr[0]);
    p_Vid->recovery_point = 0;
    return current_header;
}

static int parse_dpa(slice_t *currSlice)
{
    VideoParameters* p_Vid = currSlice->p_Vid;
    InputParameters* p_Inp = currSlice->p_Inp;
    nalu_t* nalu = p_Vid->nalu; 
    data_partition_t* dp;
    int current_header = 0;

    int slice_id_a, slice_id_b, slice_id_c;

    if (p_Vid->recovery_point_found == 0)
        return current_header;

    // read DP_A
    currSlice->dpB_NotPresent = 1;
    currSlice->dpC_NotPresent = 1;

    currSlice->idr_flag    = FALSE;
    currSlice->nal_ref_idc = nalu->nal_ref_idc;
    currSlice->dp_mode     = PAR_DP_3;
    currSlice->max_part_nr = 3;
#if MVC_EXTENSION_ENABLE
    currSlice->p_Dpb = p_Vid->p_Dpb_layer[0];
#endif
    dp = &currSlice->partArr[0];
    dp->init(nalu);
#if MVC_EXTENSION_ENABLE
    currSlice->view_id = GetBaseViewId(p_Vid, &p_Vid->active_subset_sps);
    currSlice->inter_view_flag = 1;
    currSlice->layer_id = currSlice->view_id = GetVOIdx( p_Vid, currSlice->view_id );
    currSlice->anchor_pic_flag = currSlice->idr_flag;
#endif

    slice_header(currSlice);
#if MVC_EXTENSION_ENABLE
    currSlice->p_Dpb = p_Vid->p_Dpb_layer[currSlice->view_id];
#endif

    assign_quant_params (currSlice);        


    if (is_new_picture(p_Vid->dec_picture, currSlice, p_Vid->old_slice)) {
        if (p_Vid->iSliceNumOfCurrPic==0)
            init_picture(p_Vid, currSlice, p_Inp);

        current_header = SOP;
        //check zero_byte if it is also the first NAL unit in the access unit
        p_Vid->bitstream.CheckZeroByteVCL(nalu);
    } else
        current_header = SOS;

    // From here on, p_Vid->active_sps, p_Vid->active_pps and the slice header are valid
    if (currSlice->MbaffFrameFlag)
        currSlice->current_mb_nr = currSlice->first_mb_in_slice << 1;
    else
        currSlice->current_mb_nr = currSlice->first_mb_in_slice;

    // Now I need to read the slice ID, which depends on the value of
    // redundant_pic_cnt_present_flag

    slice_id_a = dp->ue("NALU: DP_A slice_id");

    if (p_Vid->active_pps->entropy_coding_mode_flag)
        error ("received data partition with CABAC, this is not allowed", 500);

    // continue with reading next DP
    if (0 == p_Vid->bitstream.read_next_nalu(nalu))
        return current_header;

    if ( NALU_TYPE_DPB == nalu->nal_unit_type) {
        // we got a DPB
        dp = &currSlice->partArr[1];
        dp->init(nalu);

        slice_id_b = dp->ue("NALU: DP_B slice_id");

        currSlice->dpB_NotPresent = 0; 

        if ((slice_id_b != slice_id_a) || (nalu->lost_packets)) {
            printf ("Waning: got a data partition B which does not match DP_A (DP loss!)\n");
            currSlice->dpB_NotPresent = 1;
            currSlice->dpC_NotPresent = 1;
        } else {
            if (p_Vid->active_pps->redundant_pic_cnt_present_flag)
                dp->ue("NALU: DP_B redundant_pic_cnt");

            // we're finished with DP_B, so let's continue with next DP
            if (0 == p_Vid->bitstream.read_next_nalu(nalu))
                return current_header;
        }
    } else
        currSlice->dpB_NotPresent = 1;

    // check if we got DP_C
    if ( NALU_TYPE_DPC == nalu->nal_unit_type) {
        dp = &currSlice->partArr[2];
        dp->init(nalu);

        currSlice->dpC_NotPresent = 0;

        slice_id_c = dp->ue("NALU: DP_C slice_id");
        if ((slice_id_c != slice_id_a)|| (nalu->lost_packets)) {
            printf ("Warning: got a data partition C which does not match DP_A(DP loss!)\n");
            currSlice->dpC_NotPresent =1;
        }

        if (p_Vid->active_pps->redundant_pic_cnt_present_flag)
            dp->ue("NALU:SLICE_C redudand_pic_cnt");
    } else
        currSlice->dpC_NotPresent = 1;

    // check if we read anything else than the expected partitions
    if ((nalu->nal_unit_type != NALU_TYPE_DPB) && (nalu->nal_unit_type != NALU_TYPE_DPC)) {
        // we have a NALI that we can't process here, so restart processing
        return 100;
        // yes, "goto" should not be used, but it's really the best way here before we restructure the decoding loop
        // (which should be taken care of anyway)
    }

    return current_header;
}


/*!
 ************************************************************************
 * \brief
 *    Reads new slice from data_partition_t
 ************************************************************************
 */
int read_new_slice(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    InputParameters *p_Inp = currSlice->p_Inp;

    nalu_t *nalu = p_Vid->nalu; 
    int current_header = 0;
    data_partition_t* dp = &currSlice->partArr[0];

    currSlice->num_dec_mb        = 0;
    currSlice->is_reset_coeff    = 0;
    currSlice->is_reset_coeff_cr = 0;

    for (;;) {
#if (MVC_EXTENSION_ENABLE)
        currSlice->svc_extension_flag = -1;
#endif
        if (0 == p_Vid->bitstream.read_next_nalu(nalu))
            return EOS;

#if (MVC_EXTENSION_ENABLE)
        if (p_Inp->DecodeAllLayers == 1 &&
            (nalu->nal_unit_type == NALU_TYPE_PREFIX || nalu->nal_unit_type == NALU_TYPE_SLC_EXT)) {
            dp->init(nalu);

            currSlice->svc_extension_flag = dp->u(1, "svc_extension_flag");

            if (currSlice->svc_extension_flag)
                nal_unit_header_svc_extension();
            else {
                nal_unit_header_mvc_extension(&currSlice->NaluHeaderMVCExt, dp);
                currSlice->NaluHeaderMVCExt.iPrefixNALU = (nalu->nal_unit_type == NALU_TYPE_PREFIX);
            }

            if (nalu->nal_unit_type == NALU_TYPE_SLC_EXT) {        
                if (currSlice->svc_extension_flag) {
                    //to be implemented for Annex G;
                } else
                    nalu->nal_unit_type = NALU_TYPE_SLICE;
            }
        }
#endif

process_nalu:
        switch (nalu->nal_unit_type) {
        case NALU_TYPE_SLICE:
        case NALU_TYPE_IDR:
            current_header = parse_idr(currSlice);
            if (current_header != 0)
                return current_header;
            break;

        case NALU_TYPE_DPA:
            current_header = parse_dpa(currSlice);
            if (current_header == 100)
                goto process_nalu;
            if (current_header != 0)
                return current_header;
            break;

        case NALU_TYPE_DPB:
            if (p_Inp->silent == FALSE)
                printf ("found data partition B without matching DP A, discarding\n");
            break;

        case NALU_TYPE_DPC:
            if (p_Inp->silent == FALSE)
                printf ("found data partition C without matching DP A, discarding\n");
            break;

        case NALU_TYPE_SEI:
            parse_sei(nalu->buf, nalu->len, p_Vid, currSlice);
            break;

        case NALU_TYPE_PPS:
            ProcessPPS(p_Vid, nalu);
            break;

        case NALU_TYPE_SPS:
            ProcessSPS(p_Vid, nalu);
            break;

        case NALU_TYPE_AUD:
            break;

        case NALU_TYPE_EOSEQ:
            break;

        case NALU_TYPE_EOSTREAM:
            break;

        case NALU_TYPE_FILL:
            break;

#if (MVC_EXTENSION_ENABLE)
        case NALU_TYPE_VDRD:
            break;

        case NALU_TYPE_PREFIX:
            if (currSlice->svc_extension_flag==1)
                prefix_nal_unit_svc();
            break;

        case NALU_TYPE_SUB_SPS:
            if (p_Inp->DecodeAllLayers== 1)
                ProcessSubsetSPS(p_Vid, nalu);
            else {
                if (p_Inp->silent == FALSE)
                    printf ("Found Subsequence SPS NALU. Ignoring.\n");
            }
            break;

        case NALU_TYPE_SLC_EXT:
            if (p_Inp->DecodeAllLayers == 0 &&  (p_Inp->silent == FALSE))
                printf ("Found SVC extension NALU (%d). Ignoring.\n", (int) nalu->nal_unit_type);
            break;
#endif

        default:
            if (p_Inp->silent == FALSE)
                printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", (int) nalu->nal_unit_type, (int) nalu->len);
            break;
        }
    }
}
