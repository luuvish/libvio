#include "global.h"
#include "input_parameters.h"
#include "h264decoder.h"
#include "report.h"

#include "slice.h"
#include "data_partition.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"

#include "sei.h"
#include "output.h"
#include "memalloc.h"
#include "macroblock.h"
#include "neighbour.h"

using namespace vio::h264;

#include "erc_api.h"


#if (MVC_EXTENSION_ENABLE)
static int GetBaseViewId(VideoParameters* p_Vid, sub_sps_t** subset_sps)
{
    sub_sps_t* curr_subset_sps;
    int i, iBaseViewId = 0; //-1;

    *subset_sps = NULL;
    curr_subset_sps = p_Vid->SubsetSeqParSet;
    for (i = 0; i < MAXSPS; i++) {
        if (curr_subset_sps->num_views_minus1 >= 0 && curr_subset_sps->sps.Valid) {
            iBaseViewId = curr_subset_sps->view_id[0];
            break;
        }
        curr_subset_sps++;
    }

    if (i < MAXSPS)
        *subset_sps = curr_subset_sps;
    return iBaseViewId;
}

int GetVOIdx(VideoParameters *p_Vid, int iViewId)
{
    int iVOIdx = -1;
    if (p_Vid->active_subset_sps) {
        int* piViewIdMap = p_Vid->active_subset_sps->view_id;
        for (iVOIdx = p_Vid->active_subset_sps->num_views_minus1; iVOIdx >= 0; iVOIdx--) {
            if (piViewIdMap[iVOIdx] == iViewId)
                break;
        }
    } else {
        int i;
        sub_sps_t* curr_subset_sps = p_Vid->SubsetSeqParSet;
        for (i = 0; i < MAXSPS; i++) {
            if (curr_subset_sps->num_views_minus1 >= 0 && curr_subset_sps->sps.Valid)
                break;
            curr_subset_sps++;
        }

        if (i < MAXSPS) {
            p_Vid->active_subset_sps = curr_subset_sps;
            int* piViewIdMap = p_Vid->active_subset_sps->view_id;
            for (iVOIdx = p_Vid->active_subset_sps->num_views_minus1; iVOIdx >= 0; iVOIdx--)
                if (piViewIdMap[iVOIdx] == iViewId)
                    break;

            return iVOIdx;
        } else
            iVOIdx = 0;
    }

    return iVOIdx;
}
#endif


static void Error_tracking(VideoParameters *p_Vid, slice_t *currSlice)
{
    shr_t& shr = currSlice->header;

    if (shr.redundant_pic_cnt == 0)
        p_Vid->Is_primary_correct = p_Vid->Is_redundant_correct = 1;

    if (shr.redundant_pic_cnt == 0 && p_Vid->type != I_slice) {
        for (int i = 0; i < shr.num_ref_idx_l0_active_minus1 + 1 ; i++) {
            if (currSlice->ref_flag[i] == 0)  // any reference of primary slice is incorrect
                p_Vid->Is_primary_correct = 0; // primary slice is incorrect
        }
    } else if (shr.redundant_pic_cnt != 0 && p_Vid->type != I_slice) {
        int redundant_slice_ref_idx = shr.abs_diff_pic_num_minus1[0][0] + 1;
        if (currSlice->ref_flag[redundant_slice_ref_idx] == 0)  // reference of redundant slice is incorrect
            p_Vid->Is_redundant_correct = 0;  // redundant slice is incorrect
    }
}


static int parse_idr(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    nalu_t *nalu = p_Vid->nalu; 
    int current_header = 0;

    if (p_Vid->recovery_point || nalu->nal_unit_type == NALU_TYPE_IDR) {
        if (!p_Vid->recovery_point_found) {
            if (nalu->nal_unit_type != NALU_TYPE_IDR) {
                printf("Warning: Decoding does not start with an IDR picture.\n");
                p_Vid->non_conforming_stream = 1;
            } else
                p_Vid->non_conforming_stream = 0;
        }
        p_Vid->recovery_point_found = true;
    }

    if (!p_Vid->recovery_point_found)
        return current_header;

    currSlice->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
    currSlice->nal_ref_idc = nalu->nal_ref_idc;
    currSlice->parser.dp_mode = PAR_DP_1;
#if (MVC_EXTENSION_ENABLE)
    if (currSlice->svc_extension_flag != 0)
#endif
        currSlice->parser.partArr[0].init(nalu);

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
    currSlice->layer_id = currSlice->view_id = GetVOIdx(p_Vid, currSlice->view_id);
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

    currSlice->decoder.assign_quant_params(*currSlice);

    shr_t& shr = currSlice->header;

    // if primary slice is replaced with redundant slice, set the correct image type
    if (shr.redundant_pic_cnt && p_Vid->Is_primary_correct == 0 && p_Vid->Is_redundant_correct)
        p_Vid->dec_picture->slice.slice_type = p_Vid->type;

    if (!p_Vid->dec_picture || *(p_Vid->old_slice) != *currSlice) {
        if (p_Vid->iSliceNumOfCurrPic == 0)
            init_picture(currSlice);
    //if (*(p_Vid->old_slice) != *currSlice) {
        current_header = SOP;
        //check zero_byte if it is also the first NAL unit in the access unit
        p_Vid->bitstream.CheckZeroByteVCL(nalu);
    } else
        current_header = SOS;

    p_Vid->recovery_point = false;
    return current_header;
}

static int parse_dpa(slice_t *currSlice)
{
    VideoParameters* p_Vid = currSlice->p_Vid;
    nalu_t* nalu = p_Vid->nalu;
    data_partition_t* dp;
    int current_header = 0;

    int slice_id_a, slice_id_b, slice_id_c;

    if (!p_Vid->recovery_point_found)
        return current_header;

    // read DP_A
    currSlice->dpB_NotPresent = 1;
    currSlice->dpC_NotPresent = 1;

    currSlice->idr_flag    = 0;
    currSlice->nal_ref_idc = nalu->nal_ref_idc;
    currSlice->parser.dp_mode = PAR_DP_3;
#if MVC_EXTENSION_ENABLE
    currSlice->p_Dpb = p_Vid->p_Dpb_layer[0];
#endif
    dp = &currSlice->parser.partArr[0];
    dp->init(nalu);
#if MVC_EXTENSION_ENABLE
    currSlice->view_id = GetBaseViewId(p_Vid, &p_Vid->active_subset_sps);
    currSlice->inter_view_flag = 1;
    currSlice->layer_id = currSlice->view_id = GetVOIdx(p_Vid, currSlice->view_id);
    currSlice->anchor_pic_flag = currSlice->idr_flag;
#endif

    slice_header(currSlice);
#if MVC_EXTENSION_ENABLE
    currSlice->p_Dpb = p_Vid->p_Dpb_layer[currSlice->view_id];
#endif

    currSlice->decoder.assign_quant_params(*currSlice);

    if (!p_Vid->dec_picture || *(p_Vid->old_slice) != *currSlice) {
        if (p_Vid->iSliceNumOfCurrPic == 0)
            init_picture(currSlice);
    //if (*(p_Vid->old_slice) != *currSlice) {
        current_header = SOP;
        //check zero_byte if it is also the first NAL unit in the access unit
        p_Vid->bitstream.CheckZeroByteVCL(nalu);
    } else
        current_header = SOS;

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
        dp = &currSlice->parser.partArr[1];
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
        dp = &currSlice->parser.partArr[2];
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

static int read_new_slice(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    InputParameters *p_Inp = currSlice->p_Vid->p_Inp;

    nalu_t *nalu = p_Vid->nalu; 
    int current_header = 0;
    data_partition_t* dp = &currSlice->parser.partArr[0];

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
            if (!p_Inp->silent)
                printf ("found data partition B without matching DP A, discarding\n");
            break;

        case NALU_TYPE_DPC:
            if (!p_Inp->silent)
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
                if (!p_Inp->silent)
                    printf ("Found Subsequence SPS NALU. Ignoring.\n");
            }
            break;

        case NALU_TYPE_SLC_EXT:
            if (p_Inp->DecodeAllLayers == 0 && !p_Inp->silent)
                printf ("Found SVC extension NALU (%d). Ignoring.\n", (int) nalu->nal_unit_type);
            break;
#endif

        default:
            if (!p_Inp->silent)
                printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", (int) nalu->nal_unit_type, (int) nalu->len);
            break;
        }
    }
}

int DecoderParams::decode_slice_headers()
{
    VideoParameters *p_Vid = this->p_Vid;
    int current_header = 0;
    slice_t *currSlice;
    slice_t **ppSliceList = p_Vid->ppSliceList;
    //read one picture first;
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

        init_picture(currSlice);

        p_Vid->iSliceNumOfCurrPic++;
        current_header = SOS;
        //p_Vid->newframe = 0;
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

        shr_t& shr = currSlice->header;

        // error tracking of primary and redundant slices.
        Error_tracking(p_Vid, currSlice);
        // If primary and redundant are received and primary is correct, discard the redundant
        // else, primary slice will be replaced with redundant slice.
        if (shr.frame_num == p_Vid->previous_frame_num &&
            shr.redundant_pic_cnt != 0 &&
            p_Vid->Is_primary_correct != 0 && current_header != EOS)
            continue;

        if ((current_header == SOS) || (current_header == SOP && p_Vid->iSliceNumOfCurrPic == 0)) {
            currSlice->current_slice_nr = (short) p_Vid->iSliceNumOfCurrPic;
            if (p_Vid->iSliceNumOfCurrPic > 0) {
                shr.PicOrderCnt         = (*ppSliceList)->header.PicOrderCnt;
                shr.TopFieldOrderCnt    = (*ppSliceList)->header.TopFieldOrderCnt;
                shr.BottomFieldOrderCnt = (*ppSliceList)->header.BottomFieldOrderCnt;  
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
        }
        if (current_header == SOP && p_Vid->iSliceNumOfCurrPic > 0) {
            p_Vid->newframe = 1;
            currSlice->current_slice_nr = 0;
            //keep it in currentslice;
            ppSliceList[p_Vid->iSliceNumOfCurrPic] = p_Vid->pNextSlice;
            p_Vid->pNextSlice = currSlice; 
        }

        *(p_Vid->old_slice) = *currSlice;
    }
    
    return current_header;
}
