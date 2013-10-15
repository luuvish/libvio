#include "global.h"
#include "data_partition.h"
#include "dpb.h"
#include "slice.h"
#include "memalloc.h"


bool slice_t::operator!=(const slice_t& slice)
{
    const sps_t& sps = *slice.active_sps;
    const pps_t& pps = *slice.active_pps;
    const shr_t& shr = slice.header;

    bool result = false;

    result |= this->header.pic_parameter_set_id != shr.pic_parameter_set_id;
    result |= this->header.frame_num            != shr.frame_num;
    result |= this->header.field_pic_flag       != shr.field_pic_flag;

    if (shr.field_pic_flag && this->header.field_pic_flag)
        result |= this->header.bottom_field_flag != shr.bottom_field_flag;

    result |= this->nal_ref_idc != slice.nal_ref_idc && (this->nal_ref_idc == 0 || slice.nal_ref_idc == 0);
    result |= this->idr_flag    != slice.idr_flag;

    if (slice.idr_flag && this->idr_flag)
        result |= this->header.idr_pic_id != shr.idr_pic_id;

    if (sps.pic_order_cnt_type == 0) {
        result |= this->header.pic_order_cnt_lsb != shr.pic_order_cnt_lsb;
        if (pps.bottom_field_pic_order_in_frame_present_flag && !shr.field_pic_flag)
            result |= this->header.delta_pic_order_cnt_bottom != shr.delta_pic_order_cnt_bottom;
    }
    if (sps.pic_order_cnt_type == 1) {
        if (!sps.delta_pic_order_always_zero_flag) {
            result |= this->header.delta_pic_order_cnt[0] != shr.delta_pic_order_cnt[0];
            if (pps.bottom_field_pic_order_in_frame_present_flag && !shr.field_pic_flag)
                result |= this->header.delta_pic_order_cnt[1] != shr.delta_pic_order_cnt[1];
        }
    }

#if (MVC_EXTENSION_ENABLE)
    result |= this->view_id         != slice.view_id;
    result |= this->inter_view_flag != slice.inter_view_flag;
    result |= this->anchor_pic_flag != slice.anchor_pic_flag;
#endif
    result |= this->layer_id        != slice.layer_id;

    return result;
}


slice_t::slice_t()
{
    shr_t& shr = this->header;

    get_mem3Dpel(&this->mb_pred, 3, 16, 16);

#if (MVC_EXTENSION_ENABLE)
    this->view_id         = -1;
    this->inter_view_flag = 0;
    this->anchor_pic_flag = 0;
#endif
    // reference flag initialization
    for (int i = 0; i < 17; i++)
        this->ref_flag[i] = 1;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < MAX_LIST_SIZE; i++)
            this->RefPicList[j][i] = NULL;
        this->RefPicSize[j] = 0;
    }

    shr.MbToSliceGroupMap      = nullptr;
    shr.MapUnitToSliceGroupMap = nullptr;
}

slice_t::~slice_t()
{
    shr_t& shr = this->header;

    free_mem3Dpel(this->mb_pred);

    while (shr.dec_ref_pic_marking_buffer) {
        drpm_t* tmp_drpm = shr.dec_ref_pic_marking_buffer;
        shr.dec_ref_pic_marking_buffer = tmp_drpm->Next;
        free(tmp_drpm);
    }

    if (shr.MbToSliceGroupMap)
        delete []shr.MbToSliceGroupMap;
    if (shr.MapUnitToSliceGroupMap)
        delete []shr.MapUnitToSliceGroupMap;
}

void slice_header(slice_t *currSlice)
{
    shr_t& shr = currSlice->header;

    VideoParameters *p_Vid = currSlice->p_Vid;
    data_partition_t *s = &currSlice->parser.partArr[0];

    shr.first_mb_in_slice = s->ue("SH: first_mb_in_slice");
    p_Vid->type = shr.slice_type = s->ue("SH: slice_type") % 5;
    shr.pic_parameter_set_id = s->ue("SH: pic_parameter_set_id");

    assert(shr.pic_parameter_set_id >= 0 && shr.pic_parameter_set_id <= 255);

    UseParameterSet(currSlice);
    sps_t *sps = currSlice->active_sps = p_Vid->active_sps;
    pps_t *pps = currSlice->active_pps = p_Vid->active_pps;

    shr.colour_plane_id = 0;
    if (sps->separate_colour_plane_flag)
        shr.colour_plane_id = s->u(2, "SH: colour_plane_id");

    assert(shr.colour_plane_id >= 0 && shr.colour_plane_id <= 2);

    shr.frame_num = s->u(sps->log2_max_frame_num_minus4 + 4, "SH: frame_num");

    /* Tian Dong: frame_num gap processing, if found */
    if (currSlice->idr_flag) {
        p_Vid->pre_frame_num = shr.frame_num;
        // picture error concealment
        p_Vid->last_ref_pic_poc = 0;
        assert(shr.frame_num == 0);
    }

    p_Vid->structure = FRAME;
    shr.field_pic_flag    = 0;
    shr.bottom_field_flag = 0;
    if (!sps->frame_mbs_only_flag) {
        shr.field_pic_flag = s->u(1, "SH: field_pic_flag");
        if (shr.field_pic_flag) {
            shr.bottom_field_flag = s->u(1, "SH: bottom_field_flag");
            p_Vid->structure = shr.bottom_field_flag ? BOTTOM_FIELD : TOP_FIELD;
        }
    }

    shr.structure = (PictureStructure) p_Vid->structure;

    shr.MbaffFrameFlag     = (sps->mb_adaptive_frame_field_flag && !shr.field_pic_flag);
    shr.PicHeightInMbs     = sps->FrameHeightInMbs / (1 + shr.field_pic_flag);
    shr.PicHeightInSampleL = shr.PicHeightInMbs * 16;
    shr.PicHeightInSampleC = shr.PicHeightInMbs * sps->MbHeightC;
    shr.PicSizeInMbs       = sps->PicWidthInMbs * shr.PicHeightInMbs;
    if (!shr.field_pic_flag) {
        shr.MaxPicNum  = sps->MaxFrameNum;
        shr.CurrPicNum = shr.frame_num;
    } else {
        shr.MaxPicNum  = 2 * sps->MaxFrameNum;
        shr.CurrPicNum = 2 * shr.frame_num + 1;
    }

    if (currSlice->idr_flag)
        shr.idr_pic_id = s->ue("SH: idr_pic_id");
#if (MVC_EXTENSION_ENABLE)
    else if ( currSlice->svc_extension_flag == 0 && currSlice->NaluHeaderMVCExt.non_idr_flag == 0 )
        shr.idr_pic_id = s->ue("SH: idr_pic_id");
#endif

    assert(shr.idr_pic_id >= 0 && shr.idr_pic_id <= 65535);

    shr.delta_pic_order_cnt_bottom = 0;
    shr.delta_pic_order_cnt[0]     = 0;
    shr.delta_pic_order_cnt[1]     = 0;
    if (sps->pic_order_cnt_type == 0) {
        shr.pic_order_cnt_lsb = s->u(sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "SH: pic_order_cnt_lsb");
        if (pps->bottom_field_pic_order_in_frame_present_flag && !shr.field_pic_flag)
            shr.delta_pic_order_cnt_bottom = s->se("SH: delta_pic_order_cnt_bottom");
    }
    if (sps->pic_order_cnt_type == 1 && !sps->delta_pic_order_always_zero_flag) {
        shr.delta_pic_order_cnt[0] = s->se("SH: delta_pic_order_cnt[0]");
        if (pps->bottom_field_pic_order_in_frame_present_flag && !shr.field_pic_flag)
            shr.delta_pic_order_cnt[1] = s->se("SH: delta_pic_order_cnt[1]");
    }

    assert(shr.delta_pic_order_cnt_bottom >= -(1 << 31) + 1 && shr.delta_pic_order_cnt_bottom <=  (1 << 31) - 1);
    assert(shr.delta_pic_order_cnt[0] >= -(1 << 31) + 1 && shr.delta_pic_order_cnt[0] <=  (1 << 31) - 1);
    assert(shr.delta_pic_order_cnt[1] >= -(1 << 31) + 1 && shr.delta_pic_order_cnt[1] <=  (1 << 31) - 1);

    shr.redundant_pic_cnt = 0;
    if (pps->redundant_pic_cnt_present_flag)
        shr.redundant_pic_cnt = s->ue("SH: redundant_pic_cnt");

    assert(shr.redundant_pic_cnt >= 0 && shr.redundant_pic_cnt <= 127);

    if (shr.slice_type == B_slice)
        shr.direct_spatial_mv_pred_flag = s->u(1, "SH: direct_spatial_mv_pred_flag");

    shr.num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    shr.num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;
    if (shr.slice_type == P_slice || shr.slice_type == SP_slice ||
        shr.slice_type == B_slice) {
        shr.num_ref_idx_active_override_flag = s->u(1, "SH: num_ref_idx_override_flag");
        if (shr.num_ref_idx_active_override_flag) {
            shr.num_ref_idx_l0_active_minus1 = s->ue("SH: num_ref_idx_l0_active_minus1");
            if (shr.slice_type == B_slice)
                shr.num_ref_idx_l1_active_minus1 = s->ue("SH: num_ref_idx_l1_active_minus1");
        }
    }

    if ((shr.slice_type == P_slice || shr.slice_type == SP_slice || shr.slice_type == B_slice) &&
        !shr.field_pic_flag && pps->num_ref_idx_l0_default_active_minus1 > 15)
        assert(shr.num_ref_idx_active_override_flag == 1);
    if (shr.slice_type == B_slice &&
        !shr.field_pic_flag && pps->num_ref_idx_l1_default_active_minus1 > 15)
        assert(shr.num_ref_idx_active_override_flag == 1);

    if (!shr.field_pic_flag)
        assert(shr.num_ref_idx_l0_active_minus1 >= 0 && shr.num_ref_idx_l0_active_minus1 <= 15);
    else
        assert(shr.num_ref_idx_l0_active_minus1 >= 0 && shr.num_ref_idx_l0_active_minus1 <= 31);
    if (!shr.field_pic_flag)
        assert(shr.num_ref_idx_l1_active_minus1 >= 0 && shr.num_ref_idx_l1_active_minus1 <= 15);
    else
        assert(shr.num_ref_idx_l1_active_minus1 >= 0 && shr.num_ref_idx_l1_active_minus1 <= 31);

#if (MVC_EXTENSION_ENABLE)
    // if (nal_unit_type == 20 || nal_unit_type == 21)
    if (currSlice->svc_extension_flag == 0 || currSlice->svc_extension_flag == 1)
        ref_pic_list_mvc_modification(currSlice);
    else
#endif
        ref_pic_list_modification(currSlice);

    if ((pps->weighted_pred_flag && (shr.slice_type == P_slice || shr.slice_type == SP_slice)) ||
        (pps->weighted_bipred_idc == 1 && shr.slice_type == B_slice))
        pred_weight_table(currSlice);
    else {
        shr.luma_log2_weight_denom   = 5;
        shr.chroma_log2_weight_denom = 5;
    }

    if (currSlice->nal_ref_idc != 0)
        dec_ref_pic_marking(p_Vid, s, currSlice);

    shr.cabac_init_idc = 0;
    if (pps->entropy_coding_mode_flag && shr.slice_type != I_slice && shr.slice_type != SI_slice)
        shr.cabac_init_idc = s->ue("SH: cabac_init_idc");

    assert(shr.cabac_init_idc >= 0 && shr.cabac_init_idc <= 2);

    shr.slice_qp_delta = s->se("SH: slice_qp_delta");
    shr.SliceQpY = 26 + pps->pic_init_qp_minus26 + shr.slice_qp_delta;

    shr.sp_for_switch_flag = 0;
    shr.QsY = 0;
    if (shr.slice_type == SP_slice || shr.slice_type == SI_slice) {
        if (shr.slice_type == SP_slice)
            shr.sp_for_switch_flag = s->u(1, "SH: sp_for_switch_flag");
        shr.slice_qs_delta = s->se("SH: slice_qs_delta");
        shr.QsY = 26 + pps->pic_init_qs_minus26 + shr.slice_qs_delta;
    }

    assert(shr.SliceQpY >= -(sps->QpBdOffsetY) && shr.SliceQpY <= 51);
    assert(shr.QsY >= 0 && shr.QsY <= 51);

    shr.disable_deblocking_filter_idc = 0;
    shr.slice_alpha_c0_offset_div2    = 0;
    shr.slice_beta_offset_div2        = 0;
    if (pps->deblocking_filter_control_present_flag) {
        shr.disable_deblocking_filter_idc  = s->ue("SH: disable_deblocking_filter_idc");
        if (shr.disable_deblocking_filter_idc != 1) {
            shr.slice_alpha_c0_offset_div2 = s->se("SH: slice_alpha_c0_offset_div2");
            shr.slice_beta_offset_div2     = s->se("SH: slice_beta_offset_div2");
        }
    }
    shr.FilterOffsetA = shr.slice_alpha_c0_offset_div2 << 1;
    shr.FilterOffsetB = shr.slice_beta_offset_div2 << 1;

    assert(shr.disable_deblocking_filter_idc >= 0 && shr.disable_deblocking_filter_idc <= 2);
    assert(shr.slice_alpha_c0_offset_div2 >= -6 && shr.slice_alpha_c0_offset_div2 <=  6);
    assert(shr.slice_beta_offset_div2 >= -6 && shr.slice_beta_offset_div2 <=  6);

    if (pps->num_slice_groups_minus1 > 0 &&
        pps->slice_group_map_type >= 3 && pps->slice_group_map_type <= 5) {
        int len = (sps->pic_height_in_map_units_minus1+1)*(sps->pic_width_in_mbs_minus1+1)/
              (pps->slice_group_change_rate_minus1+1);
        if (((sps->pic_height_in_map_units_minus1+1)*(sps->pic_width_in_mbs_minus1+1))%
            (pps->slice_group_change_rate_minus1+1))
            len += 1;
        len = ceil(log2(len + 1));
        shr.slice_group_change_cycle = s->u(len, "SH: slice_group_change_cycle");
    }

    shr.MapUnitsInSliceGroup0 = min(shr.slice_group_change_cycle * pps->SliceGroupChangeRate, sps->PicSizeInMapUnits);
}

void ref_pic_list_modification(slice_t *currSlice)
{
    shr_t& shr = currSlice->header;
    data_partition_t *s = &currSlice->parser.partArr[0];

    shr.ref_pic_list_modification_flag_l0 = 0;
    if (shr.slice_type != I_slice && shr.slice_type != SI_slice) {
        shr.ref_pic_list_modification_flag_l0 = s->u(1, "SH: ref_pic_list_modification_flag_l0");
        if (shr.ref_pic_list_modification_flag_l0) {
            int i = 0;
            do {
                shr.modification_of_pic_nums_idc[0][i] = s->ue("SH: modification_of_pic_nums_idc_l0");
                if (shr.modification_of_pic_nums_idc[0][i] == 0 ||
                    shr.modification_of_pic_nums_idc[0][i] == 1)
                    shr.abs_diff_pic_num_minus1[0][i] = s->ue("SH: abs_diff_pic_num_minus1_l0");
                else if (shr.modification_of_pic_nums_idc[0][i] == 2)
                    shr.long_term_pic_num[0][i] = s->ue("SH: long_term_pic_idx_l0");
            } while (shr.modification_of_pic_nums_idc[0][i++] != 3);
        }
    }

    shr.ref_pic_list_modification_flag_l1 = 0;
    if (shr.slice_type == B_slice) {
        shr.ref_pic_list_modification_flag_l1 = s->u(1, "SH: ref_pic_list_reordering_flag_l1");
        if (shr.ref_pic_list_modification_flag_l1) {
            int i = 0;
            do {
                shr.modification_of_pic_nums_idc[1][i] = s->ue("SH: modification_of_pic_nums_idc_l1");
                if (shr.modification_of_pic_nums_idc[1][i] == 0 ||
                    shr.modification_of_pic_nums_idc[1][i] == 1)
                    shr.abs_diff_pic_num_minus1[1][i] = s->ue("SH: abs_diff_pic_num_minus1_l1");
                else if (shr.modification_of_pic_nums_idc[1][i] == 2)
                    shr.long_term_pic_num[1][i] = s->ue("SH: long_term_pic_idx_l1");
            } while (shr.modification_of_pic_nums_idc[1][i++] != 3);
        }
    }
}

#if (MVC_EXTENSION_ENABLE)
void ref_pic_list_mvc_modification(slice_t *currSlice)
{
    shr_t& shr = currSlice->header;
    data_partition_t *s = &currSlice->parser.partArr[0];

    shr.ref_pic_list_modification_flag_l0 = 0;
    if (shr.slice_type != I_slice && shr.slice_type != SI_slice) {
        shr.ref_pic_list_modification_flag_l0 = s->u(1, "SH: ref_pic_list_modification_flag_l0");
        if (shr.ref_pic_list_modification_flag_l0) {
            int i = 0;
            do {
                shr.modification_of_pic_nums_idc[0][i] = s->ue("SH: modification_of_pic_nums_idc_l0");
                if (shr.modification_of_pic_nums_idc[0][i] == 0 ||
                    shr.modification_of_pic_nums_idc[0][i] == 1)
                    shr.abs_diff_pic_num_minus1[0][i] = s->ue("SH: abs_diff_pic_num_minus1_l0");
                else if (shr.modification_of_pic_nums_idc[0][i] == 2)
                    shr.long_term_pic_num[0][i] = s->ue("SH: long_term_pic_idx_l0");
                else if (shr.modification_of_pic_nums_idc[0][i] == 4 ||
                         shr.modification_of_pic_nums_idc[0][i] == 5)
                    shr.abs_diff_view_idx_minus1[0][i] = s->ue("SH: abs_diff_view_idx_minus1_l0");
            } while (shr.modification_of_pic_nums_idc[0][i++] != 3);
        }
    }

    shr.ref_pic_list_modification_flag_l1 = 0;
    if (shr.slice_type == B_slice) {
        shr.ref_pic_list_modification_flag_l1 = s->u(1, "SH: ref_pic_list_reordering_flag_l1");
        if (shr.ref_pic_list_modification_flag_l1) {
            int i = 0;
            do {
                shr.modification_of_pic_nums_idc[1][i] = s->ue("SH: modification_of_pic_nums_idc_l1");
                if (shr.modification_of_pic_nums_idc[1][i] == 0 ||
                    shr.modification_of_pic_nums_idc[1][i] == 1)
                    shr.abs_diff_pic_num_minus1[1][i] = s->ue("SH: abs_diff_pic_num_minus1_l1");
                else if (shr.modification_of_pic_nums_idc[1][i] == 2)
                    shr.long_term_pic_num[1][i] = s->ue("SH: long_term_pic_idx_l1");
                else if (shr.modification_of_pic_nums_idc[1][i] == 4 ||
                         shr.modification_of_pic_nums_idc[1][i] == 5)
                    shr.abs_diff_view_idx_minus1[1][i] = s->ue("SH: abs_diff_view_idx_minus1_l1");
            } while (shr.modification_of_pic_nums_idc[1][i++] != 3);
        }
    }
}
#endif

void pred_weight_table(slice_t *currSlice)
{
    slice_t& slice = *currSlice;
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;
    data_partition_t *s = &slice.parser.partArr[0];

    shr.luma_log2_weight_denom = s->ue("SH: luma_log2_weight_denom");
    shr.chroma_log2_weight_denom = 0;
    if (sps.ChromaArrayType != 0)
        shr.chroma_log2_weight_denom = s->ue("SH: chroma_log2_weight_denom");

    assert(shr.luma_log2_weight_denom >= 0 && shr.luma_log2_weight_denom <= 7);
    assert(shr.chroma_log2_weight_denom >= 0 && shr.chroma_log2_weight_denom <= 7);

    for (int i = 0; i <= shr.num_ref_idx_l0_active_minus1; ++i) {
        shr.luma_weight_l0_flag[i] = s->u(1, "SH: luma_weight_flag_l0");
        shr.luma_weight_l0     [i] = 1 << shr.luma_log2_weight_denom;
        shr.luma_offset_l0     [i] = 0;
        if (shr.luma_weight_l0_flag[i]) {
            shr.luma_weight_l0[i] = s->se("SH: luma_weight_l0");
            shr.luma_offset_l0[i] = s->se("SH: luma_offset_l0");
        }

        assert(shr.luma_weight_l0[i] >= -128 && shr.luma_weight_l0[i] <= 127);
        assert(shr.luma_offset_l0[i] >= -128 && shr.luma_offset_l0[i] <= 127);

        if (sps.ChromaArrayType != 0) {
            shr.chroma_weight_l0_flag[i] = s->u(1, "SH: chroma_weight_flag_l0");
            for (int comp = 0; comp < 2; ++comp) {
                shr.chroma_weight_l0[i][comp] = 1 << shr.chroma_log2_weight_denom;
                shr.chroma_offset_l0[i][comp] = 0;
                if (shr.chroma_weight_l0_flag[i]) {
                    shr.chroma_weight_l0[i][comp] = s->se("SH: chroma_weight_l0");
                    shr.chroma_offset_l0[i][comp] = s->se("SH: chroma_offset_l0");
                }

                assert(shr.chroma_weight_l0[i][comp] >= -128 && shr.chroma_weight_l0[i][comp] <=  127);
                assert(shr.chroma_offset_l0[i][comp] >= -128 && shr.chroma_offset_l0[i][comp] <=  127);
            }
        }
    }

    if (shr.slice_type != B_slice)
        return;

    for (int i = 0; i <= shr.num_ref_idx_l1_active_minus1; ++i) {
        shr.luma_weight_l1_flag[i] = s->u(1, "SH: luma_weight_flag_l1");
        shr.luma_weight_l1     [i] = 1 << shr.luma_log2_weight_denom;
        shr.luma_offset_l1     [i] = 0;
        if (shr.luma_weight_l1_flag[i]) {
            shr.luma_weight_l1[i] = s->se("SH: luma_weight_l1");
            shr.luma_offset_l1[i] = s->se("SH: luma_offset_l1");
        }

        assert(shr.luma_weight_l1[i] >= -128 && shr.luma_weight_l1[i] <= 127);
        assert(shr.luma_offset_l1[i] >= -128 && shr.luma_offset_l1[i] <= 127);

        if (sps.ChromaArrayType != 0) {
            shr.chroma_weight_l1_flag[i] = s->u(1, "SH: chroma_weight_flag_l1");
            for (int comp = 0; comp < 2; ++comp) {
                shr.chroma_weight_l1[i][comp] = 1 << shr.chroma_log2_weight_denom;
                shr.chroma_offset_l1[i][comp] = 0;
                if (shr.chroma_weight_l1_flag[i]) {
                    shr.chroma_weight_l1[i][comp] = s->se("SH: chroma_weight_l1");
                    shr.chroma_offset_l1[i][comp] = s->se("SH: chroma_offset_l1");
                }

                assert(shr.chroma_weight_l1[i][comp] >= -128 && shr.chroma_weight_l1[i][comp] <=  127);
                assert(shr.chroma_offset_l1[i][comp] >= -128 && shr.chroma_offset_l1[i][comp] <=  127);
            }
        }
    }
}

void dec_ref_pic_marking(VideoParameters *p_Vid, data_partition_t *s, slice_t *currSlice)
{
    //data_partition_t *s = currSlice->parser.partArr[0].bitstream;

    sps_t& sps = *currSlice->active_sps;
    shr_t& shr = currSlice->header;
    int val;

    drpm_t* tmp_drpm, *tmp_drpm2;

    // free old buffer content
    while (shr.dec_ref_pic_marking_buffer) {
        tmp_drpm = shr.dec_ref_pic_marking_buffer;
        shr.dec_ref_pic_marking_buffer = tmp_drpm->Next;
        delete tmp_drpm;
    }

    if (currSlice->idr_flag
#if (MVC_EXTENSION_ENABLE)
        || (currSlice->svc_extension_flag == 0 && currSlice->NaluHeaderMVCExt.non_idr_flag == 0)
#endif
    ) {
        shr.no_output_of_prior_pics_flag = s->u(1, "SH: no_output_of_prior_pics_flag");
        p_Vid->no_output_of_prior_pics_flag = shr.no_output_of_prior_pics_flag;
        shr.long_term_reference_flag = s->u(1, "SH: long_term_reference_flag");

        assert(sps.max_num_ref_frames > 0 || shr.long_term_reference_flag == 1);
    } else {
        shr.adaptive_ref_pic_marking_mode_flag = s->u(1, "SH: adaptive_ref_pic_marking_mode_flag");
        if (shr.adaptive_ref_pic_marking_mode_flag) {
            do {
                tmp_drpm = new decoded_reference_picture_marking_t {};
                tmp_drpm->Next = NULL;

                val = tmp_drpm->memory_management_control_operation = s->ue("SH: memory_management_control_operation");

                if (val == 1 || val == 3)
                    tmp_drpm->difference_of_pic_nums_minus1 = s->ue("SH: difference_of_pic_nums_minus1");
                if (val == 2)
                    tmp_drpm->long_term_pic_num = s->ue("SH: long_term_pic_num");

                if (val == 3 || val == 6)
                    tmp_drpm->long_term_frame_idx = s->ue("SH: long_term_frame_idx");
                if (val == 4)
                    tmp_drpm->max_long_term_frame_idx_plus1 = s->ue("SH: max_long_term_pic_idx_plus1");

                // add command
                if (!shr.dec_ref_pic_marking_buffer)
                    shr.dec_ref_pic_marking_buffer = tmp_drpm;
                else {
                    tmp_drpm2 = shr.dec_ref_pic_marking_buffer;
                    while (tmp_drpm2->Next)
                        tmp_drpm2 = tmp_drpm2->Next;
                    tmp_drpm2->Next = tmp_drpm;
                }
            } while (val != 0);
        }
    }
}


void decode_poc(VideoParameters *p_Vid, slice_t *pSlice)
{
    slice_t& slice = *pSlice;
    sps_t& sps = *p_Vid->active_sps;
    shr_t& shr = pSlice->header;

    switch (sps.pic_order_cnt_type) {
    case 0:
        if (pSlice->idr_flag) {
            p_Vid->prevPicOrderCntMsb = 0;
            p_Vid->prevPicOrderCntLsb = 0;
        } else if (p_Vid->last_has_mmco_5) {
            p_Vid->prevPicOrderCntMsb = 0;
            p_Vid->prevPicOrderCntLsb = !p_Vid->last_pic_bottom_field ? shr.TopFieldOrderCnt : 0;
        }

        if ((shr.pic_order_cnt_lsb < p_Vid->prevPicOrderCntLsb) &&
            (p_Vid->prevPicOrderCntLsb - shr.pic_order_cnt_lsb >= sps.MaxPicOrderCntLsb / 2))
            shr.PicOrderCntMsb = p_Vid->prevPicOrderCntMsb + sps.MaxPicOrderCntLsb;
        else if ((shr.pic_order_cnt_lsb > p_Vid->prevPicOrderCntLsb) &&
                 (shr.pic_order_cnt_lsb - p_Vid->prevPicOrderCntLsb > sps.MaxPicOrderCntLsb / 2))
            shr.PicOrderCntMsb = p_Vid->prevPicOrderCntMsb - sps.MaxPicOrderCntLsb;
        else
            shr.PicOrderCntMsb = p_Vid->prevPicOrderCntMsb;

        shr.TopFieldOrderCnt    = 0;
        shr.BottomFieldOrderCnt = 0;
        if (!shr.field_pic_flag || !shr.bottom_field_flag)
            shr.TopFieldOrderCnt = shr.PicOrderCntMsb + shr.pic_order_cnt_lsb;
        if (!shr.field_pic_flag)
            shr.BottomFieldOrderCnt = shr.TopFieldOrderCnt + shr.delta_pic_order_cnt_bottom;
        else if (shr.bottom_field_flag)
            shr.BottomFieldOrderCnt = shr.PicOrderCntMsb + shr.pic_order_cnt_lsb;

        if (slice.nal_ref_idc) {
            p_Vid->prevPicOrderCntMsb = shr.PicOrderCntMsb;
            p_Vid->prevPicOrderCntLsb = shr.pic_order_cnt_lsb;
        }
        break;

    case 1:
        if (p_Vid->last_has_mmco_5) {
            p_Vid->prevFrameNum       = 0;
            p_Vid->prevFrameNumOffset = 0;
        }

        if (pSlice->idr_flag)
            shr.FrameNumOffset = 0;
        else if (p_Vid->prevFrameNum > shr.frame_num)
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset + sps.MaxFrameNum;
        else
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset;

        int32_t absFrameNum;
        if (sps.num_ref_frames_in_pic_order_cnt_cycle != 0)
            absFrameNum = shr.FrameNumOffset + shr.frame_num;
        else
            absFrameNum = 0;
        if (slice.nal_ref_idc == 0 && absFrameNum > 0)
            absFrameNum--;

        int32_t expectedPicOrderCnt;
        if (absFrameNum > 0) {
            int32_t picOrderCntCycleCnt        = (absFrameNum - 1) / sps.num_ref_frames_in_pic_order_cnt_cycle;
            int32_t frameNumInPicOrderCntCycle = (absFrameNum - 1) % sps.num_ref_frames_in_pic_order_cnt_cycle;
            expectedPicOrderCnt = picOrderCntCycleCnt * sps.ExpectedDeltaPerPicOrderCntCycle;
            for (int i = 0; i <= frameNumInPicOrderCntCycle; i++)
                expectedPicOrderCnt += sps.offset_for_ref_frame[i];
        } else
            expectedPicOrderCnt = 0;
        if (slice.nal_ref_idc == 0)
            expectedPicOrderCnt += sps.offset_for_non_ref_pic;

        shr.TopFieldOrderCnt    = 0;
        shr.BottomFieldOrderCnt = 0;
        if (!shr.field_pic_flag || !shr.bottom_field_flag)
            shr.TopFieldOrderCnt = expectedPicOrderCnt + shr.delta_pic_order_cnt[0];
        if (!shr.field_pic_flag)
            shr.BottomFieldOrderCnt = shr.TopFieldOrderCnt + sps.offset_for_top_to_bottom_field + shr.delta_pic_order_cnt[1];
        else if (shr.bottom_field_flag)
            shr.BottomFieldOrderCnt = expectedPicOrderCnt + sps.offset_for_top_to_bottom_field + shr.delta_pic_order_cnt[0];

        p_Vid->prevFrameNum       = shr.frame_num;
        p_Vid->prevFrameNumOffset = shr.FrameNumOffset;
        break;

    case 2:
        if (p_Vid->last_has_mmco_5) {
            p_Vid->prevFrameNum       = 0;
            p_Vid->prevFrameNumOffset = 0;
        }

        if (slice.idr_flag)
            shr.FrameNumOffset = 0;
        else if (p_Vid->prevFrameNum > shr.frame_num)
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset + sps.MaxFrameNum;
        else
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset;

        int tempPicOrderCnt;
        if (slice.idr_flag)
            tempPicOrderCnt = 0;
        else if (slice.nal_ref_idc == 0)
            tempPicOrderCnt = 2 * (shr.FrameNumOffset + shr.frame_num) - 1;
        else
            tempPicOrderCnt = 2 * (shr.FrameNumOffset + shr.frame_num);

        shr.TopFieldOrderCnt    = 0;
        shr.BottomFieldOrderCnt = 0;
        if (!shr.field_pic_flag || !shr.bottom_field_flag)
            shr.TopFieldOrderCnt = tempPicOrderCnt;
        if (!shr.field_pic_flag || shr.bottom_field_flag)
            shr.BottomFieldOrderCnt = tempPicOrderCnt;

        p_Vid->prevFrameNum       = shr.frame_num;
        p_Vid->prevFrameNumOffset = shr.FrameNumOffset;
        break;

    default:
        assert(false);
        break;
    }

    if (!shr.field_pic_flag)
        shr.PicOrderCnt = min(shr.TopFieldOrderCnt, shr.BottomFieldOrderCnt);
    else if (!shr.bottom_field_flag)
        shr.PicOrderCnt = shr.TopFieldOrderCnt;
    else
        shr.PicOrderCnt = shr.BottomFieldOrderCnt;
    p_Vid->PicOrderCnt  = shr.PicOrderCnt;
    p_Vid->prevFrameNum = shr.frame_num;
}
