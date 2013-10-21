#include "global.h"
#include "interpret.h"
#include "dpb.h"
#include "slice.h"
#include "memalloc.h"


using vio::h264::data_partition_t;


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
    result |= this->IdrPicFlag  != slice.IdrPicFlag;

    if (slice.IdrPicFlag && this->IdrPicFlag)
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


// 7.3.3 Slice header syntax

void data_partition_t::slice_header(slice_t& slice)
{
    shr_t& shr = slice.header;

    shr.first_mb_in_slice    = this->ue("SH: first_mb_in_slice");
    shr.slice_type           = this->ue("SH: slice_type") % 5;
    shr.pic_parameter_set_id = this->ue("SH: pic_parameter_set_id");

    VideoParameters* p_Vid = slice.p_Vid;
    const pps_t& pps = p_Vid->PicParSet[shr.pic_parameter_set_id];
    const sps_t& sps = p_Vid->SeqParSet[pps.seq_parameter_set_id];
    slice.active_pps = &p_Vid->PicParSet[shr.pic_parameter_set_id];
    slice.active_sps = &p_Vid->SeqParSet[pps.seq_parameter_set_id];

    shr.colour_plane_id = 0;
    if (sps.separate_colour_plane_flag)
        shr.colour_plane_id = this->u(2, "SH: colour_plane_id");

    assert(shr.colour_plane_id <= 2);

    shr.frame_num = this->u(sps.log2_max_frame_num_minus4 + 4, "SH: frame_num");
    assert(!slice.IdrPicFlag || shr.frame_num == 0);

    shr.field_pic_flag    = 0;
    shr.bottom_field_flag = 0;
    if (!sps.frame_mbs_only_flag) {
        shr.field_pic_flag = this->u(1, "SH: field_pic_flag");
        if (shr.field_pic_flag)
            shr.bottom_field_flag = this->u(1, "SH: bottom_field_flag");
    }

    shr.MbaffFrameFlag     = sps.mb_adaptive_frame_field_flag && !shr.field_pic_flag;
    shr.PicHeightInMbs     = sps.FrameHeightInMbs / (1 + shr.field_pic_flag);
    shr.PicHeightInSampleL = shr.PicHeightInMbs * 16;
    shr.PicHeightInSampleC = shr.PicHeightInMbs * sps.MbHeightC;
    shr.PicSizeInMbs       = sps.PicWidthInMbs * shr.PicHeightInMbs;
    shr.MaxPicNum          = sps.MaxFrameNum * (1 + shr.field_pic_flag);
    shr.CurrPicNum         = shr.frame_num * (1 + shr.field_pic_flag) + shr.field_pic_flag;

    if (slice.IdrPicFlag || (slice.mvc_extension_flag && !slice.non_idr_flag))
        shr.idr_pic_id = this->ue("SH: idr_pic_id");

    shr.delta_pic_order_cnt_bottom = 0;
    shr.delta_pic_order_cnt[0]     = 0;
    shr.delta_pic_order_cnt[1]     = 0;
    if (sps.pic_order_cnt_type == 0) {
        shr.pic_order_cnt_lsb = this->u(sps.log2_max_pic_order_cnt_lsb_minus4 + 4, "SH: pic_order_cnt_lsb");
        assert(shr.pic_order_cnt_lsb < sps.MaxPicOrderCntLsb);
        if (pps.bottom_field_pic_order_in_frame_present_flag && !shr.field_pic_flag)
            shr.delta_pic_order_cnt_bottom = this->se("SH: delta_pic_order_cnt_bottom");
    }
    if (sps.pic_order_cnt_type == 1 && !sps.delta_pic_order_always_zero_flag) {
        shr.delta_pic_order_cnt[0] = this->se("SH: delta_pic_order_cnt[0]");
        if (pps.bottom_field_pic_order_in_frame_present_flag && !shr.field_pic_flag)
            shr.delta_pic_order_cnt[1] = this->se("SH: delta_pic_order_cnt[1]");
    }

    shr.redundant_pic_cnt = 0;
    if (pps.redundant_pic_cnt_present_flag)
        shr.redundant_pic_cnt = this->ue("SH: redundant_pic_cnt");

    assert(shr.redundant_pic_cnt < 128);

    if (shr.slice_type == B_slice)
        shr.direct_spatial_mv_pred_flag = this->u(1, "SH: direct_spatial_mv_pred_flag");

    shr.num_ref_idx_l0_active_minus1 = pps.num_ref_idx_l0_default_active_minus1;
    shr.num_ref_idx_l1_active_minus1 = pps.num_ref_idx_l1_default_active_minus1;
    if (shr.slice_type == P_slice || shr.slice_type == SP_slice || shr.slice_type == B_slice) {
        shr.num_ref_idx_active_override_flag = this->u(1, "SH: num_ref_idx_override_flag");
        if (shr.num_ref_idx_active_override_flag) {
            shr.num_ref_idx_l0_active_minus1 = this->ue("SH: num_ref_idx_l0_active_minus1");
            if (shr.slice_type == B_slice)
                shr.num_ref_idx_l1_active_minus1 = this->ue("SH: num_ref_idx_l1_active_minus1");
        }
    }

    if (shr.slice_type == P_slice || shr.slice_type == SP_slice || shr.slice_type == B_slice) {
        if (!shr.field_pic_flag && pps.num_ref_idx_l0_default_active_minus1 > 15)
            assert(shr.num_ref_idx_active_override_flag);
        if (shr.slice_type == B_slice) {
            if (!shr.field_pic_flag && pps.num_ref_idx_l1_default_active_minus1 > 15)
                assert(shr.num_ref_idx_active_override_flag);
        }
    }
    assert(shr.num_ref_idx_l0_active_minus1 < MAX_NUM_REF_IDX / (1 + !shr.field_pic_flag));
    assert(shr.num_ref_idx_l1_active_minus1 < MAX_NUM_REF_IDX / (1 + !shr.field_pic_flag));

    if (slice.nal_unit_type == 20 || slice.nal_unit_type == 21)
        this->ref_pic_list_mvc_modification(slice);
    else
        this->ref_pic_list_modification(slice);

    if ((pps.weighted_pred_flag && (shr.slice_type == P_slice || shr.slice_type == SP_slice)) ||
        (pps.weighted_bipred_idc == 1 && shr.slice_type == B_slice))
        this->pred_weight_table(slice);
    else {
        shr.luma_log2_weight_denom   = 5;
        shr.chroma_log2_weight_denom = 5;
    }

    if (slice.nal_ref_idc != 0)
        this->dec_ref_pic_marking(slice);

    shr.cabac_init_idc = 0;
    if (pps.entropy_coding_mode_flag && shr.slice_type != I_slice && shr.slice_type != SI_slice)
        shr.cabac_init_idc = this->ue("SH: cabac_init_idc");

    assert(shr.cabac_init_idc <= 2);

    shr.slice_qp_delta = this->se("SH: slice_qp_delta");
    shr.SliceQpY = 26 + pps.pic_init_qp_minus26 + shr.slice_qp_delta;

    shr.sp_for_switch_flag = 0;
    shr.QsY = 0;
    if (shr.slice_type == SP_slice || shr.slice_type == SI_slice) {
        if (shr.slice_type == SP_slice)
            shr.sp_for_switch_flag = this->u(1, "SH: sp_for_switch_flag");
        shr.slice_qs_delta = this->se("SH: slice_qs_delta");
        shr.QsY = 26 + pps.pic_init_qs_minus26 + shr.slice_qs_delta;
    }

    assert(shr.SliceQpY >= -(sps.QpBdOffsetY) && shr.SliceQpY <= 51);
    assert(shr.QsY >= 0 && shr.QsY <= 51);

    shr.disable_deblocking_filter_idc = 0;
    shr.slice_alpha_c0_offset_div2    = 0;
    shr.slice_beta_offset_div2        = 0;
    if (pps.deblocking_filter_control_present_flag) {
        shr.disable_deblocking_filter_idc  = this->ue("SH: disable_deblocking_filter_idc");
        if (shr.disable_deblocking_filter_idc != 1) {
            shr.slice_alpha_c0_offset_div2 = this->se("SH: slice_alpha_c0_offset_div2");
            shr.slice_beta_offset_div2     = this->se("SH: slice_beta_offset_div2");
        }
    }
    shr.FilterOffsetA = shr.slice_alpha_c0_offset_div2 << 1;
    shr.FilterOffsetB = shr.slice_beta_offset_div2 << 1;

    assert(shr.disable_deblocking_filter_idc <= 2);
    assert(shr.slice_alpha_c0_offset_div2 >= -6 && shr.slice_alpha_c0_offset_div2 <= 6);
    assert(shr.slice_beta_offset_div2 >= -6 && shr.slice_beta_offset_div2 <= 6);

    if (pps.num_slice_groups_minus1 > 0 &&
        pps.slice_group_map_type >= 3 && pps.slice_group_map_type <= 5) {
        int len = ceil(log2(sps.PicSizeInMapUnits / pps.SliceGroupChangeRate + 1));
        shr.slice_group_change_cycle = this->u(len, "SH: slice_group_change_cycle");
        shr.MapUnitsInSliceGroup0 =
            min(shr.slice_group_change_cycle * pps.SliceGroupChangeRate, sps.PicSizeInMapUnits);
    }

    shr.structure = !shr.field_pic_flag ? FRAME :
                    !shr.bottom_field_flag ? TOP_FIELD : BOTTOM_FIELD;
}

// 7.3.3.1 Reference picture list modification syntax

void data_partition_t::ref_pic_list_modification(slice_t& slice)
{
    shr_t& shr = slice.header;

    shr.ref_pic_list_modification_flag_l0 = 0;
    if (shr.slice_type != I_slice && shr.slice_type != SI_slice) {
        shr.ref_pic_list_modification_flag_l0 = this->u(1, "SH: ref_pic_list_modification_flag_l0");
        if (shr.ref_pic_list_modification_flag_l0) {
            slice_header_t::ref_pic_list_modification_t idc;
            shr.ref_pic_list_modifications[0].clear();
            do {
                idc.modification_of_pic_nums_idc = this->ue("SH: modification_of_pic_nums_idc_l0");
                if (idc.modification_of_pic_nums_idc == 0 ||
                    idc.modification_of_pic_nums_idc == 1) {
                    idc.abs_diff_pic_num_minus1 = this->ue("SH: abs_diff_pic_num_minus1_l0");
                    assert(idc.abs_diff_pic_num_minus1 < shr.MaxPicNum);
                } else if (idc.modification_of_pic_nums_idc == 2)
                    idc.long_term_pic_num = this->ue("SH: long_term_pic_idx_l0");
                if (idc.modification_of_pic_nums_idc != 3)
                    shr.ref_pic_list_modifications[0].push_back(idc);
            } while (idc.modification_of_pic_nums_idc != 3);

            assert(shr.ref_pic_list_modifications[0].size() <= shr.num_ref_idx_l0_active_minus1 + 1);
        }
    }

    shr.ref_pic_list_modification_flag_l1 = 0;
    if (shr.slice_type == B_slice) {
        shr.ref_pic_list_modification_flag_l1 = this->u(1, "SH: ref_pic_list_reordering_flag_l1");
        if (shr.ref_pic_list_modification_flag_l1) {
            slice_header_t::ref_pic_list_modification_t idc;
            shr.ref_pic_list_modifications[1].clear();
            do {
                idc.modification_of_pic_nums_idc = this->ue("SH: modification_of_pic_nums_idc_l1");
                if (idc.modification_of_pic_nums_idc == 0 ||
                    idc.modification_of_pic_nums_idc == 1) {
                    idc.abs_diff_pic_num_minus1 = this->ue("SH: abs_diff_pic_num_minus1_l1");
                    assert(idc.abs_diff_pic_num_minus1 < shr.MaxPicNum);
                } else if (idc.modification_of_pic_nums_idc == 2)
                    idc.long_term_pic_num = this->ue("SH: long_term_pic_idx_l1");
                if (idc.modification_of_pic_nums_idc != 3)
                    shr.ref_pic_list_modifications[1].push_back(idc);
            } while (idc.modification_of_pic_nums_idc != 3);

            assert(shr.ref_pic_list_modifications[1].size() <= shr.num_ref_idx_l1_active_minus1 + 1);
        }
    }
}

// H.7.3.3.1.1 Reference picture list MVC modification syntax

void data_partition_t::ref_pic_list_mvc_modification(slice_t& slice)
{
    shr_t& shr = slice.header;

    shr.ref_pic_list_modification_flag_l0 = 0;
    if (shr.slice_type != I_slice && shr.slice_type != SI_slice) {
        shr.ref_pic_list_modification_flag_l0 = this->u(1, "SH: ref_pic_list_modification_flag_l0");
        if (shr.ref_pic_list_modification_flag_l0) {
            slice_header_t::ref_pic_list_modification_t idc;
            shr.ref_pic_list_modifications[0].clear();
            do {
                idc.modification_of_pic_nums_idc = this->ue("SH: modification_of_pic_nums_idc_l0");
                if (idc.modification_of_pic_nums_idc == 0 ||
                    idc.modification_of_pic_nums_idc == 1) {
                    idc.abs_diff_pic_num_minus1 = this->ue("SH: abs_diff_pic_num_minus1_l0");
                    assert(idc.abs_diff_pic_num_minus1 < shr.MaxPicNum);
                } else if (idc.modification_of_pic_nums_idc == 2)
                    idc.long_term_pic_num = this->ue("SH: long_term_pic_idx_l0");
                else if (idc.modification_of_pic_nums_idc == 4 ||
                         idc.modification_of_pic_nums_idc == 5)
                    idc.abs_diff_view_idx_minus1 = this->ue("SH: abs_diff_view_idx_minus1_l0");
                if (idc.modification_of_pic_nums_idc != 3)
                    shr.ref_pic_list_modifications[0].push_back(idc);
            } while (idc.modification_of_pic_nums_idc != 3);

            assert(shr.ref_pic_list_modifications[0].size() <= shr.num_ref_idx_l0_active_minus1 + 1);
        }
    }

    shr.ref_pic_list_modification_flag_l1 = 0;
    if (shr.slice_type == B_slice) {
        shr.ref_pic_list_modification_flag_l1 = this->u(1, "SH: ref_pic_list_reordering_flag_l1");
        if (shr.ref_pic_list_modification_flag_l1) {
            slice_header_t::ref_pic_list_modification_t idc;
            shr.ref_pic_list_modifications[1].clear();
            do {
                idc.modification_of_pic_nums_idc = this->ue("SH: modification_of_pic_nums_idc_l1");
                if (idc.modification_of_pic_nums_idc == 0 ||
                    idc.modification_of_pic_nums_idc == 1) {
                    idc.abs_diff_pic_num_minus1 = this->ue("SH: abs_diff_pic_num_minus1_l1");
                    assert(idc.abs_diff_pic_num_minus1 < shr.MaxPicNum);
                } else if (idc.modification_of_pic_nums_idc == 2)
                    idc.long_term_pic_num = this->ue("SH: long_term_pic_idx_l1");
                else if (idc.modification_of_pic_nums_idc == 4 ||
                         idc.modification_of_pic_nums_idc == 5)
                    idc.abs_diff_view_idx_minus1 = this->ue("SH: abs_diff_view_idx_minus1_l1");
                if (idc.modification_of_pic_nums_idc != 3)
                    shr.ref_pic_list_modifications[1].push_back(idc);
            } while (idc.modification_of_pic_nums_idc != 3);

            assert(shr.ref_pic_list_modifications[1].size() <= shr.num_ref_idx_l1_active_minus1 + 1);
        }
    }
}

// 7.3.3.2 Prediction weight table syntax

void data_partition_t::pred_weight_table(slice_t& slice)
{
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;

    shr.luma_log2_weight_denom = this->ue("SH: luma_log2_weight_denom");
    shr.chroma_log2_weight_denom = 0;
    if (sps.ChromaArrayType != 0)
        shr.chroma_log2_weight_denom = this->ue("SH: chroma_log2_weight_denom");

    assert(shr.luma_log2_weight_denom   < 8);
    assert(shr.chroma_log2_weight_denom < 8);

    shr.pred_weight_l[0][0] = slice_header_t::pred_weight_v(shr.num_ref_idx_l0_active_minus1 + 1);
    if (sps.ChromaArrayType != 0) {
        shr.pred_weight_l[0][1] = slice_header_t::pred_weight_v(shr.num_ref_idx_l0_active_minus1 + 1);
        shr.pred_weight_l[0][2] = slice_header_t::pred_weight_v(shr.num_ref_idx_l0_active_minus1 + 1);
    }

    for (int i = 0; i <= shr.num_ref_idx_l0_active_minus1; ++i) {
        auto& luma_weight = shr.pred_weight_l[0][0][i];
        luma_weight.weight_flag = this->u(1, "SH: luma_weight_flag_l0");
        luma_weight.weight      = 1 << shr.luma_log2_weight_denom;
        luma_weight.offset      = 0;
        if (luma_weight.weight_flag) {
            luma_weight.weight = this->se("SH: luma_weight_l0");
            luma_weight.offset = this->se("SH: luma_offset_l0");
        }

        if (sps.ChromaArrayType != 0) {
            bool chroma_weight_flag = this->u(1, "SH: chroma_weight_flag_l0");
            for (int j = 0; j < 2; ++j) {
                auto& chroma_weight = shr.pred_weight_l[0][j + 1][i];
                chroma_weight.weight_flag = chroma_weight_flag;;
                chroma_weight.weight      = 1 << shr.chroma_log2_weight_denom;
                chroma_weight.offset      = 0;
                if (chroma_weight.weight_flag) {
                    chroma_weight.weight = this->se("SH: chroma_weight_l0");
                    chroma_weight.offset = this->se("SH: chroma_offset_l0");
                }
            }
        }
    }

    if (shr.slice_type != B_slice)
        return;

    shr.pred_weight_l[1][0] = slice_header_t::pred_weight_v(shr.num_ref_idx_l1_active_minus1 + 1);
    if (sps.ChromaArrayType != 0) {
        shr.pred_weight_l[1][1] = slice_header_t::pred_weight_v(shr.num_ref_idx_l1_active_minus1 + 1);
        shr.pred_weight_l[1][2] = slice_header_t::pred_weight_v(shr.num_ref_idx_l1_active_minus1 + 1);
    }

    for (int i = 0; i <= shr.num_ref_idx_l1_active_minus1; ++i) {
        auto& luma_weight = shr.pred_weight_l[1][0][i];
        luma_weight.weight_flag = this->u(1, "SH: luma_weight_flag_l1");
        luma_weight.weight      = 1 << shr.luma_log2_weight_denom;
        luma_weight.offset      = 0;
        if (luma_weight.weight_flag) {
            luma_weight.weight = this->se("SH: luma_weight_l1");
            luma_weight.offset = this->se("SH: luma_offset_l1");
        }

        if (sps.ChromaArrayType != 0) {
            bool chroma_weight_flag = this->u(1, "SH: chroma_weight_flag_l1");
            for (int j = 0; j < 2; ++j) {
                auto& chroma_weight = shr.pred_weight_l[1][j + 1][i];
                chroma_weight.weight_flag = chroma_weight_flag;
                chroma_weight.weight      = 1 << shr.chroma_log2_weight_denom;
                chroma_weight.offset      = 0;
                if (chroma_weight.weight_flag) {
                    chroma_weight.weight = this->se("SH: chroma_weight_l1");
                    chroma_weight.offset = this->se("SH: chroma_offset_l1");
                }
            }
        }
    }
}

// 7.3.3.3 Decoded reference picture marking syntax

void data_partition_t::dec_ref_pic_marking(slice_t& slice)
{
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;

    shr.adaptive_ref_pic_markings.clear();

    if (slice.IdrPicFlag || (slice.mvc_extension_flag && !slice.non_idr_flag)) {
        shr.no_output_of_prior_pics_flag = this->u(1, "SH: no_output_of_prior_pics_flag");
        shr.long_term_reference_flag     = this->u(1, "SH: long_term_reference_flag");

        assert(sps.max_num_ref_frames != 0 || shr.long_term_reference_flag);
    } else {
        shr.adaptive_ref_pic_marking_mode_flag = this->u(1, "SH: adaptive_ref_pic_marking_mode_flag");
        if (shr.adaptive_ref_pic_marking_mode_flag) {
            slice_header_t::adaptive_ref_pic_marking_t mmco;
            do {
                mmco.memory_management_control_operation = this->ue("SH: memory_management_control_operation");

                if (mmco.memory_management_control_operation == 1 ||
                    mmco.memory_management_control_operation == 3)
                    mmco.difference_of_pic_nums_minus1 = this->ue("SH: difference_of_pic_nums_minus1");
                if (mmco.memory_management_control_operation == 2)
                    mmco.long_term_pic_num = this->ue("SH: long_term_pic_num");
                if (mmco.memory_management_control_operation == 3 ||
                    mmco.memory_management_control_operation == 6)
                    mmco.long_term_frame_idx = this->ue("SH: long_term_frame_idx");
                if (mmco.memory_management_control_operation == 4) {
                    mmco.max_long_term_frame_idx_plus1 = this->ue("SH: max_long_term_pic_idx_plus1");
                    assert(mmco.max_long_term_frame_idx_plus1 <= sps.max_num_ref_frames);
                }

                if (mmco.memory_management_control_operation != 0)
                    shr.adaptive_ref_pic_markings.push_back(mmco);
            } while (mmco.memory_management_control_operation != 0);
        }
    }
}

// 7.3.4 Slice data syntax

void data_partition_t::slice_data()
{
    
}
