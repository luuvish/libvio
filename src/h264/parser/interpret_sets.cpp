/*
 * =============================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * =============================================================================
 *
 *  File      : interpret_sets.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * =============================================================================
 */

#include <cassert>
#include "global.h"
#include "defines.h"
#include "sets.h"
#include "slice.h"
#include "interpret.h"


namespace vio  {
namespace h264 {


// Table A-1 Level limits

static const int MaxDpbMbs[15] = {
    396, 900, 2376, 2376, 4752, 8100, 8100, 18000, 20480,
    32768, 32768, 34816, 110400, 184320, 184320
};

// 7.3.2.1 Sequence parameter set data syntax

void data_partition_t::seq_parameter_set_rbsp(sps_t& sps)
{
    this->seq_parameter_set_data(sps);
    this->rbsp_trailing_bits();
}

// 7.3.2.1.1 Sequence parameter set data syntax

void data_partition_t::seq_parameter_set_data(sps_t& sps)
{
    #define MAX_NUM_REF_FRAMES_IN_POC_CYCLE 256

    uint8_t reserved_zero_2bits;

    sps.profile_idc          = this->u(8, "SPS: profile_idc");
    sps.constraint_set0_flag = this->u(1, "SPS: constrained_set0_flag");
    sps.constraint_set1_flag = this->u(1, "SPS: constrained_set1_flag");
    sps.constraint_set2_flag = this->u(1, "SPS: constrained_set2_flag");
    sps.constraint_set3_flag = this->u(1, "SPS: constrained_set3_flag");
    sps.constraint_set4_flag = this->u(1, "SPS: constrained_set4_flag");
    sps.constraint_set5_flag = this->u(1, "SPS: constrained_set5_flag");
    reserved_zero_2bits      = this->u(2, "SPS: reserved_zero_2bits");

    if ((sps.profile_idc != BASELINE       ) &&
        (sps.profile_idc != MAIN           ) &&
        (sps.profile_idc != EXTENDED       ) &&
        (sps.profile_idc != FREXT_HP       ) &&
        (sps.profile_idc != FREXT_Hi10P    ) &&
        (sps.profile_idc != FREXT_Hi422    ) &&
        (sps.profile_idc != FREXT_Hi444    ) &&
        (sps.profile_idc != FREXT_CAVLC444 )
        && (sps.profile_idc != MVC_HIGH)
        && (sps.profile_idc != STEREO_HIGH)) {
        printf("Invalid Profile IDC (%d) encountered. \n", sps.profile_idc);
        return;
    }

    sps.level_idc            = this->u(8, "SPS: level_idc");
    sps.seq_parameter_set_id = this->ue("SPS: seq_parameter_set_id");

    assert(sps.seq_parameter_set_id < MAX_NUM_SPS);

    sps.chroma_format_idc                    = CHROMA_FORMAT_420;
    sps.separate_colour_plane_flag           = 0;
    sps.bit_depth_luma_minus8                = 0;
    sps.bit_depth_chroma_minus8              = 0;
    sps.qpprime_y_zero_transform_bypass_flag = 0;
    sps.seq_scaling_matrix_present_flag      = 0;
    if (sps.profile_idc == FREXT_HP    || sps.profile_idc == FREXT_Hi10P ||
        sps.profile_idc == FREXT_Hi422 || sps.profile_idc == FREXT_Hi444 ||
        sps.profile_idc == FREXT_CAVLC444 || sps.profile_idc == MVC_HIGH ||
        sps.profile_idc == STEREO_HIGH) {

        sps.chroma_format_idc                    = this->ue("SPS: chroma_format_idc");
        if (sps.chroma_format_idc == CHROMA_FORMAT_444)
            sps.separate_colour_plane_flag       = this->u(1, "SPS: separate_colour_plane_flag");
        sps.bit_depth_luma_minus8                = this->ue("SPS: bit_depth_luma_minus8");
        sps.bit_depth_chroma_minus8              = this->ue("SPS: bit_depth_chroma_minus8");
        sps.qpprime_y_zero_transform_bypass_flag = this->u(1, "SPS: lossless_qpprime_y_zero_flag");
        sps.seq_scaling_matrix_present_flag      = this->u(1, "SPS: seq_scaling_matrix_present_flag");
        if (sps.seq_scaling_matrix_present_flag) {
            for (int i = 0; i < (sps.chroma_format_idc != CHROMA_FORMAT_444 ? 8 : 12); ++i) {
                sps.seq_scaling_list_present_flag[i] = this->u(1, "SPS: seq_scaling_list_present_flag");
                if (sps.seq_scaling_list_present_flag[i]) {
                    if (i < 6)
                        this->scaling_list(sps.ScalingList4x4[i], 16,
                                           sps.UseDefaultScalingMatrix4x4Flag[i]);
                    else
                        this->scaling_list(sps.ScalingList8x8[i - 6], 64,
                                           sps.UseDefaultScalingMatrix8x8Flag[i - 6]);
                }
            }
        }
    }

    assert(sps.profile_idc != 183 || sps.chroma_format_idc == CHROMA_FORMAT_400);
    assert(sps.chroma_format_idc <= CHROMA_FORMAT_444);

    sps.ChromaArrayType = !sps.separate_colour_plane_flag ? sps.chroma_format_idc : 0;
    sps.SubWidthC  = sps.ChromaArrayType != 0 && sps.ChromaArrayType < 3 ? 2 : 1;
    sps.SubHeightC = sps.ChromaArrayType != 0 && sps.ChromaArrayType < 2 ? 2 : 1;
    sps.MbWidthC   = sps.ChromaArrayType == 0 ? 0 : 16 / sps.SubWidthC;
    sps.MbHeightC  = sps.ChromaArrayType == 0 ? 0 : 16 / sps.SubHeightC;

    assert(sps.bit_depth_luma_minus8   <= 6);
    assert(sps.bit_depth_chroma_minus8 <= 6);

    sps.BitDepthY   = 8 + sps.bit_depth_luma_minus8;
    sps.QpBdOffsetY = 6 * sps.bit_depth_luma_minus8;
    sps.BitDepthC   = 8 + sps.bit_depth_chroma_minus8;
    sps.QpBdOffsetC = 6 * sps.bit_depth_chroma_minus8;
    sps.RawMbBits   = 256 * sps.BitDepthY + 2 * sps.MbWidthC * sps.MbHeightC * sps.BitDepthC;

    sps.log2_max_frame_num_minus4 = this->ue("SPS: log2_max_frame_num_minus4");

    assert(sps.log2_max_frame_num_minus4 <= 12);
    sps.MaxFrameNum = 1 << (sps.log2_max_frame_num_minus4 + 4);

    sps.pic_order_cnt_type = this->ue("SPS: pic_order_cnt_type");

    assert(sps.pic_order_cnt_type <= 2);

    if (sps.pic_order_cnt_type == 0) {
        sps.log2_max_pic_order_cnt_lsb_minus4 = this->ue("SPS: log2_max_pic_order_cnt_lsb_minus4");

        assert(sps.log2_max_pic_order_cnt_lsb_minus4 <= 12);
        sps.MaxPicOrderCntLsb = 1 << (sps.log2_max_pic_order_cnt_lsb_minus4 + 4);
    } else if (sps.pic_order_cnt_type == 1) {
        sps.delta_pic_order_always_zero_flag      = this->u(1, "SPS: delta_pic_order_always_zero_flag");
        sps.offset_for_non_ref_pic                = this->se("SPS: offset_for_non_ref_pic");
        sps.offset_for_top_to_bottom_field        = this->se("SPS: offset_for_top_to_bottom_field");
        sps.num_ref_frames_in_pic_order_cnt_cycle = this->ue("SPS: num_ref_frames_in_pic_order_cnt_cycle");

        if (sps.num_ref_frames_in_pic_order_cnt_cycle > 0)
            sps.offset_for_ref_frame = sps_t::int32_v(sps.num_ref_frames_in_pic_order_cnt_cycle);

        for (int i = 0; i < sps.num_ref_frames_in_pic_order_cnt_cycle; ++i)
            sps.offset_for_ref_frame[i] = this->se("SPS: offset_for_ref_frame[i]");

        sps.ExpectedDeltaPerPicOrderCntCycle = 0;
        for (int i = 0; i < sps.num_ref_frames_in_pic_order_cnt_cycle; ++i)
            sps.ExpectedDeltaPerPicOrderCntCycle += sps.offset_for_ref_frame[i];
    }

    sps.max_num_ref_frames = this->ue("SPS: num_ref_frames");

    sps.gaps_in_frame_num_value_allowed_flag = this->u(1, "SPS: gaps_in_frame_num_value_allowed_flag");

    sps.pic_width_in_mbs_minus1        = this->ue("SPS: pic_width_in_mbs_minus1");
    sps.pic_height_in_map_units_minus1 = this->ue("SPS: pic_height_in_map_units_minus1");

    sps.PicWidthInMbs       = sps.pic_width_in_mbs_minus1 + 1;
    sps.PicWidthInSamplesL  = sps.PicWidthInMbs * 16;
    sps.PicWidthInSamplesC  = sps.PicWidthInMbs * sps.MbWidthC;
    sps.PicHeightInMapUnits = sps.pic_height_in_map_units_minus1 + 1;
    sps.PicSizeInMapUnits   = sps.PicWidthInMbs * sps.PicHeightInMapUnits;

    sps.frame_mbs_only_flag = this->u(1, "SPS: frame_mbs_only_flag");

    sps.FrameHeightInMbs = (2 - sps.frame_mbs_only_flag) * sps.PicHeightInMapUnits;

    sps.mb_adaptive_frame_field_flag = 0;
    if (!sps.frame_mbs_only_flag)
        sps.mb_adaptive_frame_field_flag = this->u(1, "SPS: mb_adaptive_frame_field_flag");

    sps.direct_8x8_inference_flag = this->u(1, "SPS: direct_8x8_inference_flag");

    assert(sps.frame_mbs_only_flag || sps.direct_8x8_inference_flag);

    sps.frame_cropping_flag = this->u(1, "SPS: frame_cropping_flag");
    if (sps.frame_cropping_flag) {
        sps.frame_crop_left_offset   = this->ue("SPS: frame_crop_left_offset");
        sps.frame_crop_right_offset  = this->ue("SPS: frame_crop_right_offset");
        sps.frame_crop_top_offset    = this->ue("SPS: frame_crop_top_offset");
        sps.frame_crop_bottom_offset = this->ue("SPS: frame_crop_bottom_offset");

        sps.CropUnitX = sps.SubWidthC;
        sps.CropUnitY = sps.SubHeightC * (2 - sps.frame_mbs_only_flag);

        assert(sps.frame_crop_left_offset <=
               (sps.PicWidthInSamplesL / sps.CropUnitX) - (sps.frame_crop_right_offset + 1));
        assert(sps.frame_crop_top_offset <=
               (16 * sps.FrameHeightInMbs / sps.CropUnitY) - (sps.frame_crop_bottom_offset + 1));
    }

    sps.vui_parameters_present_flag = this->u(1, "SPS: vui_parameters_present_flag");
    sps.vui_parameters.matrix_coefficients = 2;
    if (sps.vui_parameters_present_flag)
        this->vui_parameters(sps.vui_parameters);

    if (sps.vui_parameters_present_flag) {
        if (!(sps.BitDepthC == sps.BitDepthY && sps.chroma_format_idc == 3))
            assert(sps.vui_parameters.matrix_coefficients != 0);
        if (!(sps.BitDepthC == sps.BitDepthY ||
             (sps.BitDepthC == sps.BitDepthY + 1 && sps.chroma_format_idc == 3)))
            assert(sps.vui_parameters.matrix_coefficients != 8);
        if (sps.chroma_format_idc != 1)
            assert(!sps.vui_parameters.chroma_loc_info_present_flag);
    }

    int level_idc;
    bool level_1b = (sps.profile_idc == 66 || sps.profile_idc == 77 || sps.profile_idc == 88) &&
                    (sps.level_idc == 11 && sps.constraint_set3_flag);
    level_idc = level_1b ? 10 : sps.level_idc;
    level_idc = 3 * (level_idc / 10) + (level_idc % 10) - 3;
    level_idc = level_idc < 0 ? 0 : level_idc;
    assert(level_idc <= 14);

    sps.MaxDpbFrames =
        min<int>(MaxDpbMbs[level_idc] / (sps.PicWidthInMbs * sps.FrameHeightInMbs), 16);

    if (!sps.vui_parameters_present_flag || !sps.vui_parameters.bitstream_restriction_flag) {
        sps.vui_parameters.max_num_reorder_frames =
            (sps.profile_idc == 44 || sps.profile_idc == 86 || sps.profile_idc == 100 ||
             sps.profile_idc == 110 || sps.profile_idc == 122 || sps.profile_idc == 244) &&
            sps.constraint_set3_flag ? 0 : sps.MaxDpbFrames;
        sps.vui_parameters.max_dec_frame_buffering =
            (sps.profile_idc == 44 || sps.profile_idc == 86 || sps.profile_idc == 100 ||
             sps.profile_idc == 110 || sps.profile_idc == 122 || sps.profile_idc == 244) &&
            sps.constraint_set3_flag ? 0 : sps.MaxDpbFrames;
    }

    assert(sps.vui_parameters.max_num_reorder_frames <= sps.MaxDpbFrames);
    assert(sps.vui_parameters.max_dec_frame_buffering <= sps.MaxDpbFrames);
    assert(sps.max_num_ref_frames <= sps.MaxDpbFrames);

    sps.Valid = true;
}

static const uint8_t ZZ_SCAN_4x4[16] = {
     0,  1,  4,  8,  5,  2,  3,  6,
     9, 12, 13, 10,  7, 11, 14, 15
};

static const uint8_t ZZ_SCAN_8x8[64] = {
     0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};

// 7.3.2.1.1.1 Scaling list syntax

void data_partition_t::scaling_list(int* scalingList, int sizeOfScalingList, bool& useDefaultScalingMatrixFlag)
{
    int lastScale = 8;
    int nextScale = 8;

    for (int j = 0; j < sizeOfScalingList; j++) {
        int scanj = (sizeOfScalingList == 16) ? ZZ_SCAN_4x4[j] : ZZ_SCAN_8x8[j];
        if (nextScale != 0) {
            int8_t delta_scale = this->se("   : delta_sl");
            nextScale = (lastScale + delta_scale + 256) % 256;
            useDefaultScalingMatrixFlag = (scanj == 0 && nextScale == 0);
        }
        scalingList[scanj] = (nextScale == 0) ? lastScale : nextScale;
        lastScale = scalingList[scanj];
    }
}

// 7.3.2.1.2 Sequence parameter set extension RBSP syntax

void data_partition_t::seq_parameter_set_extension_rbsp(sps_ext_t& sps_ext)
{
    sps_ext.seq_parameter_set_id = this->ue("SPS_EXT: seq_parameter_set_id");
    sps_ext.aux_format_idc       = this->ue("SPS_EXT: aux_format_idc");

    assert(sps_ext.seq_parameter_set_id < MAX_NUM_SPS);
    assert(sps_ext.aux_format_idc <= 3);

    if (sps_ext.aux_format_idc != 0) {
        sps_ext.bit_depth_aux_minus8 = this->ue("SPS_EXT: bit_depth_aux_minus8");

        assert(sps_ext.bit_depth_aux_minus8 <= 4);

        uint8_t bit_depth_aux = sps_ext.bit_depth_aux_minus8 + 9;

        sps_ext.alpha_incr_flag         = this->u(1, "SPS_EXT: alpha_incr_flag");
        sps_ext.alpha_opaque_value      = this->u(bit_depth_aux, "SPS_EXT: alpha_opaque_value");
        sps_ext.alpha_transparent_value = this->u(bit_depth_aux, "SPS_EXT: alpha_transparent_value");
    }
    sps_ext.additional_extension_flag = this->u(1, "SPS_EXT: additional_extension_flag");

    this->rbsp_trailing_bits();
}

// 7.3.2.1.3 Subset sequence parameter set RBSP syntax

void data_partition_t::subset_seq_parameter_set_rbsp(sub_sps_t& sub_sps)
{
    this->seq_parameter_set_data(sub_sps.sps);

    if (sub_sps.sps.profile_idc == 83 || sub_sps.sps.profile_idc == 86) {
        this->seq_parameter_set_svc_extension(sub_sps.sps_svc);
        sub_sps.svc_vui_parameters_present_flag = this->u(1, "svc_vui_parameters_present_flag");
        if (sub_sps.svc_vui_parameters_present_flag)
            this->svc_vui_parameters_extension(sub_sps.svc_vui_parameters);
    } else if (sub_sps.sps.profile_idc == 118 || sub_sps.sps.profile_idc == 128) {
        sub_sps.bit_equal_to_one = this->u(1, "bit_equal_to_one");
        assert(sub_sps.bit_equal_to_one == 1);
        this->seq_parameter_set_mvc_extension(sub_sps.sps_mvc);
        sub_sps.mvc_vui_parameters_present_flag = this->u(1, "mvc_vui_parameters_present_flag");
        if (sub_sps.mvc_vui_parameters_present_flag)
            this->mvc_vui_parameters_extension(sub_sps.mvc_vui_parameters);
    } else if (sub_sps.sps.profile_idc == 138) {
        sub_sps.bit_equal_to_one = this->u(1, "bit_equal_to_one");
        assert(sub_sps.bit_equal_to_one == 1);
        this->seq_parameter_set_mvcd_extension(sub_sps.sps_mvcd);
    }

    bool additional_extension2_flag = this->u(1, "additional_extension2_flag");
    if (additional_extension2_flag) {
        bool additional_extension2_data_flag;
        while (this->more_rbsp_data())
            additional_extension2_data_flag = this->u(1, "additional_extension2_flag");
    }

    this->rbsp_trailing_bits();

    sub_sps.Valid = sub_sps.sps.Valid;
}

// 7.3.2.2 Picture parameter set RBSP syntax

void data_partition_t::pic_parameter_set_rbsp(VideoParameters* p_Vid, pps_t& pps)
{
    pps.pic_parameter_set_id                         = this->ue("PPS: pic_parameter_set_id");
    pps.seq_parameter_set_id                         = this->ue("PPS: seq_parameter_set_id");
    pps.entropy_coding_mode_flag                     = this->u(1, "PPS: entropy_coding_mode_flag");
    pps.bottom_field_pic_order_in_frame_present_flag = this->u(1, "PPS: bottom_field_pic_order_in_frame_present_flag");

    assert(pps.seq_parameter_set_id < MAX_NUM_SPS);

    const sps_t& sps = p_Vid->SeqParSet[pps.seq_parameter_set_id];

    if (!sps.Valid) {
        pps.Valid = false;
        return;
    }

    pps.num_slice_groups_minus1 = this->ue("PPS: num_slice_groups_minus1");
    assert(pps.num_slice_groups_minus1 < MAX_NUM_SLICE_GROUPS);
    if (pps.num_slice_groups_minus1 > 0) {
        pps.slice_group_map_type = this->ue("PPS: slice_group_map_type");

        assert(pps.slice_group_map_type <= 6);
        assert((pps.num_slice_groups_minus1 == 1) ||
               (pps.slice_group_map_type != 3 && pps.slice_group_map_type != 4 && pps.slice_group_map_type != 5));

        if (pps.slice_group_map_type == 0) {
            pps.slice_groups = pps_t::slice_group_v(pps.num_slice_groups_minus1 + 1);
            for (int iGroup = 0; iGroup <= pps.num_slice_groups_minus1; ++iGroup) {
                auto& slice_group = pps.slice_groups[iGroup];
                slice_group.run_length_minus1 = this->ue("PPS: run_length_minus1 [i]");

                assert(slice_group.run_length_minus1 < sps.PicSizeInMapUnits);
            }
        } else if (pps.slice_group_map_type == 2) {
            pps.slice_groups = pps_t::slice_group_v(pps.num_slice_groups_minus1 + 1);
            for (int iGroup = 0; iGroup < pps.num_slice_groups_minus1; ++iGroup) {
                auto& slice_group = pps.slice_groups[iGroup];
                slice_group.top_left     = this->ue("PPS: top_left [i]");
                slice_group.bottom_right = this->ue("PPS: bottom_right [i]");

                assert(slice_group.top_left <= slice_group.bottom_right);
                assert(slice_group.bottom_right < sps.PicSizeInMapUnits);
                assert(slice_group.top_left % sps.PicWidthInMbs <= slice_group.bottom_right % sps.PicWidthInMbs);
            }
        } else if (pps.slice_group_map_type == 3 ||
                   pps.slice_group_map_type == 4 ||
                   pps.slice_group_map_type == 5) {
            pps.slice_group_change_direction_flag = this->u(1, "PPS: slice_group_change_direction_flag");
            pps.slice_group_change_rate_minus1    = this->ue("PPS: slice_group_change_rate_minus1");

            assert(pps.slice_group_change_rate_minus1 < sps.PicSizeInMapUnits);

            pps.SliceGroupChangeRate = pps.slice_group_change_rate_minus1 + 1;
        } else if (pps.slice_group_map_type == 6) {
            pps.pic_size_in_map_units_minus1 = this->ue("PPS: pic_size_in_map_units_minus1");

            assert(pps.pic_size_in_map_units_minus1 < sps.PicSizeInMapUnits);

            pps.slice_group_id = pps_t::uint8_v(pps.pic_size_in_map_units_minus1 + 1);
            int bitsSliceGroupId = ceil(log2(pps.num_slice_groups_minus1 + 1));
            for (int i = 0; i <= pps.pic_size_in_map_units_minus1; ++i) {
                pps.slice_group_id[i] = this->u(bitsSliceGroupId, "slice_group_id[i]");

                assert(pps.slice_group_id[i] <= pps.num_slice_groups_minus1);
            }
        }
    }

    pps.num_ref_idx_l0_default_active_minus1 = this->ue("PPS: num_ref_idx_l0_default_active_minus1");
    pps.num_ref_idx_l1_default_active_minus1 = this->ue("PPS: num_ref_idx_l1_default_active_minus1");

    assert(pps.num_ref_idx_l0_default_active_minus1 < MAX_NUM_REF_IDX);
    assert(pps.num_ref_idx_l1_default_active_minus1 < MAX_NUM_REF_IDX);

    pps.weighted_pred_flag  = this->u(1, "PPS: weighted_pred_flag");
    pps.weighted_bipred_idc = this->u(2, "PPS: weighted_bipred_idc");

    assert(pps.weighted_bipred_idc <= 2);

    pps.pic_init_qp_minus26    = this->se("PPS: pic_init_qp_minus26");
    pps.pic_init_qs_minus26    = this->se("PPS: pic_init_qs_minus26");
    pps.chroma_qp_index_offset = this->se("PPS: chroma_qp_index_offset");

    assert(pps.pic_init_qp_minus26 >= -(26 + sps.QpBdOffsetY) && pps.pic_init_qp_minus26 <= 25);
    assert(pps.pic_init_qs_minus26 >= -26 && pps.pic_init_qs_minus26 <= 25);
    assert(pps.chroma_qp_index_offset >= -12 && pps.chroma_qp_index_offset <= 12);

    pps.deblocking_filter_control_present_flag = this->u(1, "PPS: deblocking_filter_control_present_flag");
    pps.constrained_intra_pred_flag            = this->u(1, "PPS: constrained_intra_pred_flag");
    pps.redundant_pic_cnt_present_flag         = this->u(1, "PPS: redundant_pic_cnt_present_flag");

    pps.transform_8x8_mode_flag         = 0;
    pps.pic_scaling_matrix_present_flag = 0;
    pps.second_chroma_qp_index_offset   = pps.chroma_qp_index_offset;
    if (this->more_rbsp_data()) {
        pps.transform_8x8_mode_flag         = this->u(1, "PPS: transform_8x8_mode_flag");
        pps.pic_scaling_matrix_present_flag = this->u(1, "PPS: pic_scaling_matrix_present_flag");
        if (pps.pic_scaling_matrix_present_flag) {
            for (int i = 0; i < 6 + ((sps.chroma_format_idc != CHROMA_FORMAT_444) ? 2 : 6) * pps.transform_8x8_mode_flag; ++i) {
                pps.pic_scaling_list_present_flag[i] = this->u(1, "PPS: pic_scaling_list_present_flag");
                if (pps.pic_scaling_list_present_flag[i]) {
                    if (i < 6)
                        this->scaling_list(pps.ScalingList4x4[i], 16,
                                           pps.UseDefaultScalingMatrix4x4Flag[i]);
                    else
                        this->scaling_list(pps.ScalingList8x8[i - 6], 64,
                                           pps.UseDefaultScalingMatrix8x8Flag[i - 6]);
                }
            }
        }
        pps.second_chroma_qp_index_offset = this->se("PPS: second_chroma_qp_index_offset");
    }
    assert(pps.second_chroma_qp_index_offset >= -12 && pps.second_chroma_qp_index_offset <= 12);

    this->rbsp_trailing_bits();

    pps.Valid = true;
}

// 7.3.2.4 Access unit delimiter RBSP syntax

void data_partition_t::access_unit_delimiter_rbsp(void)
{
    uint8_t primary_pic_type;
    primary_pic_type = this->u(3);

    this->rbsp_trailing_bits();
}

// 7.3.2.5 End of sequence RBSP syntax

void data_partition_t::end_of_seq_rbsp(void)
{
}

// 7.3.2.6 End of stream RBSP syntax

void data_partition_t::end_of_stream_rbsp(void)
{
}

// 7.3.2.7 Filler data RBSP syntax

void data_partition_t::filler_data_rbsp(void)
{
    uint8_t ff_byte;
    while (this->next_bits(8) == 0xFF)
        ff_byte = this->f(8);

    this->rbsp_trailing_bits();
}

// 7.3.2.8 Slice layer without partitioning RBSP syntax

void data_partition_t::slice_layer_without_partitioning_rbsp(void)
{
    slice_t slice;
    this->slice_header(slice);
    this->slice_data();
    this->rbsp_trailing_bits();
}

// 7.3.2.9.1 Slice data partition A RBSP syntax

void data_partition_t::slice_data_partition_a_layer_rbsp(void)
{
    slice_t slice;
    this->slice_header(slice);
//    uint32_t slice_id = this->ue();
    this->slice_data();
    this->rbsp_trailing_bits();
}

// 7.3.2.9.2 Slice data partition B RBSP syntax

void data_partition_t::slice_data_partition_b_layer_rbsp(void)
{
//    uint32_t slice_id = this->ue();
//    if (slice.seperate_colour_plane_flag)
//        colour_plane_id = this->u(2);
//    if (slice.redundant_pic_cnt_present_flag)
//        redundant_pic_cnt = this->ue();
    this->slice_data();
    this->rbsp_trailing_bits();
}

// 7.3.2.9.3 Slice data partition C RBSP syntax

void data_partition_t::slice_data_partition_c_layer_rbsp(void)
{
//    uint32_t slice_id = this->ue();
//    if (slice.seperate_colour_plane_flag)
//        colour_plane_id = this->u(2);
//    if (slice.redundant_pic_cnt_present_flag)
//        redundant_pic_cnt = this->ue();
    this->slice_data();
    this->rbsp_trailing_bits();
}

// 7.3.2.10 RBSP slice trailing bits syntax

void data_partition_t::rbsp_slice_trailing_bits(void)
{
    this->rbsp_trailing_bits();

//    if (slice.entropy_coding_mode_flag) {
//        while (this->more_rbsp_trailing_data())
//            cabac_zero_word = this->f(16);
//    }
}

// 7.3.2.11 RBSP trailing bits syntax

void data_partition_t::rbsp_trailing_bits(void)
{
    bool rbsp_stop_one_bit;
    bool rbsp_alignment_zero_bit;

    rbsp_stop_one_bit = this->f(1);
    while (!this->byte_aligned())
        rbsp_alignment_zero_bit = this->f(1);
}

// 7.3.2.12 Prefix NAL unit RBSP syntax

void data_partition_t::prefix_nal_unit_rbsp(void)
{
    //to be implemented for Annex G;
}


// E.1.1 VUI parameter syntax

void data_partition_t::vui_parameters(vui_t& vui)
{
    #define Extended_SAR 255

    vui.aspect_ratio_idc = 0;
    vui.aspect_ratio_info_present_flag = this->u(1, "VUI: aspect_ratio_info_present_flag");
    if (vui.aspect_ratio_info_present_flag) {
        vui.aspect_ratio_idc = this->u(8, "VUI: aspect_ratio_idc");
        if (vui.aspect_ratio_idc == Extended_SAR) {
            vui.sar_width  = this->u(16, "VUI: sar_width");
            vui.sar_height = this->u(16, "VUI: sar_height");
        }
    }

    vui.overscan_info_present_flag = this->u(1, "VUI: overscan_info_present_flag");
    if (vui.overscan_info_present_flag)
        vui.overscan_appropriate_flag  = this->u(1, "VUI: overscan_appropriate_flag");

    vui.video_format             = 5; // Unspecified video format
    vui.video_full_range_flag    = 0;
    vui.colour_primaries         = 2; // Unspecified
    vui.transfer_characteristics = 2; // Unspecified
    vui.matrix_coefficients      = 2; // Unspecified
    vui.video_signal_type_present_flag = this->u(1, "VUI: video_signal_type_present_flag");
    if (vui.video_signal_type_present_flag) {
        vui.video_format          = this->u(3, "VUI: video_format");
        vui.video_full_range_flag = this->u(1, "VUI: video_full_range_flag");
        vui.colour_description_present_flag = this->u(1, "VUI: color_description_present_flag");
        if (vui.colour_description_present_flag) {
            vui.colour_primaries         = this->u(8, "VUI: colour_primaries");
            vui.transfer_characteristics = this->u(8, "VUI: transfer_characteristics");
            vui.matrix_coefficients      = this->u(8, "VUI: matrix_coefficients");
        }
    }

    vui.chroma_sample_loc_type_top_field    = 0;
    vui.chroma_sample_loc_type_bottom_field = 0;
    vui.chroma_loc_info_present_flag = this->u(1, "VUI: chroma_loc_info_present_flag");
    if (vui.chroma_loc_info_present_flag) {
        vui.chroma_sample_loc_type_top_field    = this->ue("VUI: chroma_sample_loc_type_top_field");
        vui.chroma_sample_loc_type_bottom_field = this->ue("VUI: chroma_sample_loc_type_bottom_field");
    }

    assert(vui.chroma_sample_loc_type_top_field    <= 5);
    assert(vui.chroma_sample_loc_type_bottom_field <= 5);

    vui.fixed_frame_rate_flag = 0;
    vui.timing_info_present_flag = this->u(1, "VUI: timing_info_present_flag");
    if (vui.timing_info_present_flag) {
        vui.num_units_in_tick     = this->u(32, "VUI: num_units_in_tick");
        vui.time_scale            = this->u(32, "VUI: time_scale");
        vui.fixed_frame_rate_flag = this->u(1, "VUI: fixed_frame_rate_flag");
    }

    assert(!vui.timing_info_present_flag || vui.num_units_in_tick > 0);
    assert(!vui.timing_info_present_flag || vui.time_scale > 0);

    vui.nal_hrd_parameters_present_flag = this->u(1, "VUI: nal_hrd_parameters_present_flag");
    if (vui.nal_hrd_parameters_present_flag)
        this->hrd_parameters(vui.nal_hrd_parameters);

    vui.vcl_hrd_parameters_present_flag = this->u(1, "VUI: vcl_hrd_parameters_present_flag");
    if (vui.vcl_hrd_parameters_present_flag)
        this->hrd_parameters(vui.vcl_hrd_parameters);

    vui.low_delay_hrd_flag = 1 - vui.fixed_frame_rate_flag;
    if (vui.nal_hrd_parameters_present_flag || vui.vcl_hrd_parameters_present_flag)
        vui.low_delay_hrd_flag = this->u(1, "VUI: low_delay_hrd_flag");

    assert(!vui.fixed_frame_rate_flag || !vui.low_delay_hrd_flag);

    vui.motion_vectors_over_pic_boundaries_flag = 1;
    vui.max_bytes_per_pic_denom                 = 2;
    vui.max_bits_per_mb_denom                   = 1;
    vui.log2_max_mv_length_horizontal           = 16;
    vui.log2_max_mv_length_vertical             = 16;

    vui.pic_struct_present_flag    = this->u(1, "VUI: pic_struct_present_flag");
    vui.bitstream_restriction_flag = this->u(1, "VUI: bitstream_restriction_flag");
    if (vui.bitstream_restriction_flag) {
        vui.motion_vectors_over_pic_boundaries_flag = this->u(1, "VUI: motion_vectors_over_pic_boundaries_flag");
        vui.max_bytes_per_pic_denom                 = this->ue("VUI: max_bytes_per_pic_denom");
        vui.max_bits_per_mb_denom                   = this->ue("VUI: max_bits_per_mb_denom");
        vui.log2_max_mv_length_horizontal           = this->ue("VUI: log2_max_mv_length_horizontal");
        vui.log2_max_mv_length_vertical             = this->ue("VUI: log2_max_mv_length_vertical");
        vui.max_num_reorder_frames                  = this->ue("VUI: num_reorder_frames");
        vui.max_dec_frame_buffering                 = this->ue("VUI: max_dec_frame_buffering");
    }

    assert(vui.max_bytes_per_pic_denom       <= 16);
    assert(vui.max_bits_per_mb_denom         <= 16);
    assert(vui.log2_max_mv_length_horizontal <= 16);
    assert(vui.log2_max_mv_length_vertical   <= 16);
    assert(vui.max_num_reorder_frames        <= vui.max_dec_frame_buffering);
}

// E.1.2 HRD parameters syntax

void data_partition_t::hrd_parameters(hrd_t& hrd)
{
    #define MAX_NUM_HRD_CPB 32

    hrd.cpb_cnt_minus1 = this->ue("VUI: cpb_cnt_minus1");
    hrd.bit_rate_scale = this->u(4, "VUI: bit_rate_scale");
    hrd.cpb_size_scale = this->u(4, "VUI: cpb_size_scale");

    assert(hrd.cpb_cnt_minus1 < MAX_NUM_HRD_CPB);

    // if (profile_idc == 66, 77, 88)
    //     BitRate[SchedSelIdx] = 1200 * MaxBR;
    //     CpbSize[SchedSelIdx] = 1000 * MaxCPB;
    // else
    //     BitRate[SchedSelIdx] = cpbBrVcl(Nal)Factor * MaxBR
    //     CpbSize[SchedSelIdx] = cpbBrVcl(Nal)Factor * MaxCPB

    hrd.cpbs = hrd_t::cpb_v(hrd.cpb_cnt_minus1 + 1);
    for (int SchedSelIdx = 0; SchedSelIdx <= hrd.cpb_cnt_minus1; ++SchedSelIdx) {
        hrd_t::cpb_t& cpb = hrd.cpbs[SchedSelIdx];
        cpb.bit_rate_value_minus1 = this->ue("VUI: bit_rate_value_minus1");
        cpb.cpb_size_value_minus1 = this->ue("VUI: cpb_size_value_minus1");
        cpb.cbr_flag              = this->u(1, "VUI: cbr_flag");
        // BitRate[SchedSelIdx] = (bit_rate_value_minus1[SchedSelIdx] + 1) * (1 << (6 + bit_rate_scale));
        // CpbSize[SchedSelIdx] = (cpb_size_value_minus1[SchedSelIdx] + 1) * (1 << (4 + cpb_size_scale));
    }

    hrd.initial_cpb_removal_delay_length_minus1 = 23;
    hrd.cpb_removal_delay_length_minus1         = 23;
    hrd.dpb_output_delay_length_minus1          = 23;
    hrd.time_offset_length                      = 24;
    hrd.initial_cpb_removal_delay_length_minus1 = this->u(5, "VUI: initial_cpb_removal_delay_length_minus1");
    hrd.cpb_removal_delay_length_minus1         = this->u(5, "VUI: cpb_removal_delay_length_minus1");
    hrd.dpb_output_delay_length_minus1          = this->u(5, "VUI: dpb_output_delay_length_minus1");
    hrd.time_offset_length                      = this->u(5, "VUI: time_offset_length");
}

// G.7.3.1.4 Sequence parameter set SVC extension syntax

void data_partition_t::seq_parameter_set_svc_extension(sps_svc_t& sps_svc)
{
    //to be implemented for Annex G;
}

// G.14.1 SVC VUI parameters extension syntax

void data_partition_t::svc_vui_parameters_extension(svc_vui_t& svc_vui)
{
    //to be implemented for Annex G;
    #define MAX_NUM_SVC_VUI_ENTRIES 1024
}

// H.7.3.2.1.4 Sequence parameter set MVC extension syntax

void data_partition_t::seq_parameter_set_mvc_extension(sps_mvc_t& sps_mvc)
{
    #define MAX_NUM_SPS_MVC_VIEWS                     1024
    #define MAX_NUM_SPS_MVC_ANCHOR                      16
    #define MAX_NUM_SPS_MVC_LEVEL_VALUES_SIGNALLED      64
    #define MAX_NUM_SPS_MVC_APPLICABLE_OPS            1024
    #define MAX_NUM_SPS_MVC_APPLICABLE_OP_TARGET_VIEW 1024

    sps_mvc.num_views_minus1 = this->ue("num_views_minus1");
    assert(sps_mvc.num_views_minus1 < MAX_NUM_SPS_MVC_VIEWS);
    sps_mvc.views = sps_mvc_t::view_v(sps_mvc.num_views_minus1 + 1);

    for (int i = 0; i <= sps_mvc.num_views_minus1; ++i) {
        sps_mvc.views[i].view_id = this->ue("view_id");
        assert(sps_mvc.views[i].view_id < MAX_NUM_SPS_MVC_VIEWS);
    }

    for (int i = 1; i <= sps_mvc.num_views_minus1; ++i) {
        sps_mvc_t::view_t& view = sps_mvc.views[i];

        view.num_anchor_refs_l0 = this->ue("num_anchor_refs_l0");
        assert(view.num_anchor_refs_l0 <= min<int>(15, sps_mvc.num_views_minus1));
        if (view.num_anchor_refs_l0 > 0) {
            view.anchor_ref_l0 = sps_mvc_t::uint16_v(view.num_anchor_refs_l0);
            for (int j = 0; j < view.num_anchor_refs_l0; ++j) {
                view.anchor_ref_l0[j] = this->ue("anchor_ref_l0");
                assert(view.anchor_ref_l0[j] < MAX_NUM_SPS_MVC_VIEWS);
            }
        }

        view.num_anchor_refs_l1 = this->ue("num_anchor_refs_l1");
        assert(view.num_anchor_refs_l1 <= min<int>(15, sps_mvc.num_views_minus1));
        if (view.num_anchor_refs_l1 > 0) {
            view.anchor_ref_l1 = sps_mvc_t::uint16_v(view.num_anchor_refs_l1);
            for (int j = 0; j < view.num_anchor_refs_l1; ++j) {
                view.anchor_ref_l1[j] = this->ue("anchor_ref_l1");
                assert(view.anchor_ref_l1[j] < MAX_NUM_SPS_MVC_VIEWS);
            }
        }
    }

    for (int i = 1; i <= sps_mvc.num_views_minus1; ++i) {
        sps_mvc_t::view_t& view = sps_mvc.views[i];

        view.num_non_anchor_refs_l0 = this->ue("num_non_anchor_refs_l0");
        assert(view.num_non_anchor_refs_l0 <= min<int>(15, sps_mvc.num_views_minus1));
        if (view.num_non_anchor_refs_l0 > 0) {
            view.non_anchor_ref_l0 = sps_mvc_t::uint16_v(view.num_non_anchor_refs_l0);
            for (int j = 0; j < view.num_non_anchor_refs_l0; ++j) {
                view.non_anchor_ref_l0[j] = this->ue("non_anchor_ref_l0");
                assert(view.non_anchor_ref_l0[j] < MAX_NUM_SPS_MVC_VIEWS);
            }
        }

        view.num_non_anchor_refs_l1 = this->ue("num_non_anchor_refs_l1");
        assert(view.num_non_anchor_refs_l1 <= min<int>(15, sps_mvc.num_views_minus1));
        if (view.num_non_anchor_refs_l1 > 0) {
            view.non_anchor_ref_l1 = sps_mvc_t::uint16_v(view.num_non_anchor_refs_l1);
            for (int j = 0; j < view.num_non_anchor_refs_l1; ++j) {
                view.non_anchor_ref_l1[j] = this->ue("non_anchor_ref_l1");
                assert(view.non_anchor_ref_l1[j] < MAX_NUM_SPS_MVC_VIEWS);
            }
        }
    }

    sps_mvc.num_level_values_signalled_minus1 = this->ue("num_level_values_signalled_minus1");
    assert(sps_mvc.num_level_values_signalled_minus1 < MAX_NUM_SPS_MVC_LEVEL_VALUES_SIGNALLED);
    sps_mvc.level_values_signalled =
        sps_mvc_t::level_value_v(sps_mvc.num_level_values_signalled_minus1 + 1);

    for (int i = 0; i <= sps_mvc.num_level_values_signalled_minus1; ++i) {
        sps_mvc_t::level_value_t& level_value = sps_mvc.level_values_signalled[i];
        level_value.level_idc                 = this->u(8, "level_idc");
        level_value.num_applicable_ops_minus1 = this->ue("num_applicable_ops_minus1");
        assert(level_value.num_applicable_ops_minus1 < MAX_NUM_SPS_MVC_APPLICABLE_OPS);
        level_value.applicable_ops =
            sps_mvc_t::level_value_t::applicable_op_v(level_value.num_applicable_ops_minus1 + 1);

        for (int j = 0; j <= level_value.num_applicable_ops_minus1; ++j) {
            sps_mvc_t::level_value_t::applicable_op_t& applicable_op = level_value.applicable_ops[j];
            applicable_op.applicable_op_temporal_id = this->u(3, "applicable_op_temporal_id");
            applicable_op.applicable_op_num_target_views_minus1 = this->ue("applicable_op_num_target_views_minus1");
            assert(applicable_op.applicable_op_num_target_views_minus1 < MAX_NUM_SPS_MVC_APPLICABLE_OP_TARGET_VIEW);
            applicable_op.applicable_op_target_view_id =
                sps_mvc_t::uint16_v(applicable_op.applicable_op_num_target_views_minus1 + 1);
            for (int k = 0; k <= applicable_op.applicable_op_num_target_views_minus1; ++k) {
                applicable_op.applicable_op_target_view_id[k] = this->ue("applicable_op_target_view_id");
                assert(applicable_op.applicable_op_target_view_id[k] < MAX_NUM_SPS_MVC_APPLICABLE_OP_TARGET_VIEW);
            }
            applicable_op.applicable_op_num_views_minus1 = this->ue("applicable_op_num_views_minus1");
            assert(applicable_op.applicable_op_num_views_minus1 < MAX_NUM_SPS_MVC_APPLICABLE_OP_TARGET_VIEW);
        }
    }
}

// H.14.1 MVC VUI parameters extension syntax

void data_partition_t::mvc_vui_parameters_extension(mvc_vui_t& mvc_vui)
{
    #define MAX_NUM_MVC_VUI_OPS   1024
    #define MAX_NUM_MVC_VUI_VIEWS 1024

    mvc_vui.vui_mvc_num_ops_minus1 = this->ue("vui_mvc_num_ops_minus1");
    assert(mvc_vui.vui_mvc_num_ops_minus1 < MAX_NUM_MVC_VUI_OPS);
    mvc_vui.vui_mvc_ops = mvc_vui_t::vui_mvc_op_v(mvc_vui.vui_mvc_num_ops_minus1 + 1);

    for (int i = 0; i <= mvc_vui.vui_mvc_num_ops_minus1; ++i) {
        mvc_vui_t::vui_mvc_op_t& vui_mvc_op = mvc_vui.vui_mvc_ops[i];

        vui_mvc_op.vui_mvc_temporal_id = this->u(3, "vui_mvc_temporal_id");
        vui_mvc_op.vui_mvc_num_target_output_views_minus1 = this->ue("vui_mvc_num_target_output_views_minus1");
        assert(vui_mvc_op.vui_mvc_num_target_output_views_minus1 < MAX_NUM_MVC_VUI_VIEWS);
        vui_mvc_op.vui_mvc_view_id = mvc_vui_t::uint16_v(vui_mvc_op.vui_mvc_num_target_output_views_minus1 + 1);

        for (int j = 0; j <= vui_mvc_op.vui_mvc_num_target_output_views_minus1; ++j) {
            vui_mvc_op.vui_mvc_view_id[j] = this->ue("vui_mvc_view_id");
            assert(vui_mvc_op.vui_mvc_view_id[j] < MAX_NUM_MVC_VUI_VIEWS);
        }

        vui_mvc_op.vui_mvc_timing_info_present_flag = this->u(1, "vui_mvc_timing_info_present_flag");
        if (vui_mvc_op.vui_mvc_timing_info_present_flag) {
            vui_mvc_op.vui_mvc_num_units_in_tick     = this->u(32, "vui_mvc_num_units_in_tick");
            vui_mvc_op.vui_mvc_time_scale            = this->u(32, "vui_mvc_time_scale");
            vui_mvc_op.vui_mvc_fixed_frame_rate_flag = this->u(1, "vui_mvc_fixed_frame_rate_flag");
        }

        vui_mvc_op.vui_mvc_nal_hrd_parameters_present_flag = this->u(1, "vui_mvc_nal_hrd_parameters_present_flag");
        if (vui_mvc_op.vui_mvc_nal_hrd_parameters_present_flag)
            this->hrd_parameters(vui_mvc_op.vui_mvc_nal_hrd_parameters);
        vui_mvc_op.vui_mvc_vcl_hrd_parameters_present_flag = this->u(1, "vui_mvc_vcl_hrd_parameters_present_flag");
        if (vui_mvc_op.vui_mvc_vcl_hrd_parameters_present_flag)
            this->hrd_parameters(vui_mvc_op.vui_mvc_vcl_hrd_parameters);
        if (vui_mvc_op.vui_mvc_nal_hrd_parameters_present_flag ||
            vui_mvc_op.vui_mvc_vcl_hrd_parameters_present_flag)
            vui_mvc_op.vui_mvc_low_delay_hrd_flag = this->u(1, "vui_mvc_low_delay_hrd_flag");
        vui_mvc_op.vui_mvc_pic_struct_present_flag = this->u(1, "vui_mvc_pic_struct_present_flag");
    }
}

// I.7.3.2.1.5 Sequence parameter set MVCD extension syntax

void data_partition_t::seq_parameter_set_mvcd_extension(sps_mvcd_t& sps_mvcd)
{
    //to be implemented for Annex I;
}


bool operator==(const sps_t& l, const sps_t& r)
{
    if (!l.Valid || !r.Valid)
        return false;

    bool equal = true;

    equal &= (l.profile_idc               == r.profile_idc);
    equal &= (l.constraint_set0_flag      == r.constraint_set0_flag);
    equal &= (l.constraint_set1_flag      == r.constraint_set1_flag);
    equal &= (l.constraint_set2_flag      == r.constraint_set2_flag);
    equal &= (l.level_idc                 == r.level_idc);
    equal &= (l.seq_parameter_set_id      == r.seq_parameter_set_id);
    equal &= (l.log2_max_frame_num_minus4 == r.log2_max_frame_num_minus4);
    equal &= (l.pic_order_cnt_type        == r.pic_order_cnt_type);
    if (!equal)
        return false;

    if (l.pic_order_cnt_type == 0)
        equal &= (l.log2_max_pic_order_cnt_lsb_minus4 == r.log2_max_pic_order_cnt_lsb_minus4);
    else if (l.pic_order_cnt_type == 1) {
        equal &= (l.delta_pic_order_always_zero_flag == r.delta_pic_order_always_zero_flag);
        equal &= (l.offset_for_non_ref_pic == r.offset_for_non_ref_pic);
        equal &= (l.offset_for_top_to_bottom_field == r.offset_for_top_to_bottom_field);
        equal &= (l.num_ref_frames_in_pic_order_cnt_cycle == r.num_ref_frames_in_pic_order_cnt_cycle);
        if (!equal)
            return false;

        for (int i = 0; i < l.num_ref_frames_in_pic_order_cnt_cycle; ++i)
            equal &= (l.offset_for_ref_frame[i] == r.offset_for_ref_frame[i]);
    }

    equal &= (l.max_num_ref_frames                   == r.max_num_ref_frames);
    equal &= (l.gaps_in_frame_num_value_allowed_flag == r.gaps_in_frame_num_value_allowed_flag);
    equal &= (l.pic_width_in_mbs_minus1              == r.pic_width_in_mbs_minus1);
    equal &= (l.pic_height_in_map_units_minus1       == r.pic_height_in_map_units_minus1);
    equal &= (l.frame_mbs_only_flag                  == r.frame_mbs_only_flag);
    if (!equal)
        return false;
  
    if (!l.frame_mbs_only_flag)
        equal &= (l.mb_adaptive_frame_field_flag == r.mb_adaptive_frame_field_flag);

    equal &= (l.direct_8x8_inference_flag == r.direct_8x8_inference_flag);
    equal &= (l.frame_cropping_flag       == r.frame_cropping_flag);
    if (!equal)
        return false;

    if (l.frame_cropping_flag) {
        equal &= (l.frame_crop_left_offset   == r.frame_crop_left_offset);
        equal &= (l.frame_crop_right_offset  == r.frame_crop_right_offset);
        equal &= (l.frame_crop_top_offset    == r.frame_crop_top_offset);
        equal &= (l.frame_crop_bottom_offset == r.frame_crop_bottom_offset);
    }
    equal &= (l.vui_parameters_present_flag == r.vui_parameters_present_flag);

    return equal;
}

bool operator==(const pps_t& l, const pps_t& r)
{
    if (!l.Valid || !r.Valid)
        return false;

    bool equal = true;

    equal &= (l.pic_parameter_set_id     == r.pic_parameter_set_id);
    equal &= (l.seq_parameter_set_id     == r.seq_parameter_set_id);
    equal &= (l.entropy_coding_mode_flag == r.entropy_coding_mode_flag);
    equal &= (l.bottom_field_pic_order_in_frame_present_flag == r.bottom_field_pic_order_in_frame_present_flag);
    equal &= (l.num_slice_groups_minus1  == r.num_slice_groups_minus1);

    if (!equal)
        return false;

    if (l.num_slice_groups_minus1 > 0) {
        equal &= (l.slice_group_map_type == r.slice_group_map_type);
        if (!equal)
            return false;

        if (l.slice_group_map_type == 0) {
            for (int i = 0; i <= l.num_slice_groups_minus1; ++i)
                equal &= (l.slice_groups[i].run_length_minus1 == r.slice_groups[i].run_length_minus1);
        } else if (l.slice_group_map_type == 2) {
            for (int i = 0; i < l.num_slice_groups_minus1; ++i) {
                equal &= (l.slice_groups[i].top_left     == r.slice_groups[i].top_left);
                equal &= (l.slice_groups[i].bottom_right == r.slice_groups[i].bottom_right);
            }
        } else if (l.slice_group_map_type == 3 || l.slice_group_map_type == 4 || l.slice_group_map_type == 5) {
            equal &= (l.slice_group_change_direction_flag == r.slice_group_change_direction_flag);
            equal &= (l.slice_group_change_rate_minus1 == r.slice_group_change_rate_minus1);
        } else if (l.slice_group_map_type == 6) {
            equal &= (l.pic_size_in_map_units_minus1 == r.pic_size_in_map_units_minus1);
            if (!equal)
                return false;

            for (int i = 0; i <= l.pic_size_in_map_units_minus1; ++i)
                equal &= (l.slice_group_id[i] == r.slice_group_id[i]);
        }
    }

    equal &= (l.num_ref_idx_l0_default_active_minus1 == r.num_ref_idx_l0_default_active_minus1);
    equal &= (l.num_ref_idx_l1_default_active_minus1 == r.num_ref_idx_l1_default_active_minus1);
    equal &= (l.weighted_pred_flag     == r.weighted_pred_flag);
    equal &= (l.weighted_bipred_idc    == r.weighted_bipred_idc);
    equal &= (l.pic_init_qp_minus26    == r.pic_init_qp_minus26);
    equal &= (l.pic_init_qs_minus26    == r.pic_init_qs_minus26);
    equal &= (l.chroma_qp_index_offset == r.chroma_qp_index_offset);
    equal &= (l.deblocking_filter_control_present_flag == r.deblocking_filter_control_present_flag);
    equal &= (l.constrained_intra_pred_flag == r.constrained_intra_pred_flag);
    equal &= (l.redundant_pic_cnt_present_flag == r.redundant_pic_cnt_present_flag);
    if (!equal)
        return false;

    //Fidelity Range Extensions Stuff
    //It is initialized to zero, so should be ok to check all the time.
    equal &= (l.transform_8x8_mode_flag == r.transform_8x8_mode_flag);
    equal &= (l.pic_scaling_matrix_present_flag == r.pic_scaling_matrix_present_flag);
    if (l.pic_scaling_matrix_present_flag) {
        for (int i = 0; i < 6 + (l.transform_8x8_mode_flag << 1); ++i) {
            equal &= (l.pic_scaling_list_present_flag[i] == r.pic_scaling_list_present_flag[i]);
            if (l.pic_scaling_list_present_flag[i]) {
                if (i < 6) {
                    for (int j = 0; j < 16; ++j)
                        equal &= (l.ScalingList4x4[i][j] == r.ScalingList4x4[i][j]);
                } else {
                    for (int j = 0; j < 64; ++j)
                        equal &= (l.ScalingList8x8[i - 6][j] == r.ScalingList8x8[i - 6][j]);
                }
            }
        }
    }
    equal &= (l.second_chroma_qp_index_offset == r.second_chroma_qp_index_offset);
    return equal;
}


}
}
