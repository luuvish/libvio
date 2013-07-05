
/*!
 ************************************************************************
 *  \file
 *     parset.c
 *  \brief
 *     Parameter Sets
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Stephan Wenger          <stewe@cs.tu-berlin.de>
 *
 ***********************************************************************
 */

#include "global.h"
#include "slice.h"
#include "bitstream_elements.h"
#include "bitstream_nal.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "image.h"
#include "parset.h"
#include "memalloc.h"
#include "fmo.h"
#include "dpb.h"
#include "erc_api.h"
#include "macroblock.h"

#define MAX_QP          51

#define BASE_VIEW_IDX             0

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


// E.1.1 VUI parameter syntax
void vui_parameters(DataPartition *p, vui_t *vui)
{
    Bitstream *s = p->bitstream;

    vui->aspect_ratio_idc = 0;
    vui->aspect_ratio_info_present_flag = s->u(1, "VUI: aspect_ratio_info_present_flag");
    if (vui->aspect_ratio_info_present_flag) {
        vui->aspect_ratio_idc           = s->u(8, "VUI: aspect_ratio_idc");
        if (vui->aspect_ratio_idc == Extended_SAR) {
            vui->sar_width              = s->u(16, "VUI: sar_width");
            vui->sar_height             = s->u(16, "VUI: sar_height");
        }
    }

    vui->overscan_info_present_flag     = s->u(1, "VUI: overscan_info_present_flag");
    if (vui->overscan_info_present_flag)
        vui->overscan_appropriate_flag  = s->u(1, "VUI: overscan_appropriate_flag");

    vui->video_format             = 5; // Unspecified video format
    vui->video_full_range_flag    = 0;
    vui->colour_primaries         = 2; // Unspecified
    vui->transfer_characteristics = 2; // Unspecified
    vui->matrix_coefficients      = 2; // Unspecified
    vui->video_signal_type_present_flag = s->u(1, "VUI: video_signal_type_present_flag");
    if (vui->video_signal_type_present_flag) {
        vui->video_format                    = s->u(3, "VUI: video_format");
        vui->video_full_range_flag           = s->u(1, "VUI: video_full_range_flag");
        vui->colour_description_present_flag = s->u(1, "VUI: color_description_present_flag");
        if (vui->colour_description_present_flag) {
            vui->colour_primaries            = s->u(8, "VUI: colour_primaries");
            vui->transfer_characteristics    = s->u(8, "VUI: transfer_characteristics");
            vui->matrix_coefficients         = s->u(8, "VUI: matrix_coefficients");
        }
    }

    // if (!(BitDepthC == BitDepthY && chroma_format_idc == 3))
    //     assert(vui->matrix_coefficients != 0);
    // if (!(BitDepthC == BitDepthY || (BitDepthC == BitDepthY + 1 && chroma_format_idc == 3)))
    //     assert(vui->matrix_coefficients != 8);

    vui->chroma_sample_loc_type_top_field    = 0;
    vui->chroma_sample_loc_type_bottom_field = 0;
    vui->chroma_loc_info_present_flag = s->u(1, "VUI: chroma_loc_info_present_flag");
    if (vui->chroma_loc_info_present_flag) {
        vui->chroma_sample_loc_type_top_field    = s->ue("VUI: chroma_sample_loc_type_top_field");
        vui->chroma_sample_loc_type_bottom_field = s->ue("VUI: chroma_sample_loc_type_bottom_field");
    }

    // if (chroma_format_idc != 1)
    //     assert(vui->chroma_loc_info_present_flag == 0);

    assert(vui->chroma_sample_loc_type_top_field >= 0 && vui->chroma_sample_loc_type_top_field <= 5);
    assert(vui->chroma_sample_loc_type_bottom_field >= 0 && vui->chroma_sample_loc_type_bottom_field <= 5);

    vui->fixed_frame_rate_flag = 0;
    vui->timing_info_present_flag            = s->u(1, "VUI: timing_info_present_flag");
    if (vui->timing_info_present_flag) {
        vui->num_units_in_tick               = s->u(32, "VUI: num_units_in_tick");
        vui->time_scale                      = s->u(32, "VUI: time_scale");
        vui->fixed_frame_rate_flag           = s->u(1, "VUI: fixed_frame_rate_flag");
    }

    if (vui->timing_info_present_flag) {
        assert(vui->num_units_in_tick > 0);
        assert(vui->time_scale > 0);
    }

    vui->nal_hrd_parameters_present_flag   = s->u(1, "VUI: nal_hrd_parameters_present_flag");
    if (vui->nal_hrd_parameters_present_flag)
        hrd_parameters(p, &vui->nal_hrd_parameters);

    vui->vcl_hrd_parameters_present_flag   = s->u(1, "VUI: vcl_hrd_parameters_present_flag");
    if (vui->vcl_hrd_parameters_present_flag)
        hrd_parameters(p, &vui->vcl_hrd_parameters);

    vui->low_delay_hrd_flag = 1 - vui->fixed_frame_rate_flag;
    if (vui->nal_hrd_parameters_present_flag || vui->vcl_hrd_parameters_present_flag)
        vui->low_delay_hrd_flag            = s->u(1, "VUI: low_delay_hrd_flag");

    if (vui->fixed_frame_rate_flag)
        assert(vui->low_delay_hrd_flag == 0);

    vui->motion_vectors_over_pic_boundaries_flag = 1;
    vui->max_bytes_per_pic_denom                 = 2;
    vui->max_bits_per_mb_denom                   = 1;
    vui->log2_max_mv_length_horizontal           = 16;
    vui->log2_max_mv_length_vertical             = 16;
    //if (profile_idc == 44, 86, 100, 110, 122, 244 && constraint_set3_flag)
    //    vui->max_num_reorder_frames = 0;
    //else
    //    vui->max_num_reorder_frames = MaxDpbFrames;
    //if (profile_idc == 44, 86, 100, 110, 122, 244 && constraint_set3_flag)
    //    vui->max_dec_frame_buffering = 0;
    //else
    //    vui->max_dec_frame_buffering = MaxDpbFrames;
    vui->pic_struct_present_flag          = s->u(1, "VUI: pic_struct_present_flag");
    vui->bitstream_restriction_flag       = s->u(1, "VUI: bitstream_restriction_flag");
    if (vui->bitstream_restriction_flag) {
        vui->motion_vectors_over_pic_boundaries_flag = s->u(1, "VUI: motion_vectors_over_pic_boundaries_flag");
        vui->max_bytes_per_pic_denom                 = s->ue("VUI: max_bytes_per_pic_denom");
        vui->max_bits_per_mb_denom                   = s->ue("VUI: max_bits_per_mb_denom");
        vui->log2_max_mv_length_horizontal           = s->ue("VUI: log2_max_mv_length_horizontal");
        vui->log2_max_mv_length_vertical             = s->ue("VUI: log2_max_mv_length_vertical");
        vui->max_num_reorder_frames                  = s->ue("VUI: num_reorder_frames");
        vui->max_dec_frame_buffering                 = s->ue("VUI: max_dec_frame_buffering");
    }

    assert(vui->max_bytes_per_pic_denom >= 0 && vui->max_bytes_per_pic_denom <= 16);
    assert(vui->max_bits_per_mb_denom >= 0 && vui->max_bits_per_mb_denom <= 16);
    assert(vui->log2_max_mv_length_horizontal >= 0 && vui->log2_max_mv_length_horizontal <= 16);
    assert(vui->log2_max_mv_length_vertical >= 0 && vui->log2_max_mv_length_vertical <= 16);
    assert(vui->max_num_reorder_frames >= 0 &&
           vui->max_num_reorder_frames <= vui->max_dec_frame_buffering);
}

// E.1.2 HRD parameters syntax
void hrd_parameters(DataPartition *p, hrd_t *hrd)
{
    Bitstream *s = p->bitstream;
    int SchedSelIdx;

    hrd->cpb_cnt_minus1                                      = s->ue("VUI: cpb_cnt_minus1");
    hrd->bit_rate_scale                                      = s->u(4, "VUI: bit_rate_scale");
    hrd->cpb_size_scale                                      = s->u(4, "VUI: cpb_size_scale");

    assert(hrd->cpb_cnt_minus1 >= 0 && hrd->cpb_cnt_minus1 <= 31);

    // if (profile_idc == 66, 77, 88)
    //     BitRate[SchedSelIdx] = 1200 * MaxBR;
    //     CpbSize[SchedSelIdx] = 1000 * MaxCPB;
    // else
    //     BitRate[SchedSelIdx] = cpbBrVcl(Nal)Factor * MaxBR
    //     CpbSize[SchedSelIdx] = cpbBrVcl(Nal)Factor * MaxCPB

    for (SchedSelIdx = 0; SchedSelIdx <= hrd->cpb_cnt_minus1; SchedSelIdx++) {
        hrd->bit_rate_value_minus1[SchedSelIdx]             = s->ue("VUI: bit_rate_value_minus1");
        hrd->cpb_size_value_minus1[SchedSelIdx]             = s->ue("VUI: cpb_size_value_minus1");
        hrd->cbr_flag             [SchedSelIdx]             = s->u(1, "VUI: cbr_flag");
        // BitRate[SchedSelIdx] = (bit_rate_value_minus1[SchedSelIdx] + 1) * (1 << (6 + bit_rate_scale));
        // CpbSize[SchedSelIdx] = (cpb_size_value_minus1[SchedSelIdx] + 1) * (1 << (4 + cpb_size_scale));
    }

    hrd->initial_cpb_removal_delay_length_minus1 = 23;
    hrd->cpb_removal_delay_length_minus1         = 23;
    hrd->dpb_output_delay_length_minus1          = 23;
    hrd->time_offset_length                      = 24;
    hrd->initial_cpb_removal_delay_length_minus1 = s->u(5, "VUI: initial_cpb_removal_delay_length_minus1");
    hrd->cpb_removal_delay_length_minus1         = s->u(5, "VUI: cpb_removal_delay_length_minus1");
    hrd->dpb_output_delay_length_minus1          = s->u(5, "VUI: dpb_output_delay_length_minus1");
    hrd->time_offset_length                      = s->u(5, "VUI: time_offset_length");
}

// 7.3.2.1 Sequence parameter set data syntax
void seq_parameter_set_rbsp(DataPartition *p, sps_t *sps)
{
    unsigned i;
    Bitstream *s = p->bitstream;

    uint8_t reserved_zero_2bits;

    sps->profile_idc                           = s->u(8, "SPS: profile_idc");
    sps->constraint_set0_flag                  = s->u(1, "SPS: constrained_set0_flag");
    sps->constraint_set1_flag                  = s->u(1, "SPS: constrained_set1_flag");
    sps->constraint_set2_flag                  = s->u(1, "SPS: constrained_set2_flag");
    sps->constraint_set3_flag                  = s->u(1, "SPS: constrained_set3_flag");
    sps->constraint_set4_flag                  = s->u(1, "SPS: constrained_set4_flag");
    sps->constraint_set5_flag                  = s->u(1, "SPS: constrained_set5_flag");
    reserved_zero_2bits                        = s->u(2, "SPS: reserved_zero_2bits");

    if ((sps->profile_idc!=BASELINE       ) &&
        (sps->profile_idc!=MAIN           ) &&
        (sps->profile_idc!=EXTENDED       ) &&
        (sps->profile_idc!=FREXT_HP       ) &&
        (sps->profile_idc!=FREXT_Hi10P    ) &&
        (sps->profile_idc!=FREXT_Hi422    ) &&
        (sps->profile_idc!=FREXT_Hi444    ) &&
        (sps->profile_idc!=FREXT_CAVLC444 )
        && (sps->profile_idc!=MVC_HIGH)
        && (sps->profile_idc!=STEREO_HIGH)) {
        printf("Invalid Profile IDC (%d) encountered. \n", sps->profile_idc);
        return;
    }

    // assert(reserved_zero_2bits == 0);

    sps->level_idc                            = s->u(8, "SPS: level_idc");
    sps->seq_parameter_set_id                 = s->ue("SPS: seq_parameter_set_id");

    assert(sps->seq_parameter_set_id >= 0 && sps->seq_parameter_set_id <= 31);

    sps->chroma_format_idc                    = YUV420;
    sps->separate_colour_plane_flag           = 0;
    sps->bit_depth_luma_minus8                = 0;
    sps->bit_depth_chroma_minus8              = 0;
    sps->qpprime_y_zero_transform_bypass_flag = 0;
    sps->seq_scaling_matrix_present_flag      = 0;
    if (sps->profile_idc == FREXT_HP    || sps->profile_idc == FREXT_Hi10P ||
        sps->profile_idc == FREXT_Hi422 || sps->profile_idc == FREXT_Hi444 ||
        sps->profile_idc == FREXT_CAVLC444 || sps->profile_idc == MVC_HIGH ||
        sps->profile_idc == STEREO_HIGH) {

        sps->chroma_format_idc                      = s->ue("SPS: chroma_format_idc");
        if (sps->chroma_format_idc == YUV444)
            sps->separate_colour_plane_flag         = s->u(1, "SPS: separate_colour_plane_flag");
        sps->bit_depth_luma_minus8                  = s->ue("SPS: bit_depth_luma_minus8");
        sps->bit_depth_chroma_minus8                = s->ue("SPS: bit_depth_chroma_minus8");
        sps->qpprime_y_zero_transform_bypass_flag   = s->u(1, "SPS: lossless_qpprime_y_zero_flag");
        sps->seq_scaling_matrix_present_flag        = s->u(1, "SPS: seq_scaling_matrix_present_flag");
        if (sps->seq_scaling_matrix_present_flag) {
            for (i = 0; i < ((sps->chroma_format_idc != YUV444) ? 8 : 12); i++) {
                sps->seq_scaling_list_present_flag[i] = s->u(1, "SPS: seq_scaling_list_present_flag");
                if (sps->seq_scaling_list_present_flag[i]) {
                    if (i < 6)
                        scaling_list(sps->ScalingList4x4[i], 16,
                                     &sps->UseDefaultScalingMatrix4x4Flag[i], s);
                    else
                        scaling_list(sps->ScalingList8x8[i - 6], 64,
                                     &sps->UseDefaultScalingMatrix8x8Flag[i - 6], s);
                }
            }
        }
    }

    if (!sps->separate_colour_plane_flag)
        sps->ChromaArrayType = sps->chroma_format_idc;
    else
        sps->ChromaArrayType = 0;
    sps->SubWidthC = sps->chroma_format_idc < 3 ? 2 : 1;
    sps->SubHeightC = sps->chroma_format_idc < 2 ? 2 : 1;
    if (sps->chroma_format_idc == 0 || sps->separate_colour_plane_flag == 1) {
        sps->MbWidthC = 0;
        sps->MbHeightC = 0;
    } else {
        sps->MbWidthC = 16 / sps->SubWidthC;
        sps->MbHeightC = 16 / sps->SubHeightC;
    }

    assert(sps->bit_depth_luma_minus8 >= 0 && sps->bit_depth_luma_minus8 <= 6);
    assert(sps->bit_depth_chroma_minus8 >= 0 && sps->bit_depth_chroma_minus8 <= 6);

    sps->BitDepthY   = 8 + sps->bit_depth_luma_minus8;
    sps->QpBdOffsetY = 6 * sps->bit_depth_luma_minus8;
    sps->BitDepthC   = 8 + sps->bit_depth_chroma_minus8;
    sps->QpBdOffsetC = 6 * sps->bit_depth_chroma_minus8;
    sps->RawMbBits   = 256 * sps->BitDepthY + 2 * sps->MbWidthC * sps->MbHeightC * sps->BitDepthC;

    sps->log2_max_frame_num_minus4              = s->ue("SPS: log2_max_frame_num_minus4");

    assert(sps->log2_max_frame_num_minus4 >= 0 && sps->log2_max_frame_num_minus4 <= 12);
    sps->MaxFrameNum = 1 << (sps->log2_max_frame_num_minus4 + 4);

    sps->pic_order_cnt_type                     = s->ue("SPS: pic_order_cnt_type");

    assert(sps->pic_order_cnt_type >= 0 && sps->pic_order_cnt_type <= 2);

    if (sps->pic_order_cnt_type == 0)
        sps->log2_max_pic_order_cnt_lsb_minus4     = s->ue("SPS: log2_max_pic_order_cnt_lsb_minus4");
    else if (sps->pic_order_cnt_type == 1) {
        sps->delta_pic_order_always_zero_flag      = s->u(1, "SPS: delta_pic_order_always_zero_flag");
        sps->offset_for_non_ref_pic                = s->se("SPS: offset_for_non_ref_pic");
        sps->offset_for_top_to_bottom_field        = s->se("SPS: offset_for_top_to_bottom_field");
        sps->num_ref_frames_in_pic_order_cnt_cycle = s->ue("SPS: num_ref_frames_in_pic_order_cnt_cycle");
        for (i = 0; i < sps->num_ref_frames_in_pic_order_cnt_cycle; i++)
            sps->offset_for_ref_frame[i]           = s->se("SPS: offset_for_ref_frame[i]");
    }

    assert(sps->log2_max_pic_order_cnt_lsb_minus4 >= 0 && sps->log2_max_pic_order_cnt_lsb_minus4 <= 12);
    sps->MaxPicOrderCntLsb = 1 << (sps->log2_max_pic_order_cnt_lsb_minus4 + 4);
    assert(sps->offset_for_non_ref_pic >= -(1 << 31) + 1 && sps->offset_for_non_ref_pic <= (1 << 31) - 1);
    assert(sps->offset_for_top_to_bottom_field >= -(1 << 31) + 1 && sps->offset_for_top_to_bottom_field <= (1 << 31) - 1);
    assert(sps->num_ref_frames_in_pic_order_cnt_cycle >= 0 && sps->num_ref_frames_in_pic_order_cnt_cycle <= 255);

    // ExpectedDeltaPerPicOrderCntCycle = 0;
    // for (i = 0; i < num_ref_frames_in_pic_order_cnt_cycle; i++)
    //     ExpectedDeltaPerPicOrderCntCycle += offset_for_ref_frame[i];

    sps->max_num_ref_frames                        = s->ue("SPS: num_ref_frames");

    // assert(sps->max_num_ref_frames >= 0 && sps->max_num_ref_frames <= MaxDpbFrames);

    sps->gaps_in_frame_num_value_allowed_flag  = s->u(1, "SPS: gaps_in_frame_num_value_allowed_flag");

    sps->pic_width_in_mbs_minus1               = s->ue("SPS: pic_width_in_mbs_minus1");
    sps->pic_height_in_map_units_minus1        = s->ue("SPS: pic_height_in_map_units_minus1");

    sps->PicWidthInMbs = sps->pic_width_in_mbs_minus1 + 1;
    sps->PicWidthInSamplesL = sps->PicWidthInMbs * 16;
    sps->PicWidthInSamplesC = sps->PicWidthInMbs * sps->MbWidthC;
    sps->PicHeightInMapUnits = sps->pic_height_in_map_units_minus1 + 1;
    sps->PicSizeInMapUnits = sps->PicWidthInMbs * sps->PicHeightInMapUnits;

    sps->frame_mbs_only_flag                   = s->u(1, "SPS: frame_mbs_only_flag");

    sps->FrameHeightInMbs = (2 - sps->frame_mbs_only_flag) * sps->PicHeightInMapUnits;

    sps->mb_adaptive_frame_field_flag = 0;
    if (!sps->frame_mbs_only_flag)
        sps->mb_adaptive_frame_field_flag      = s->u(1, "SPS: mb_adaptive_frame_field_flag");

    sps->direct_8x8_inference_flag             = s->u(1, "SPS: direct_8x8_inference_flag");

    if (!sps->frame_mbs_only_flag)
        assert(sps->direct_8x8_inference_flag == 1);

    sps->frame_cropping_flag                   = s->u(1, "SPS: frame_cropping_flag");
    if (sps->frame_cropping_flag) {
        sps->frame_crop_left_offset      = s->ue("SPS: frame_crop_left_offset");
        sps->frame_crop_right_offset     = s->ue("SPS: frame_crop_right_offset");
        sps->frame_crop_top_offset       = s->ue("SPS: frame_crop_top_offset");
        sps->frame_crop_bottom_offset    = s->ue("SPS: frame_crop_bottom_offset");
    }

    if (sps->frame_cropping_flag) {
        if (sps->ChromaArrayType == 0) {
            sps->CropUnitX = 1;
            sps->CropUnitY = 2 - sps->frame_mbs_only_flag;
        } else {
            sps->CropUnitX = sps->SubWidthC;
            sps->CropUnitY = sps->SubHeightC * (2 - sps->frame_mbs_only_flag);
        }
        assert(sps->frame_crop_left_offset >= 0 &&
               sps->frame_crop_left_offset <= (sps->PicWidthInSamplesL / sps->CropUnitX) -
                                              (sps->frame_crop_right_offset + 1));
        assert(sps->frame_crop_top_offset >= 0 &&
               sps->frame_crop_top_offset <= (16 * sps->FrameHeightInMbs / sps->CropUnitY) -
                                             (sps->frame_crop_bottom_offset + 1));
    }

    sps->vui_parameters_present_flag           = s->u(1, "SPS: vui_parameters_present_flag");
    sps->vui_parameters.matrix_coefficients = 2;
    if (sps->vui_parameters_present_flag)
        vui_parameters(p, &sps->vui_parameters);

    sps->Valid = TRUE;
}

static const byte ZZ_SCAN[16] = {
     0,  1,  4,  8,  5,  2,  3,  6,
     9, 12, 13, 10,  7, 11, 14, 15
};

static const byte ZZ_SCAN8[64] = {
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
void scaling_list(int *scalingList, int sizeOfScalingList, bool *useDefaultScalingMatrixFlag, Bitstream *s)
{
    int lastScale = 8;
    int nextScale = 8;
    int j;

    for (j = 0; j < sizeOfScalingList; j++) {
        int scanj = (sizeOfScalingList == 16) ? ZZ_SCAN[j] : ZZ_SCAN8[j];
        if (nextScale != 0) {
            int delta_scale = s->se("   : delta_sl");
            nextScale = (lastScale + delta_scale + 256) % 256;
            *useDefaultScalingMatrixFlag = (scanj == 0 && nextScale == 0);
        }
        scalingList[scanj] = (nextScale == 0) ? lastScale : nextScale;
        lastScale = scalingList[scanj];
    }
}

// 7.3.2.2 Picture parameter set RBSP syntax
void pic_parameter_set_rbsp(VideoParameters *p_Vid, DataPartition *p, pps_t *pps)
{
    int iGroup;
    int i;
    int chroma_format_idc;
    Bitstream *s = p->bitstream;

    pps->pic_parameter_set_id                         = s->ue("PPS: pic_parameter_set_id");
    pps->seq_parameter_set_id                         = s->ue("PPS: seq_parameter_set_id");
    pps->entropy_coding_mode_flag                     = s->u(1, "PPS: entropy_coding_mode_flag");
    pps->bottom_field_pic_order_in_frame_present_flag = s->u(1, "PPS: bottom_field_pic_order_in_frame_present_flag");

    assert(pps->pic_parameter_set_id >= 0 && pps->pic_parameter_set_id <= 255);
    assert(pps->seq_parameter_set_id >= 0 && pps->seq_parameter_set_id <= 31);

    chroma_format_idc = p_Vid->SeqParSet[pps->seq_parameter_set_id].chroma_format_idc;

    pps->num_slice_groups_minus1                = s->ue("PPS: num_slice_groups_minus1");
    if (pps->num_slice_groups_minus1 > 0) {
        pps->slice_group_map_type               = s->ue("PPS: slice_group_map_type");
        if (pps->slice_group_map_type == 0) {
            for (iGroup = 0; iGroup <= pps->num_slice_groups_minus1; iGroup++)
                pps->run_length_minus1[iGroup]                 = s->ue("PPS: run_length_minus1 [i]");
        } else if (pps->slice_group_map_type == 2) {
            for (iGroup = 0; iGroup < pps->num_slice_groups_minus1; iGroup++) {
                pps->top_left    [iGroup]                      = s->ue("PPS: top_left [i]");
                pps->bottom_right[iGroup]                      = s->ue("PPS: bottom_right [i]");
            }
        } else if (pps->slice_group_map_type == 3 ||
                   pps->slice_group_map_type == 4 ||
                   pps->slice_group_map_type == 5) {
            pps->slice_group_change_direction_flag     = s->u(1, "PPS: slice_group_change_direction_flag");
            pps->slice_group_change_rate_minus1        = s->ue("PPS: slice_group_change_rate_minus1");
        } else if (pps->slice_group_map_type == 6) {
            const int bitsSliceGroupId[8] = { 1, 1, 2, 2, 3, 3, 3, 3 };
            pps->pic_size_in_map_units_minus1      = s->ue("PPS: pic_size_in_map_units_minus1");
            if ((pps->slice_group_id = (byte *)calloc (pps->pic_size_in_map_units_minus1+1, 1)) == NULL)
                no_mem_exit ("InterpretPPS: slice_group_id");
            for (i = 0; i <= pps->pic_size_in_map_units_minus1; i++)
                pps->slice_group_id[i] = (byte) s->u(bitsSliceGroupId[pps->num_slice_groups_minus1], "slice_group_id[i]");
        }
    }

    assert(pps->num_slice_groups_minus1 >= 0 && pps->num_slice_groups_minus1 <= 7);
    assert(pps->slice_group_map_type >= 0 && pps->slice_group_map_type <= 6);
    if (pps->num_slice_groups_minus1 != 1)
        assert(pps->slice_group_map_type != 3 && pps->slice_group_map_type != 4 && pps->slice_group_map_type != 5);
    // assert(pps->run_length_minus1[i] <= PicSizeInMapUnits - 1);
    // assert(pps->top_left[i] <= pps->bottom_right[i] && pps->bottom_right[i] < PicSizeInMapUnits);
    // assert(pps->top_left[i] % PicWidthInMbs <= pps->bottom_right[i] % PicWidthInMbs);

    // assert(pps->slice_group_change_rate_minus1 >= 0 &&
    //        pps->slice_group_change_rate_minus1 <= PicSizeInMapUnits - 1);
    pps->SliceGroupChangeRate = pps->slice_group_change_rate_minus1 + 1;
    //assert(pps->pic_size_in_map_units_minus1 == PicSizeInMapUnits - 1);
    //assert(pps->slice_group_id[i] >= 0 && pps->slice_group_id[i] <= pps->num_slice_groups_minus1);

    pps->num_ref_idx_l0_default_active_minus1  = s->ue("PPS: num_ref_idx_l0_default_active_minus1");
    pps->num_ref_idx_l1_default_active_minus1  = s->ue("PPS: num_ref_idx_l1_default_active_minus1");

    assert(pps->num_ref_idx_l0_default_active_minus1 >= 0 && pps->num_ref_idx_l0_default_active_minus1 <= 31);
    assert(pps->num_ref_idx_l1_default_active_minus1 >= 0 && pps->num_ref_idx_l1_default_active_minus1 <= 31);

    pps->weighted_pred_flag                    = s->u(1, "PPS: weighted_pred_flag");
    pps->weighted_bipred_idc                   = s->u(2, "PPS: weighted_bipred_idc");

    assert(pps->weighted_bipred_idc >= 0 && pps->weighted_bipred_idc <= 2);

    pps->pic_init_qp_minus26                   = s->se("PPS: pic_init_qp_minus26");
    pps->pic_init_qs_minus26                   = s->se("PPS: pic_init_qs_minus26");
    pps->chroma_qp_index_offset                = s->se("PPS: chroma_qp_index_offset");

    //assert(pps->pic_init_qp_minus26 >= -(26 + QpBdOffsetY) && pps->pic_init_qp_minus26 <= 25);
    assert(pps->pic_init_qs_minus26 >= -26 && pps->pic_init_qs_minus26 <= 25);
    assert(pps->chroma_qp_index_offset >= -12 && pps->chroma_qp_index_offset <= 12);

    pps->deblocking_filter_control_present_flag = s->u(1, "PPS: deblocking_filter_control_present_flag");
    pps->constrained_intra_pred_flag            = s->u(1, "PPS: constrained_intra_pred_flag");
    pps->redundant_pic_cnt_present_flag         = s->u(1, "PPS: redundant_pic_cnt_present_flag");

    pps->transform_8x8_mode_flag         = 0;
    pps->pic_scaling_matrix_present_flag = 0;
    pps->second_chroma_qp_index_offset   = pps->chroma_qp_index_offset;
    if (s->more_rbsp_data()) {
        pps->transform_8x8_mode_flag           = s->u(1, "PPS: transform_8x8_mode_flag");
        pps->pic_scaling_matrix_present_flag   = s->u(1, "PPS: pic_scaling_matrix_present_flag");
        if (pps->pic_scaling_matrix_present_flag) {
            for (i = 0; i < 6 + ((chroma_format_idc != YUV444) ? 2 : 6) * pps->transform_8x8_mode_flag; i++) {
                pps->pic_scaling_list_present_flag[i]= s->u(1, "PPS: pic_scaling_list_present_flag");
                if (pps->pic_scaling_list_present_flag[i]) {
                    if (i < 6)
                        scaling_list(pps->ScalingList4x4[i], 16,
                                     &pps->UseDefaultScalingMatrix4x4Flag[i], s);
                    else
                        scaling_list(pps->ScalingList8x8[i - 6], 64,
                                     &pps->UseDefaultScalingMatrix8x8Flag[i - 6], s);
                }
            }
        }
        pps->second_chroma_qp_index_offset      = s->se("PPS: second_chroma_qp_index_offset");
    }
    assert(pps->second_chroma_qp_index_offset >= -12 && pps->second_chroma_qp_index_offset <= 12);

    pps->Valid = TRUE;
}






void PPSConsistencyCheck (pps_t *pps);
void SPSConsistencyCheck (sps_t *sps);

void MakeSPSavailable (VideoParameters *p_Vid, int id, sps_t *sps);

#if (MVC_EXTENSION_ENABLE)
void SubsetSPSConsistencyCheck (sub_sps_t *subset_sps);
#endif



// fill subset_sps with content of p
#if (MVC_EXTENSION_ENABLE)

static void get_max_dec_frame_buf_size(sps_t *sps)
{
  int pic_size = (sps->pic_width_in_mbs_minus1 + 1) * (sps->pic_height_in_map_units_minus1 + 1) * (sps->frame_mbs_only_flag?1:2) * 384;

  int size = 0;

    switch (sps->level_idc) {
    case 9:
        size = 152064;
        break;
    case 10:
        size = 152064;
        break;
    case 11:
        if (!is_FREXT_profile(sps->profile_idc) && (sps->constraint_set3_flag == 1))
            size = 152064;
        else
            size = 345600;
        break;
    case 12:
        size = 912384;
        break;
    case 13:
        size = 912384;
        break;
    case 20:
        size = 912384;
        break;
    case 21:
        size = 1824768;
        break;
    case 22:
        size = 3110400;
        break;
    case 30:
        size = 3110400;
        break;
    case 31:
        size = 6912000;
        break;
    case 32:
        size = 7864320;
        break;
    case 40:
        size = 12582912;
        break;
    case 41:
        size = 12582912;
        break;
    case 42:
        size = 13369344;
        break;
    case 50:
        size = 42393600;
        break;
    case 51:
        size = 70778880;
        break;
    case 52:
        size = 70778880;
        break;
    default:
        error("undefined level", 500);
        break;
    }

    size /= pic_size;
    size = imin( size, 16);
    sps->max_dec_frame_buffering = size;
}

static void seq_parameter_set_mvc_extension(sub_sps_t *subset_sps, Bitstream *s)
{
  int i, j, num_views;

  subset_sps->num_views_minus1 = s->ue("num_views_minus1");
  num_views = 1+subset_sps->num_views_minus1;
  if( num_views >0)
  {
    if ((subset_sps->view_id = (int*) calloc(num_views, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->view_id");
    if ((subset_sps->num_anchor_refs_l0 = (int*) calloc(num_views, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->num_anchor_refs_l0");
    if ((subset_sps->num_anchor_refs_l1 = (int*) calloc(num_views, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->num_anchor_refs_l1");
    if ((subset_sps->anchor_ref_l0 = (int**) calloc(num_views, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->anchor_ref_l0");
    if ((subset_sps->anchor_ref_l1 = (int**) calloc(num_views, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->anchor_ref_l1");
    if ((subset_sps->num_non_anchor_refs_l0 = (int*) calloc(num_views, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->num_non_anchor_refs_l0");
    if ((subset_sps->num_non_anchor_refs_l1 = (int*) calloc(num_views, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->num_non_anchor_refs_l1");
    if ((subset_sps->non_anchor_ref_l0 = (int**) calloc(num_views, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->non_anchor_ref_l0");
    if ((subset_sps->non_anchor_ref_l1 = (int**) calloc(num_views, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->non_anchor_ref_l1");
  }
  for(i=0; i<num_views; i++)
  {
    subset_sps->view_id[i] = s->ue("view_id");
  }
  for(i=1; i<num_views; i++)
  {
    subset_sps->num_anchor_refs_l0[i] = s->ue("num_anchor_refs_l0");
    if(subset_sps->num_anchor_refs_l0[i]>0)
    {
      if ((subset_sps->anchor_ref_l0[i] = (int*) calloc(subset_sps->num_anchor_refs_l0[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->anchor_ref_l0[i]");
      for(j=0; j<subset_sps->num_anchor_refs_l0[i]; j++)
        subset_sps->anchor_ref_l0[i][j] = s->ue("anchor_ref_l0");
    }

    subset_sps->num_anchor_refs_l1[i] = s->ue("num_anchor_refs_l1");
    if(subset_sps->num_anchor_refs_l1[i]>0)
    {
      if ((subset_sps->anchor_ref_l1[i] = (int*) calloc(subset_sps->num_anchor_refs_l1[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->anchor_ref_l1[i]");
      for(j=0; j<subset_sps->num_anchor_refs_l1[i]; j++)
        subset_sps->anchor_ref_l1[i][j] = s->ue("anchor_ref_l1");
    }
  }
  for(i=1; i<num_views; i++)
  {
    subset_sps->num_non_anchor_refs_l0[i] = s->ue("num_non_anchor_refs_l0");
    if(subset_sps->num_non_anchor_refs_l0[i]>0)
    {
      if ((subset_sps->non_anchor_ref_l0[i] = (int*) calloc(subset_sps->num_non_anchor_refs_l0[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->non_anchor_ref_l0[i]");
      for(j=0; j<subset_sps->num_non_anchor_refs_l0[i]; j++)
        subset_sps->non_anchor_ref_l0[i][j] = s->ue("non_anchor_ref_l0");
    }
    subset_sps->num_non_anchor_refs_l1[i] = s->ue("num_non_anchor_refs_l1");
    if(subset_sps->num_non_anchor_refs_l1[i]>0)
    {
      if ((subset_sps->non_anchor_ref_l1[i] = (int*) calloc(subset_sps->num_non_anchor_refs_l1[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->non_anchor_ref_l1[i]");
      for(j=0; j<subset_sps->num_non_anchor_refs_l1[i]; j++)
        subset_sps->non_anchor_ref_l1[i][j] = s->ue("non_anchor_ref_l1");
    }
  }
  subset_sps->num_level_values_signalled_minus1 = s->ue("num_level_values_signalled_minus1");
  if(subset_sps->num_level_values_signalled_minus1 >=0)
  {
    i = 1+ subset_sps->num_level_values_signalled_minus1;
    if ((subset_sps->level_idc = (int*) calloc(i, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->level_idc");
    if ((subset_sps->num_applicable_ops_minus1 = (int*) calloc(i, sizeof(int))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->num_applicable_ops_minus1");
    if ((subset_sps->applicable_op_temporal_id = (int**) calloc(i, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_temporal_id");
    if ((subset_sps->applicable_op_num_target_views_minus1 = (int**) calloc(i, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_num_target_views_minus1");
    if ((subset_sps->applicable_op_target_view_id = (int***) calloc(i, sizeof(int**))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_target_view_id");
    if ((subset_sps->applicable_op_num_views_minus1 = (int**) calloc(i, sizeof(int*))) == NULL)
      no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_num_views_minus1");
  }
  for(i=0; i<=subset_sps->num_level_values_signalled_minus1; i++)
  {
    subset_sps->level_idc[i] = s->u(8, "level_idc");
    subset_sps->num_applicable_ops_minus1[i] = s->ue("num_applicable_ops_minus1");
    if(subset_sps->num_applicable_ops_minus1[i]>=0)
    {
      if ((subset_sps->applicable_op_temporal_id[i] = (int*) calloc(1+subset_sps->num_applicable_ops_minus1[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_temporal_id[i]");
      if ((subset_sps->applicable_op_num_target_views_minus1[i] = (int*) calloc(1+subset_sps->num_applicable_ops_minus1[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_num_target_views_minus1[i]");
      if ((subset_sps->applicable_op_target_view_id[i] = (int**) calloc(1+subset_sps->num_applicable_ops_minus1[i], sizeof(int *))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_target_view_id[i]");
      if ((subset_sps->applicable_op_num_views_minus1[i] = (int*) calloc(1+subset_sps->num_applicable_ops_minus1[i], sizeof(int))) == NULL)
        no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_num_views_minus1[i]");

      for(j=0; j<=subset_sps->num_applicable_ops_minus1[i]; j++)
      {
        int k;
        subset_sps->applicable_op_temporal_id[i][j] = s->u(3, "applicable_op_temporal_id");
        subset_sps->applicable_op_num_target_views_minus1[i][j] = s->ue("applicable_op_num_target_views_minus1");
        if(subset_sps->applicable_op_num_target_views_minus1[i][j]>=0)
        {
          if ((subset_sps->applicable_op_target_view_id[i][j] = (int*) calloc(1+subset_sps->applicable_op_num_target_views_minus1[i][j], sizeof(int))) == NULL)
            no_mem_exit("init_subset_seq_parameter_set: subset_sps->applicable_op_target_view_id[i][j]");
          for(k = 0; k <= subset_sps->applicable_op_num_target_views_minus1[i][j]; k++)
            subset_sps->applicable_op_target_view_id[i][j][k] = s->ue("applicable_op_target_view_id");
        }
        subset_sps->applicable_op_num_views_minus1[i][j] = s->ue("applicable_op_num_views_minus1");
      }
    }
  }
}

static int MemAlloc1D(void** ppBuf, int iEleSize, int iNum)
{
  if(iEleSize*iNum <=0)
    return 1;

  *ppBuf = calloc(iNum, iEleSize);
  return (*ppBuf == NULL);
}

static void mvc_hrd_parameters(mvc_vui_t *pMVCVUI, Bitstream *s)
{
  int i;

  pMVCVUI->cpb_cnt_minus1 = (char) s->ue("cpb_cnt_minus1");
  assert(pMVCVUI->cpb_cnt_minus1<=31);
  pMVCVUI->bit_rate_scale = (char) s->u(4, "bit_rate_scale");
  pMVCVUI->cpb_size_scale = (char) s->u(4, "cpb_size_scale");
  for(i=0; i<=pMVCVUI->cpb_cnt_minus1; i++)
  {
    pMVCVUI->bit_rate_value_minus1[i] = s->ue("bit_rate_value_minus1");
    pMVCVUI->cpb_size_value_minus1[i] = s->ue("cpb_size_value_minus1");
    pMVCVUI->cbr_flag[i]              = s->u(1, "cbr_flag");
  }
  pMVCVUI->initial_cpb_removal_delay_length_minus1 = (char) s->u(5, "initial_cpb_removal_delay_length_minus1");
  pMVCVUI->cpb_removal_delay_length_minus1         = (char) s->u(5, "cpb_removal_delay_length_minus1");
  pMVCVUI->dpb_output_delay_length_minus1          = (char) s->u(5, "dpb_output_delay_length_minus1");
  pMVCVUI->time_offset_length                      = (char) s->u(5, "time_offset_length");

}

static void mvc_vui_parameters_extension(mvc_vui_t *pMVCVUI, Bitstream *s)
{
  int i, j, iNumOps;

  pMVCVUI->num_ops_minus1 = s->ue("vui_mvc_num_ops_minus1");
  iNumOps = 1+ pMVCVUI->num_ops_minus1;
  if(iNumOps > 0)
  {
    MemAlloc1D((void **)&(pMVCVUI->temporal_id), sizeof(pMVCVUI->temporal_id[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->num_target_output_views_minus1), sizeof(pMVCVUI->num_target_output_views_minus1[0]), iNumOps);
    if ((pMVCVUI->view_id = (int**) calloc(iNumOps, sizeof(int*))) == NULL)
      no_mem_exit("mvc_vui_parameters_extension: pMVCVUI->view_id");
    MemAlloc1D((void **)&(pMVCVUI->timing_info_present_flag), sizeof(pMVCVUI->timing_info_present_flag[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->num_units_in_tick), sizeof(pMVCVUI->num_units_in_tick[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->time_scale), sizeof(pMVCVUI->time_scale[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->fixed_frame_rate_flag), sizeof(pMVCVUI->fixed_frame_rate_flag[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->nal_hrd_parameters_present_flag), sizeof(pMVCVUI->nal_hrd_parameters_present_flag[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->vcl_hrd_parameters_present_flag), sizeof(pMVCVUI->vcl_hrd_parameters_present_flag[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->low_delay_hrd_flag), sizeof(pMVCVUI->low_delay_hrd_flag[0]), iNumOps);
    MemAlloc1D((void **)&(pMVCVUI->pic_struct_present_flag), sizeof(pMVCVUI->pic_struct_present_flag[0]), iNumOps);

    for(i=0; i<iNumOps; i++)
    {
      pMVCVUI->temporal_id[i] = (char) s->u(3, "vui_mvc_temporal_id");
      pMVCVUI->num_target_output_views_minus1[i] = s->ue("vui_mvc_num_target_output_views_minus1");
      if(pMVCVUI->num_target_output_views_minus1[i] >= 0)
        MemAlloc1D((void **)&(pMVCVUI->view_id[i]), sizeof(pMVCVUI->view_id[0][0]), pMVCVUI->num_target_output_views_minus1[i]+1);
      for(j=0; j<=pMVCVUI->num_target_output_views_minus1[i]; j++)
        pMVCVUI->view_id[i][j] = s->ue("vui_mvc_view_id");
      pMVCVUI->timing_info_present_flag[i] = (char) s->u(1, "vui_mvc_timing_info_present_flag");
      if(pMVCVUI->timing_info_present_flag[i])
      {
        pMVCVUI->num_units_in_tick[i]     = s->u(32, "vui_mvc_num_units_in_tick");
        pMVCVUI->time_scale[i]            = s->u(32, "vui_mvc_time_scale");
        pMVCVUI->fixed_frame_rate_flag[i] = s->u(1, "vui_mvc_fixed_frame_rate_flag");
      }
      pMVCVUI->nal_hrd_parameters_present_flag[i] = (char) s->u(1, "vui_mvc_nal_hrd_parameters_present_flag");
      if(pMVCVUI->nal_hrd_parameters_present_flag[i])
        mvc_hrd_parameters(pMVCVUI, s);
      pMVCVUI->vcl_hrd_parameters_present_flag[i] = (char) s->u(1, "vcl_hrd_parameters_present_flag");
      if(pMVCVUI->vcl_hrd_parameters_present_flag[i])
        mvc_hrd_parameters(pMVCVUI, s);
      if(pMVCVUI->nal_hrd_parameters_present_flag[i]||pMVCVUI->vcl_hrd_parameters_present_flag[i])
        pMVCVUI->low_delay_hrd_flag[i]    = (char) s->u(1, "vui_mvc_low_delay_hrd_flag");
      pMVCVUI->pic_struct_present_flag[i] = (char) s->u(1, "vui_mvc_pic_struct_present_flag");
    }
  }
}

// 7.3.2.1.3 Subset sequence parameter set RBSP syntax
static int subset_seq_parameter_set_rbsp(VideoParameters *p_Vid, DataPartition *p, int *curr_seq_set_id)
{
    sub_sps_t *subset_sps;
    bool additional_extension2_flag;
    bool additional_extension2_data_flag;
    Bitstream *s = p->bitstream;
    sps_t *sps = AllocSPS();

  assert (p != NULL);
  assert (p->bitstream != NULL);
  assert (p->bitstream->streamBuffer != 0);

  seq_parameter_set_rbsp(p, sps);
  get_max_dec_frame_buf_size(sps);

  *curr_seq_set_id = sps->seq_parameter_set_id;
  subset_sps = p_Vid->SubsetSeqParSet + sps->seq_parameter_set_id;
  if(subset_sps->Valid || subset_sps->num_views_minus1>=0)
  {
    if(memcmp(&subset_sps->sps, sps, sizeof (sps_t)-sizeof(int)))
      assert(0);
    reset_subset_sps(subset_sps);
  }
  memcpy (&subset_sps->sps, sps, sizeof (sps_t));

  assert (subset_sps != NULL);
  subset_sps->Valid = FALSE;

    if (subset_sps->sps.profile_idc == 83 || subset_sps->sps.profile_idc == 86) {

    } else if (subset_sps->sps.profile_idc == MVC_HIGH ||
               subset_sps->sps.profile_idc == STEREO_HIGH) {
        subset_sps->bit_equal_to_one = s->u(1, "bit_equal_to_one");

        if (subset_sps->bit_equal_to_one != 1) {
            printf("\nbit_equal_to_one is not equal to 1!\n");
            return 0;
        }

        seq_parameter_set_mvc_extension(subset_sps, s);

        subset_sps->mvc_vui_parameters_present_flag = s->u(1, "mvc_vui_parameters_present_flag");
        if (subset_sps->mvc_vui_parameters_present_flag)
            mvc_vui_parameters_extension(&subset_sps->MVCVUIParams, s);
    }

    additional_extension2_flag = s->u(1, "additional_extension2_flag");
    if (additional_extension2_flag) {
        while (s->more_rbsp_data())
            additional_extension2_data_flag = s->u(1, "additional_extension2_flag");
    }

    if (subset_sps->sps.Valid)
        subset_sps->Valid = TRUE;

    FreeSPS (sps);
    return 0;
}
#endif


#if (MVC_EXTENSION_ENABLE)
void nal_unit_header_mvc_extension(NALUnitHeaderMVCExt_t *NaluHeaderMVCExt, Bitstream *s)
{  
    //to be implemented;  
    NaluHeaderMVCExt->non_idr_flag     = s->u(1, "non_idr_flag");
    NaluHeaderMVCExt->priority_id      = s->u(6, "priority_id");
    NaluHeaderMVCExt->view_id          = s->u(10, "view_id");
    NaluHeaderMVCExt->temporal_id      = s->u(3, "temporal_id");
    NaluHeaderMVCExt->anchor_pic_flag  = s->u(1, "anchor_pic_flag");
    NaluHeaderMVCExt->inter_view_flag  = s->u(1, "inter_view_flag");
    NaluHeaderMVCExt->reserved_one_bit = s->u(1, "reserved_one_bit");
    if (NaluHeaderMVCExt->reserved_one_bit != 1)
        printf("Nalu Header MVC Extension: reserved_one_bit is not 1!\n");
}

void nal_unit_header_svc_extension(void)
{
    //to be implemented for Annex G;
}

void prefix_nal_unit_svc(void)
{
    //to be implemented for Annex G;
}
#endif


static int sps_is_equal(sps_t *sps1, sps_t *sps2)
{
  unsigned i;
  int equal = 1;

  if ((!sps1->Valid) || (!sps2->Valid))
    return 0;

  equal &= (sps1->profile_idc == sps2->profile_idc);
  equal &= (sps1->constraint_set0_flag == sps2->constraint_set0_flag);
  equal &= (sps1->constraint_set1_flag == sps2->constraint_set1_flag);
  equal &= (sps1->constraint_set2_flag == sps2->constraint_set2_flag);
  equal &= (sps1->level_idc == sps2->level_idc);
  equal &= (sps1->seq_parameter_set_id == sps2->seq_parameter_set_id);
  equal &= (sps1->log2_max_frame_num_minus4 == sps2->log2_max_frame_num_minus4);
  equal &= (sps1->pic_order_cnt_type == sps2->pic_order_cnt_type);

  if (!equal) return equal;

  if( sps1->pic_order_cnt_type == 0 )
  {
    equal &= (sps1->log2_max_pic_order_cnt_lsb_minus4 == sps2->log2_max_pic_order_cnt_lsb_minus4);
  }

  else if( sps1->pic_order_cnt_type == 1 )
  {
    equal &= (sps1->delta_pic_order_always_zero_flag == sps2->delta_pic_order_always_zero_flag);
    equal &= (sps1->offset_for_non_ref_pic == sps2->offset_for_non_ref_pic);
    equal &= (sps1->offset_for_top_to_bottom_field == sps2->offset_for_top_to_bottom_field);
    equal &= (sps1->num_ref_frames_in_pic_order_cnt_cycle == sps2->num_ref_frames_in_pic_order_cnt_cycle);
    if (!equal) return equal;

    for ( i = 0 ; i< sps1->num_ref_frames_in_pic_order_cnt_cycle ;i ++)
      equal &= (sps1->offset_for_ref_frame[i] == sps2->offset_for_ref_frame[i]);
  }

  equal &= (sps1->max_num_ref_frames == sps2->max_num_ref_frames);
  equal &= (sps1->gaps_in_frame_num_value_allowed_flag == sps2->gaps_in_frame_num_value_allowed_flag);
  equal &= (sps1->pic_width_in_mbs_minus1 == sps2->pic_width_in_mbs_minus1);
  equal &= (sps1->pic_height_in_map_units_minus1 == sps2->pic_height_in_map_units_minus1);
  equal &= (sps1->frame_mbs_only_flag == sps2->frame_mbs_only_flag);

  if (!equal) return equal;
  if( !sps1->frame_mbs_only_flag )
    equal &= (sps1->mb_adaptive_frame_field_flag == sps2->mb_adaptive_frame_field_flag);

  equal &= (sps1->direct_8x8_inference_flag == sps2->direct_8x8_inference_flag);
  equal &= (sps1->frame_cropping_flag == sps2->frame_cropping_flag);
  if (!equal) return equal;
  if (sps1->frame_cropping_flag)
  {
    equal &= (sps1->frame_crop_left_offset == sps2->frame_crop_left_offset);
    equal &= (sps1->frame_crop_right_offset == sps2->frame_crop_right_offset);
    equal &= (sps1->frame_crop_top_offset == sps2->frame_crop_top_offset);
    equal &= (sps1->frame_crop_bottom_offset == sps2->frame_crop_bottom_offset);
  }
  equal &= (sps1->vui_parameters_present_flag == sps2->vui_parameters_present_flag);

  return equal;
}

static int pps_is_equal(pps_t *pps1, pps_t *pps2)
{
  unsigned i, j;
  int equal = 1;

  if ((!pps1->Valid) || (!pps2->Valid))
    return 0;

  equal &= (pps1->pic_parameter_set_id == pps2->pic_parameter_set_id);
  equal &= (pps1->seq_parameter_set_id == pps2->seq_parameter_set_id);
  equal &= (pps1->entropy_coding_mode_flag == pps2->entropy_coding_mode_flag);
  equal &= (pps1->bottom_field_pic_order_in_frame_present_flag == pps2->bottom_field_pic_order_in_frame_present_flag);
  equal &= (pps1->num_slice_groups_minus1 == pps2->num_slice_groups_minus1);

  if (!equal) return equal;

  if (pps1->num_slice_groups_minus1>0)
  {
      equal &= (pps1->slice_group_map_type == pps2->slice_group_map_type);
      if (!equal) return equal;
      if (pps1->slice_group_map_type == 0)
      {
        for (i=0; i<=pps1->num_slice_groups_minus1; i++)
          equal &= (pps1->run_length_minus1[i] == pps2->run_length_minus1[i]);
      }
      else if( pps1->slice_group_map_type == 2 )
      {
        for (i=0; i<pps1->num_slice_groups_minus1; i++)
        {
          equal &= (pps1->top_left[i] == pps2->top_left[i]);
          equal &= (pps1->bottom_right[i] == pps2->bottom_right[i]);
        }
      }
      else if( pps1->slice_group_map_type == 3 || pps1->slice_group_map_type==4 || pps1->slice_group_map_type==5 )
      {
        equal &= (pps1->slice_group_change_direction_flag == pps2->slice_group_change_direction_flag);
        equal &= (pps1->slice_group_change_rate_minus1 == pps2->slice_group_change_rate_minus1);
      }
      else if( pps1->slice_group_map_type == 6 )
      {
        equal &= (pps1->pic_size_in_map_units_minus1 == pps2->pic_size_in_map_units_minus1);
        if (!equal) return equal;
        for (i=0; i<=pps1->pic_size_in_map_units_minus1; i++)
          equal &= (pps1->slice_group_id[i] == pps2->slice_group_id[i]);
      }
  }

  equal &= (pps1->num_ref_idx_l0_default_active_minus1 == pps2->num_ref_idx_l0_default_active_minus1);
  equal &= (pps1->num_ref_idx_l1_default_active_minus1 == pps2->num_ref_idx_l1_default_active_minus1);
  equal &= (pps1->weighted_pred_flag == pps2->weighted_pred_flag);
  equal &= (pps1->weighted_bipred_idc == pps2->weighted_bipred_idc);
  equal &= (pps1->pic_init_qp_minus26 == pps2->pic_init_qp_minus26);
  equal &= (pps1->pic_init_qs_minus26 == pps2->pic_init_qs_minus26);
  equal &= (pps1->chroma_qp_index_offset == pps2->chroma_qp_index_offset);
  equal &= (pps1->deblocking_filter_control_present_flag == pps2->deblocking_filter_control_present_flag);
  equal &= (pps1->constrained_intra_pred_flag == pps2->constrained_intra_pred_flag);
  equal &= (pps1->redundant_pic_cnt_present_flag == pps2->redundant_pic_cnt_present_flag);

  if (!equal) return equal;

  //Fidelity Range Extensions Stuff
  //It is initialized to zero, so should be ok to check all the time.
  equal &= (pps1->transform_8x8_mode_flag == pps2->transform_8x8_mode_flag);
  equal &= (pps1->pic_scaling_matrix_present_flag == pps2->pic_scaling_matrix_present_flag);
  if(pps1->pic_scaling_matrix_present_flag)
  {
    for(i = 0; i < (6 + ((unsigned)pps1->transform_8x8_mode_flag << 1)); i++)
    {
      equal &= (pps1->pic_scaling_list_present_flag[i] == pps2->pic_scaling_list_present_flag[i]);
      if(pps1->pic_scaling_list_present_flag[i])
      {
        if(i < 6)
        {
          for (j = 0; j < 16; j++)
            equal &= (pps1->ScalingList4x4[i][j] == pps2->ScalingList4x4[i][j]);
        }
        else
        {
          for (j = 0; j < 64; j++)
            equal &= (pps1->ScalingList8x8[i-6][j] == pps2->ScalingList8x8[i-6][j]);
        }
      }
    }
  }
  equal &= (pps1->second_chroma_qp_index_offset == pps2->second_chroma_qp_index_offset);

  return equal;
}


void PPSConsistencyCheck (pps_t *pps)
{
  printf ("Consistency checking a picture parset, to be implemented\n");
//  if (pps->seq_parameter_set_id invalid then do something)
}

void SPSConsistencyCheck (sps_t *sps)
{
  printf ("Consistency checking a sequence parset, to be implemented\n");
}

#if (MVC_EXTENSION_ENABLE)
void SubsetSPSConsistencyCheck (sub_sps_t *subset_sps)
{
  printf ("Consistency checking a subset sequence parset, to be implemented\n");
}
#endif

void MakePPSavailable (VideoParameters *p_Vid, int id, pps_t *pps)
{
  assert (pps->Valid == TRUE);

  if (p_Vid->PicParSet[id].Valid == TRUE && p_Vid->PicParSet[id].slice_group_id != NULL)
    free (p_Vid->PicParSet[id].slice_group_id);

  memcpy (&p_Vid->PicParSet[id], pps, sizeof (pps_t));

  // we can simply use the memory provided with the pps. the PPS is destroyed after this function
  // call and will not try to free if pps->slice_group_id == NULL
  p_Vid->PicParSet[id].slice_group_id = pps->slice_group_id;
  pps->slice_group_id          = NULL;
}

void CleanUpPPS(VideoParameters *p_Vid)
{
  int i;

  for (i=0; i<MAXPPS; i++)
  {
    if (p_Vid->PicParSet[i].Valid == TRUE && p_Vid->PicParSet[i].slice_group_id != NULL)
      free (p_Vid->PicParSet[i].slice_group_id);

    p_Vid->PicParSet[i].Valid = FALSE;
  }
}


void MakeSPSavailable (VideoParameters *p_Vid, int id, sps_t *sps)
{
  assert (sps->Valid == TRUE);
  memcpy (&p_Vid->SeqParSet[id], sps, sizeof (sps_t));
}


void ProcessSPS(VideoParameters *p_Vid, NALU_t *nalu)
{  
    DataPartition *dp = AllocPartition(1);
    sps_t *sps = AllocSPS();

    InitPartition(dp, nalu);
    seq_parameter_set_rbsp(dp, sps);
#if (MVC_EXTENSION_ENABLE)
    get_max_dec_frame_buf_size(sps);
#endif

    if (sps->Valid) {
        if (p_Vid->active_sps) {
            if (sps->seq_parameter_set_id == p_Vid->active_sps->seq_parameter_set_id) {
                if (!sps_is_equal(sps, p_Vid->active_sps)) {
                    if (p_Vid->dec_picture)
                        // this may only happen on slice loss
                        exit_picture(p_Vid, &p_Vid->dec_picture);
                    p_Vid->active_sps = NULL;
                }
            }
        }
        // SPSConsistencyCheck (pps);
        MakeSPSavailable (p_Vid, sps->seq_parameter_set_id, sps);

#if (MVC_EXTENSION_ENABLE)
        if (p_Vid->profile_idc < (int) sps->profile_idc)
#endif
            p_Vid->profile_idc = sps->profile_idc;
    }

    FreePartition(dp, 1);
    FreeSPS(sps);
}

#if (MVC_EXTENSION_ENABLE)
void ProcessSubsetSPS(VideoParameters *p_Vid, NALU_t *nalu)
{
    DataPartition *dp = AllocPartition(1);
    sub_sps_t *subset_sps;
    int curr_seq_set_id;

    InitPartition(dp, nalu);
    subset_seq_parameter_set_rbsp(p_Vid, dp, &curr_seq_set_id);

    subset_sps = p_Vid->SubsetSeqParSet + curr_seq_set_id;
    get_max_dec_frame_buf_size(&(subset_sps->sps));
    //check capability;
    if (subset_sps->num_views_minus1 > 1) {
        printf("Warning: num_views:%d is greater than 2, only decode baselayer!\n", subset_sps->num_views_minus1+1);
        subset_sps->Valid = 0;
        subset_sps->sps.Valid = 0;
        p_Vid->p_Inp->DecodeAllLayers = 0;
    } else if (subset_sps->num_views_minus1==1 && (subset_sps->view_id[0]!=0 || subset_sps->view_id[1]!=1))
        OpenOutputFiles(p_Vid, subset_sps->view_id[0], subset_sps->view_id[1]);

    if (subset_sps->Valid) {
        // SubsetSPSConsistencyCheck (subset_sps);
        p_Vid->profile_idc = subset_sps->sps.profile_idc;
    }

    FreePartition(dp, 1);
}
#endif

void ProcessPPS(VideoParameters *p_Vid, NALU_t *nalu)
{
    DataPartition *dp = AllocPartition(1);
    pps_t *pps = AllocPPS();

    InitPartition(dp, nalu);
    pic_parameter_set_rbsp(p_Vid, dp, pps);
    // PPSConsistencyCheck (pps);
    if (p_Vid->active_pps) {
        if (pps->pic_parameter_set_id == p_Vid->active_pps->pic_parameter_set_id) {
            if (!pps_is_equal(pps, p_Vid->active_pps)) {
                //copy to next PPS;
                memcpy(p_Vid->pNextPPS, p_Vid->active_pps, sizeof (pps_t));
                if (p_Vid->dec_picture)
                    // this may only happen on slice loss
                    exit_picture(p_Vid, &p_Vid->dec_picture);
                p_Vid->active_pps = NULL;
            }
        }
    }
    MakePPSavailable(p_Vid, pps->pic_parameter_set_id, pps);
    FreePartition(dp, 1);
    FreePPS(pps);
}

static void init_qp_process(CodingParameters *cps)
{
  int bitdepth_qp_scale = imax(cps->bitdepth_luma_qp_scale, cps->bitdepth_chroma_qp_scale);
  int i;

  // We should allocate memory outside of this process since maybe we will have a change of SPS 
  // and we may need to recreate these. Currently should only support same bitdepth
  if (cps->qp_per_matrix == NULL)
    if ((cps->qp_per_matrix = (int*)malloc((MAX_QP + 1 +  bitdepth_qp_scale)*sizeof(int))) == NULL)
      no_mem_exit("init_qp_process: cps->qp_per_matrix");

  if (cps->qp_rem_matrix == NULL)
    if ((cps->qp_rem_matrix = (int*)malloc((MAX_QP + 1 +  bitdepth_qp_scale)*sizeof(int))) == NULL)
      no_mem_exit("init_qp_process: cps->qp_rem_matrix");

  for (i = 0; i < MAX_QP + bitdepth_qp_scale + 1; i++)
  {
    cps->qp_per_matrix[i] = i / 6;
    cps->qp_rem_matrix[i] = i % 6;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 *
 *  \par Input:
 *    Input Parameters VideoParameters *p_Vid
 *
 *  \par Output:
 *     Number of allocated bytes
 ***********************************************************************
 */
static int init_global_buffers(VideoParameters *p_Vid, int layer_id)
{
  int memory_size=0;
  int i;
  CodingParameters *cps = p_Vid->p_EncodePar[layer_id];
  BlockPos* PicPos;

  if (p_Vid->global_init_done[layer_id])
  {
    free_layer_buffers(p_Vid, layer_id);
  }

  // allocate memory for reference frame in find_snr
  memory_size += get_mem2Dpel(&cps->imgY_ref, cps->height, cps->width);
  if (cps->yuv_format != YUV400)
  {
    memory_size += get_mem3Dpel(&cps->imgUV_ref, 2, cps->height_cr, cps->width_cr);
  }
  else
    cps->imgUV_ref = NULL;

  // allocate memory in structure p_Vid
  if( (cps->separate_colour_plane_flag != 0) )
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      if(((cps->mb_data_JV[i]) = (Macroblock *) calloc(cps->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
        no_mem_exit("init_global_buffers: cps->mb_data_JV");
    }
    cps->mb_data = NULL;
  }
  else
  {
    if(((cps->mb_data) = (Macroblock *) calloc(cps->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
      no_mem_exit("init_global_buffers: cps->mb_data");
  }
  if( (cps->separate_colour_plane_flag != 0) )
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      if(((cps->intra_block_JV[i]) = (char*) calloc(cps->FrameSizeInMbs, sizeof(char))) == NULL)
        no_mem_exit("init_global_buffers: cps->intra_block_JV");
    }
    cps->intra_block = NULL;
  }
  else
  {
    if(((cps->intra_block) = (char*) calloc(cps->FrameSizeInMbs, sizeof(char))) == NULL)
      no_mem_exit("init_global_buffers: cps->intra_block");
  }


  if(((cps->PicPos) = (BlockPos*) calloc(cps->FrameSizeInMbs + 1, sizeof(BlockPos))) == NULL)
    no_mem_exit("init_global_buffers: PicPos");

  PicPos = cps->PicPos;
  for (i = 0; i < (int) cps->FrameSizeInMbs + 1;++i)
  {
    PicPos[i].x = (short) (i % cps->PicWidthInMbs);
    PicPos[i].y = (short) (i / cps->PicWidthInMbs);
  }

  if( (cps->separate_colour_plane_flag != 0) )
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      get_mem2D(&(cps->ipredmode_JV[i]), 4*cps->FrameHeightInMbs, 4*cps->PicWidthInMbs);
    }
    cps->ipredmode = NULL;
  }
  else
   memory_size += get_mem2D(&(cps->ipredmode), 4*cps->FrameHeightInMbs, 4*cps->PicWidthInMbs);

  // CAVLC mem
  memory_size += get_mem4D(&(cps->nz_coeff), cps->FrameSizeInMbs, 3, BLOCK_SIZE, BLOCK_SIZE);
  if( (cps->separate_colour_plane_flag != 0) )
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      get_mem2Dint(&(cps->siblock_JV[i]), cps->FrameHeightInMbs, cps->PicWidthInMbs);
      if(cps->siblock_JV[i]== NULL)
        no_mem_exit("init_global_buffers: p_Vid->siblock_JV");
    }
    cps->siblock = NULL;
  }
  else
  {
    memory_size += get_mem2Dint(&(cps->siblock), cps->FrameHeightInMbs, cps->PicWidthInMbs);
  }
  init_qp_process(cps);
  cps->oldFrameSizeInMbs = cps->FrameSizeInMbs;

  if(layer_id == 0 )
    init_output(cps, ((cps->pic_unit_bitsize_on_disk+7) >> 3));
  else
    cps->img2buf = p_Vid->p_EncodePar[0]->img2buf;
  p_Vid->global_init_done[layer_id] = 1;

  return (memory_size);
}

/*!
 ************************************************************************
 * \brief
 *    Updates images max values
 *
 ************************************************************************
 */
static void updateMaxValue(FrameFormat *format)
{
  format->max_value[0] = (1 << format->bit_depth[0]) - 1;
  format->max_value_sq[0] = format->max_value[0] * format->max_value[0];
  format->max_value[1] = (1 << format->bit_depth[1]) - 1;
  format->max_value_sq[1] = format->max_value[1] * format->max_value[1];
  format->max_value[2] = (1 << format->bit_depth[2]) - 1;
  format->max_value_sq[2] = format->max_value[2] * format->max_value[2];
}

/*!
 ************************************************************************
 * \brief
 *    Reset format information
 *
 ************************************************************************
 */
static void reset_format_info(sps_t *sps, VideoParameters *p_Vid, FrameFormat *source, FrameFormat *output)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  static const int SubWidthC  [4]= { 1, 2, 2, 1};
  static const int SubHeightC [4]= { 1, 2, 1, 1};

  int crop_left, crop_right;
  int crop_top, crop_bottom;

  // cropping for luma
  if (sps->frame_cropping_flag)
  {
    crop_left   = SubWidthC [sps->chroma_format_idc] * sps->frame_crop_left_offset;
    crop_right  = SubWidthC [sps->chroma_format_idc] * sps->frame_crop_right_offset;
    crop_top    = SubHeightC[sps->chroma_format_idc] * ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_top_offset;
    crop_bottom = SubHeightC[sps->chroma_format_idc] * ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  source->width[0] = p_Vid->width - crop_left - crop_right;
  source->height[0] = p_Vid->height - crop_top - crop_bottom;

  // cropping for chroma
  if (sps->frame_cropping_flag)
  {
    crop_left   = sps->frame_crop_left_offset;
    crop_right  = sps->frame_crop_right_offset;
    crop_top    = ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_top_offset;
    crop_bottom = ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  if ((sps->chroma_format_idc==YUV400) && p_Inp->write_uv)
  {
    source->width[1]  = (source->width[0] >> 1);
    source->width[2]  = source->width[1];
    source->height[1] = (source->height[0] >> 1);
    source->height[2] = source->height[1];
  }
  else
  {
    source->width[1]  = p_Vid->width_cr - crop_left - crop_right;
    source->width[2]  = source->width[1];
    source->height[1] = p_Vid->height_cr - crop_top - crop_bottom;
    source->height[2] = source->height[1];
  }

  output->width[0]  = p_Vid->width;
  source->width[1]  = p_Vid->width_cr;
  source->width[2]  = p_Vid->width_cr;
  output->height[0] = p_Vid->height;
  output->height[1] = p_Vid->height_cr;
  output->height[2] = p_Vid->height_cr;

  source->size_cmp[0] = source->width[0] * source->height[0];
  source->size_cmp[1] = source->width[1] * source->height[1];
  source->size_cmp[2] = source->size_cmp[1];
  source->size        = source->size_cmp[0] + source->size_cmp[1] + source->size_cmp[2];
  source->mb_width    = source->width[0]  / MB_BLOCK_SIZE;
  source->mb_height   = source->height[0] / MB_BLOCK_SIZE;

  // output size (excluding padding)
  output->size_cmp[0] = output->width[0] * output->height[0];
  output->size_cmp[1] = output->width[1] * output->height[1];
  output->size_cmp[2] = output->size_cmp[1];
  output->size        = output->size_cmp[0] + output->size_cmp[1] + output->size_cmp[2];
  output->mb_width    = output->width[0]  / MB_BLOCK_SIZE;
  output->mb_height   = output->height[0] / MB_BLOCK_SIZE;


  output->bit_depth[0] = source->bit_depth[0] = p_Vid->bitdepth_luma;
  output->bit_depth[1] = source->bit_depth[1] = p_Vid->bitdepth_chroma;
  output->bit_depth[2] = source->bit_depth[2] = p_Vid->bitdepth_chroma;  
  output->pic_unit_size_on_disk = (imax(output->bit_depth[0], output->bit_depth[1]) > 8) ? 16 : 8;
  output->pic_unit_size_shift3 = output->pic_unit_size_on_disk >> 3;

  output->frame_rate  = source->frame_rate;
  output->color_model = source->color_model;
  output->yuv_format  = source->yuv_format = (ColorFormat) sps->chroma_format_idc;

  output->auto_crop_bottom    = crop_bottom;
  output->auto_crop_right     = crop_right;
  output->auto_crop_bottom_cr = (crop_bottom * p_Vid->mb_cr_size_y) / MB_BLOCK_SIZE;
  output->auto_crop_right_cr  = (crop_right * p_Vid->mb_cr_size_x) / MB_BLOCK_SIZE;

  source->auto_crop_bottom    = output->auto_crop_bottom;
  source->auto_crop_right     = output->auto_crop_right;
  source->auto_crop_bottom_cr = output->auto_crop_bottom_cr;
  source->auto_crop_right_cr  = output->auto_crop_right_cr;

  updateMaxValue(source);
  updateMaxValue(output);

  if (p_Vid->first_sps) {
    p_Vid->first_sps = 0;
    if(!p_Inp->bDisplayDecParams) {
      fprintf(stdout,"Profile IDC  : %d\n", sps->profile_idc);
      fprintf(stdout,"Image Format : %dx%d (%dx%d)\n", source->width[0], source->height[0], p_Vid->width, p_Vid->height);
      if (p_Vid->yuv_format == YUV400)
        fprintf(stdout,"Color Format : 4:0:0 ");
      else if (p_Vid->yuv_format == YUV420)
        fprintf(stdout,"Color Format : 4:2:0 ");
      else if (p_Vid->yuv_format == YUV422)
        fprintf(stdout,"Color Format : 4:2:2 ");
      else
        fprintf(stdout,"Color Format : 4:4:4 ");

      fprintf(stdout,"(%d:%d:%d)\n", source->bit_depth[0], source->bit_depth[1], source->bit_depth[2]);
      fprintf(stdout,"--------------------------------------------------------------------------\n");
    }
    if (!p_Inp->silent)
    {
      fprintf(stdout,"POC must = frame# or field# for SNRs to be correct\n");
      fprintf(stdout,"--------------------------------------------------------------------------\n");
      fprintf(stdout,"  Frame          POC  Pic#   QP    SnrY     SnrU     SnrV   Y:U:V Time(ms)\n");
      fprintf(stdout,"--------------------------------------------------------------------------\n");
    }
  }
}

static void setup_layer_info(VideoParameters *p_Vid, sps_t *sps, LayerParameters *p_Lps)
{
  int layer_id = p_Lps->layer_id;
  p_Lps->p_Vid = p_Vid;
  p_Lps->p_Cps = p_Vid->p_EncodePar[layer_id];
  p_Lps->p_SPS = sps;
  p_Lps->p_Dpb = p_Vid->p_Dpb_layer[layer_id];
}

static void set_coding_par(sps_t *sps, CodingParameters *cps)
{
  // maximum vertical motion vector range in luma quarter pixel units
  cps->profile_idc = sps->profile_idc;
  cps->lossless_qpprime_flag   = sps->qpprime_y_zero_transform_bypass_flag;
  if (sps->level_idc <= 10)
  {
    cps->max_vmv_r = 64 * 4;
  }
  else if (sps->level_idc <= 20)
  {
    cps->max_vmv_r = 128 * 4;
  }
  else if (sps->level_idc <= 30)
  {
    cps->max_vmv_r = 256 * 4;
  }
  else
  {
    cps->max_vmv_r = 512 * 4; // 512 pixels in quarter pixels
  }

  // Fidelity Range Extensions stuff (part 1)
  cps->bitdepth_chroma = 0;
  cps->width_cr        = 0;
  cps->height_cr       = 0;
  cps->bitdepth_luma       = (short) (sps->bit_depth_luma_minus8 + 8);
  cps->bitdepth_scale[0]   = 1 << sps->bit_depth_luma_minus8;
  if (sps->chroma_format_idc != YUV400)
  {
    cps->bitdepth_chroma   = (short) (sps->bit_depth_chroma_minus8 + 8);
    cps->bitdepth_scale[1] = 1 << sps->bit_depth_chroma_minus8;
  }

  cps->max_frame_num = 1<<(sps->log2_max_frame_num_minus4+4);
  cps->PicWidthInMbs = (sps->pic_width_in_mbs_minus1 +1);
  cps->PicHeightInMapUnits = (sps->pic_height_in_map_units_minus1 +1);
  cps->FrameHeightInMbs = ( 2 - sps->frame_mbs_only_flag ) * cps->PicHeightInMapUnits;
  cps->FrameSizeInMbs = cps->PicWidthInMbs * cps->FrameHeightInMbs;

  cps->yuv_format=sps->chroma_format_idc;
  cps->separate_colour_plane_flag = sps->separate_colour_plane_flag;

  cps->width = cps->PicWidthInMbs * MB_BLOCK_SIZE;
  cps->height = cps->FrameHeightInMbs * MB_BLOCK_SIZE;  

  cps->iLumaPadX = MCBUF_LUMA_PAD_X;
  cps->iLumaPadY = MCBUF_LUMA_PAD_Y;
  cps->iChromaPadX = MCBUF_CHROMA_PAD_X;
  cps->iChromaPadY = MCBUF_CHROMA_PAD_Y;
  if (sps->chroma_format_idc == YUV420)
  {
    cps->width_cr  = (cps->width  >> 1);
    cps->height_cr = (cps->height >> 1);
  }
  else if (sps->chroma_format_idc == YUV422)
  {
    cps->width_cr  = (cps->width >> 1);
    cps->height_cr = cps->height;
    cps->iChromaPadY = MCBUF_CHROMA_PAD_Y*2;
  }
  else if (sps->chroma_format_idc == YUV444)
  {
    //YUV444
    cps->width_cr = cps->width;
    cps->height_cr = cps->height;
    cps->iChromaPadX = cps->iLumaPadX;
    cps->iChromaPadY = cps->iLumaPadY;
  }
  //pel bitdepth init
  cps->bitdepth_luma_qp_scale   = 6 * (cps->bitdepth_luma - 8);

  if(cps->bitdepth_luma > cps->bitdepth_chroma || sps->chroma_format_idc == YUV400)
    cps->pic_unit_bitsize_on_disk = (cps->bitdepth_luma > 8)? 16:8;
  else
    cps->pic_unit_bitsize_on_disk = (cps->bitdepth_chroma > 8)? 16:8;
  cps->dc_pred_value_comp[0]    = 1<<(cps->bitdepth_luma - 1);
  cps->max_pel_value_comp[0] = (1<<cps->bitdepth_luma) - 1;
  cps->mb_size[0][0] = cps->mb_size[0][1] = MB_BLOCK_SIZE;

  if (sps->chroma_format_idc != YUV400)
  {
    //for chrominance part
    cps->bitdepth_chroma_qp_scale = 6 * (cps->bitdepth_chroma - 8);
    cps->dc_pred_value_comp[1]    = (1 << (cps->bitdepth_chroma - 1));
    cps->dc_pred_value_comp[2]    = cps->dc_pred_value_comp[1];
    cps->max_pel_value_comp[1]    = (1 << cps->bitdepth_chroma) - 1;
    cps->max_pel_value_comp[2]    = (1 << cps->bitdepth_chroma) - 1;
    cps->num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    cps->num_uv_blocks = (cps->num_blk8x8_uv >> 1);
    cps->num_cdc_coeff = (cps->num_blk8x8_uv << 1);
    cps->mb_size[1][0] = cps->mb_size[2][0] = cps->mb_cr_size_x  = (sps->chroma_format_idc==YUV420 || sps->chroma_format_idc==YUV422)?  8 : 16;
    cps->mb_size[1][1] = cps->mb_size[2][1] = cps->mb_cr_size_y  = (sps->chroma_format_idc==YUV444 || sps->chroma_format_idc==YUV422)? 16 :  8;

    cps->subpel_x    = cps->mb_cr_size_x == 8 ? 7 : 3;
    cps->subpel_y    = cps->mb_cr_size_y == 8 ? 7 : 3;
    cps->shiftpel_x  = cps->mb_cr_size_x == 8 ? 3 : 2;
    cps->shiftpel_y  = cps->mb_cr_size_y == 8 ? 3 : 2;
    cps->total_scale = cps->shiftpel_x + cps->shiftpel_y;
  }
  else
  {
    cps->bitdepth_chroma_qp_scale = 0;
    cps->max_pel_value_comp[1] = 0;
    cps->max_pel_value_comp[2] = 0;
    cps->num_blk8x8_uv = 0;
    cps->num_uv_blocks = 0;
    cps->num_cdc_coeff = 0;
    cps->mb_size[1][0] = cps->mb_size[2][0] = cps->mb_cr_size_x  = 0;
    cps->mb_size[1][1] = cps->mb_size[2][1] = cps->mb_cr_size_y  = 0;
    cps->subpel_x      = 0;
    cps->subpel_y      = 0;
    cps->shiftpel_x    = 0;
    cps->shiftpel_y    = 0;
    cps->total_scale   = 0;
  }

  cps->mb_cr_size = cps->mb_cr_size_x * cps->mb_cr_size_y;
  cps->mb_size_blk[0][0] = cps->mb_size_blk[0][1] = cps->mb_size[0][0] >> 2;
  cps->mb_size_blk[1][0] = cps->mb_size_blk[2][0] = cps->mb_size[1][0] >> 2;
  cps->mb_size_blk[1][1] = cps->mb_size_blk[2][1] = cps->mb_size[1][1] >> 2;

  cps->mb_size_shift[0][0] = cps->mb_size_shift[0][1] = CeilLog2_sf (cps->mb_size[0][0]);
  cps->mb_size_shift[1][0] = cps->mb_size_shift[2][0] = CeilLog2_sf (cps->mb_size[1][0]);
  cps->mb_size_shift[1][1] = cps->mb_size_shift[2][1] = CeilLog2_sf (cps->mb_size[1][1]);

  cps->rgb_output =  (sps->vui_parameters.matrix_coefficients==0);
}

/*!
 ***********************************************************************
 * \brief
 *    Initialize FREXT variables
 ***********************************************************************
 */
void init_frext(VideoParameters *p_Vid)  //!< video parameters
{
  //pel bitdepth init
  p_Vid->bitdepth_luma_qp_scale   = 6 * (p_Vid->bitdepth_luma - 8);

  if(p_Vid->bitdepth_luma > p_Vid->bitdepth_chroma || p_Vid->active_sps->chroma_format_idc == YUV400)
    p_Vid->pic_unit_bitsize_on_disk = (p_Vid->bitdepth_luma > 8)? 16:8;
  else
    p_Vid->pic_unit_bitsize_on_disk = (p_Vid->bitdepth_chroma > 8)? 16:8;
  p_Vid->dc_pred_value_comp[0]    = 1<<(p_Vid->bitdepth_luma - 1);
  p_Vid->max_pel_value_comp[0] = (1<<p_Vid->bitdepth_luma) - 1;
  p_Vid->mb_size[0][0] = p_Vid->mb_size[0][1] = MB_BLOCK_SIZE;

  if (p_Vid->active_sps->chroma_format_idc != YUV400)
  {
    //for chrominance part
    p_Vid->bitdepth_chroma_qp_scale = 6 * (p_Vid->bitdepth_chroma - 8);
    p_Vid->dc_pred_value_comp[1]    = (1 << (p_Vid->bitdepth_chroma - 1));
    p_Vid->dc_pred_value_comp[2]    = p_Vid->dc_pred_value_comp[1];
    p_Vid->max_pel_value_comp[1]    = (1 << p_Vid->bitdepth_chroma) - 1;
    p_Vid->max_pel_value_comp[2]    = (1 << p_Vid->bitdepth_chroma) - 1;
    p_Vid->num_blk8x8_uv = (1 << p_Vid->active_sps->chroma_format_idc) & (~(0x1));
    p_Vid->num_uv_blocks = (p_Vid->num_blk8x8_uv >> 1);
    p_Vid->num_cdc_coeff = (p_Vid->num_blk8x8_uv << 1);
    p_Vid->mb_size[1][0] = p_Vid->mb_size[2][0] = p_Vid->mb_cr_size_x  = (p_Vid->active_sps->chroma_format_idc==YUV420 || p_Vid->active_sps->chroma_format_idc==YUV422)?  8 : 16;
    p_Vid->mb_size[1][1] = p_Vid->mb_size[2][1] = p_Vid->mb_cr_size_y  = (p_Vid->active_sps->chroma_format_idc==YUV444 || p_Vid->active_sps->chroma_format_idc==YUV422)? 16 :  8;

    p_Vid->subpel_x    = p_Vid->mb_cr_size_x == 8 ? 7 : 3;
    p_Vid->subpel_y    = p_Vid->mb_cr_size_y == 8 ? 7 : 3;
    p_Vid->shiftpel_x  = p_Vid->mb_cr_size_x == 8 ? 3 : 2;
    p_Vid->shiftpel_y  = p_Vid->mb_cr_size_y == 8 ? 3 : 2;
    p_Vid->total_scale = p_Vid->shiftpel_x + p_Vid->shiftpel_y;
  }
  else
  {
    p_Vid->bitdepth_chroma_qp_scale = 0;
    p_Vid->max_pel_value_comp[1] = 0;
    p_Vid->max_pel_value_comp[2] = 0;
    p_Vid->num_blk8x8_uv = 0;
    p_Vid->num_uv_blocks = 0;
    p_Vid->num_cdc_coeff = 0;
    p_Vid->mb_size[1][0] = p_Vid->mb_size[2][0] = p_Vid->mb_cr_size_x  = 0;
    p_Vid->mb_size[1][1] = p_Vid->mb_size[2][1] = p_Vid->mb_cr_size_y  = 0;
    p_Vid->subpel_x      = 0;
    p_Vid->subpel_y      = 0;
    p_Vid->shiftpel_x    = 0;
    p_Vid->shiftpel_y    = 0;
    p_Vid->total_scale   = 0;
  }

  p_Vid->mb_cr_size = p_Vid->mb_cr_size_x * p_Vid->mb_cr_size_y;
  p_Vid->mb_size_blk[0][0] = p_Vid->mb_size_blk[0][1] = p_Vid->mb_size[0][0] >> 2;
  p_Vid->mb_size_blk[1][0] = p_Vid->mb_size_blk[2][0] = p_Vid->mb_size[1][0] >> 2;
  p_Vid->mb_size_blk[1][1] = p_Vid->mb_size_blk[2][1] = p_Vid->mb_size[1][1] >> 2;

  p_Vid->mb_size_shift[0][0] = p_Vid->mb_size_shift[0][1] = CeilLog2_sf (p_Vid->mb_size[0][0]);
  p_Vid->mb_size_shift[1][0] = p_Vid->mb_size_shift[2][0] = CeilLog2_sf (p_Vid->mb_size[1][0]);
  p_Vid->mb_size_shift[1][1] = p_Vid->mb_size_shift[2][1] = CeilLog2_sf (p_Vid->mb_size[1][1]);
}

static void set_global_coding_par(VideoParameters *p_Vid, CodingParameters *cps)
{
    p_Vid->bitdepth_chroma = 0;
    p_Vid->width_cr        = 0;
    p_Vid->height_cr       = 0;
    p_Vid->lossless_qpprime_flag   = cps->lossless_qpprime_flag;
    p_Vid->max_vmv_r = cps->max_vmv_r;

    // Fidelity Range Extensions stuff (part 1)
    p_Vid->bitdepth_luma       = cps->bitdepth_luma;
    p_Vid->bitdepth_scale[0]   = cps->bitdepth_scale[0];
    p_Vid->bitdepth_chroma = cps->bitdepth_chroma;
    p_Vid->bitdepth_scale[1] = cps->bitdepth_scale[1];

    p_Vid->max_frame_num = cps->max_frame_num;
    p_Vid->PicWidthInMbs = cps->PicWidthInMbs;
    p_Vid->PicHeightInMapUnits = cps->PicHeightInMapUnits;
    p_Vid->FrameHeightInMbs = cps->FrameHeightInMbs;
    p_Vid->FrameSizeInMbs = cps->FrameSizeInMbs;

    p_Vid->yuv_format = cps->yuv_format;

    p_Vid->width = cps->width;
    p_Vid->height = cps->height;
    p_Vid->iLumaPadX = MCBUF_LUMA_PAD_X;
    p_Vid->iLumaPadY = MCBUF_LUMA_PAD_Y;
    p_Vid->iChromaPadX = MCBUF_CHROMA_PAD_X;
    p_Vid->iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (p_Vid->yuv_format == YUV420)
    {
      p_Vid->width_cr  = (p_Vid->width  >> 1);
      p_Vid->height_cr = (p_Vid->height >> 1);
    }
    else if (p_Vid->yuv_format == YUV422)
    {
      p_Vid->width_cr  = (p_Vid->width >> 1);
      p_Vid->height_cr = p_Vid->height;
      p_Vid->iChromaPadY = MCBUF_CHROMA_PAD_Y*2;
    }
    else if (p_Vid->yuv_format == YUV444)
    {
      //YUV444
      p_Vid->width_cr = p_Vid->width;
      p_Vid->height_cr = p_Vid->height;
      p_Vid->iChromaPadX = p_Vid->iLumaPadX;
      p_Vid->iChromaPadY = p_Vid->iLumaPadY;
    }

    init_frext(p_Vid);
}

/*!
 ************************************************************************
 * \brief
 *    Activate Sequence Parameter Sets
 *
 ************************************************************************
 */
void activate_sps (VideoParameters *p_Vid, sps_t *sps)
{
  InputParameters *p_Inp = p_Vid->p_Inp;  

  if (p_Vid->active_sps != sps)
  {
    if (p_Vid->dec_picture)
    {
      // this may only happen on slice loss
      exit_picture(p_Vid, &p_Vid->dec_picture);
    }
    p_Vid->active_sps = sps;

    if(p_Vid->dpb_layer_id==0 && is_BL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[0]->init_done)
    {
      set_coding_par(sps, p_Vid->p_EncodePar[0]);
      setup_layer_info( p_Vid, sps, p_Vid->p_LayerPar[0]);
    }
    else if(p_Vid->dpb_layer_id==1 && is_EL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[1]->init_done)
    {
      set_coding_par(sps, p_Vid->p_EncodePar[1]);
      setup_layer_info(p_Vid, sps, p_Vid->p_LayerPar[1]);
    }

//to be removed in future;
    set_global_coding_par(p_Vid, p_Vid->p_EncodePar[p_Vid->dpb_layer_id]);
//end;

#if (MVC_EXTENSION_ENABLE)
    if (/*p_Vid->last_pic_width_in_mbs_minus1 != p_Vid->active_sps->pic_width_in_mbs_minus1
        || p_Vid->last_pic_height_in_map_units_minus1 != p_Vid->active_sps->pic_height_in_map_units_minus1
        || p_Vid->last_max_dec_frame_buffering != GetMaxDecFrameBuffering(p_Vid)
        || */(p_Vid->last_profile_idc != p_Vid->active_sps->profile_idc && is_BL_profile(p_Vid->active_sps->profile_idc) && !p_Vid->p_Dpb_layer[0]->init_done /*&& is_BL_profile(p_Vid->last_profile_idc)*/))
    {
      init_global_buffers(p_Vid, 0);

      if (!p_Vid->no_output_of_prior_pics_flag)
      {
        flush_dpb(p_Vid->p_Dpb_layer[0]);
        flush_dpb(p_Vid->p_Dpb_layer[1]);
      }
      init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);
    }
    else if(p_Vid->last_profile_idc != p_Vid->active_sps->profile_idc && (
            is_MVC_profile(p_Vid->last_profile_idc) || is_MVC_profile(p_Vid->active_sps->profile_idc)
            )&& (!p_Vid->p_Dpb_layer[1]->init_done))
    {
      assert(p_Vid->p_Dpb_layer[0]->init_done);
      if(p_Vid->p_Dpb_layer[0]->init_done)
      {
        free_dpb(p_Vid->p_Dpb_layer[0]);
        init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);
      }
      init_global_buffers(p_Vid, 1);
      // for now lets re_init both buffers. Later, we should only re_init appropriate one
      // Note that we seem to be doing this for every frame which seems not good.
      //re_init_dpb(p_Vid, p_Vid->p_Dpb_layer[1], 2);
#if MVC_EXTENSION_ENABLE
      init_dpb(p_Vid, p_Vid->p_Dpb_layer[1], 2);
#endif
      //p_Vid->last_profile_idc = p_Vid->active_sps->profile_idc;
    }
    //p_Vid->p_Dpb_layer[0]->num_ref_frames = p_Vid->active_sps->num_ref_frames;
    //p_Vid->p_Dpb_layer[1]->num_ref_frames = p_Vid->active_sps->num_ref_frames;
    p_Vid->last_pic_width_in_mbs_minus1 = p_Vid->active_sps->pic_width_in_mbs_minus1;  
    p_Vid->last_pic_height_in_map_units_minus1 = p_Vid->active_sps->pic_height_in_map_units_minus1;
    p_Vid->last_max_dec_frame_buffering = GetMaxDecFrameBuffering(p_Vid);
    p_Vid->last_profile_idc = p_Vid->active_sps->profile_idc;

#else
    init_global_buffers(p_Vid, 0);

    if (!p_Vid->no_output_of_prior_pics_flag)
    {
      flush_dpb(p_Vid->p_Dpb_layer[0]);
    }
    init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 0);
    // for now lets init both buffers. Later, we should only re_init appropriate one
    //init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);
    // obviously this is not needed her but just adding it for completeness
    //init_dpb(p_Vid, p_Vid->p_Dpb_layer[1], 2);
#endif

#if (DISABLE_ERC == 0)
    ercInit(p_Vid, p_Vid->width, p_Vid->height, 1);
    if(p_Vid->dec_picture)
    {
      ercReset(p_Vid->erc_errorVar, p_Vid->PicSizeInMbs, p_Vid->PicSizeInMbs, p_Vid->dec_picture->size_x);
      p_Vid->erc_mvperMB = 0;
    }
#endif
  }
  
  reset_format_info(sps, p_Vid, &p_Inp->source, &p_Inp->output);
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



void UseParameterSet (Slice *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  int PicParsetId = currSlice->pic_parameter_set_id;  
  pps_t *pps = &p_Vid->PicParSet[PicParsetId];
  sps_t *sps = &p_Vid->SeqParSet[pps->seq_parameter_set_id];
  int i;

  if (pps->Valid != TRUE)
    printf ("Trying to use an invalid (uninitialized) Picture Parameter Set with ID %d, expect the unexpected...\n", PicParsetId);
#if (MVC_EXTENSION_ENABLE)
  if (currSlice->svc_extension_flag == -1)
  {
    if (sps->Valid != TRUE)
      printf ("PicParset %d references an invalid (uninitialized) Sequence Parameter Set with ID %d, expect the unexpected...\n", 
      PicParsetId, (int) pps->seq_parameter_set_id);
  }
  else
  {
    // Set SPS to the subset SPS parameters
    p_Vid->active_subset_sps = p_Vid->SubsetSeqParSet + pps->seq_parameter_set_id;
    sps = &(p_Vid->active_subset_sps->sps);
    if (p_Vid->SubsetSeqParSet[pps->seq_parameter_set_id].Valid != TRUE)
      printf ("PicParset %d references an invalid (uninitialized) Subset Sequence Parameter Set with ID %d, expect the unexpected...\n", 
      PicParsetId, (int) pps->seq_parameter_set_id);
  }
#else
  if (sps->Valid != TRUE)
    printf ("PicParset %d references an invalid (uninitialized) Sequence Parameter Set with ID %d, expect the unexpected...\n", 
    PicParsetId, (int) pps->seq_parameter_set_id);
#endif

  // In theory, and with a well-designed software, the lines above
  // are everything necessary.  In practice, we need to patch many values
  // in p_Vid-> (but no more in p_Inp-> -- these have been taken care of)

  // Set Sequence Parameter Stuff first
  //  printf ("Using Picture Parameter set %d and associated Sequence Parameter Set %d\n", PicParsetId, pps->seq_parameter_set_id);
  if ((int) sps->pic_order_cnt_type < 0 || sps->pic_order_cnt_type > 2)  // != 1
  {
    printf("invalid sps->pic_order_cnt_type = %d\n", (int) sps->pic_order_cnt_type);
    error("pic_order_cnt_type != 1", -1000);
  }

  if (sps->pic_order_cnt_type == 1)
  {
    if(sps->num_ref_frames_in_pic_order_cnt_cycle >= MAX_NUM_REF_FRAMES)
    {
      error("num_ref_frames_in_pic_order_cnt_cycle too large",-1011);
    }
  }
  p_Vid->dpb_layer_id = currSlice->layer_id;
  activate_sps(p_Vid, sps);
  activate_pps(p_Vid, pps);

  // currSlice->dp_mode is set by read_new_slice (NALU first byte available there)
  if (pps->entropy_coding_mode_flag == (Boolean) CAVLC)
  {
    for (i=0; i<3; i++)
    {
      currSlice->partArr[i].readSyntaxElement = readSyntaxElement_UVLC;      
    }
  }
  else
  {
    for (i=0; i<3; i++)
    {
      currSlice->partArr[i].readSyntaxElement = readSyntaxElement_CABAC;
    }
  }
  p_Vid->type = currSlice->slice_type;
}

#if (MVC_EXTENSION_ENABLE)

void init_subset_sps_list(sub_sps_t *subset_sps_list, int iSize)
{
  int i;
  memset(subset_sps_list, 0, iSize*sizeof(subset_sps_list[0]));
  for(i=0; i<iSize; i++)
  {
    subset_sps_list[i].sps.seq_parameter_set_id = -1;
    subset_sps_list[i].num_views_minus1 = -1;
    subset_sps_list[i].num_level_values_signalled_minus1 = -1;
    subset_sps_list[i].MVCVUIParams.num_ops_minus1 = -1;
  }
}

void reset_subset_sps(sub_sps_t *subset_sps)
{
  int i, j;

  if(subset_sps && subset_sps->num_views_minus1>=0)
  {
    subset_sps->sps.seq_parameter_set_id = -1;

    free_pointer(subset_sps->view_id);
    for(i=0; i<=subset_sps->num_views_minus1; i++)
    {
      free_pointer(subset_sps->anchor_ref_l0[i]);
      free_pointer(subset_sps->anchor_ref_l1[i]);
    }
    free_pointer(subset_sps->anchor_ref_l0);
    free_pointer(subset_sps->anchor_ref_l1);
    free_pointer(subset_sps->num_anchor_refs_l0);
    free_pointer(subset_sps->num_anchor_refs_l1);

    for(i=0; i<=subset_sps->num_views_minus1; i++)
    {
      free_pointer(subset_sps->non_anchor_ref_l0[i]);
      free_pointer(subset_sps->non_anchor_ref_l1[i]);
    }
    free_pointer(subset_sps->non_anchor_ref_l0);
    free_pointer(subset_sps->non_anchor_ref_l1);
    free_pointer(subset_sps->num_non_anchor_refs_l0);
    free_pointer(subset_sps->num_non_anchor_refs_l1);

    if(subset_sps->num_level_values_signalled_minus1 >= 0)
    {
      free_pointer(subset_sps->level_idc);
      for(i=0; i<=subset_sps->num_level_values_signalled_minus1; i++)
      {
        for(j=0; j<=subset_sps->num_applicable_ops_minus1[i]; j++)
        {
          free_pointer(subset_sps->applicable_op_target_view_id[i][j]);
        }
        free_pointer(subset_sps->applicable_op_target_view_id[i]);
        free_pointer(subset_sps->applicable_op_temporal_id[i]);
        free_pointer(subset_sps->applicable_op_num_target_views_minus1[i]);
        free_pointer(subset_sps->applicable_op_num_views_minus1[i]);
      }
      free_pointer(subset_sps->applicable_op_target_view_id);
      free_pointer(subset_sps->applicable_op_temporal_id);
      free_pointer(subset_sps->applicable_op_num_target_views_minus1);
      free_pointer(subset_sps->applicable_op_num_views_minus1);      
      free_pointer(subset_sps->num_applicable_ops_minus1);

      subset_sps->num_level_values_signalled_minus1 = -1;
    }

    //end;
    subset_sps->num_views_minus1 = -1;
  }

  if(subset_sps && subset_sps->mvc_vui_parameters_present_flag)
  {
    mvc_vui_t *pMVCVUI = &(subset_sps->MVCVUIParams);
    if(pMVCVUI->num_ops_minus1 >=0)
    {
      free_pointer(pMVCVUI->temporal_id);
      free_pointer(pMVCVUI->num_target_output_views_minus1);
      for(i=0; i<=pMVCVUI->num_ops_minus1; i++)
        free_pointer(pMVCVUI->view_id[i]);
      free_pointer(pMVCVUI->view_id);
      free_pointer(pMVCVUI->timing_info_present_flag);
      free_pointer(pMVCVUI->num_units_in_tick);
      free_pointer(pMVCVUI->time_scale);
      free_pointer(pMVCVUI->fixed_frame_rate_flag);
      free_pointer(pMVCVUI->nal_hrd_parameters_present_flag);
      free_pointer(pMVCVUI->vcl_hrd_parameters_present_flag);
      free_pointer(pMVCVUI->low_delay_hrd_flag);
      free_pointer(pMVCVUI->pic_struct_present_flag);

      pMVCVUI->num_ops_minus1 = -1;
    }
    subset_sps->mvc_vui_parameters_present_flag = 0;
  }
}

int GetBaseViewId(VideoParameters *p_Vid, sub_sps_t **subset_sps)
{
  sub_sps_t *curr_subset_sps;
  int i, iBaseViewId=0; //-1;

  *subset_sps = NULL;
  curr_subset_sps = p_Vid->SubsetSeqParSet;
  for(i=0; i<MAXSPS; i++)
  {
    if(curr_subset_sps->num_views_minus1>=0 && curr_subset_sps->sps.Valid) // && curr_subset_sps->sps.seq_parameter_set_id < MAXSPS)
    {
      iBaseViewId = curr_subset_sps->view_id[BASE_VIEW_IDX];
      break;
    }
    curr_subset_sps++;
  }

  if(i<MAXSPS)
    *subset_sps = curr_subset_sps;
  return iBaseViewId;
}
#endif

/*!
 *************************************************************************************
 * \brief
 *    Allocates memory for a picture paramater set
 *
 * \return
 *    pointer to a pps
 *************************************************************************************
 */

pps_t *AllocPPS ()
 {
   pps_t *p;

   if ((p=(pps_t *)calloc (1, sizeof (pps_t))) == NULL)
     no_mem_exit ("AllocPPS: PPS");
   p->slice_group_id = NULL;
   return p;
 }


/*!
 *************************************************************************************
 * \brief
 *    Allocates memory for am sequence paramater set
 *
 * \return
 *    pointer to a sps
 *************************************************************************************
 */

sps_t *AllocSPS ()
 {
   sps_t *p;

   if ((p=(sps_t *)calloc (1, sizeof (sps_t))) == NULL)
     no_mem_exit ("AllocSPS: SPS");
   return p;
 }


/*!
 *************************************************************************************
 * \brief
 *    Frees a picture parameter set
 *
 * \param pps to be freed
 *   Picture parameter set to be freed
 *************************************************************************************
 */

 void FreePPS (pps_t *pps)
 {
   assert (pps != NULL);
   if (pps->slice_group_id != NULL) 
     free (pps->slice_group_id);
   free (pps);
 }


 /*!
 *************************************************************************************
 * \brief
 *    Frees a sps
 *
 * \param sps
 *   Sequence parameter set to be freed
 *************************************************************************************
 */

 void FreeSPS (sps_t *sps)
 {
   assert (sps != NULL);
   free (sps);
 }

