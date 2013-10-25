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
 *  File      : interpret_sei.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * =============================================================================
 */

#include "interpret.h"
#include "global.h"
#include "sei.h"
#include "slice.h"
#include "sets.h"


using namespace vio::h264;
using vio::h264::InterpreterRbsp;
using vio::h264::vui_t;
using vio::h264::hrd_t;

struct tone_mapping_t {
    uint32_t tone_map_id;
    bool     tone_map_cancel_flag;
    uint32_t tone_map_repetition_period;
    uint8_t  coded_data_bit_depth;
    uint8_t  target_bit_depth;
    uint32_t tone_map_model_id;
    // variables for model 0
    uint32_t min_value;
    uint32_t max_value;
    // variables for model 1
    uint32_t sigmoid_midpoint;
    uint32_t sigmoid_width;
    // variables for model 2
    uint32_t start_of_coded_interval[1<<MAX_SEI_BIT_DEPTH];
    // variables for model 3
    uint16_t num_pivots;
    uint32_t coded_pivot_value[MAX_NUM_PIVOTS];
    uint32_t target_pivot_value[MAX_NUM_PIVOTS];
};

enum {
    SEI_BUFFERING_PERIOD = 0,
    SEI_PIC_TIMING,
    SEI_PAN_SCAN_RECT,
    SEI_FILLER_PAYLOAD,
    SEI_USER_DATA_REGISTERED_ITU_T_T35,
    SEI_USER_DATA_UNREGISTERED,
    SEI_RECOVERY_POINT,
    SEI_DEC_REF_PIC_MARKING_REPETITION,
    SEI_SPARE_PIC,
    SEI_SCENE_INFO,
    SEI_SUB_SEQ_INFO,
    SEI_SUB_SEQ_LAYER_CHARACTERISTICS,
    SEI_SUB_SEQ_CHARACTERISTICS,
    SEI_FULL_FRAME_FREEZE,
    SEI_FULL_FRAME_FREEZE_RELEASE,
    SEI_FULL_FRAME_SNAPSHOT,
    SEI_PROGRESSIVE_REFINEMENT_SEGMENT_START,
    SEI_PROGRESSIVE_REFINEMENT_SEGMENT_END,
    SEI_MOTION_CONSTRAINED_SLICE_GROUP_SET,
    SEI_FILM_GRAIN_CHARACTERISTICS,
    SEI_DEBLOCKING_FILTER_DISPLAY_PREFERENCE,
    SEI_STEREO_VIDEO_INFO,
    SEI_POST_FILTER_HINTS,
    SEI_TONE_MAPPING_INFO,
    SEI_SCALABILITY_INFO,
    SEI_SUB_PIC_SCALABLE_LAYER,
    SEI_NON_REQUIRED_LAYER_REP,
    SEI_PRIORITY_LAYER_INFO,
    SEI_LAYERS_NOT_PRESENT,
    SEI_LAYER_DEPENDENCY_CHANGE,
    SEI_SCALABLE_NESTING,
    SEI_BASE_LAYER_TEMPORAL_HRD,
    SEI_QUALITY_LAYER_INTEGRITY_CHECK,
    SEI_REDUNDANT_PIC_PROPERTY,
    SEI_TL0_DEP_REP_INDEX,
    SEI_TL_SWITCHING_POINT,
    SEI_PARALLEL_DECODING_INFO,
    SEI_MVC_SCALABLE_NESTING,
    SEI_VIEW_SCALABILITY_INFO,
    SEI_MULTIVIEW_SCENE_INFO,
    SEI_MULTIVIEW_ACQUISITION_INFO,
    SEI_NON_REQUIRED_VIEW_COMPONENT,
    SEI_VIEW_DEPENDENCY_CHANGE,
    SEI_OPERATION_POINTS_NOT_PRESENT,
    SEI_BASE_VIEW_TEMPORAL_HRD,
    SEI_FRAME_PACKING_ARRANGEMENT,
    SEI_MULTIVIEW_VIEW_POSITION,
    SEI_DISPLAY_ORIENTATION,
    SEI_MVCD_SCALABLE_NESTING,
    SEI_MVCD_VIEW_SCALABILITY_INFO,
    SEI_DEPTH_REPRESENTATION_INFO,
    SEI_THREE_DIMENSIONAL_REFERENCE_DISPLAYS_INFO,
    SEI_DEPTH_TIMING,
    SEI_DEPTH_SAMPLING_INFO
};


// D.1 SEI payload syntax

void InterpreterSEI::sei_payload(uint32_t payloadType, uint32_t payloadSize)
{
    this->num_bytes_in_rbsp = this->frame_bitoffset / 8 + payloadSize;

    switch (payloadType) {
    case SEI_BUFFERING_PERIOD:
        this->buffering_period(payloadSize);
        break;
    case SEI_PIC_TIMING:
        this->pic_timing(payloadSize);
        break;
    case SEI_PAN_SCAN_RECT:
        this->pan_scan_rect(payloadSize);
        break;
    case SEI_FILLER_PAYLOAD:
        this->filler_payload(payloadSize);
        break;
    case SEI_USER_DATA_REGISTERED_ITU_T_T35:
        this->user_data_registered_itu_t_t35(payloadSize);
        break;
    case SEI_USER_DATA_UNREGISTERED:
        this->user_data_unregistered(payloadSize);
        break;
    case SEI_RECOVERY_POINT:
        this->recovery_point(payloadSize);
        break;
    case SEI_DEC_REF_PIC_MARKING_REPETITION:
        this->dec_ref_pic_marking_repetition(payloadSize);
        break;
    case SEI_SPARE_PIC:
        this->spare_pic(payloadSize);
        break;
    case SEI_SCENE_INFO:
        this->scene_info(payloadSize);
        break;
    case SEI_SUB_SEQ_INFO:
        this->sub_seq_info(payloadSize);
        break;
    case SEI_SUB_SEQ_LAYER_CHARACTERISTICS:
        this->sub_seq_layer_characteristics(payloadSize);
        break;
    case SEI_SUB_SEQ_CHARACTERISTICS:
        this->sub_seq_characteristics(payloadSize);
        break;
    case SEI_FULL_FRAME_FREEZE:
        this->full_frame_freeze(payloadSize);
        break;
    case SEI_FULL_FRAME_FREEZE_RELEASE:
        this->full_frame_freeze_release(payloadSize);
        break;
    case SEI_FULL_FRAME_SNAPSHOT:
        this->full_frame_snapshot(payloadSize);
        break;
    case SEI_PROGRESSIVE_REFINEMENT_SEGMENT_START:
        this->progressive_refinement_segment_start(payloadSize);
        break;
    case SEI_PROGRESSIVE_REFINEMENT_SEGMENT_END:
        this->progressive_refinement_segment_end(payloadSize);
        break;
    case SEI_MOTION_CONSTRAINED_SLICE_GROUP_SET:
        this->motion_constrained_slice_group_set(payloadSize);
        break;
    case SEI_FILM_GRAIN_CHARACTERISTICS:
        this->film_grain_characteristics(payloadSize);
        break;
    case SEI_DEBLOCKING_FILTER_DISPLAY_PREFERENCE:
        this->deblocking_filter_display_preference(payloadSize);
        break;
    case SEI_STEREO_VIDEO_INFO:
        this->stereo_video_info(payloadSize);
        break;
    case SEI_POST_FILTER_HINTS:
        this->post_filter_hint(payloadSize);
        break;
    case SEI_TONE_MAPPING_INFO:
        this->tone_mapping_info(payloadSize);
        break;
    case SEI_SCALABILITY_INFO:
    case SEI_SUB_PIC_SCALABLE_LAYER:
    case SEI_NON_REQUIRED_LAYER_REP:
    case SEI_PRIORITY_LAYER_INFO:
    case SEI_LAYERS_NOT_PRESENT:
    case SEI_LAYER_DEPENDENCY_CHANGE:
    case SEI_SCALABLE_NESTING:
    case SEI_BASE_LAYER_TEMPORAL_HRD:
    case SEI_QUALITY_LAYER_INTEGRITY_CHECK:
    case SEI_REDUNDANT_PIC_PROPERTY:
    case SEI_TL0_DEP_REP_INDEX:
    case SEI_TL_SWITCHING_POINT:
    case SEI_PARALLEL_DECODING_INFO:
    case SEI_MVC_SCALABLE_NESTING:
    case SEI_VIEW_SCALABILITY_INFO:
    case SEI_MULTIVIEW_SCENE_INFO:
    case SEI_MULTIVIEW_ACQUISITION_INFO:
    case SEI_NON_REQUIRED_VIEW_COMPONENT:
    case SEI_VIEW_DEPENDENCY_CHANGE:
    case SEI_OPERATION_POINTS_NOT_PRESENT:
    case SEI_BASE_VIEW_TEMPORAL_HRD:
        this->reserved_sei_message(payloadSize);
        break;
    case SEI_FRAME_PACKING_ARRANGEMENT:
        this->frame_packing_arrangement(payloadSize);
        break;
    case SEI_MULTIVIEW_VIEW_POSITION:
        this->reserved_sei_message(payloadSize);
        break;
    case SEI_DISPLAY_ORIENTATION:
        this->display_orientation(payloadSize);
        break;
    case SEI_MVCD_SCALABLE_NESTING:
    case SEI_MVCD_VIEW_SCALABILITY_INFO:
    case SEI_DEPTH_REPRESENTATION_INFO:
    case SEI_THREE_DIMENSIONAL_REFERENCE_DISPLAYS_INFO:
    case SEI_DEPTH_TIMING:
    case SEI_DEPTH_SAMPLING_INFO:
    default:
        this->reserved_sei_message(payloadSize);
        break;
    }

    bool bit_equal_to_one;
    bool bit_equal_to_zero;

    if (!this->byte_aligned()) {
        bit_equal_to_one = this->f(1);
        while (!this->byte_aligned())
            bit_equal_to_zero = this->f(1);
    }
}

void InterpreterSEI::buffering_period(uint32_t payloadSize)
{
    uint8_t seq_parameter_set_id = this->ue("SEI: seq_parameter_set_id");

    sps_t& sps = p_Vid->SeqParSet[seq_parameter_set_id];
    vui_t& vui = sps.vui_parameters;
    hrd_t& nal = vui.nal_hrd_parameters;
    hrd_t& vcl = vui.vcl_hrd_parameters;

    if (sps.vui_parameters_present_flag) {
        uint8_t initial_cpb_removal_delay_length;
        uint8_t initial_cpb_removal_delay, initial_cpb_removal_delay_offset;
        if (vui.nal_hrd_parameters_present_flag) {
            initial_cpb_removal_delay_length = nal.initial_cpb_removal_delay_length_minus1 + 1;
            for (int k = 0; k <= nal.cpb_cnt_minus1; ++k) {
                initial_cpb_removal_delay        = this->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay");
                initial_cpb_removal_delay_offset = this->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay_offset");
            }
        }
        if (vui.vcl_hrd_parameters_present_flag) {
            initial_cpb_removal_delay_length = vcl.initial_cpb_removal_delay_length_minus1 + 1;
            for (int k = 0; k <= sps.vui_parameters.vcl_hrd_parameters.cpb_cnt_minus1; ++k) {
                initial_cpb_removal_delay        = this->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay");
                initial_cpb_removal_delay_offset = this->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay_offset");
            }
        }
    }
}

void InterpreterSEI::pic_timing(uint32_t payloadSize)
{
    if (!p_Vid->active_sps) {
        fprintf(stderr, "Warning: no active SPS, timing SEI cannot be parsed\n");
        return;
    }

    sps_t& sps = *p_Vid->active_sps;
    vui_t& vui = sps.vui_parameters;
    hrd_t& nal = vui.nal_hrd_parameters;
    hrd_t& vcl = vui.vcl_hrd_parameters;

    bool CpbDpbDelaysPresentFlag = sps.vui_parameters_present_flag &&
        (vui.nal_hrd_parameters_present_flag || vui.vcl_hrd_parameters_present_flag);
    if (CpbDpbDelaysPresentFlag) {
        uint8_t cpb_removal_len = 24;
        uint8_t dpb_output_len  = 24;
        if (sps.vui_parameters_present_flag) {
            if (vui.nal_hrd_parameters_present_flag) {
                cpb_removal_len = nal.cpb_removal_delay_length_minus1 + 1;
                dpb_output_len  = nal.dpb_output_delay_length_minus1  + 1;
            } else if (vui.vcl_hrd_parameters_present_flag) {
                cpb_removal_len = vcl.cpb_removal_delay_length_minus1 + 1;
                dpb_output_len  = vcl.dpb_output_delay_length_minus1  + 1;
            }
        }

        uint8_t cpb_removal_delay, dpb_output_delay;
        if (vui.nal_hrd_parameters_present_flag || vui.vcl_hrd_parameters_present_flag) {
            cpb_removal_delay = this->u(cpb_removal_len, "SEI: cpb_removal_delay");
            dpb_output_delay  = this->u(dpb_output_len,  "SEI: dpb_output_delay");
        }
    }

    bool pic_struct_present_flag = sps.vui_parameters_present_flag ?
                                   vui.pic_struct_present_flag : 0;
    if (pic_struct_present_flag) {
        static const uint8_t NumClockTS[16] = {1, 1, 1, 2, 2, 3, 3, 2, 3, 0, };
        uint8_t pic_struct = this->u(4, "SEI: pic_struct");

        for (int i = 0; i < NumClockTS[pic_struct]; ++i) {
            bool clock_timestamp_flag = this->u(1, "SEI: clock_timestamp_flag");
            if (clock_timestamp_flag) {
                uint8_t ct_type, counting_type, n_frames;
                bool    nuit_field_based_flag, full_timestamp_flag;
                bool    discontinuity_flag, cnt_dropped_flag;
                ct_type               = this->u(2, "SEI: ct_type");
                nuit_field_based_flag = this->u(1, "SEI: nuit_field_based_flag");
                counting_type         = this->u(5, "SEI: counting_type");
                full_timestamp_flag   = this->u(1, "SEI: full_timestamp_flag");
                discontinuity_flag    = this->u(1, "SEI: discontinuity_flag");
                cnt_dropped_flag      = this->u(1, "SEI: cnt_dropped_flag");
                n_frames              = this->u(8, "SEI: nframes");

                uint8_t seconds_value, minutes_value, hours_value;
                if (full_timestamp_flag) {
                    seconds_value = this->u(6, "SEI: seconds_value");
                    minutes_value = this->u(6, "SEI: minutes_value");
                    hours_value   = this->u(5, "SEI: hours_value");
                } else {
                    bool seconds_flag = this->u(1, "SEI: seconds_flag");
                    if (seconds_flag) {
                        seconds_value = this->u(6, "SEI: seconds_value");
                        bool minutes_flag = this->u(1, "SEI: minutes_flag");
                        if (minutes_flag) {
                            minutes_value = this->u(6, "SEI: minutes_value");
                            bool hours_flag = this->u(1, "SEI: hours_flag");
                            if (hours_flag)
                                hours_value = this->u(5, "SEI: hours_value");
                        }
                    }
                }

                uint8_t time_offset_length;
                int32_t time_offset;
                time_offset_length =
                    vui.vcl_hrd_parameters_present_flag ? vcl.time_offset_length :
                    vui.nal_hrd_parameters_present_flag ? nal.time_offset_length : 24;
                time_offset = time_offset_length > 0 ?
                    this->i(time_offset_length, "SEI: time_offset") : 0;
            }
        }
    }
}

void InterpreterSEI::pan_scan_rect(uint32_t payloadSize)
{
    uint32_t pan_scan_rect_id;
    pan_scan_rect_id = this->ue("SEI: pan_scan_rect_id");

    bool pan_scan_rect_cancel_flag = this->u(1, "SEI: pan_scan_rect_cancel_flag");
    if (!pan_scan_rect_cancel_flag) {
        int32_t  pan_scan_rect_left_offset, pan_scan_rect_right_offset;
        int32_t  pan_scan_rect_top_offset, pan_scan_rect_bottom_offset;
        uint32_t pan_scan_rect_repetition_period;
        uint32_t pan_scan_cnt_minus1 = this->ue("SEI: pan_scan_cnt_minus1");
        for (int i = 0; i <= pan_scan_cnt_minus1; ++i) {
            pan_scan_rect_left_offset   = this->se("SEI: pan_scan_rect_left_offset");
            pan_scan_rect_right_offset  = this->se("SEI: pan_scan_rect_right_offset");
            pan_scan_rect_top_offset    = this->se("SEI: pan_scan_rect_top_offset");
            pan_scan_rect_bottom_offset = this->se("SEI: pan_scan_rect_bottom_offset");
        }
        pan_scan_rect_repetition_period = this->ue("SEI: pan_scan_rect_repetition_period");
    }
}

void InterpreterSEI::filler_payload(uint32_t payloadSize)
{
    uint8_t ff_byte;
    for (int k = 0; k < payloadSize; ++k)
        ff_byte = this->f(8); // equal to 0xFF
}

void InterpreterSEI::user_data_registered_itu_t_t35(uint32_t payloadSize)
{
    uint8_t itu_t_t35_country_code;
    uint8_t itu_t_t35_country_code_extension_byte;
    uint8_t itu_t_t35_payload_byte;
    int i = 0;

    itu_t_t35_country_code = this->b(8);
    if (itu_t_t35_country_code == 0xFF)
        i = 1;
    else {
        itu_t_t35_country_code_extension_byte = this->b(8);
        i = 2;
    }

    do {
        itu_t_t35_payload_byte = this->b(8);
        i++;
    } while (i < payloadSize);
}

void InterpreterSEI::user_data_unregistered(uint32_t payloadSize)
{
    uint8_t uuid_iso_iec_11578[16];
    for (int i = 0; i < 16; ++i)
        uuid_iso_iec_11578[i] = this->u(8);

    uint8_t user_data_payload_byte;
    for (int i = 16; i < payloadSize; ++i)
        user_data_payload_byte = this->b(8);
}

void InterpreterSEI::recovery_point(uint32_t payloadSize)
{
    uint32_t recovery_frame_cnt;
    bool     exact_match_flag, broken_link_flag;
    uint8_t  changing_slice_group_idc;

    recovery_frame_cnt       = this->ue("SEI: recovery_frame_cnt");
    exact_match_flag         = this->u(1, "SEI: exact_match_flag");
    broken_link_flag         = this->u(1, "SEI: broken_link_flag");
    changing_slice_group_idc = this->u(2, "SEI: changing_slice_group_idc");

    p_Vid->recovery_point     = true;
    p_Vid->recovery_frame_cnt = recovery_frame_cnt;
}

void InterpreterSEI::dec_ref_pic_marking_repetition(uint32_t payloadSize)
{
    bool     original_idr_flag;
    uint32_t original_frame_num;

    original_idr_flag  = this->u(1, "SEI: original_idr_flag");
    original_frame_num = this->ue("SEI: original_frame_num");

    if (!p_Vid->active_sps->frame_mbs_only_flag) {
        bool original_field_pic_flag = this->u(1, "SEI: original_field_pic_flag");
        bool original_bottom_field_flag;
        if (original_field_pic_flag)
            original_bottom_field_flag = this->u(1, "SEI: original_bottom_field_flag");
    }

    shr_t& shr = slice->header;

    // we need to save everything that is probably overwritten in dec_ref_pic_marking()
    bool old_idr_flag                           = slice->IdrPicFlag;
    bool old_no_output_of_prior_pics_flag       = shr.no_output_of_prior_pics_flag;
    bool old_long_term_reference_flag           = shr.long_term_reference_flag;
    bool old_adaptive_ref_pic_marking_mode_flag = shr.adaptive_ref_pic_marking_mode_flag;
    auto old_adaptive_ref_pic_markings          = shr.adaptive_ref_pic_markings;

    // set new initial values
    slice->IdrPicFlag = original_idr_flag;

    this->dec_ref_pic_marking(*slice);

    // restore old values in p_Vid
    slice->IdrPicFlag                      = old_idr_flag;
    shr.no_output_of_prior_pics_flag       = old_no_output_of_prior_pics_flag;
    shr.long_term_reference_flag           = old_long_term_reference_flag;
    shr.adaptive_ref_pic_marking_mode_flag = old_adaptive_ref_pic_marking_mode_flag;
    shr.adaptive_ref_pic_markings          = old_adaptive_ref_pic_markings;

    p_Vid->no_output_of_prior_pics_flag = shr.no_output_of_prior_pics_flag;
}

void InterpreterSEI::spare_pic(uint32_t payloadSize)
{
    sps_t& sps = *p_Vid->active_sps;

    uint32_t target_frame_num;
    bool     spare_field_flag;
    bool     target_bottom_field_flag;
    target_frame_num         = this->ue("SEI: target_frame_num");
    spare_field_flag         = this->u(1, "SEI: spare_field_flag");
    target_bottom_field_flag = spare_field_flag ? this->u(1, "SEI: target_bottom_field_flag") : 0;

    uint32_t num_spare_pics_minus1 = this->ue("SEI: num_spare_pics_minus1");
    uint8_t** map = new uint8_t*[num_spare_pics_minus1 + 1];

    for (int i = 0; i <= num_spare_pics_minus1; ++i) {
        map[i] = new uint8_t[sps.PicSizeInMapUnits];

        uint32_t delta_spare_frame_num;
        bool     spare_bottom_field_flag;
        delta_spare_frame_num = this->ue("SEI: delta_spare_frame_num");
        spare_bottom_field_flag = spare_field_flag ? this->u(1, "SEI: spare_bottom_field_flag") : 0;

        uint32_t spare_area_idc = this->ue("SEI: ref_area_indicator");
        if (spare_area_idc == 0) {
            for (int j = 0; j < sps.PicSizeInMapUnits; ++j)
                map[i][j] = 0;
        } else if (spare_area_idc == 1) {
            for (int j = 0; j < sps.PicSizeInMapUnits; ++j)
                map[i][j] = this->u(1, "SEI: spare_unit_flag");
        } else if (spare_area_idc == 2) {
            int bit0 = 0;
            int bit1 = 1;
            int no_bit0 = -1;

            int x = (sps.PicWidthInMbs    - 1) / 2;
            int y = (sps.FrameHeightInMbs - 1) / 2;
            int left = x, right = x;
            int top = y, bottom = y;
            int directx = 0;
            int directy = 1;

            for (int m = 0; m < sps.FrameHeightInMbs; ++m) {
                for (int n = 0; n < sps.PicWidthInMbs; ++n) {
                    if (no_bit0 < 0)
                        no_bit0 = this->ue("SEI: zero_run_length");
                    if (no_bit0 > 0)
                        map[i][y * sps.PicWidthInMbs + x] = (uint8_t)bit0;
                    else
                        map[i][y * sps.PicWidthInMbs + x] = (uint8_t)bit1;
                    no_bit0--;

                    // go to the next mb:
                    if (directx == -1 && directy == 0) {
                        if (x > left)
                            x--;
                        else if (x == 0) {
                            y = bottom + 1;
                            bottom++;
                            directx = 1;
                            directy = 0;
                        } else if (x == left) {
                            x--;
                            left--;
                            directx = 0;
                            directy = 1;
                        }
                    } else if (directx == 1 && directy == 0) {
                        if (x < right)
                            x++;
                        else if (x == sps.PicWidthInMbs - 1) {
                            y = top - 1;
                            top--;
                            directx = -1;
                            directy = 0;
                        } else if (x == right) {
                            x++;
                            right++;
                            directx = 0;
                            directy = -1;
                        }
                    } else if (directx == 0 && directy == -1) {
                        if (y > top)
                            y--;
                        else if (y == 0) {
                            x = left - 1;
                            left--;
                            directx = 0;
                            directy = 1;
                        } else if (y == top) {
                            y--;
                            top--;
                            directx = -1;
                            directy = 0;
                        }
                    } else if (directx == 0 && directy == 1) {
                        if (y < bottom)
                            y++;
                        else if (y == sps.FrameHeightInMbs - 1) {
                            x = right + 1;
                            right++;
                            directx = 0;
                            directy = -1;
                        } else if (y == bottom) {
                            y++;
                            bottom++;
                            directx = 1;
                            directy = 0;
                        }
                    }
                }
            }
        }

        delete []map[i];
    }

    delete []map;
}

void InterpreterSEI::scene_info(uint32_t payloadSize)
{
    bool scene_info_present_flag = this->u(1, "SEI: scene_info_present_flag");
    if (scene_info_present_flag) {
        uint32_t scene_id;
        uint32_t scene_transition_type;
        uint32_t second_scene_id;
        scene_id              = this->ue("SEI: scene_id");
        scene_transition_type = this->ue("SEI: scene_transition_type");
        if (scene_transition_type > 3)
            second_scene_id = this->ue("SEI: second_scene_id");
    }
}

void InterpreterSEI::sub_seq_info(uint32_t payloadSize)
{
    uint32_t sub_seq_layer_num, sub_seq_id;
    bool     first_ref_pic_flag, leading_non_ref_pic_flag;
    bool     last_pic_flag, sub_seq_frame_num_flag;
    uint32_t sub_seq_frame_num;

    sub_seq_layer_num        = this->ue("SEI: sub_seq_layer_num");
    sub_seq_id               = this->ue("SEI: sub_seq_id");
    first_ref_pic_flag       = this->u(1, "SEI: first_ref_pic_flag");
    leading_non_ref_pic_flag = this->u(1, "SEI: leading_non_ref_pic_flag");
    last_pic_flag            = this->u(1, "SEI: last_pic_flag");
    sub_seq_frame_num_flag   = this->u(1, "SEI: sub_seq_frame_num_flag");
    if (sub_seq_frame_num_flag)
        sub_seq_frame_num    = this->ue("SEI: sub_seq_frame_num");
}

void InterpreterSEI::sub_seq_layer_characteristics(uint32_t payloadSize)
{
    bool     accurate_statistics_flag;
    uint16_t average_bit_rate;
    uint16_t average_frame_rate;
    uint32_t num_sub_seq_layers_minus1 = this->ue("SEI: num_sub_layers_minus1");
    for (int i = 0; i <= num_sub_seq_layers_minus1; ++i) {
        accurate_statistics_flag = this->u(1, "SEI: accurate_statistics_flag");
        average_bit_rate         = this->u(16, "SEI: average_bit_rate");
        average_frame_rate       = this->u(16, "SEI: average_frame_rate");
    }
}

void InterpreterSEI::sub_seq_characteristics(uint32_t payloadSize)
{
    uint32_t sub_seq_layer_num;
    uint32_t sub_seq_id;
    bool     duration_flag;
    uint32_t sub_seq_duration;
    sub_seq_layer_num = this->ue("SEI: sub_seq_layer_num");
    sub_seq_id        = this->ue("SEI: sub_seq_id");
    duration_flag     = this->u(1, "SEI: duration_flag");
    if (duration_flag)
        sub_seq_duration = this->u(32, "SEI: duration_flag");

    bool     accurate_statistics_flag;
    uint16_t average_bit_rate;
    uint16_t average_frame_rate;
    bool average_rate_flag = this->u(1, "SEI: average_rate_flag");
    if (average_rate_flag) {
        accurate_statistics_flag = this->u(1, "SEI: accurate_statistics_flag");
        average_bit_rate         = this->u(16, "SEI: average_bit_rate");
        average_frame_rate       = this->u(16, "SEI: average_frame_rate");
    }

    uint32_t ref_sub_seq_layer_num;
    uint32_t ref_sub_seq_id;
    bool     ref_sub_seq_direction;
    uint32_t num_referenced_subseqs = this->ue("SEI: num_referenced_subseqs");
    for (int i = 0; i < num_referenced_subseqs; ++i) {
        ref_sub_seq_layer_num  = this->ue("SEI: ref_sub_seq_layer_num");
        ref_sub_seq_id         = this->ue("SEI: ref_sub_seq_id");
        ref_sub_seq_direction  = this->u(1, "SEI: ref_sub_seq_direction");
    }
}

void InterpreterSEI::full_frame_freeze(uint32_t payloadSize)
{
    uint32_t full_frame_freeze_repetition_period;
    full_frame_freeze_repetition_period = this->ue("SEI: full_frame_freeze_repetition_period");
}

void InterpreterSEI::full_frame_freeze_release(uint32_t payloadSize)
{
}

void InterpreterSEI::full_frame_snapshot(uint32_t payloadSize)
{
    uint32_t snapshot_id;
    snapshot_id = this->ue("SEI: snapshot_id");
}

void InterpreterSEI::progressive_refinement_segment_start(uint32_t payloadSize)
{
    uint32_t progressive_refinement_id;
    uint32_t num_refinement_steps_minus1;
    progressive_refinement_id   = this->ue("SEI: progressive_refinement_id");
    num_refinement_steps_minus1 = this->ue("SEI: num_refinement_steps_minus1");
}

void InterpreterSEI::progressive_refinement_segment_end(uint32_t payloadSize)
{
    uint32_t progressive_refinement_id;
    progressive_refinement_id = this->ue("SEI: progressive_refinement_id");
}

void InterpreterSEI::motion_constrained_slice_group_set(uint32_t payloadSize)
{
    uint32_t num_slice_groups_minus1 = this->ue("SEI: num_slice_groups_minus1");
    if (num_slice_groups_minus1 > 0) {
        int      sliceGroupSize = ceil(log2(num_slice_groups_minus1 + 1));
        uint32_t slice_group_id;
        for (int i = 0; i <= num_slice_groups_minus1; ++i)
            slice_group_id = this->u(sliceGroupSize, "SEI: slice_group_id");
    }

    bool     exact_match_flag, pan_scan_rect_flag;
    uint32_t pan_scan_rect_id;
    exact_match_flag   = this->u(1, "SEI: exact_match_flag");
    pan_scan_rect_flag = this->u(1, "SEI: pan_scan_rect_flag");
    if (pan_scan_rect_flag)
        pan_scan_rect_id = this->ue("SEI: pan_scan_rect_id");
}

void InterpreterSEI::film_grain_characteristics(uint32_t payloadSize)
{
    bool film_grain_characteristics_cancel_flag = this->u(1, "SEI: film_grain_characteristics_cancel_flag");
    if (!film_grain_characteristics_cancel_flag) {
        uint8_t  film_grain_model_id;
        bool     separate_colour_description_present_flag;
        uint8_t  film_grain_bit_depth_luma_minus8;
        uint8_t  film_grain_bit_depth_chroma_minus8;
        bool     film_grain_full_range_flag;
        uint8_t  film_grain_colour_primaries;
        uint8_t  film_grain_transfer_characteristics;
        uint8_t  film_grain_matrix_coefficients;
        uint8_t  blending_mode_id;
        uint8_t  log2_scale_factor;
        bool     comp_model_present_flag[3];
        uint8_t  num_intensity_intervals_minus1;
        uint8_t  num_model_values_minus1;
        uint8_t  intensity_interval_lower_bound;
        uint8_t  intensity_interval_upper_bound;
        int32_t  comp_model_value;
        uint32_t film_grain_characteristics_repetition_period;

        film_grain_model_id                      = this->u(2, "SEI: model_id");
        separate_colour_description_present_flag = this->u(1, "SEI: separate_colour_description_present_flag");
        if (separate_colour_description_present_flag) {
            film_grain_bit_depth_luma_minus8     = this->u(3, "SEI: film_grain_bit_depth_luma_minus8");
            film_grain_bit_depth_chroma_minus8   = this->u(3, "SEI: film_grain_bit_depth_chroma_minus8");
            film_grain_full_range_flag           = this->u(1, "SEI: film_grain_full_range_flag");
            film_grain_colour_primaries          = this->u(8, "SEI: film_grain_colour_primaries");
            film_grain_transfer_characteristics  = this->u(8, "SEI: film_grain_transfer_characteristics");
            film_grain_matrix_coefficients       = this->u(8, "SEI: film_grain_matrix_coefficients");
        }
        blending_mode_id  = this->u(2, "SEI: blending_mode_id");
        log2_scale_factor = this->u(4, "SEI: log2_scale_factor");
        for (int c = 0; c < 3; ++c)
            comp_model_present_flag[c] = this->u(1, "SEI: comp_model_present_flag");
        for (int c = 0; c < 3; ++c) {
            if (comp_model_present_flag[c]) {
                num_intensity_intervals_minus1 = this->u(8, "SEI: num_intensity_intervals_minus1");
                num_model_values_minus1        = this->u(3, "SEI: num_model_values_minus1");
                for (int i = 0; i <= num_intensity_intervals_minus1; ++i) {
                    intensity_interval_lower_bound = this->u(8, "SEI: intensity_interval_lower_bound");
                    intensity_interval_upper_bound = this->u(8, "SEI: intensity_interval_upper_bound");
                    for (int j = 0; j <= num_model_values_minus1; ++j)
                        comp_model_value = this->se("SEI: comp_model_value");
                }
            }
        }
        film_grain_characteristics_repetition_period = this->ue("SEI: film_grain_characteristics_repetition_period");
    }
}

void InterpreterSEI::deblocking_filter_display_preference(uint32_t payloadSize)
{
    bool     display_prior_to_deblocking_preferred_flag;
    bool     dec_frame_buffering_constraint_flag;
    uint32_t deblocking_display_preference_repetition_period;
    bool deblocking_display_preference_cancel_flag = this->u(1, "SEI: deblocking_display_preference_cancel_flag");
    if (!deblocking_display_preference_cancel_flag) {
        display_prior_to_deblocking_preferred_flag      = this->u(1, "SEI: display_prior_to_deblocking_preferred_flag");
        dec_frame_buffering_constraint_flag             = this->u(1, "SEI: dec_frame_buffering_constraint_flag");
        deblocking_display_preference_repetition_period = this->ue("SEI: deblocking_display_preference_repetition_period");
    }
}

void InterpreterSEI::stereo_video_info(uint32_t payloadSize)
{
    bool field_views_flags;
    bool top_field_is_left_view_flag;
    bool current_frame_is_left_view_flag;
    bool next_frame_is_second_view_flag;
    bool left_view_self_contained_flag;
    bool right_view_self_contained_flag;

    field_views_flags = this->u(1, "SEI: field_views_flags");
    if (field_views_flags)
        top_field_is_left_view_flag     = this->u(1, "SEI: top_field_is_left_view_flag");
    else {
        current_frame_is_left_view_flag = this->u(1, "SEI: current_frame_is_left_view_flag");
        next_frame_is_second_view_flag  = this->u(1, "SEI: next_frame_is_second_view_flag");
    }

    left_view_self_contained_flag  = this->u(1, "SEI: left_view_self_contained_flag");
    right_view_self_contained_flag = this->u(1, "SEI: right_view_self_contained_flag");
}

void InterpreterSEI::post_filter_hint(uint32_t payloadSize)
{
    uint32_t filter_hint_size_y;
    uint32_t filter_hint_size_x;
    uint8_t  filter_hint_type;
    int32_t  filter_hint;
    bool     additional_extension_flag;

    filter_hint_size_y = this->ue("SEI: filter_hint_size_y");
    filter_hint_size_x = this->ue("SEI: filter_hint_size_x");
    filter_hint_type   = this->u(2, "SEI: filter_hint_type");
    for (int colour_component = 0; colour_component < 3; ++colour_component) {
        for (int cy = 0; cy < filter_hint_size_y; ++cy) {
            for (int cx = 0; cx < filter_hint_size_x; ++cx)
                filter_hint = this->se("SEI: filter_hint");
        }
    }

    additional_extension_flag = this->u(1, "SEI: additional_extension_flag");
}

void InterpreterSEI::tone_mapping_info(uint32_t payloadSize)
{
    tone_mapping_t seiToneMappingTmp {};
    int i = 0, max_coded_num, max_output_num;

    seiToneMappingTmp.tone_map_id          = this->ue("SEI: tone_map_id");
    seiToneMappingTmp.tone_map_cancel_flag = this->u(1, "SEI: tone_map_cancel_flag");

    if (!seiToneMappingTmp.tone_map_cancel_flag)  {
        seiToneMappingTmp.tone_map_repetition_period = this->ue("SEI: tone_map_repetition_period");
        seiToneMappingTmp.coded_data_bit_depth       = this->u(8, "SEI: coded_data_bit_depth");
        seiToneMappingTmp.target_bit_depth           = this->u(8, "SEI: sei_bit_depth");
        seiToneMappingTmp.tone_map_model_id          = this->ue("SEI: model_id");

        max_coded_num  = 1 << seiToneMappingTmp.coded_data_bit_depth;
        max_output_num = 1 << seiToneMappingTmp.target_bit_depth;

        if (seiToneMappingTmp.tone_map_model_id == 0) { // linear mapping with clipping
            seiToneMappingTmp.min_value = this->u(32, "SEI: min_value");
            seiToneMappingTmp.max_value = this->u(32, "SEI: min_value");
        } else if (seiToneMappingTmp.tone_map_model_id == 1) { // sigmoidal mapping
            seiToneMappingTmp.sigmoid_midpoint = this->u(32, "SEI: sigmoid_midpoint");
            seiToneMappingTmp.sigmoid_width    = this->u(32, "SEI: sigmoid_width");
        } else if (seiToneMappingTmp.tone_map_model_id == 2) { // user defined table mapping
            for (int i = 0; i < max_output_num; ++i)
                seiToneMappingTmp.start_of_coded_interval[i] = this->u((((seiToneMappingTmp.coded_data_bit_depth+7)>>3)<<3), "SEI: start_of_coded_interval");
        } else if (seiToneMappingTmp.tone_map_model_id == 3) { // piece-wise linear mapping
            seiToneMappingTmp.num_pivots = this->u(16, "SEI: num_pivots");
            seiToneMappingTmp.coded_pivot_value[0] = 0;
            seiToneMappingTmp.target_pivot_value[0] = 0;
            seiToneMappingTmp.coded_pivot_value[seiToneMappingTmp.num_pivots+1] = max_coded_num-1;
            seiToneMappingTmp.target_pivot_value[seiToneMappingTmp.num_pivots+1] = max_output_num-1;

            for (int i = 1; i <= seiToneMappingTmp.num_pivots; ++i) {
                seiToneMappingTmp.coded_pivot_value[i] = this->u( (((seiToneMappingTmp.coded_data_bit_depth+7)>>3)<<3), "SEI: coded_pivot_value");
                seiToneMappingTmp.target_pivot_value[i] = this->u( (((seiToneMappingTmp.target_bit_depth+7)>>3)<<3), "SEI: sei_pivot_value");
            }
        }

        // Currently, only when the map_id == 0, the tone-mapping is actually applied.
        if (seiToneMappingTmp.tone_map_id == 0) {
            p_Vid->seiToneMapping->seiHasTone_mapping = 1;
            p_Vid->seiToneMapping->tone_map_repetition_period = seiToneMappingTmp.tone_map_repetition_period;
            p_Vid->seiToneMapping->coded_data_bit_depth = seiToneMappingTmp.coded_data_bit_depth;
            p_Vid->seiToneMapping->sei_bit_depth = seiToneMappingTmp.target_bit_depth;
            p_Vid->seiToneMapping->model_id = seiToneMappingTmp.tone_map_model_id;
            p_Vid->seiToneMapping->count = 0;

            // generate the look up table of tone mapping
            switch (seiToneMappingTmp.tone_map_model_id) {
            case 0:            // linear mapping with clipping
                for (int i = 0; i <= seiToneMappingTmp.min_value; ++i)
                    p_Vid->seiToneMapping->lut[i] = 0;

                for (int i = seiToneMappingTmp.min_value + 1; i < seiToneMappingTmp.max_value; ++i)
                    p_Vid->seiToneMapping->lut[i] = (px_t) ((i-seiToneMappingTmp.min_value) * (max_output_num-1)/(seiToneMappingTmp.max_value- seiToneMappingTmp.min_value));

                for (int i = seiToneMappingTmp.max_value; i < max_coded_num; ++i)
                    p_Vid->seiToneMapping->lut[i] = (px_t) (max_output_num - 1);
                break;
            case 1: // sigmoid mapping
                for (int i = 0; i < max_coded_num; ++i) {
                    double tmp = 1.0 + exp( -6*(double)(i-seiToneMappingTmp.sigmoid_midpoint)/seiToneMappingTmp.sigmoid_width);
                    p_Vid->seiToneMapping->lut[i] = (px_t)( (double)(max_output_num-1)/ tmp + 0.5);
                }
                break;
            case 2: // user defined table
                if (0 < max_output_num - 1) {
                    for (int j = 0; j < max_output_num - 1; ++j) {
                        for (int i = seiToneMappingTmp.start_of_coded_interval[j]; i<seiToneMappingTmp.start_of_coded_interval[j+1]; i++) 
                            p_Vid->seiToneMapping->lut[i] = (px_t) j;
                    }
                    p_Vid->seiToneMapping->lut[i] = (px_t) (max_output_num - 1);
                }
                break;
            case 3: // piecewise linear mapping
                for (int j = 0; j < seiToneMappingTmp.num_pivots + 1; ++j) {
                    double slope = (double)(seiToneMappingTmp.target_pivot_value[j+1] - seiToneMappingTmp.target_pivot_value[j])/(seiToneMappingTmp.coded_pivot_value[j+1]-seiToneMappingTmp.coded_pivot_value[j]);
                    for (int i = seiToneMappingTmp.coded_pivot_value[j]; i <= seiToneMappingTmp.coded_pivot_value[j+1]; i++) 
                        p_Vid->seiToneMapping->lut[i] = (px_t) (seiToneMappingTmp.target_pivot_value[j] + (int)(( (i - seiToneMappingTmp.coded_pivot_value[j]) * slope)));
                }
                break;
            default:
                break;
            }
        }
    }
}

void InterpreterSEI::frame_packing_arrangement(uint32_t payloadSize)
{
    uint32_t frame_packing_arrangement_id;
    bool     frame_packing_arrangement_cancel_flag;
    uint8_t  frame_packing_arrangement_type;
    bool     quincunx_sampling_flag;
    uint8_t  content_interpretation_type;
    bool     spatial_flipping_flag;
    bool     frame0_flipped_flag;
    bool     field_views_flag;
    bool     current_frame_is_frame0_flag;
    bool     frame0_self_contained_flag;
    bool     frame1_self_contained_flag;
    uint8_t  frame0_grid_position_x;
    uint8_t  frame0_grid_position_y;
    uint8_t  frame1_grid_position_x;
    uint8_t  frame1_grid_position_y;
    uint8_t  frame_packing_arrangement_reserved_byte;
    uint32_t frame_packing_arrangement_repetition_period;
    bool     frame_packing_arrangement_extension_flag;
  
    frame_packing_arrangement_id          = this->ue( "SEI: frame_packing_arrangement_id");
    frame_packing_arrangement_cancel_flag = this->u(1, "SEI: frame_packing_arrangement_cancel_flag");
    if (!frame_packing_arrangement_cancel_flag) {
        frame_packing_arrangement_type = this->u(7, "SEI: frame_packing_arrangement_type");
        quincunx_sampling_flag         = this->u(1, "SEI: quincunx_sampling_flag");
        content_interpretation_type    = this->u(6, "SEI: content_interpretation_type");
        spatial_flipping_flag          = this->u(1, "SEI: spatial_flipping_flag");
        frame0_flipped_flag            = this->u(1, "SEI: frame0_flipped_flag");
        field_views_flag               = this->u(1, "SEI: field_views_flag");
        current_frame_is_frame0_flag   = this->u(1, "SEI: current_frame_is_frame0_flag");
        frame0_self_contained_flag     = this->u(1, "SEI: frame0_self_contained_flag");
        frame1_self_contained_flag     = this->u(1, "SEI: frame1_self_contained_flag");
        if (!quincunx_sampling_flag && frame_packing_arrangement_type != 5) {
            frame0_grid_position_x = this->u(4, "SEI: frame0_grid_position_x");
            frame0_grid_position_y = this->u(4, "SEI: frame0_grid_position_y");
            frame1_grid_position_x = this->u(4, "SEI: frame1_grid_position_x");
            frame1_grid_position_y = this->u(4, "SEI: frame1_grid_position_y");
        }
        frame_packing_arrangement_reserved_byte     = this->u(8, "SEI: frame_packing_arrangement_reserved_byte");
        frame_packing_arrangement_repetition_period = this->ue("SEI: frame_packing_arrangement_repetition_period");
    }
    frame_packing_arrangement_extension_flag = this->u(1, "SEI: frame_packing_arrangement_extension_flag");
}

void InterpreterSEI::display_orientation(uint32_t payloadSize)
{
    bool display_orientation_cancel_flag = this->u(1, "SEI: display_orientation_cancel_flag");
    if (!display_orientation_cancel_flag) {
        bool     hor_flip, ver_flip;
        uint16_t anticlockwise_rotation;
        uint32_t display_orientation_repetition_period;
        bool     display_orientation_extension_flag;
        hor_flip                              = this->u(1, "SEI: hor_flip");
        ver_flip                              = this->u(1, "SEI: ver_flip");
        anticlockwise_rotation                = this->u(16, "SEI: anticlockwise_rotation");
        display_orientation_repetition_period = this->ue("SEI: display_orientation_repetition_period");
        display_orientation_extension_flag    = this->u(1, "SEI: display_orientation_extension_flag");
    }
}

// D.1.27 Reserved SEI message syntax

void InterpreterSEI::reserved_sei_message(uint32_t payloadSize)
{
    uint8_t reserved_sei_message_payload_byte;
    for (int i = 0; i < payloadSize; ++i)
        reserved_sei_message_payload_byte = this->b(8);
}






// tone map using the look-up-table generated according to SEI tone mapping message
void tone_map (px_t** imgX, px_t* lut, int size_x, int size_y)
{
    for (int i = 0; i < size_y; ++i) {
        for (int j = 0; j < size_x; ++j)
            imgX[i][j] = (px_t)lut[imgX[i][j]];
    }
}

void init_tone_mapping_sei(ToneMappingSEI *seiToneMapping) 
{
    seiToneMapping->seiHasTone_mapping = 0;
    seiToneMapping->count = 0;
}

void update_tone_mapping_sei(ToneMappingSEI *seiToneMapping) 
{
    if (seiToneMapping->tone_map_repetition_period == 0) {
        seiToneMapping->seiHasTone_mapping = 0;
        seiToneMapping->count = 0;
    } else if (seiToneMapping->tone_map_repetition_period > 1) {
        seiToneMapping->count++;
        if (seiToneMapping->count>=seiToneMapping->tone_map_repetition_period) {
            seiToneMapping->seiHasTone_mapping = 0;
            seiToneMapping->count = 0;
        }
    }
}
