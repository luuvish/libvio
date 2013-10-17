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


using vio::h264::data_partition_t;
using vio::h264::vui_t;
using vio::h264::hrd_t;

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


void buffering_period(byte* payload, int size, VideoParameters* p_Vid)
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint8_t seq_parameter_set_id = buf->ue("SEI: seq_parameter_set_id");
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
                initial_cpb_removal_delay        = buf->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay");
                initial_cpb_removal_delay_offset = buf->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay_offset");
            }
        }
        if (vui.vcl_hrd_parameters_present_flag) {
            initial_cpb_removal_delay_length = vcl.initial_cpb_removal_delay_length_minus1 + 1;
            for (int k = 0; k <= sps.vui_parameters.vcl_hrd_parameters.cpb_cnt_minus1; ++k) {
                initial_cpb_removal_delay        = buf->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay");
                initial_cpb_removal_delay_offset = buf->u(initial_cpb_removal_delay_length, "SEI: initial_cpb_removal_delay_offset");
            }
        }
    }

    free (buf);
}

void pic_timing(byte* payload, int size, VideoParameters* p_Vid)
{
    if (!p_Vid->active_sps) {
        fprintf(stderr, "Warning: no active SPS, timing SEI cannot be parsed\n");
        return;
    }

    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

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
            cpb_removal_delay = buf->u(cpb_removal_len, "SEI: cpb_removal_delay");
            dpb_output_delay  = buf->u(dpb_output_len,  "SEI: dpb_output_delay");
        }
    }

    bool pic_struct_present_flag = sps.vui_parameters_present_flag ?
                                   vui.pic_struct_present_flag : 0;
    if (pic_struct_present_flag) {
        static const uint8_t NumClockTS[16] = {1, 1, 1, 2, 2, 3, 3, 2, 3, 0, };
        uint8_t pic_struct = buf->u(4, "SEI: pic_struct");

        for (int i = 0; i < NumClockTS[pic_struct]; ++i) {
            bool clock_timestamp_flag = buf->u(1, "SEI: clock_timestamp_flag");
            if (clock_timestamp_flag) {
                uint8_t ct_type, counting_type, n_frames;
                bool    nuit_field_based_flag, full_timestamp_flag;
                bool    discontinuity_flag, cnt_dropped_flag;
                ct_type               = buf->u(2, "SEI: ct_type");
                nuit_field_based_flag = buf->u(1, "SEI: nuit_field_based_flag");
                counting_type         = buf->u(5, "SEI: counting_type");
                full_timestamp_flag   = buf->u(1, "SEI: full_timestamp_flag");
                discontinuity_flag    = buf->u(1, "SEI: discontinuity_flag");
                cnt_dropped_flag      = buf->u(1, "SEI: cnt_dropped_flag");
                n_frames              = buf->u(8, "SEI: nframes");

                uint8_t seconds_value, minutes_value, hours_value;
                if (full_timestamp_flag) {
                    seconds_value = buf->u(6, "SEI: seconds_value");
                    minutes_value = buf->u(6, "SEI: minutes_value");
                    hours_value   = buf->u(5, "SEI: hours_value");
                } else {
                    bool seconds_flag = buf->u(1, "SEI: seconds_flag");
                    if (seconds_flag) {
                        seconds_value = buf->u(6, "SEI: seconds_value");
                        bool minutes_flag = buf->u(1, "SEI: minutes_flag");
                        if (minutes_flag) {
                            minutes_value = buf->u(6, "SEI: minutes_value");
                            bool hours_flag = buf->u(1, "SEI: hours_flag");
                            if (hours_flag)
                                hours_value = buf->u(5, "SEI: hours_value");
                        }
                    }
                }

                uint8_t time_offset_length;
                int32_t time_offset;
                time_offset_length =
                    vui.vcl_hrd_parameters_present_flag ? vcl.time_offset_length :
                    vui.nal_hrd_parameters_present_flag ? nal.time_offset_length : 24;
                time_offset = time_offset_length > 0 ?
                    buf->i(time_offset_length, "SEI: time_offset") : 0;
            }
        }
    }

    free (buf);
}

void pan_scan_rect(byte* payload, int size, VideoParameters *p_Vid)
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t pan_scan_rect_id;
    pan_scan_rect_id = buf->ue("SEI: pan_scan_rect_id");

    bool pan_scan_rect_cancel_flag = buf->u(1, "SEI: pan_scan_rect_cancel_flag");
    if (!pan_scan_rect_cancel_flag) {
        int32_t  pan_scan_rect_left_offset, pan_scan_rect_right_offset;
        int32_t  pan_scan_rect_top_offset, pan_scan_rect_bottom_offset;
        uint32_t pan_scan_rect_repetition_period;
        uint32_t pan_scan_cnt_minus1 = buf->ue("SEI: pan_scan_cnt_minus1");
        for (int i = 0; i <= pan_scan_cnt_minus1; ++i) {
            pan_scan_rect_left_offset   = buf->se("SEI: pan_scan_rect_left_offset");
            pan_scan_rect_right_offset  = buf->se("SEI: pan_scan_rect_right_offset");
            pan_scan_rect_top_offset    = buf->se("SEI: pan_scan_rect_top_offset");
            pan_scan_rect_bottom_offset = buf->se("SEI: pan_scan_rect_bottom_offset");
        }
        pan_scan_rect_repetition_period = buf->ue("SEI: pan_scan_rect_repetition_period");
    }

    free (buf);
}

void filler_payload( byte* payload, int payloadSize, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = payloadSize;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint8_t ff_byte;
    for (int k = 0; k < payloadSize; ++k)
        ff_byte = buf->f(8); // equal to 0xFF

    free (buf);
}

void user_data_registered_itu_t_t35( byte* payload, int payloadSize, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = payloadSize;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint8_t itu_t_t35_country_code;
    uint8_t itu_t_t35_country_code_extension_byte;
    uint8_t itu_t_t35_payload_byte;
    int i = 0;

    itu_t_t35_country_code = buf->b(8);
    if (itu_t_t35_country_code == 0xFF)
        i = 1;
    else {
        itu_t_t35_country_code_extension_byte = buf->b(8);
        i = 2;
    }

    do {
        itu_t_t35_payload_byte = buf->b(8);
        i++;
    } while (i < payloadSize);

    free (buf);
}

void user_data_unregistered( byte* payload, int payloadSize, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = payloadSize;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint8_t uuid_iso_iec_11578[16];
    for (int i = 0; i < 16; ++i)
        uuid_iso_iec_11578[i] = buf->u(8);

    uint8_t user_data_payload_byte;
    for (int i = 16; i < payloadSize; ++i)
        user_data_payload_byte = buf->b(8);

    free (buf);
}

void recovery_point( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t recovery_frame_cnt;
    bool     exact_match_flag, broken_link_flag;
    uint8_t  changing_slice_group_idc;

    recovery_frame_cnt       = buf->ue("SEI: recovery_frame_cnt");
    exact_match_flag         = buf->u(1, "SEI: exact_match_flag");
    broken_link_flag         = buf->u(1, "SEI: broken_link_flag");
    changing_slice_group_idc = buf->u(2, "SEI: changing_slice_group_idc");

    p_Vid->recovery_point     = true;
    p_Vid->recovery_frame_cnt = recovery_frame_cnt;

    free(buf);
}

void dec_ref_pic_marking_repetition( byte* payload, int size, VideoParameters *p_Vid, slice_t *pSlice )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    bool     original_idr_flag;
    uint32_t original_frame_num;

    original_idr_flag  = buf->u(1, "SEI: original_idr_flag");
    original_frame_num = buf->ue("SEI: original_frame_num");

    if (!p_Vid->active_sps->frame_mbs_only_flag) {
        bool original_field_pic_flag = buf->u(1, "SEI: original_field_pic_flag");
        bool original_bottom_field_flag;
        if (original_field_pic_flag)
            original_bottom_field_flag = buf->u(1, "SEI: original_bottom_field_flag");
    }

    shr_t& shr = pSlice->header;

    // we need to save everything that is probably overwritten in dec_ref_pic_marking()
    drpm_t* old_drpm = shr.dec_ref_pic_marking_buffer;
    bool old_idr_flag                        = pSlice->idr_flag;
    bool old_no_output_of_prior_pics_flag    = shr.no_output_of_prior_pics_flag;
    bool old_long_term_reference_flag        = shr.long_term_reference_flag;
    bool old_adaptive_ref_pic_buffering_flag = shr.adaptive_ref_pic_marking_mode_flag;

    // set new initial values
    pSlice->idr_flag = original_idr_flag;
    shr.dec_ref_pic_marking_buffer = NULL;

    buf->dec_ref_pic_marking(p_Vid, *pSlice);

    while (shr.dec_ref_pic_marking_buffer) {
        drpm_t* tmp_drpm = shr.dec_ref_pic_marking_buffer;
        shr.dec_ref_pic_marking_buffer = tmp_drpm->Next;
        delete tmp_drpm;
    }

    // restore old values in p_Vid
    shr.dec_ref_pic_marking_buffer = old_drpm;
    pSlice->idr_flag = old_idr_flag;
    shr.no_output_of_prior_pics_flag = old_no_output_of_prior_pics_flag;
    p_Vid->no_output_of_prior_pics_flag = shr.no_output_of_prior_pics_flag;
    shr.long_term_reference_flag = old_long_term_reference_flag;
    shr.adaptive_ref_pic_marking_mode_flag = old_adaptive_ref_pic_buffering_flag;

    free (buf);
}

void spare_pic( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    sps_t& sps = *p_Vid->active_sps;

    uint32_t target_frame_num;
    bool     spare_field_flag;
    bool     target_bottom_field_flag;
    target_frame_num         = buf->ue("SEI: target_frame_num");
    spare_field_flag         = buf->u(1, "SEI: spare_field_flag");
    target_bottom_field_flag = spare_field_flag ? buf->u(1, "SEI: target_bottom_field_flag") : 0;

    uint32_t num_spare_pics_minus1    = buf->ue("SEI: num_spare_pics_minus1");
    uint8_t** map = new uint8_t*[num_spare_pics_minus1 + 1];

    for (int i = 0; i <= num_spare_pics_minus1; ++i) {
        map[i] = new uint8_t[sps.PicSizeInMapUnits];

        uint32_t delta_spare_frame_num;
        bool     spare_bottom_field_flag;
        delta_spare_frame_num = buf->ue("SEI: delta_spare_frame_num");
        spare_bottom_field_flag = spare_field_flag ? buf->u(1, "SEI: spare_bottom_field_flag") : 0;

        uint32_t spare_area_idc = buf->ue("SEI: ref_area_indicator");
        if (spare_area_idc == 0) {
            for (int j = 0; j < sps.PicSizeInMapUnits; ++j)
                map[i][j] = 0;
        } else if (spare_area_idc == 1) {
            for (int j = 0; j < sps.PicSizeInMapUnits; ++j)
                map[i][j] = buf->u(1, "SEI: spare_unit_flag");
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
                        no_bit0 = buf->ue("SEI: zero_run_length");
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
    } // end of num_spare_pics

    delete []map;
    free (buf);
}

void scene_info( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    bool scene_info_present_flag = buf->u(1, "SEI: scene_info_present_flag");
    if (scene_info_present_flag) {
        uint32_t scene_id;
        uint32_t scene_transition_type;
        uint32_t second_scene_id;
        scene_id              = buf->ue("SEI: scene_id");
        scene_transition_type = buf->ue("SEI: scene_transition_type");
        if (scene_transition_type > 3)
            second_scene_id = buf->ue("SEI: second_scene_id");
    }

    free (buf);
}

void sub_seq_info( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t sub_seq_layer_num, sub_seq_id;
    bool     first_ref_pic_flag, leading_non_ref_pic_flag;
    bool     last_pic_flag, sub_seq_frame_num_flag;
    uint32_t sub_seq_frame_num;

    sub_seq_layer_num        = buf->ue("SEI: sub_seq_layer_num");
    sub_seq_id               = buf->ue("SEI: sub_seq_id");
    first_ref_pic_flag       = buf->u(1, "SEI: first_ref_pic_flag");
    leading_non_ref_pic_flag = buf->u(1, "SEI: leading_non_ref_pic_flag");
    last_pic_flag            = buf->u(1, "SEI: last_pic_flag");
    sub_seq_frame_num_flag   = buf->u(1, "SEI: sub_seq_frame_num_flag");
    if (sub_seq_frame_num_flag)
        sub_seq_frame_num    = buf->ue("SEI: sub_seq_frame_num");

    free (buf);
}

void sub_seq_layer_characteristics( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    bool     accurate_statistics_flag;
    uint16_t average_bit_rate;
    uint16_t average_frame_rate;
    uint32_t num_sub_seq_layers_minus1 = buf->ue("SEI: num_sub_layers_minus1");
    for (int i = 0; i <= num_sub_seq_layers_minus1; ++i) {
        accurate_statistics_flag = buf->u(1, "SEI: accurate_statistics_flag");
        average_bit_rate         = buf->u(16, "SEI: average_bit_rate");
        average_frame_rate       = buf->u(16, "SEI: average_frame_rate");
    }

    free (buf);
}

void sub_seq_characteristics( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t sub_seq_layer_num;
    uint32_t sub_seq_id;
    bool     duration_flag;
    uint32_t sub_seq_duration;
    sub_seq_layer_num = buf->ue("SEI: sub_seq_layer_num");
    sub_seq_id        = buf->ue("SEI: sub_seq_id");
    duration_flag     = buf->u(1, "SEI: duration_flag");
    if (duration_flag)
        sub_seq_duration = buf->u(32, "SEI: duration_flag");

    bool     accurate_statistics_flag;
    uint16_t average_bit_rate;
    uint16_t average_frame_rate;
    bool average_rate_flag = buf->u(1, "SEI: average_rate_flag");
    if (average_rate_flag) {
        accurate_statistics_flag = buf->u(1, "SEI: accurate_statistics_flag");
        average_bit_rate         = buf->u(16, "SEI: average_bit_rate");
        average_frame_rate       = buf->u(16, "SEI: average_frame_rate");
    }

    uint32_t ref_sub_seq_layer_num;
    uint32_t ref_sub_seq_id;
    bool     ref_sub_seq_direction;
    uint32_t num_referenced_subseqs = buf->ue("SEI: num_referenced_subseqs");
    for (int i = 0; i < num_referenced_subseqs; ++i) {
        ref_sub_seq_layer_num  = buf->ue("SEI: ref_sub_seq_layer_num");
        ref_sub_seq_id         = buf->ue("SEI: ref_sub_seq_id");
        ref_sub_seq_direction  = buf->u(1, "SEI: ref_sub_seq_direction");
    }

    free (buf);
}

void full_frame_freeze( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t full_frame_freeze_repetition_period;
    full_frame_freeze_repetition_period = buf->ue("SEI: full_frame_freeze_repetition_period");

    free (buf);
}

void full_frame_freeze_release( byte* payload, int size, VideoParameters *p_Vid )
{
}

void full_frame_snapshot( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t snapshot_id;
    snapshot_id = buf->ue("SEI: snapshot_id");

    free (buf);
}

void progressive_refinement_segment_start( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t progressive_refinement_id;
    uint32_t num_refinement_steps_minus1;
    progressive_refinement_id   = buf->ue("SEI: progressive_refinement_id");
    num_refinement_steps_minus1 = buf->ue("SEI: num_refinement_steps_minus1");

    free (buf);
}

void progressive_refinement_segment_end( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t progressive_refinement_id;
    progressive_refinement_id = buf->ue("SEI: progressive_refinement_id");

    free (buf);
}

void motion_constrained_slice_group_set( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t num_slice_groups_minus1 = buf->ue("SEI: num_slice_groups_minus1");
    if (num_slice_groups_minus1 > 0) {
        int      sliceGroupSize = ceil(log2(num_slice_groups_minus1 + 1));
        uint32_t slice_group_id;
        for (int i = 0; i <= num_slice_groups_minus1; ++i)
            slice_group_id = buf->u(sliceGroupSize, "SEI: slice_group_id");
    }

    bool     exact_match_flag, pan_scan_rect_flag;
    uint32_t pan_scan_rect_id;
    exact_match_flag   = buf->u(1, "SEI: exact_match_flag");
    pan_scan_rect_flag = buf->u(1, "SEI: pan_scan_rect_flag");
    if (pan_scan_rect_flag)
        pan_scan_rect_id = buf->ue("SEI: pan_scan_rect_id");

    free (buf);
}

void film_grain_characteristics( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    bool film_grain_characteristics_cancel_flag = buf->u(1, "SEI: film_grain_characteristics_cancel_flag");
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

        film_grain_model_id                      = buf->u(2, "SEI: model_id");
        separate_colour_description_present_flag = buf->u(1, "SEI: separate_colour_description_present_flag");
        if (separate_colour_description_present_flag) {
            film_grain_bit_depth_luma_minus8     = buf->u(3, "SEI: film_grain_bit_depth_luma_minus8");
            film_grain_bit_depth_chroma_minus8   = buf->u(3, "SEI: film_grain_bit_depth_chroma_minus8");
            film_grain_full_range_flag           = buf->u(1, "SEI: film_grain_full_range_flag");
            film_grain_colour_primaries          = buf->u(8, "SEI: film_grain_colour_primaries");
            film_grain_transfer_characteristics  = buf->u(8, "SEI: film_grain_transfer_characteristics");
            film_grain_matrix_coefficients       = buf->u(8, "SEI: film_grain_matrix_coefficients");
        }
        blending_mode_id  = buf->u(2, "SEI: blending_mode_id");
        log2_scale_factor = buf->u(4, "SEI: log2_scale_factor");
        for (int c = 0; c < 3; ++c)
            comp_model_present_flag[c] = buf->u(1, "SEI: comp_model_present_flag");
        for (int c = 0; c < 3; ++c) {
            if (comp_model_present_flag[c]) {
                num_intensity_intervals_minus1 = buf->u(8, "SEI: num_intensity_intervals_minus1");
                num_model_values_minus1        = buf->u(3, "SEI: num_model_values_minus1");
                for (int i = 0; i <= num_intensity_intervals_minus1; ++i) {
                    intensity_interval_lower_bound = buf->u(8, "SEI: intensity_interval_lower_bound");
                    intensity_interval_upper_bound = buf->u(8, "SEI: intensity_interval_upper_bound");
                    for (int j = 0; j <= num_model_values_minus1; ++j)
                        comp_model_value = buf->se("SEI: comp_model_value");
                }
            }
        }
        film_grain_characteristics_repetition_period = buf->ue("SEI: film_grain_characteristics_repetition_period");
    }

    free (buf);
}

void deblocking_filter_display_preference( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    bool     display_prior_to_deblocking_preferred_flag;
    bool     dec_frame_buffering_constraint_flag;
    uint32_t deblocking_display_preference_repetition_period;
    bool deblocking_display_preference_cancel_flag = buf->u(1, "SEI: deblocking_display_preference_cancel_flag");
    if (!deblocking_display_preference_cancel_flag) {
        display_prior_to_deblocking_preferred_flag      = buf->u(1, "SEI: display_prior_to_deblocking_preferred_flag");
        dec_frame_buffering_constraint_flag             = buf->u(1, "SEI: dec_frame_buffering_constraint_flag");
        deblocking_display_preference_repetition_period = buf->ue("SEI: deblocking_display_preference_repetition_period");
    }

    free (buf);
}

void stereo_video_info( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    bool field_views_flags;
    bool top_field_is_left_view_flag;
    bool current_frame_is_left_view_flag;
    bool next_frame_is_second_view_flag;
    bool left_view_self_contained_flag;
    bool right_view_self_contained_flag;

    field_views_flags = buf->u(1, "SEI: field_views_flags");
    if (field_views_flags)
        top_field_is_left_view_flag     = buf->u(1, "SEI: top_field_is_left_view_flag");
    else {
        current_frame_is_left_view_flag = buf->u(1, "SEI: current_frame_is_left_view_flag");
        next_frame_is_second_view_flag  = buf->u(1, "SEI: next_frame_is_second_view_flag");
    }

    left_view_self_contained_flag  = buf->u(1, "SEI: left_view_self_contained_flag");
    right_view_self_contained_flag = buf->u(1, "SEI: right_view_self_contained_flag");

    free (buf);
}

void post_filter_hints( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    uint32_t filter_hint_size_y;
    uint32_t filter_hint_size_x;
    uint8_t  filter_hint_type;
    int32_t  filter_hint;
    bool     additional_extension_flag;

    filter_hint_size_y = buf->ue("SEI: filter_hint_size_y");
    filter_hint_size_x = buf->ue("SEI: filter_hint_size_x");
    filter_hint_type   = buf->u(2, "SEI: filter_hint_type");
    for (int colour_component = 0; colour_component < 3; ++colour_component) {
        for (int cy = 0; cy < filter_hint_size_y; ++cy) {
            for (int cx = 0; cx < filter_hint_size_x; ++cx)
                filter_hint = buf->se("SEI: filter_hint");
        }
    }

    additional_extension_flag = buf->u(1, "SEI: additional_extension_flag");

    free (buf);
}

void tone_mapping_info( byte* payload, int size, VideoParameters *p_Vid )
{
    tone_mapping_t seiToneMappingTmp {};
    int i = 0, max_coded_num, max_output_num;

    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

    seiToneMappingTmp.tone_map_id          = buf->ue("SEI: tone_map_id");
    seiToneMappingTmp.tone_map_cancel_flag = buf->u(1, "SEI: tone_map_cancel_flag");

    if (!seiToneMappingTmp.tone_map_cancel_flag)  {
        seiToneMappingTmp.tone_map_repetition_period = buf->ue("SEI: tone_map_repetition_period");
        seiToneMappingTmp.coded_data_bit_depth       = buf->u(8, "SEI: coded_data_bit_depth");
        seiToneMappingTmp.target_bit_depth           = buf->u(8, "SEI: sei_bit_depth");
        seiToneMappingTmp.tone_map_model_id          = buf->ue("SEI: model_id");

        max_coded_num  = 1 << seiToneMappingTmp.coded_data_bit_depth;
        max_output_num = 1 << seiToneMappingTmp.target_bit_depth;

        if (seiToneMappingTmp.tone_map_model_id == 0) { // linear mapping with clipping
            seiToneMappingTmp.min_value = buf->u(32, "SEI: min_value");
            seiToneMappingTmp.max_value = buf->u(32, "SEI: min_value");
        } else if (seiToneMappingTmp.tone_map_model_id == 1) { // sigmoidal mapping
            seiToneMappingTmp.sigmoid_midpoint = buf->u(32, "SEI: sigmoid_midpoint");
            seiToneMappingTmp.sigmoid_width    = buf->u(32, "SEI: sigmoid_width");
        } else if (seiToneMappingTmp.tone_map_model_id == 2) { // user defined table mapping
            for (int i = 0; i < max_output_num; ++i)
                seiToneMappingTmp.start_of_coded_interval[i] = buf->u((((seiToneMappingTmp.coded_data_bit_depth+7)>>3)<<3), "SEI: start_of_coded_interval");
        } else if (seiToneMappingTmp.tone_map_model_id == 3) { // piece-wise linear mapping
            seiToneMappingTmp.num_pivots = buf->u(16, "SEI: num_pivots");
            seiToneMappingTmp.coded_pivot_value[0] = 0;
            seiToneMappingTmp.target_pivot_value[0] = 0;
            seiToneMappingTmp.coded_pivot_value[seiToneMappingTmp.num_pivots+1] = max_coded_num-1;
            seiToneMappingTmp.target_pivot_value[seiToneMappingTmp.num_pivots+1] = max_output_num-1;

            for (int i = 1; i <= seiToneMappingTmp.num_pivots; ++i) {
                seiToneMappingTmp.coded_pivot_value[i] = buf->u( (((seiToneMappingTmp.coded_data_bit_depth+7)>>3)<<3), "SEI: coded_pivot_value");
                seiToneMappingTmp.target_pivot_value[i] = buf->u( (((seiToneMappingTmp.target_bit_depth+7)>>3)<<3), "SEI: sei_pivot_value");
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

    free (buf);
}

void frame_packing_arrangement( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;

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
  
    frame_packing_arrangement_id          = buf->ue( "SEI: frame_packing_arrangement_id");
    frame_packing_arrangement_cancel_flag = buf->u(1, "SEI: frame_packing_arrangement_cancel_flag");
    if (!frame_packing_arrangement_cancel_flag) {
        frame_packing_arrangement_type = buf->u(7, "SEI: frame_packing_arrangement_type");
        quincunx_sampling_flag         = buf->u(1, "SEI: quincunx_sampling_flag");
        content_interpretation_type    = buf->u(6, "SEI: content_interpretation_type");
        spatial_flipping_flag          = buf->u(1, "SEI: spatial_flipping_flag");
        frame0_flipped_flag            = buf->u(1, "SEI: frame0_flipped_flag");
        field_views_flag               = buf->u(1, "SEI: field_views_flag");
        current_frame_is_frame0_flag   = buf->u(1, "SEI: current_frame_is_frame0_flag");
        frame0_self_contained_flag     = buf->u(1, "SEI: frame0_self_contained_flag");
        frame1_self_contained_flag     = buf->u(1, "SEI: frame1_self_contained_flag");
        if (!quincunx_sampling_flag && frame_packing_arrangement_type != 5) {
            frame0_grid_position_x = buf->u(4, "SEI: frame0_grid_position_x");
            frame0_grid_position_y = buf->u(4, "SEI: frame0_grid_position_y");
            frame1_grid_position_x = buf->u(4, "SEI: frame1_grid_position_x");
            frame1_grid_position_y = buf->u(4, "SEI: frame1_grid_position_y");
        }
        frame_packing_arrangement_reserved_byte     = buf->u(8, "SEI: frame_packing_arrangement_reserved_byte");
        frame_packing_arrangement_repetition_period = buf->ue("SEI: frame_packing_arrangement_repetition_period");
    }
    frame_packing_arrangement_extension_flag = buf->u(1, "SEI: frame_packing_arrangement_extension_flag");

    free (buf);
}

void display_orientation( byte* payload, int size, VideoParameters *p_Vid )
{
    data_partition_t* buf = new data_partition_t;
    buf->num_bytes_in_rbsp = size;
    buf->rbsp_byte = payload;
    buf->frame_bitoffset = 0;
  
    bool display_orientation_cancel_flag = buf->u(1, "SEI: display_orientation_cancel_flag");
    if (!display_orientation_cancel_flag) {
        bool     hor_flip, ver_flip;
        uint16_t anticlockwise_rotation;
        uint32_t display_orientation_repetition_period;
        bool     display_orientation_extension_flag;
        hor_flip                              = buf->u(1, "SEI: hor_flip");
        ver_flip                              = buf->u(1, "SEI: ver_flip");
        anticlockwise_rotation                = buf->u(16, "SEI: anticlockwise_rotation");
        display_orientation_repetition_period = buf->ue("SEI: display_orientation_repetition_period");
        display_orientation_extension_flag    = buf->u(1, "SEI: display_orientation_extension_flag");
    }

    free (buf);
}

void reserved_sei_message( byte* payload, int size, VideoParameters *p_Vid )
{
    int offset = 0;
    byte payload_byte;

    while (offset < size) {
        payload_byte = payload[offset];
        offset ++;
    }
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


void parse_sei(byte *msg, int size, VideoParameters *p_Vid, slice_t *pSlice)
{
    int payload_type = 0;
    int payload_size = 0;
    int offset = 1;
    byte tmp_byte;

    do {
        // sei_message();
        payload_type = 0;
        tmp_byte = msg[offset++];
        while (tmp_byte == 0xFF) {
            payload_type += 255;
            tmp_byte = msg[offset++];
        }
        payload_type += tmp_byte;   // this is the last byte

        payload_size = 0;
        tmp_byte = msg[offset++];
        while (tmp_byte == 0xFF) {
            payload_size += 255;
            tmp_byte = msg[offset++];
        }
        payload_size += tmp_byte;   // this is the last byte

        switch (payload_type) {    // sei_payload( type, size );
        case SEI_BUFFERING_PERIOD:
            buffering_period( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PIC_TIMING:
            pic_timing( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PAN_SCAN_RECT:
            pan_scan_rect( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FILLER_PAYLOAD:
            filler_payload( msg+offset, payload_size, p_Vid );
            break;
        case SEI_USER_DATA_REGISTERED_ITU_T_T35:
            user_data_registered_itu_t_t35( msg+offset, payload_size, p_Vid );
            break;
        case SEI_USER_DATA_UNREGISTERED:
            user_data_unregistered( msg+offset, payload_size, p_Vid );
            break;
        case SEI_RECOVERY_POINT:
            recovery_point( msg+offset, payload_size, p_Vid );
            break;
        case SEI_DEC_REF_PIC_MARKING_REPETITION:
            dec_ref_pic_marking_repetition( msg+offset, payload_size, p_Vid, pSlice );
            break;
        case SEI_SPARE_PIC:
            spare_pic( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SCENE_INFO:
            scene_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SUB_SEQ_INFO:
            sub_seq_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SUB_SEQ_LAYER_CHARACTERISTICS:
            sub_seq_layer_characteristics( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SUB_SEQ_CHARACTERISTICS:
            sub_seq_characteristics( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FULL_FRAME_FREEZE:
            full_frame_freeze( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FULL_FRAME_FREEZE_RELEASE:
            full_frame_freeze_release( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FULL_FRAME_SNAPSHOT:
            full_frame_snapshot( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PROGRESSIVE_REFINEMENT_SEGMENT_START:
            progressive_refinement_segment_start( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PROGRESSIVE_REFINEMENT_SEGMENT_END:
            progressive_refinement_segment_end( msg+offset, payload_size, p_Vid );
            break;
        case SEI_MOTION_CONSTRAINED_SLICE_GROUP_SET:
            motion_constrained_slice_group_set( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FILM_GRAIN_CHARACTERISTICS:
            film_grain_characteristics( msg+offset, payload_size, p_Vid );
            break;
        case SEI_DEBLOCKING_FILTER_DISPLAY_PREFERENCE:
            deblocking_filter_display_preference( msg+offset, payload_size, p_Vid );
            break;
        case SEI_STEREO_VIDEO_INFO:
            stereo_video_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_POST_FILTER_HINTS:
            post_filter_hints( msg+offset, payload_size, p_Vid );
            break;
        case SEI_TONE_MAPPING_INFO:
            tone_mapping_info( msg+offset, payload_size, p_Vid );
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
            reserved_sei_message( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FRAME_PACKING_ARRANGEMENT:
            frame_packing_arrangement( msg+offset, payload_size, p_Vid );
            break;
        case SEI_MULTIVIEW_VIEW_POSITION:
            reserved_sei_message( msg+offset, payload_size, p_Vid );
            break;
        case SEI_DISPLAY_ORIENTATION:
            display_orientation(msg+offset, payload_size, p_Vid);
            break;
        case SEI_MVCD_SCALABLE_NESTING:
        case SEI_MVCD_VIEW_SCALABILITY_INFO:
        case SEI_DEPTH_REPRESENTATION_INFO:
        case SEI_THREE_DIMENSIONAL_REFERENCE_DISPLAYS_INFO:
        case SEI_DEPTH_TIMING:
        case SEI_DEPTH_SAMPLING_INFO:
        default:
            reserved_sei_message( msg+offset, payload_size, p_Vid );
            break;
        }
        offset += payload_size;

    } while (msg[offset] != 0x80);    // more_rbsp_data()  msg[offset] != 0x80

    // ignore the trailing bits rbsp_trailing_bits();
    assert(msg[offset] == 0x80);      // this is the trailing bits
    assert(offset + 1 == size);
}


// 7.3.2.3 Supplemental enhancement information RBSP syntax

void data_partition_t::sei_rbsp(void)
{
    do {
        this->sei_message();
    } while (this->more_rbsp_data());

    this->rbsp_trailing_bits();
}

// 7.3.2.3.1 Supplemental enhancement information message syntax

void data_partition_t::sei_message()
{

}
