#include "global.h"
#include "data_partition.h"
#include "fmo.h"
#include "dpb.h"
#include "slice.h"
#include "memalloc.h"


slice_backup_t& slice_backup_t::operator=(const slice_t& slice)
{
    sps_t& sps = *slice.active_sps;

    this->pic_parameter_set_id = slice.pic_parameter_set_id;
    this->frame_num            = slice.frame_num;
    this->field_pic_flag       = slice.field_pic_flag;

    if (slice.field_pic_flag)
        this->bottom_field_flag = slice.bottom_field_flag;

    this->nal_ref_idc = slice.nal_ref_idc;
    this->idr_flag    = slice.idr_flag;

    if (slice.idr_flag)
        this->idr_pic_id = slice.idr_pic_id;

    if (sps.pic_order_cnt_type == 0) {
        this->pic_order_cnt_lsb          = slice.pic_order_cnt_lsb;
        this->delta_pic_order_cnt_bottom = slice.delta_pic_order_cnt_bottom;
    }
    if (sps.pic_order_cnt_type == 1) {
        this->delta_pic_order_cnt[0] = slice.delta_pic_order_cnt[0];
        this->delta_pic_order_cnt[1] = slice.delta_pic_order_cnt[1];
    }

#if (MVC_EXTENSION_ENABLE)
    this->view_id         = slice.view_id;
    this->inter_view_flag = slice.inter_view_flag; 
    this->anchor_pic_flag = slice.anchor_pic_flag;
#endif
    this->layer_id        = slice.layer_id;

    return *this;
}

bool slice_backup_t::operator!=(const slice_t& slice)
{
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;

    bool result = false;

    result |= this->pic_parameter_set_id != slice.pic_parameter_set_id;
    result |= this->frame_num            != slice.frame_num;
    result |= this->field_pic_flag       != slice.field_pic_flag;

    if (slice.field_pic_flag && this->field_pic_flag)
        result |= this->bottom_field_flag != slice.bottom_field_flag;

    result |= this->nal_ref_idc != slice.nal_ref_idc && (this->nal_ref_idc == 0 || slice.nal_ref_idc == 0);
    result |= this->idr_flag    != slice.idr_flag;

    if (slice.idr_flag && this->idr_flag)
        result |= this->idr_pic_id != slice.idr_pic_id;

    if (sps.pic_order_cnt_type == 0) {
        result |= this->pic_order_cnt_lsb != slice.pic_order_cnt_lsb;
        if (pps.bottom_field_pic_order_in_frame_present_flag && !slice.field_pic_flag)
            result |= this->delta_pic_order_cnt_bottom != slice.delta_pic_order_cnt_bottom;
    }
    if (sps.pic_order_cnt_type == 1) {
        if (!sps.delta_pic_order_always_zero_flag) {
            result |= this->delta_pic_order_cnt[0] != slice.delta_pic_order_cnt[0];
            if (pps.bottom_field_pic_order_in_frame_present_flag && !slice.field_pic_flag)
                result |= this->delta_pic_order_cnt[1] != slice.delta_pic_order_cnt[1];
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
}

slice_t::~slice_t()
{
    free_mem3Dpel(this->mb_pred);

    while (this->dec_ref_pic_marking_buffer) {
        drpm_t* tmp_drpm = this->dec_ref_pic_marking_buffer;
        this->dec_ref_pic_marking_buffer=tmp_drpm->Next;
        free(tmp_drpm);
    }
}

void slice_header(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    data_partition_t *s = &currSlice->parser.partArr[0];

    currSlice->first_mb_in_slice = s->ue("SH: first_mb_in_slice");
    p_Vid->type = currSlice->slice_type = s->ue("SH: slice_type") % 5;
    currSlice->pic_parameter_set_id = s->ue("SH: pic_parameter_set_id");

    assert(currSlice->pic_parameter_set_id >= 0 && currSlice->pic_parameter_set_id <= 255);

    UseParameterSet(currSlice);
    sps_t *sps = currSlice->active_sps = p_Vid->active_sps;
    pps_t *pps = currSlice->active_pps = p_Vid->active_pps;

    currSlice->colour_plane_id = 0;
    if (sps->separate_colour_plane_flag)
        currSlice->colour_plane_id = s->u(2, "SH: colour_plane_id");

    assert(currSlice->colour_plane_id >= 0 && currSlice->colour_plane_id <= 2);

    currSlice->frame_num = s->u(sps->log2_max_frame_num_minus4 + 4, "SH: frame_num");

    /* Tian Dong: frame_num gap processing, if found */
    if (currSlice->idr_flag) {
        p_Vid->pre_frame_num = currSlice->frame_num;
        // picture error concealment
        p_Vid->last_ref_pic_poc = 0;
        assert(currSlice->frame_num == 0);
    }

    p_Vid->structure = FRAME;
    currSlice->field_pic_flag    = 0;
    currSlice->bottom_field_flag = 0;
    if (!sps->frame_mbs_only_flag) {
        currSlice->field_pic_flag = s->u(1, "SH: field_pic_flag");
        if (currSlice->field_pic_flag) {
            currSlice->bottom_field_flag = s->u(1, "SH: bottom_field_flag");
            p_Vid->structure = currSlice->bottom_field_flag ? BOTTOM_FIELD : TOP_FIELD;
        }
    }

    currSlice->structure = (PictureStructure) p_Vid->structure;

    currSlice->MbaffFrameFlag     = (sps->mb_adaptive_frame_field_flag && !currSlice->field_pic_flag);
    currSlice->PicHeightInMbs     = sps->FrameHeightInMbs / (1 + currSlice->field_pic_flag);
    currSlice->PicHeightInSampleL = currSlice->PicHeightInMbs * 16;
    currSlice->PicHeightInSampleC = currSlice->PicHeightInMbs * sps->MbHeightC;
    currSlice->PicSizeInMbs       = sps->PicWidthInMbs * currSlice->PicHeightInMbs;
    if (!currSlice->field_pic_flag) {
        currSlice->MaxPicNum  = sps->MaxFrameNum;
        currSlice->CurrPicNum = currSlice->frame_num;
    } else {
        currSlice->MaxPicNum  = 2 * sps->MaxFrameNum;
        currSlice->CurrPicNum = 2 * currSlice->frame_num + 1;
    }

    if (currSlice->idr_flag)
        currSlice->idr_pic_id = s->ue("SH: idr_pic_id");
#if (MVC_EXTENSION_ENABLE)
    else if ( currSlice->svc_extension_flag == 0 && currSlice->NaluHeaderMVCExt.non_idr_flag == 0 )
        currSlice->idr_pic_id = s->ue("SH: idr_pic_id");
#endif

    assert(currSlice->idr_pic_id >= 0 && currSlice->idr_pic_id <= 65535);

    currSlice->delta_pic_order_cnt_bottom = 0;
    currSlice->delta_pic_order_cnt[0]     = 0;
    currSlice->delta_pic_order_cnt[1]     = 0;
    if (sps->pic_order_cnt_type == 0) {
        currSlice->pic_order_cnt_lsb = s->u(sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "SH: pic_order_cnt_lsb");
        if (pps->bottom_field_pic_order_in_frame_present_flag && !currSlice->field_pic_flag)
            currSlice->delta_pic_order_cnt_bottom = s->se("SH: delta_pic_order_cnt_bottom");
    }
    if (sps->pic_order_cnt_type == 1 && !sps->delta_pic_order_always_zero_flag) {
        currSlice->delta_pic_order_cnt[0] = s->se("SH: delta_pic_order_cnt[0]");
        if (pps->bottom_field_pic_order_in_frame_present_flag && !currSlice->field_pic_flag)
            currSlice->delta_pic_order_cnt[1] = s->se("SH: delta_pic_order_cnt[1]");
    }

    assert(currSlice->delta_pic_order_cnt_bottom >= -(1 << 31) + 1 &&
           currSlice->delta_pic_order_cnt_bottom <=  (1 << 31) - 1);
    assert(currSlice->delta_pic_order_cnt[0] >= -(1 << 31) + 1 &&
           currSlice->delta_pic_order_cnt[0] <=  (1 << 31) - 1);
    assert(currSlice->delta_pic_order_cnt[1] >= -(1 << 31) + 1 &&
           currSlice->delta_pic_order_cnt[1] <=  (1 << 31) - 1);

    currSlice->redundant_pic_cnt = 0;
    if (pps->redundant_pic_cnt_present_flag)
        currSlice->redundant_pic_cnt = s->ue("SH: redundant_pic_cnt");

    assert(currSlice->redundant_pic_cnt >= 0 && currSlice->redundant_pic_cnt <= 127);

    if (currSlice->slice_type == B_slice)
        currSlice->direct_spatial_mv_pred_flag = s->u(1, "SH: direct_spatial_mv_pred_flag");

    currSlice->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    currSlice->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;
    if (currSlice->slice_type == P_slice || currSlice->slice_type == SP_slice ||
        currSlice->slice_type == B_slice) {
        currSlice->num_ref_idx_active_override_flag = s->u(1, "SH: num_ref_idx_override_flag");
        if (currSlice->num_ref_idx_active_override_flag) {
            currSlice->num_ref_idx_l0_active_minus1 = s->ue("SH: num_ref_idx_l0_active_minus1");
            if (currSlice->slice_type == B_slice)
                currSlice->num_ref_idx_l1_active_minus1 = s->ue("SH: num_ref_idx_l1_active_minus1");
        }
    }

    if ((currSlice->slice_type == P_slice || currSlice->slice_type == SP_slice ||
         currSlice->slice_type == B_slice) &&
        !currSlice->field_pic_flag && pps->num_ref_idx_l0_default_active_minus1 > 15)
        assert(currSlice->num_ref_idx_active_override_flag == 1);
    if (currSlice->slice_type == B_slice &&
        !currSlice->field_pic_flag && pps->num_ref_idx_l1_default_active_minus1 > 15)
        assert(currSlice->num_ref_idx_active_override_flag == 1);

    if (!currSlice->field_pic_flag)
        assert(currSlice->num_ref_idx_l0_active_minus1 >= 0 &&
               currSlice->num_ref_idx_l0_active_minus1 <= 15);
    else
        assert(currSlice->num_ref_idx_l0_active_minus1 >= 0 &&
               currSlice->num_ref_idx_l0_active_minus1 <= 31);
    if (!currSlice->field_pic_flag)
        assert(currSlice->num_ref_idx_l1_active_minus1 >= 0 &&
               currSlice->num_ref_idx_l1_active_minus1 <= 15);
    else
        assert(currSlice->num_ref_idx_l1_active_minus1 >= 0 &&
               currSlice->num_ref_idx_l1_active_minus1 <= 31);

#if (MVC_EXTENSION_ENABLE)
    // if (nal_unit_type == 20 || nal_unit_type == 21)
    if (currSlice->svc_extension_flag == 0 || currSlice->svc_extension_flag == 1)
        ref_pic_list_mvc_modification(currSlice);
    else
#endif
        ref_pic_list_modification(currSlice);

    if ((pps->weighted_pred_flag && (currSlice->slice_type == P_slice || currSlice->slice_type == SP_slice)) ||
        (pps->weighted_bipred_idc == 1 && currSlice->slice_type == B_slice))
        pred_weight_table(currSlice);
    else {
        currSlice->luma_log2_weight_denom   = 5;
        currSlice->chroma_log2_weight_denom = 5;
    }

    if (currSlice->nal_ref_idc != 0)
        dec_ref_pic_marking(p_Vid, s, currSlice);

    currSlice->cabac_init_idc = 0;
    if (pps->entropy_coding_mode_flag && currSlice->slice_type != I_slice && currSlice->slice_type != SI_slice)
        currSlice->cabac_init_idc = s->ue("SH: cabac_init_idc");

    assert(currSlice->cabac_init_idc >= 0 && currSlice->cabac_init_idc <= 2);

    currSlice->slice_qp_delta = s->se("SH: slice_qp_delta");
    currSlice->SliceQpY = 26 + pps->pic_init_qp_minus26 + currSlice->slice_qp_delta;

    currSlice->sp_for_switch_flag = 0;
    currSlice->QsY = 0;
    if (currSlice->slice_type == SP_slice || currSlice->slice_type == SI_slice) {
        if (currSlice->slice_type == SP_slice)
            currSlice->sp_for_switch_flag = s->u(1, "SH: sp_for_switch_flag");
        currSlice->slice_qs_delta = s->se("SH: slice_qs_delta");
        currSlice->QsY = 26 + pps->pic_init_qs_minus26 + currSlice->slice_qs_delta;
    }

    assert(currSlice->SliceQpY >= -(sps->QpBdOffsetY) &&
           currSlice->SliceQpY <= 51);
    assert(currSlice->QsY >= 0 && currSlice->QsY <= 51);

    currSlice->disable_deblocking_filter_idc = 0;
    currSlice->slice_alpha_c0_offset_div2    = 0;
    currSlice->slice_beta_offset_div2        = 0;
    if (pps->deblocking_filter_control_present_flag) {
        currSlice->disable_deblocking_filter_idc  = s->ue("SH: disable_deblocking_filter_idc");
        if (currSlice->disable_deblocking_filter_idc != 1) {
            currSlice->slice_alpha_c0_offset_div2 = s->se("SH: slice_alpha_c0_offset_div2");
            currSlice->slice_beta_offset_div2     = s->se("SH: slice_beta_offset_div2");
        }
    }
    currSlice->FilterOffsetA = currSlice->slice_alpha_c0_offset_div2 << 1;
    currSlice->FilterOffsetB = currSlice->slice_beta_offset_div2 << 1;

    assert(currSlice->disable_deblocking_filter_idc >= 0 &&
           currSlice->disable_deblocking_filter_idc <= 2);
    assert(currSlice->slice_alpha_c0_offset_div2 >= -6 &&
           currSlice->slice_alpha_c0_offset_div2 <=  6);
    assert(currSlice->slice_beta_offset_div2 >= -6 &&
           currSlice->slice_beta_offset_div2 <=  6);

    if (pps->num_slice_groups_minus1 > 0 &&
        pps->slice_group_map_type >= 3 && pps->slice_group_map_type <= 5) {
        int len = (sps->pic_height_in_map_units_minus1+1)*(sps->pic_width_in_mbs_minus1+1)/
              (pps->slice_group_change_rate_minus1+1);
        if (((sps->pic_height_in_map_units_minus1+1)*(sps->pic_width_in_mbs_minus1+1))%
            (pps->slice_group_change_rate_minus1+1))
            len += 1;
        len = ceil(log2(len + 1));
        currSlice->slice_group_change_cycle = s->u(len, "SH: slice_group_change_cycle");
    }

    currSlice->MapUnitsInSliceGroup0 = min(currSlice->slice_group_change_cycle * pps->SliceGroupChangeRate, sps->PicSizeInMapUnits);
}

void ref_pic_list_modification(slice_t *currSlice)
{
    data_partition_t *s = &currSlice->parser.partArr[0];

    currSlice->ref_pic_list_modification_flag_l0 = 0;
    if (currSlice->slice_type != I_slice && currSlice->slice_type != SI_slice) {
        currSlice->ref_pic_list_modification_flag_l0 = s->u(1, "SH: ref_pic_list_modification_flag_l0");
        if (currSlice->ref_pic_list_modification_flag_l0) {
            int i = 0;
            do {
                currSlice->modification_of_pic_nums_idc[0][i] = s->ue("SH: modification_of_pic_nums_idc_l0");
                if (currSlice->modification_of_pic_nums_idc[0][i] == 0 ||
                    currSlice->modification_of_pic_nums_idc[0][i] == 1)
                    currSlice->abs_diff_pic_num_minus1[0][i] = s->ue("SH: abs_diff_pic_num_minus1_l0");
                else if (currSlice->modification_of_pic_nums_idc[0][i] == 2)
                    currSlice->long_term_pic_num[0][i] = s->ue("SH: long_term_pic_idx_l0");
            } while (currSlice->modification_of_pic_nums_idc[0][i++] != 3);
        }
    }

    currSlice->ref_pic_list_modification_flag_l1 = 0;
    if (currSlice->slice_type == B_slice) {
        currSlice->ref_pic_list_modification_flag_l1 = s->u(1, "SH: ref_pic_list_reordering_flag_l1");
        if (currSlice->ref_pic_list_modification_flag_l1) {
            int i = 0;
            do {
                currSlice->modification_of_pic_nums_idc[1][i] = s->ue("SH: modification_of_pic_nums_idc_l1");
                if (currSlice->modification_of_pic_nums_idc[1][i] == 0 ||
                    currSlice->modification_of_pic_nums_idc[1][i] == 1)
                    currSlice->abs_diff_pic_num_minus1[1][i] = s->ue("SH: abs_diff_pic_num_minus1_l1");
                else if (currSlice->modification_of_pic_nums_idc[1][i] == 2)
                    currSlice->long_term_pic_num[1][i] = s->ue("SH: long_term_pic_idx_l1");
            } while (currSlice->modification_of_pic_nums_idc[1][i++] != 3);
        }
    }
}

#if (MVC_EXTENSION_ENABLE)
void ref_pic_list_mvc_modification(slice_t *currSlice)
{
    data_partition_t *s = &currSlice->parser.partArr[0];

    currSlice->ref_pic_list_modification_flag_l0 = 0;
    if (currSlice->slice_type != I_slice && currSlice->slice_type != SI_slice) {
        currSlice->ref_pic_list_modification_flag_l0 = s->u(1, "SH: ref_pic_list_modification_flag_l0");
        if (currSlice->ref_pic_list_modification_flag_l0) {
            int i = 0;
            do {
                currSlice->modification_of_pic_nums_idc[0][i] = s->ue("SH: modification_of_pic_nums_idc_l0");
                if (currSlice->modification_of_pic_nums_idc[0][i] == 0 ||
                    currSlice->modification_of_pic_nums_idc[0][i] == 1)
                    currSlice->abs_diff_pic_num_minus1[0][i] = s->ue("SH: abs_diff_pic_num_minus1_l0");
                else if (currSlice->modification_of_pic_nums_idc[0][i] == 2)
                    currSlice->long_term_pic_num[0][i] = s->ue("SH: long_term_pic_idx_l0");
                else if (currSlice->modification_of_pic_nums_idc[0][i] == 4 ||
                         currSlice->modification_of_pic_nums_idc[0][i] == 5)
                    currSlice->abs_diff_view_idx_minus1[0][i] = s->ue("SH: abs_diff_view_idx_minus1_l0");
            } while (currSlice->modification_of_pic_nums_idc[0][i++] != 3);
        }
    }

    currSlice->ref_pic_list_modification_flag_l1 = 0;
    if (currSlice->slice_type == B_slice) {
        currSlice->ref_pic_list_modification_flag_l1 = s->u(1, "SH: ref_pic_list_reordering_flag_l1");
        if (currSlice->ref_pic_list_modification_flag_l1) {
            int i = 0;
            do {
                currSlice->modification_of_pic_nums_idc[1][i] = s->ue("SH: modification_of_pic_nums_idc_l1");
                if (currSlice->modification_of_pic_nums_idc[1][i] == 0 ||
                    currSlice->modification_of_pic_nums_idc[1][i] == 1)
                    currSlice->abs_diff_pic_num_minus1[1][i] = s->ue("SH: abs_diff_pic_num_minus1_l1");
                else if (currSlice->modification_of_pic_nums_idc[1][i] == 2)
                    currSlice->long_term_pic_num[1][i] = s->ue("SH: long_term_pic_idx_l1");
                else if (currSlice->modification_of_pic_nums_idc[1][i] == 4 ||
                         currSlice->modification_of_pic_nums_idc[1][i] == 5)
                    currSlice->abs_diff_view_idx_minus1[1][i] = s->ue("SH: abs_diff_view_idx_minus1_l1");
            } while (currSlice->modification_of_pic_nums_idc[1][i++] != 3);
        }
    }
}
#endif

void pred_weight_table(slice_t *currSlice)
{
    slice_t& slice = *currSlice;
    sps_t& sps = *slice.active_sps;
    data_partition_t *s = &slice.parser.partArr[0];

    slice.luma_log2_weight_denom = s->ue("SH: luma_log2_weight_denom");
    slice.chroma_log2_weight_denom = 0;
    if (sps.ChromaArrayType != 0)
        slice.chroma_log2_weight_denom = s->ue("SH: chroma_log2_weight_denom");

    assert(slice.luma_log2_weight_denom >= 0 && slice.luma_log2_weight_denom <= 7);
    assert(slice.chroma_log2_weight_denom >= 0 && slice.chroma_log2_weight_denom <= 7);

    for (int i = 0; i <= slice.num_ref_idx_l0_active_minus1; ++i) {
        slice.luma_weight_l0_flag[i] = s->u(1, "SH: luma_weight_flag_l0");
        slice.luma_weight_l0     [i] = 1 << slice.luma_log2_weight_denom;
        slice.luma_offset_l0     [i] = 0;
        if (slice.luma_weight_l0_flag[i]) {
            slice.luma_weight_l0[i] = s->se("SH: luma_weight_l0");
            slice.luma_offset_l0[i] = s->se("SH: luma_offset_l0");
        }

        assert(slice.luma_weight_l0[i] >= -128 && slice.luma_weight_l0[i] <= 127);
        assert(slice.luma_offset_l0[i] >= -128 && slice.luma_offset_l0[i] <= 127);

        if (sps.ChromaArrayType != 0) {
            slice.chroma_weight_l0_flag[i] = s->u(1, "SH: chroma_weight_flag_l0");
            for (int comp = 0; comp < 2; ++comp) {
                slice.chroma_weight_l0[i][comp] = 1 << slice.chroma_log2_weight_denom;
                slice.chroma_offset_l0[i][comp] = 0;
                if (slice.chroma_weight_l0_flag[i]) {
                    slice.chroma_weight_l0[i][comp] = s->se("SH: chroma_weight_l0");
                    slice.chroma_offset_l0[i][comp] = s->se("SH: chroma_offset_l0");
                }

                assert(slice.chroma_weight_l0[i][comp] >= -128 && slice.chroma_weight_l0[i][comp] <=  127);
                assert(slice.chroma_offset_l0[i][comp] >= -128 && slice.chroma_offset_l0[i][comp] <=  127);
            }
        }
    }

    if (slice.slice_type != B_slice)
        return;

    for (int i = 0; i <= slice.num_ref_idx_l1_active_minus1; ++i) {
        slice.luma_weight_l1_flag[i] = s->u(1, "SH: luma_weight_flag_l1");
        slice.luma_weight_l1     [i] = 1 << slice.luma_log2_weight_denom;
        slice.luma_offset_l1     [i] = 0;
        if (slice.luma_weight_l1_flag[i]) {
            slice.luma_weight_l1[i] = s->se("SH: luma_weight_l1");
            slice.luma_offset_l1[i] = s->se("SH: luma_offset_l1");
        }

        assert(slice.luma_weight_l1[i] >= -128 && slice.luma_weight_l1[i] <= 127);
        assert(slice.luma_offset_l1[i] >= -128 && slice.luma_offset_l1[i] <= 127);

        if (sps.ChromaArrayType != 0) {
            slice.chroma_weight_l1_flag[i] = s->u(1, "SH: chroma_weight_flag_l1");
            for (int comp = 0; comp < 2; ++comp) {
                slice.chroma_weight_l1[i][comp] = 1 << slice.chroma_log2_weight_denom;
                slice.chroma_offset_l1[i][comp] = 0;
                if (slice.chroma_weight_l1_flag[i]) {
                    slice.chroma_weight_l1[i][comp] = s->se("SH: chroma_weight_l1");
                    slice.chroma_offset_l1[i][comp] = s->se("SH: chroma_offset_l1");
                }

                assert(slice.chroma_weight_l1[i][comp] >= -128 && slice.chroma_weight_l1[i][comp] <=  127);
                assert(slice.chroma_offset_l1[i][comp] >= -128 && slice.chroma_offset_l1[i][comp] <=  127);
            }
        }
    }
}

void dec_ref_pic_marking(VideoParameters *p_Vid, data_partition_t *s, slice_t *currSlice)
{
    //data_partition_t *s = currSlice->parser.partArr[0].bitstream;

    sps_t *sps = currSlice->active_sps;
    int val;

    drpm_t* tmp_drpm, *tmp_drpm2;

    // free old buffer content
    while (currSlice->dec_ref_pic_marking_buffer) {
        tmp_drpm = currSlice->dec_ref_pic_marking_buffer;
        currSlice->dec_ref_pic_marking_buffer = tmp_drpm->Next;
        free(tmp_drpm);
    }

    if (currSlice->idr_flag
#if (MVC_EXTENSION_ENABLE)
        || (currSlice->svc_extension_flag == 0 && currSlice->NaluHeaderMVCExt.non_idr_flag == 0)
#endif
    ) {
        currSlice->no_output_of_prior_pics_flag = s->u(1, "SH: no_output_of_prior_pics_flag");
        p_Vid->no_output_of_prior_pics_flag = currSlice->no_output_of_prior_pics_flag;
        currSlice->long_term_reference_flag = s->u(1, "SH: long_term_reference_flag");

        assert(sps->max_num_ref_frames > 0 || currSlice->long_term_reference_flag == 1);
    } else {
        currSlice->adaptive_ref_pic_marking_mode_flag = s->u(1, "SH: adaptive_ref_pic_marking_mode_flag");
        if (currSlice->adaptive_ref_pic_marking_mode_flag) {
            do {
                tmp_drpm = (drpm_t*)calloc (1,sizeof (decoded_reference_picture_marking_t));
                tmp_drpm->Next = NULL;

                val = tmp_drpm->memory_management_control_operation = s->ue("SH: memory_management_control_operation");

                if ((val==1)||(val==3))
                    tmp_drpm->difference_of_pic_nums_minus1 = s->ue("SH: difference_of_pic_nums_minus1");
                if (val==2)
                    tmp_drpm->long_term_pic_num = s->ue("SH: long_term_pic_num");

                if ((val==3)||(val==6))
                    tmp_drpm->long_term_frame_idx = s->ue("SH: long_term_frame_idx");
                if (val==4)
                    tmp_drpm->max_long_term_frame_idx_plus1 = s->ue("SH: max_long_term_pic_idx_plus1");

                // add command
                if (currSlice->dec_ref_pic_marking_buffer==NULL)
                    currSlice->dec_ref_pic_marking_buffer=tmp_drpm;
                else {
                    tmp_drpm2=currSlice->dec_ref_pic_marking_buffer;
                    while (tmp_drpm2->Next!=NULL) tmp_drpm2=tmp_drpm2->Next;
                        tmp_drpm2->Next=tmp_drpm;
                }
            } while (val != 0);
        }
    }
}


void decode_poc(VideoParameters *p_Vid, slice_t *pSlice)
{
    sps_t *sps = p_Vid->active_sps;

    switch (sps->pic_order_cnt_type) {
    case 0:
        if (pSlice->idr_flag) {
            p_Vid->prevPicOrderCntMsb = 0;
            p_Vid->prevPicOrderCntLsb = 0;
        } else if (p_Vid->last_has_mmco_5) {
            p_Vid->prevPicOrderCntMsb = 0;
            p_Vid->prevPicOrderCntLsb = !p_Vid->last_pic_bottom_field ? pSlice->TopFieldOrderCnt : 0;
        }

        if ((pSlice->pic_order_cnt_lsb < p_Vid->prevPicOrderCntLsb) &&
            (p_Vid->prevPicOrderCntLsb - pSlice->pic_order_cnt_lsb >= sps->MaxPicOrderCntLsb / 2))
            pSlice->PicOrderCntMsb = p_Vid->prevPicOrderCntMsb + sps->MaxPicOrderCntLsb;
        else if ((pSlice->pic_order_cnt_lsb > p_Vid->prevPicOrderCntLsb) &&
                 (pSlice->pic_order_cnt_lsb - p_Vid->prevPicOrderCntLsb > sps->MaxPicOrderCntLsb / 2))
            pSlice->PicOrderCntMsb = p_Vid->prevPicOrderCntMsb - sps->MaxPicOrderCntLsb;
        else
            pSlice->PicOrderCntMsb = p_Vid->prevPicOrderCntMsb;

        pSlice->TopFieldOrderCnt    = 0;
        pSlice->BottomFieldOrderCnt = 0;
        if (!pSlice->field_pic_flag || !pSlice->bottom_field_flag)
            pSlice->TopFieldOrderCnt = pSlice->PicOrderCntMsb + pSlice->pic_order_cnt_lsb;
        if (!pSlice->field_pic_flag)
            pSlice->BottomFieldOrderCnt = pSlice->TopFieldOrderCnt + pSlice->delta_pic_order_cnt_bottom;
        else if (pSlice->bottom_field_flag)
            pSlice->BottomFieldOrderCnt = pSlice->PicOrderCntMsb + pSlice->pic_order_cnt_lsb;

        if (pSlice->nal_ref_idc) {
            p_Vid->prevPicOrderCntMsb = pSlice->PicOrderCntMsb;
            p_Vid->prevPicOrderCntLsb = pSlice->pic_order_cnt_lsb;
        }
        break;

    case 1:
        if (p_Vid->last_has_mmco_5) {
            p_Vid->prevFrameNum       = 0;
            p_Vid->prevFrameNumOffset = 0;
        }

        if (pSlice->idr_flag)
            pSlice->FrameNumOffset = 0;
        else if (p_Vid->prevFrameNum > pSlice->frame_num)
            pSlice->FrameNumOffset = p_Vid->prevFrameNumOffset + sps->MaxFrameNum;
        else
            pSlice->FrameNumOffset = p_Vid->prevFrameNumOffset;

        int32_t absFrameNum;
        if (sps->num_ref_frames_in_pic_order_cnt_cycle != 0)
            absFrameNum = pSlice->FrameNumOffset + pSlice->frame_num;
        else
            absFrameNum = 0;
        if (pSlice->nal_ref_idc == 0 && absFrameNum > 0)
            absFrameNum--;

        int32_t expectedPicOrderCnt;
        if (absFrameNum > 0) {
            int32_t picOrderCntCycleCnt        = (absFrameNum - 1) / sps->num_ref_frames_in_pic_order_cnt_cycle;
            int32_t frameNumInPicOrderCntCycle = (absFrameNum - 1) % sps->num_ref_frames_in_pic_order_cnt_cycle;
            expectedPicOrderCnt = picOrderCntCycleCnt * sps->ExpectedDeltaPerPicOrderCntCycle;
            for (int i = 0; i <= frameNumInPicOrderCntCycle; i++)
                expectedPicOrderCnt += sps->offset_for_ref_frame[i];
        } else
            expectedPicOrderCnt = 0;
        if (pSlice->nal_ref_idc == 0)
            expectedPicOrderCnt += sps->offset_for_non_ref_pic;

        pSlice->TopFieldOrderCnt    = 0;
        pSlice->BottomFieldOrderCnt = 0;
        if (!pSlice->field_pic_flag || !pSlice->bottom_field_flag)
            pSlice->TopFieldOrderCnt = expectedPicOrderCnt + pSlice->delta_pic_order_cnt[0];
        if (!pSlice->field_pic_flag)
            pSlice->BottomFieldOrderCnt = pSlice->TopFieldOrderCnt + sps->offset_for_top_to_bottom_field + pSlice->delta_pic_order_cnt[1];
        else if (pSlice->bottom_field_flag)
            pSlice->BottomFieldOrderCnt = expectedPicOrderCnt + sps->offset_for_top_to_bottom_field + pSlice->delta_pic_order_cnt[0];

        p_Vid->prevFrameNum       = pSlice->frame_num;
        p_Vid->prevFrameNumOffset = pSlice->FrameNumOffset;
        break;

    case 2:
        if (p_Vid->last_has_mmco_5) {
            p_Vid->prevFrameNum       = 0;
            p_Vid->prevFrameNumOffset = 0;
        }

        if (pSlice->idr_flag)
            pSlice->FrameNumOffset = 0;
        else if (p_Vid->prevFrameNum > pSlice->frame_num)
            pSlice->FrameNumOffset = p_Vid->prevFrameNumOffset + sps->MaxFrameNum;
        else
            pSlice->FrameNumOffset = p_Vid->prevFrameNumOffset;

        int tempPicOrderCnt;
        if (pSlice->idr_flag)
            tempPicOrderCnt = 0;
        else if (pSlice->nal_ref_idc == 0)
            tempPicOrderCnt = 2 * (pSlice->FrameNumOffset + pSlice->frame_num) - 1;
        else
            tempPicOrderCnt = 2 * (pSlice->FrameNumOffset + pSlice->frame_num);

        pSlice->TopFieldOrderCnt    = 0;
        pSlice->BottomFieldOrderCnt = 0;
        if (!pSlice->field_pic_flag || !pSlice->bottom_field_flag)
            pSlice->TopFieldOrderCnt = tempPicOrderCnt;
        if (!pSlice->field_pic_flag || pSlice->bottom_field_flag)
            pSlice->BottomFieldOrderCnt = tempPicOrderCnt;

        p_Vid->prevFrameNum       = pSlice->frame_num;
        p_Vid->prevFrameNumOffset = pSlice->FrameNumOffset;
        break;

    default:
        assert(false);
        break;
    }

    if (!pSlice->field_pic_flag)
        pSlice->PicOrderCnt = min(pSlice->TopFieldOrderCnt, pSlice->BottomFieldOrderCnt);
    else if (!pSlice->bottom_field_flag)
        pSlice->PicOrderCnt = pSlice->TopFieldOrderCnt;
    else
        pSlice->PicOrderCnt = pSlice->BottomFieldOrderCnt;
    p_Vid->PicOrderCnt  = pSlice->PicOrderCnt;
    p_Vid->prevFrameNum = pSlice->frame_num;
}
