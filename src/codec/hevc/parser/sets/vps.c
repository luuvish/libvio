#include "rbsp.h"
#include "vps.h"

#define MAX_OP_LAYERS (64)

void video_parameter_set_rbsp(vps_t *vps, rbsp_t *rbsp) {
	vps->video_parameter_set_id       = rbsp->u(4, "video_parameter_set_id");
	assert(vps->video_parameter_set_id < MAX_VPS_ID);
	vps->vps_temporal_id_nesting_flag = rbsp->u(1, "vps_temporal_id_nesting_flag");
	uint32_t vps_reserved_zero_2bits  = rbsp->u(2, "vps_reserved_zero_2bits");
	assert(vps_reserved_zero_2bits == 0);
	uint32_t vps_reserved_zero_6bits  = rbsp->u(6, "vps_reserved_zero_6bits");
	assert(vps_reserved_zero_6bits == 0);
	vps->vps_max_sub_layers_minus1    = rbsp->u(3, "vps_max_sub_layers_minus1");
	assert(vps->vps_max_sub_layers_minus1 < MAX_SUB_LAYERS);

	profile_tier_level(&vps->ptl, 1, vps->vps_max_sub_layers_minus1, rbsp);

	uint32_t vps_reserved_zero_12bits = rbsp->u(12, "vps_reserved_zero_12bits");
	assert(vps_reserved_zero_12bits == 0);

	for (int i = 0; i <= vps->vps_max_sub_layers_minus1; i++) {
		vps->vps_max_dec_pic_buffering[i] = rbsp->ue("vps_max_dec_pic_buffering[i]");
		vps->vps_num_reorder_pics     [i] = rbsp->ue("vps_num_reorder_pics[i]");
		vps->vps_max_latency_increase [i] = rbsp->ue("vps_max_latency_increase[i]");
		assert(vps->vps_max_dec_pic_buffering[i] <= MaxDpbSize);
		assert(vps->vps_num_reorder_pics     [i] <= vps->vps_max_dec_pic_buffering[i]);
		assert(vps->vps_max_latency_increase [i] <= ((1<<32)-2));
		if (i > 0) {
			assert(vps->vps_max_dec_pic_buffering[i-1] <= vps->vps_max_dec_pic_buffering[i]);
			assert(vps->vps_num_reorder_pics     [i-1] <= vps->vps_num_reorder_pics     [i]);
		}
	}

	vps->vps_num_hrd_parameters = rbsp->ue("vps_num_hrd_parameters");
	assert(vps->vps_num_hrd_parameters <= 1);
	assert(vps->vps_num_hrd_parameters < 1024);
	for (int i = 0; i < vps->vps_num_hrd_parameters; i++) {
		if (i > 0)
			operation_point(vps, i);
		hrd_parameters(i == 0, vps->vps_max_sub_layers_minus1);
	}

	uint32_t vps_extension_flag = rbsp->u(1, "vps_extension_flag");
	if (vps_extension_flag) {
		while (rbsp->more_rbsp_data())
			uint32_t vps_extension_data_flag = rbsp->u(1);
	}

	rbsp->rbsp_trailing_bits();
}

void profile_tier_level(ptl_t *ptl, ProfilePresentFlag, MaxNumSubLayersMinus1, rbsp_t *rbsp) {
	if (ProfilePresentFlag) {
		ptl->general.profile_space = rbsp->u(2, "XXX_profile_space[]");
		ptl->general.tier_flag     = rbsp->u(1, "XXX_tier_flag[]");
		ptl->general.profile_idc   = rbsp->u(5, "XXX_profile_idc[]");
		for (int i = 0; i < 32; i++)
			ptl->general.profile_compatability_flag[i] = rbsp->u(1, "XXX_profile_compatability_flag[][j]");

		uint32_t general_reverved_zero_16bits = rbsp->u(16, "XXX_reserved_zero_16bits[]");
		assert(general_reverved_zero_16bits == 0);
	}
	ptl->general.level_idc = rbsp->u(8, "general_level_idc");

	for (int i = 0; i < MaxNumSubLayersMinus1; i++) {
		ptl->sub_layer[i].profile_present_flag = rbsp->u(1, "sub_layer_profile_present_flag[i]");
		ptl->sub_layer[i].level_present_flag   = rbsp->u(1, "sub_layer_level_present_flag[i]");
		if (ProfilePresentFlag && ptl->sub_layer_profile_present_flag[i]) {
			ptl->sub_layer[i].profile_space    = rbsp->u(2, "XXX_profile_space[]");
			ptl->sub_layer[i].tier_flag        = rbsp->u(1, "XXX_tier_flag[]");
			ptl->sub_layer[i].profile_idc      = rbsp->u(5, "XXX_profile_idc[]");
			for (int j = 0; j < 32; j++)
				ptl->sub_layer[i].profile_compatability_flag[j] = rbsp->u(1, "XXX_profile_compatability_flag[][j]");
			uint32_t sub_layer_reverved_zero_16bits = rbsp->u(16, "XXX_reserved_zero_16bits[]");
			assert(sub_layer_reverved_zero_16bits == 0);
		}
		if (ptl->sub_layer[i].level_present_flag)
			ptl->sub_layer[i].level_idc = rbsp->u(8, "sub_layer_level_idc[i]");
	}
}

void operation_point(vps_t *vps, optIdx) {
	vps->op_num_layer_id_values_minus1[optIdx] = rbsp->ue("op_num_layer_id_values_minus1");
	assert(vps->op_num_layer_id_values_minus1[optIdx] < MAX_OP_LAYERS);
	assert(vps->op_num_layer_id_values_minus1[optIdx] == 0);
	for (int i = 0; i <= vps->op_num_layer_id_values_minus1; i++)
		vps->op_layer_id[optIdx][i] = rbsp->u(6, "op_layer_id");
}

void hrd_parameters() {

}
