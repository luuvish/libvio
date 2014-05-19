#ifndef __VPS_H__
#define __VPS_H__

#define MAX_VPS_ID     (16)
#define MAX_SUB_LAYERS ( 7)

typedef struct profile_tier_level_t {
	struct {
		uint32_t profile_space                         : 2; // u(2)
		uint32_t tier_flag                             : 1; // u(1)
		uint32_t profile_idc                           : 5; // u(5)
		uint32_t profile_compatability_flag[32]           ; // u(1)
		uint32_t level_idc                             : 8; // u(8)
	} general; 

	struct {
		uint32_t profile_present_flag                  : 1; // u(1)
		uint32_t level_present_flag                    : 1; // u(1)
		uint32_t profile_space                         : 2; // u(2)
		uint32_t tier_flag                             : 1; // u(1)
		uint32_t profile_idc                           : 5; // u(5)
		uint32_t profile_compatability_flag[32]           ; // u(1)
		uint32_t level_idc                             : 8; // u(8)
	} sub_layer[MAX_SUB_LAYERS];
} ptl_t;

typedef struct video_parameter_set_t {
	uint32_t video_parameter_set_id                    : 4; // u(4)
	uint32_t vps_temporal_id_nesting_flag              : 1; // u(1)
	uint32_t vps_max_sub_layers_minus1                 : 3; // u(3)

	ptl_t    ptl;

	uint32_t vps_max_dec_pic_buffering[MAX_SUB_LAYERS]    ; // ue(v)
	uint32_t vps_num_reorder_pics     [MAX_SUB_LAYERS]    ; // ue(v)
	uint32_t vps_max_latency_increase [MAX_SUB_LAYERS]    ; // ue(v)

	uint32_t vps_num_hrd_parameters                       ; // ue(v)
} vps_t;

#endif /* __VPS_H__ */
