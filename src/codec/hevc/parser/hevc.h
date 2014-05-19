#ifndef __HEVC_H__
#define __HEVC_H__

#include "vps.h"
#include "sps.h"
#include "pps.h"

typedef struct hevc_t {
	vps_t vps[MAX_VPS_ID];
	sps_t sps[MAX_SPS_ID];
	pps_t pps[MAX_PPS_ID];	
} hevc_t;

#endif /* __HEVC_H__ */
