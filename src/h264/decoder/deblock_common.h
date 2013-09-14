#ifndef _LOOP_FILTER_H_
#define _LOOP_FILTER_H_


#include "global.h"

void get_strength_ver    (mb_t *MbQ, int edge, int mvlimit, storable_picture *p);
void get_strength_hor    (mb_t *MbQ, int edge, int mvlimit, storable_picture *p);
void get_strength_ver_MBAff(byte *Strength, mb_t *MbQ, int edge, int mvlimit, storable_picture *p);
void get_strength_hor_MBAff(byte *Strength, mb_t *MbQ, int edge, int mvlimit, storable_picture *p);

void set_loop_filter_functions_mbaff(VideoParameters *p_Vid);
void set_loop_filter_functions_normal(VideoParameters *p_Vid);


#endif
