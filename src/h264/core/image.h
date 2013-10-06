#ifndef _IMAGE_H_
#define _IMAGE_H_


#include "dpb.h"

struct slice_t;

void init_picture  (VideoParameters *p_Vid, slice_t *currSlice, InputParameters *p_Inp);

int  read_new_slice(slice_t *currSlice);
void exit_picture  (VideoParameters *p_Vid, storable_picture **dec_picture);

void decode_picture(VideoParameters *p_Vid);

#if (MVC_EXTENSION_ENABLE)
extern int GetVOIdx(VideoParameters *p_Vid, int iViewId);
#endif


#endif
