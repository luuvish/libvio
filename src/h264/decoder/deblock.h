#ifndef _DEBLOCK_H_
#define _DEBLOCK_H_


#include "global.h"
#include "dpb.h"

struct slice_t;

void init_Deblock(VideoParameters *p_Vid, int mb_aff_frame_flag);
void pic_deblock(VideoParameters *p_Vid, storable_picture *p);
// For 4:4:4 independent mode
void change_plane_JV(VideoParameters *p_Vid, int nplane, struct slice_t *pSlice);


struct deblock_t {
	void init();
	void deblock();
};


#endif /* _DEBLOCK_H_ */
