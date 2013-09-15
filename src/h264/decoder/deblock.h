#ifndef _DEBLOCK_H_
#define _DEBLOCK_H_


#include "global.h"
#include "dpb.h"

struct slice_t;

struct deblock_t {
	void init();
	void deblock(VideoParameters *p_Vid, storable_picture *p);

private:
	void init_neighbors(VideoParameters *p_Vid);
	void make_frame_picture_JV(VideoParameters *p_Vid);

	void DeblockMb(VideoParameters *p_Vid, storable_picture *p, int MbQAddr);
	void DeblockPicture(VideoParameters *p_Vid, storable_picture *p);

  	void edge_loop_luma_ver  (ColorPlane pl, imgpel** Img, byte* Strength, macroblock_t* MbQ, int edge, storable_picture* p);
  	void edge_loop_luma_hor  (ColorPlane pl, imgpel** Img, byte* Strength, macroblock_t* MbQ, int edge, storable_picture* p);
  	void edge_loop_chroma_ver(ColorPlane pl, imgpel** Img, byte* Strength, macroblock_t* MbQ, int edge, storable_picture* p);
  	void edge_loop_chroma_hor(ColorPlane pl, imgpel** Img, byte* Strength, macroblock_t* MbQ, int edge, storable_picture* p);
};

void get_strength_ver    (mb_t *MbQ, int edge, int mvlimit, storable_picture *p);
void get_strength_hor    (mb_t *MbQ, int edge, int mvlimit, storable_picture *p);
void get_strength_ver_MBAff(byte *Strength, mb_t *MbQ, int edge, int mvlimit, storable_picture *p);
void get_strength_hor_MBAff(byte *Strength, mb_t *MbQ, int edge, int mvlimit, storable_picture *p);


extern deblock_t deblock;


#endif /* _DEBLOCK_H_ */
