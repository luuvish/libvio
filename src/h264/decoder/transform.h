#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_


#include "global.h"
#include "transform.h"

struct macroblock_t;

static const byte QP_SCALE_CR[52] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35,
    35, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 39
};


void copy_image_data_4x4  (imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2);
void copy_image_data_8x8  (imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2);
void copy_image_data_16x16(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2);

void Inv_Residual_trans_4x4   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
void Inv_Residual_trans_8x8   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
void Inv_Residual_trans_16x16 (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
void Inv_Residual_trans_Chroma(macroblock_t* mb, int uv);

void itrans4x4   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
void itrans8x8   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
void itrans16x16 (macroblock_t* mb, ColorPlane pl);
void iTransform  (macroblock_t* mb, ColorPlane pl, int smb);


struct transform_t {
	void inverse_luma_dc        (macroblock_t* mb, ColorPlane pl);
	void inverse_chroma_dc      (macroblock_t* mb, ColorPlane pl);

	void inverse_transform_4x4  (macroblock_t* mb);
	void inverse_transform_8x8  (macroblock_t* mb);
	void inverse_transform_16x16(macroblock_t* mb);
};

extern transform_t transform;


#endif /* _TRANSFORM_H_ */
