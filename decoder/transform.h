
/*!
 ************************************************************************
 * \file block.h
 *
 * \brief
 *    definitions for block decoding functions
 *
 * \author
 *  Inge Lille-Langoy               <inge.lille-langoy@telenor.com>    \n
 *  Telenor Satellite Services                                         \n
 *  P.O.Box 6914 St.Olavs plass                                        \n
 *  N-0130 Oslo, Norway
 *
 ************************************************************************
 */

#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#ifdef __cplusplus
extern "C" {
#endif

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

struct macroblock_t;

void Inv_Residual_trans_4x4   (struct macroblock_t *currMB, ColorPlane pl, int ioff, int joff);
void Inv_Residual_trans_8x8   (struct macroblock_t *currMB, ColorPlane pl, int ioff,int joff);
void Inv_Residual_trans_16x16 (struct macroblock_t *currMB, ColorPlane pl, int ioff,int joff);
void Inv_Residual_trans_Chroma(struct macroblock_t *currMB, int uv);

void itrans4x4   (struct macroblock_t *currMB, ColorPlane pl, int ioff, int joff);
void itrans8x8   (struct macroblock_t *currMB, ColorPlane pl, int ioff, int joff);
void itrans16x16 (struct macroblock_t *currMB, ColorPlane pl);
void itrans_2    (struct macroblock_t *currMB, ColorPlane pl);
void itrans_420  (struct macroblock_t *currMB, ColorPlane pl);
void itrans_422  (struct macroblock_t *currMB, ColorPlane pl);
void iTransform  (struct macroblock_t *currMB, ColorPlane pl, int smb);


#ifdef __cplusplus
}
#endif

#endif /* _TRANSFORM_H_ */
