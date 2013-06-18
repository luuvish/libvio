
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

#ifndef _BLOCK_H_
#define _BLOCK_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "transform.h"

struct macroblock_dec;

static const byte QP_SCALE_CR[52]=
{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
   12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
   28,29,29,30,31,32,32,33,34,34,35,35,36,36,37,37,
   37,38,38,38,39,39,39,39

};

//! look up tables for FRExt_chroma support
static const unsigned char subblk_offset_x[3][8][4] =
{
  {
    {0, 4, 0, 4},
    {0, 4, 0, 4},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}, 
  },
  { 
    {0, 4, 0, 4},
    {0, 4, 0, 4},
    {0, 4, 0, 4},
    {0, 4, 0, 4},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}, 
  },
  {
    {0, 4, 0, 4},
    {8,12, 8,12},
    {0, 4, 0, 4},
    {8,12, 8,12},
    {0, 4, 0, 4},
    {8,12, 8,12},
    {0, 4, 0, 4},
    {8,12, 8,12}  
  }
};


static const unsigned char subblk_offset_y[3][8][4] =
{
  {
    {0, 0, 4, 4},
    {0, 0, 4, 4},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}
  },
  { 
    {0, 0, 4, 4},
    {8, 8,12,12},
    {0, 0, 4, 4},
    {8, 8,12,12},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}
  },
  { 
    {0, 0, 4, 4},
    {0, 0, 4, 4},
    {8, 8,12,12},
    {8, 8,12,12},
    {0, 0, 4, 4},
    {0, 0, 4, 4},
    {8, 8,12,12},
    {8, 8,12,12}
  }
};

static const byte decode_block_scan[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

void iMBtrans4x4(struct macroblock_dec *currMB, ColorPlane pl, int smb);
void iMBtrans8x8(struct macroblock_dec *currMB, ColorPlane pl);

void itrans_sp_cr(struct macroblock_dec *currMB, int uv);

void Inv_Residual_trans_4x4(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
void Inv_Residual_trans_8x8(struct macroblock_dec *currMB, ColorPlane pl, int ioff,int joff);
void Inv_Residual_trans_16x16 (struct macroblock_dec *currMB, ColorPlane pl);
void Inv_Residual_trans_Chroma(struct macroblock_dec *currMB, int uv);

void itrans4x4   (struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
void itrans4x4_ls(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
void itrans_sp   (struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
void itrans_2    (struct macroblock_dec *currMB, ColorPlane pl);
void iTransform  (struct macroblock_dec *currMB, ColorPlane pl, int smb);

void copy_image_data       (imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2, int width, int height);
void copy_image_data_16x16 (imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2);
void copy_image_data_8x8   (imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2);
void copy_image_data_4x4   (imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2);

#ifdef __cplusplus
}
#endif

#endif
