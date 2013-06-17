
/*!
 *************************************************************************************
 * \file mb_access.h
 *
 * \brief
 *    Functions for macroblock neighborhoods
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Karsten Suehring
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>  
 *************************************************************************************
 */

#ifndef _NEIGHBOUR_H_
#define _NEIGHBOUR_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void CheckAvailabilityOfNeighbors(Macroblock *currMB);
extern void CheckAvailabilityOfNeighborsMBAFF(Macroblock *currMB);
extern void CheckAvailabilityOfNeighborsNormal(Macroblock *currMB);

extern void getAffNeighbour         (Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
extern void getNonAffNeighbour      (Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
extern void get4x4Neighbour         (Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
extern void get4x4NeighbourBase     (Macroblock *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix);
extern Boolean mb_is_available      (int mbAddr, Macroblock *currMB);
extern void get_mb_pos              (VideoParameters *p_Vid, int mb_addr, int mb_size[2], short *x, short *y);
extern void get_mb_block_pos_normal (BlockPos *PicPos, int mb_addr, short *x, short *y);
extern void get_mb_block_pos_mbaff  (BlockPos *PicPos, int mb_addr, short *x, short *y);

void get_neighbors(Macroblock *currMB, PixelPos *block, int mb_x, int mb_y, int blockshape_x);
void check_dp_neighbors(Macroblock *currMB);

int predict_nnz(Macroblock *currMB, int block_type, int i,int j);
int predict_nnz_chroma(Macroblock *currMB, int i,int j);

#ifdef __cplusplus
}
#endif

#endif /* _NEIGHBOUR_H_ */
