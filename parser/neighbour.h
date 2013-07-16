
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

struct macroblock_t;

void CheckAvailabilityOfNeighbors(struct macroblock_t *currMB);
void CheckAvailabilityOfNeighborsMBAFF(struct macroblock_t *currMB);
void CheckAvailabilityOfNeighborsNormal(struct macroblock_t *currMB);
void CheckAvailabilityOfNeighborsCABAC(struct macroblock_t *currMB);

void getAffNeighbour         (struct macroblock_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
void getNonAffNeighbour      (struct macroblock_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
void get4x4Neighbour         (struct macroblock_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
void get4x4NeighbourBase     (struct macroblock_t *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix);
Boolean mb_is_available      (int mbAddr, struct macroblock_t *currMB);
void get_mb_pos              (VideoParameters *p_Vid, int mb_addr, int mb_size[2], short *x, short *y);
void get_mb_block_pos_normal (BlockPos *PicPos, int mb_addr, short *x, short *y);
void get_mb_block_pos_mbaff  (BlockPos *PicPos, int mb_addr, short *x, short *y);

void get_neighbors(struct macroblock_t *currMB, PixelPos *block, int mb_x, int mb_y, int blockshape_x);
void check_dp_neighbors(struct macroblock_t *currMB);

int predict_nnz(struct macroblock_t *currMB, int block_type, int i,int j);
int predict_nnz_chroma(struct macroblock_t *currMB, int i,int j);

#ifdef __cplusplus
}
#endif

#endif /* _NEIGHBOUR_H_ */
