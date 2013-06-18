
/*!
 ***************************************************************************
 *
 * \file transform.h
 *
 * \brief
 *    prototypes of transform functions
 *
 * \date
 *    10 July 2007
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    Alexis Michael Tourapis
 **************************************************************************/

#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#ifdef __cplusplus
extern "C" {
#endif

void forward4x4   (int **block , int **tblock, int pos_y, int pos_x);
void inverse4x4   (int **tblock, int **block , int pos_y, int pos_x);
void inverse8x8   (int **tblock, int **block , int pos_x);
void ihadamard4x4 (int **tblock, int **block);
void ihadamard4x2 (int **tblock, int **block);
void ihadamard2x2 (int block[4], int tblock[4]);

struct macroblock_dec;

void itrans8x8   (struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
void icopy8x8    (struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);

#ifdef __cplusplus
}
#endif

#endif //_TRANSFORM_H_
