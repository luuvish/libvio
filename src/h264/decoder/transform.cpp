#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "transform.h"
#include "image.h"
#include "neighbour.h"
#include "transform.h"
#include "quantization.h"
#include "memalloc.h"
#include "intra_prediction.h"


// Macro defines
#define Q_BITS          15
#define DQ_BITS          6
#define DQ_BITS_8        6 

static inline int isign(int x)
{
    return ( (x > 0) - (x < 0));
}

static inline int isignab(int a, int b)
{
    return ((b) < 0) ? -abs(a) : abs(a);
}

static inline int rshift_rnd(int x, int a)
{
    return (a > 0) ? ((x + (1 << (a-1) )) >> a) : (x << (-a));
}

static inline int rshift_rnd_sf(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}

// SP decoding parameter (EQ. 8-425)
static const int A[4][4] = {
    { 16, 20, 16, 20},
    { 20, 25, 20, 25},
    { 16, 20, 16, 20},
    { 20, 25, 20, 25}
};



void copy_image_data_4x4(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2)
{
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1   + off1), (*imgBuf2   + off2), BLOCK_SIZE * sizeof (imgpel));
}

void copy_image_data_8x8(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2)
{  
    for (int j = 0; j < BLOCK_SIZE_8x8; j += 4) {
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    }
}

void copy_image_data_16x16(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2)
{
    for (int j = 0; j < MB_BLOCK_SIZE; j += 4) { 
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
    }
}

static void copy_image_data(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2, int width, int height)
{
    for (int j = 0; j < height; ++j)
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), width * sizeof (imgpel));
}


static void ihadamard2x2(int tblock[4], int block[4])
{
    int t0 = tblock[0] + tblock[1];
    int t1 = tblock[0] - tblock[1];
    int t2 = tblock[2] + tblock[3];
    int t3 = tblock[2] - tblock[3];

    block[0] = t0 + t2;
    block[1] = t1 + t3;
    block[2] = t0 - t2;
    block[3] = t1 - t3;
}

static void ihadamard2x4(int **tblock, int **block)
{
    int tmp[8];
    int *pTmp = tmp;

    // Horizontal
    for (int i = 0; i < BLOCK_SIZE; i++) {
        int *pblock = tblock[i];

        int t0 = *(pblock++);
        int t1 = *(pblock++);

        *(pTmp++) = t0 + t1;
        *(pTmp++) = t0 - t1;
    }

    // Vertical 
    for (int i = 0; i < 2; i++) {
        pTmp = tmp + i;

        int t0 = *pTmp;
        int t1 = *(pTmp += BLOCK_SIZE);
        int t2 = *(pTmp += BLOCK_SIZE);
        int t3 = *(pTmp += BLOCK_SIZE);

        int p0 = t0 + t2;
        int p1 = t0 - t2;
        int p2 = t1 - t3;
        int p3 = t1 + t3;
        
        block[0][i] = p0 + p3;
        block[1][i] = p1 + p2;
        block[2][i] = p1 - p2;
        block[3][i] = p0 - p3;
    }
}

static void ihadamard4x4(int **tblock, int **block)
{
    int tmp[16];
    int *pTmp = tmp;

    // Horizontal
    for (int i = 0; i < BLOCK_SIZE; i++) {
        int *pblock = tblock[i];

        int t0 = *(pblock++);
        int t1 = *(pblock++);
        int t2 = *(pblock++);
        int t3 = *(pblock  );

        int p0 = t0 + t2;
        int p1 = t0 - t2;
        int p2 = t1 - t3;
        int p3 = t1 + t3;

        *(pTmp++) = p0 + p3;
        *(pTmp++) = p1 + p2;
        *(pTmp++) = p1 - p2;
        *(pTmp++) = p0 - p3;
    }

    // Vertical 
    for (int i = 0; i < BLOCK_SIZE; i++) {
        pTmp = tmp + i;

        int t0 = *pTmp;
        int t1 = *(pTmp += BLOCK_SIZE);
        int t2 = *(pTmp += BLOCK_SIZE);
        int t3 = *(pTmp += BLOCK_SIZE);

        int p0 = t0 + t2;
        int p1 = t0 - t2;
        int p2 = t1 - t3;
        int p3 = t1 + t3;
        
        block[0][i] = p0 + p3;
        block[1][i] = p1 + p2;
        block[2][i] = p1 - p2;
        block[3][i] = p0 - p3;
    }
}

static void inverse4x4(int **tblock, int **block, int pos_y, int pos_x)
{
    int tmp[16];
    int *pTmp = tmp;

    // Horizontal
    for (int i = pos_y; i < pos_y + BLOCK_SIZE; i++) {
        int *pblock = &tblock[i][pos_x];

        int t0 = *(pblock++);
        int t1 = *(pblock++);
        int t2 = *(pblock++);
        int t3 = *(pblock  );

        int p0 =  t0 + t2;
        int p1 =  t0 - t2;
        int p2 = (t1 >> 1) - t3;
        int p3 =  t1 + (t3 >> 1);

        *(pTmp++) = p0 + p3;
        *(pTmp++) = p1 + p2;
        *(pTmp++) = p1 - p2;
        *(pTmp++) = p0 - p3;
    }

    // Vertical 
    for (int i = 0; i < BLOCK_SIZE; i++) {
        pTmp = tmp + i;

        int t0 = *pTmp;
        int t1 = *(pTmp += BLOCK_SIZE);
        int t2 = *(pTmp += BLOCK_SIZE);
        int t3 = *(pTmp += BLOCK_SIZE);

        int p0 = t0 + t2;
        int p1 = t0 - t2;
        int p2 = (t1 >> 1) - t3;
        int p3 = t1 + (t3 >> 1);

        block[pos_y    ][pos_x + i] = p0 + p3;
        block[pos_y + 1][pos_x + i] = p1 + p2;
        block[pos_y + 2][pos_x + i] = p1 - p2;
        block[pos_y + 3][pos_x + i] = p0 - p3;
    }
}

static void inverse8x8(int **tblock, int **block, int pos_x)
{
    int tmp[64];
    int *pTmp = tmp;

    // Horizontal  
    for (int i = 0; i < BLOCK_SIZE_8x8; i++) {
        int *pblock = &tblock[i][pos_x];

        int p0 = *(pblock++);
        int p1 = *(pblock++);
        int p2 = *(pblock++);
        int p3 = *(pblock++);
        int p4 = *(pblock++);
        int p5 = *(pblock++);
        int p6 = *(pblock++);
        int p7 = *(pblock  );

        int a0 = p0 + p4;
        int a1 = p0 - p4;
        int a2 = p6 - (p2 >> 1);
        int a3 = p2 + (p6 >> 1);

        int b0 = a0 + a3;
        int b2 = a1 - a2;
        int b4 = a1 + a2;
        int b6 = a0 - a3;

        a0 = -p3 + p5 - p7 - (p7 >> 1);    
        a1 =  p1 + p7 - p3 - (p3 >> 1);    
        a2 = -p1 + p7 + p5 + (p5 >> 1);    
        a3 =  p3 + p5 + p1 + (p1 >> 1);

        int b1 = a0 + (a3 >> 2);    
        int b3 = a1 + (a2 >> 2);    
        int b5 = a2 - (a1 >> 2);
        int b7 = a3 - (a0 >> 2);                

        *(pTmp++) = b0 + b7;
        *(pTmp++) = b2 - b5;
        *(pTmp++) = b4 + b3;
        *(pTmp++) = b6 + b1;
        *(pTmp++) = b6 - b1;
        *(pTmp++) = b4 - b3;
        *(pTmp++) = b2 + b5;
        *(pTmp++) = b0 - b7;
    }

    // Vertical 
    for (int i = 0; i < BLOCK_SIZE_8x8; i++) {
        pTmp = tmp + i;

        int p0 = *pTmp;
        int p1 = *(pTmp += BLOCK_SIZE_8x8);
        int p2 = *(pTmp += BLOCK_SIZE_8x8);
        int p3 = *(pTmp += BLOCK_SIZE_8x8);
        int p4 = *(pTmp += BLOCK_SIZE_8x8);
        int p5 = *(pTmp += BLOCK_SIZE_8x8);
        int p6 = *(pTmp += BLOCK_SIZE_8x8);
        int p7 = *(pTmp += BLOCK_SIZE_8x8);

        int a0 = p0 + p4;
        int a1 = p0 - p4;
        int a2 = p6 - (p2>>1);
        int a3 = p2 + (p6>>1);

        int b0 = a0 + a3;
        int b2 = a1 - a2;
        int b4 = a1 + a2;
        int b6 = a0 - a3;

        a0 = -p3 + p5 - p7 - (p7 >> 1);
        a1 =  p1 + p7 - p3 - (p3 >> 1);
        a2 = -p1 + p7 + p5 + (p5 >> 1);
        a3 =  p3 + p5 + p1 + (p1 >> 1);

        int b1 = a0 + (a3 >> 2);
        int b7 = a3 - (a0 >> 2);
        int b3 = a1 + (a2 >> 2);
        int b5 = a2 - (a1 >> 2);

        block[0][pos_x + i] = b0 + b7;
        block[1][pos_x + i] = b2 - b5;
        block[2][pos_x + i] = b4 + b3;
        block[3][pos_x + i] = b6 + b1;
        block[4][pos_x + i] = b6 - b1;
        block[5][pos_x + i] = b4 - b3;
        block[6][pos_x + i] = b2 + b5;
        block[7][pos_x + i] = b0 - b7;
    }
}

static void forward4x4(int **block, int **tblock, int pos_y, int pos_x)
{
    int tmp[16];
    int *pTmp = tmp;

    // Horizontal
    for (int i = pos_y; i < pos_y + BLOCK_SIZE; i++) {
        int *pblock = &block[i][pos_x];

        int p0 = *(pblock++);
        int p1 = *(pblock++);
        int p2 = *(pblock++);
        int p3 = *(pblock  );

        int t0 = p0 + p3;
        int t1 = p1 + p2;
        int t2 = p1 - p2;
        int t3 = p0 - p3;

        *(pTmp++) =  t0 + t1;
        *(pTmp++) = (t3 << 1) + t2;
        *(pTmp++) =  t0 - t1;    
        *(pTmp++) =  t3 - (t2 << 1);
    }

    // Vertical 
    for (int i = 0; i < BLOCK_SIZE; i++) {
        pTmp = tmp + i;

        int p0 = *pTmp;
        int p1 = *(pTmp += BLOCK_SIZE);
        int p2 = *(pTmp += BLOCK_SIZE);
        int p3 = *(pTmp += BLOCK_SIZE);

        int t0 = p0 + p3;
        int t1 = p1 + p2;
        int t2 = p1 - p2;
        int t3 = p0 - p3;

        tblock[pos_y    ][pos_x + i] = t0 +  t1;
        tblock[pos_y + 1][pos_x + i] = t2 + (t3 << 1);
        tblock[pos_y + 2][pos_x + i] = t0 -  t1;
        tblock[pos_y + 3][pos_x + i] = t3 - (t2 << 1);
    }
}


void Inv_Residual_trans_4x4(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    int    **cof     = currSlice->cof    [pl];
    int    **mb_rres = currSlice->mb_rres[pl];
    imgpel **mb_pred = currSlice->mb_pred[pl];
    imgpel **mb_rec  = currSlice->mb_rec [pl];
    int      temp[4][4];

    if (currMB->ipmode_DPCM == Intra_4x4_Vertical) {
        for (int i = 0; i < 4; ++i) {
            temp[0][i] = cof[joff + 0][ioff + i];
            temp[1][i] = cof[joff + 1][ioff + i] + temp[0][i];
            temp[2][i] = cof[joff + 2][ioff + i] + temp[1][i];
            temp[3][i] = cof[joff + 3][ioff + i] + temp[2][i];
        }
        for (int i = 0; i < 4; ++i) {
            mb_rres[joff    ][ioff + i] = temp[0][i];
            mb_rres[joff + 1][ioff + i] = temp[1][i];
            mb_rres[joff + 2][ioff + i] = temp[2][i];
            mb_rres[joff + 3][ioff + i] = temp[3][i];
        }
    } else if (currMB->ipmode_DPCM == Intra_4x4_Horizontal) {
        for (int j = 0; j < 4; ++j) {
            temp[j][0] = cof[joff + j][ioff    ];
            temp[j][1] = cof[joff + j][ioff + 1] + temp[j][0];
            temp[j][2] = cof[joff + j][ioff + 2] + temp[j][1];
            temp[j][3] = cof[joff + j][ioff + 3] + temp[j][2];
        }
        for (int j = 0; j < 4; ++j) {
            mb_rres[joff + j][ioff    ] = temp[j][0];
            mb_rres[joff + j][ioff + 1] = temp[j][1];
            mb_rres[joff + j][ioff + 2] = temp[j][2];
            mb_rres[joff + j][ioff + 3] = temp[j][3];
        }
    } else {
        for (int j = joff; j < joff + BLOCK_SIZE; ++j) {
            for (int i = ioff; i < ioff + BLOCK_SIZE; ++i)
                mb_rres[j][i] = cof[j][i];
        }
    }

    for (int j = joff; j < joff + BLOCK_SIZE; ++j) {
        for (int i = ioff; i < ioff + BLOCK_SIZE; ++i)
            mb_rec[j][i] = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
}

void Inv_Residual_trans_8x8(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    int    **cof     = currSlice->cof    [pl];
    int    **mb_rres = currSlice->mb_rres[pl];
    imgpel **mb_pred = currSlice->mb_pred[pl];
    imgpel **mb_rec  = currSlice->mb_rec [pl];
    int      temp[8][8];

    if (currMB->ipmode_DPCM == Intra_4x4_Vertical) {
        for (int i = 0; i < 8; ++i) {
            temp[0][i] = cof[joff + 0][ioff + i];
            temp[1][i] = cof[joff + 1][ioff + i] + temp[0][i];
            temp[2][i] = cof[joff + 2][ioff + i] + temp[1][i];
            temp[3][i] = cof[joff + 3][ioff + i] + temp[2][i];
            temp[4][i] = cof[joff + 4][ioff + i] + temp[3][i];
            temp[5][i] = cof[joff + 5][ioff + i] + temp[4][i];
            temp[6][i] = cof[joff + 6][ioff + i] + temp[5][i];
            temp[7][i] = cof[joff + 7][ioff + i] + temp[6][i];
        }
        for (int i = 0; i < 8; ++i) {
            mb_rres[joff    ][ioff + i] = temp[0][i];
            mb_rres[joff + 1][ioff + i] = temp[1][i];
            mb_rres[joff + 2][ioff + i] = temp[2][i];
            mb_rres[joff + 3][ioff + i] = temp[3][i];
            mb_rres[joff + 4][ioff + i] = temp[4][i];
            mb_rres[joff + 5][ioff + i] = temp[5][i];
            mb_rres[joff + 6][ioff + i] = temp[6][i];
            mb_rres[joff + 7][ioff + i] = temp[7][i];
        }
    } else if (currMB->ipmode_DPCM == Intra_4x4_Horizontal) {
        for (int i = 0; i < 8; ++i) {
            temp[i][0] = mb_rres[joff + i][ioff + 0];
            temp[i][1] = mb_rres[joff + i][ioff + 1] + temp[i][0];
            temp[i][2] = mb_rres[joff + i][ioff + 2] + temp[i][1];
            temp[i][3] = mb_rres[joff + i][ioff + 3] + temp[i][2];
            temp[i][4] = mb_rres[joff + i][ioff + 4] + temp[i][3];
            temp[i][5] = mb_rres[joff + i][ioff + 5] + temp[i][4];
            temp[i][6] = mb_rres[joff + i][ioff + 6] + temp[i][5];
            temp[i][7] = mb_rres[joff + i][ioff + 7] + temp[i][6];
        }
        for (int i = 0; i < 8; ++i) {
            mb_rres[joff + i][ioff + 0] = temp[i][0];
            mb_rres[joff + i][ioff + 1] = temp[i][1];
            mb_rres[joff + i][ioff + 2] = temp[i][2];
            mb_rres[joff + i][ioff + 3] = temp[i][3];
            mb_rres[joff + i][ioff + 4] = temp[i][4];
            mb_rres[joff + i][ioff + 5] = temp[i][5];
            mb_rres[joff + i][ioff + 6] = temp[i][6];
            mb_rres[joff + i][ioff + 7] = temp[i][7];
        }
    }

    for (int j = joff; j < joff + BLOCK_SIZE * 2; ++j) {
        for (int i = ioff; i < ioff + BLOCK_SIZE * 2; ++i)
            mb_rec[j][i] = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
}

void Inv_Residual_trans_16x16(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    int    **cof     = currSlice->cof    [pl];
    int    **mb_rres = currSlice->mb_rres[pl];
    imgpel **mb_pred = currSlice->mb_pred[pl];
    imgpel **mb_rec  = currSlice->mb_rec [pl];
    int      temp[16][16];

    if (currMB->ipmode_DPCM == Intra_16x16_Vertical) {
        for (int i = 0; i < MB_BLOCK_SIZE; ++i) {
            temp[0][i] = cof[0][i];
            for (int j = 1; j < MB_BLOCK_SIZE; j++)
                temp[j][i] = cof[j][i] + temp[j-1][i];
        }
        for (int i = 0; i < MB_BLOCK_SIZE; ++i) {
            for (int j = 0; j < MB_BLOCK_SIZE; j++)
                mb_rres[j][i] = temp[j][i];
        }
    } else if (currMB->ipmode_DPCM == Intra_16x16_Horizontal) {
        for (int j = 0; j < MB_BLOCK_SIZE; ++j) {
            temp[j][ 0] = cof[j][ 0  ];
            for (int i = 1; i < MB_BLOCK_SIZE; i++)
                temp[j][i] = cof[j][i] + temp[j][i-1];
        }
        for (int j = 0; j < MB_BLOCK_SIZE; ++j) {
            for (int i = 0; i < MB_BLOCK_SIZE; ++i)
                mb_rres[j][i] = temp[j][i];
        }
    } else {
        for (int j = 0; j < MB_BLOCK_SIZE; ++j) {
            for (int i = 0; i < MB_BLOCK_SIZE; ++i)
                mb_rres[j][i] = cof[j][i];
        }
    }

    for (int j = 0; j < MB_BLOCK_SIZE; ++j) {
        for (int i = 0; i < MB_BLOCK_SIZE; ++i)
            mb_rec[j][i] = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
}

void Inv_Residual_trans_Chroma(mb_t *currMB, int uv)  
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int    **cof     = currSlice->cof    [uv + 1];
    int    **mb_rres = currSlice->mb_rres[uv + 1];
    imgpel **mb_pred = currSlice->mb_pred[uv + 1];
    imgpel **mb_rec  = currSlice->mb_rec [uv + 1];
    int      temp[16][16];

    if (currMB->intra_chroma_pred_mode == Intra_Chroma_Vertical) {
        for (int i = 0; i < sps->MbWidthC; i++) {
            temp[0][i] = cof[0][i];
            for (int j = 1; j < sps->MbHeightC; j++)
                temp[j][i] = temp[j-1][i] + cof[j][i];
        }
        for (int i = 0; i < sps->MbWidthC; i++) {
            for (int j = 0; j < sps->MbHeightC; j++)
                mb_rres[j][i] = temp[j][i];
        }
    } else if (currMB->intra_chroma_pred_mode == Intra_Chroma_Horizontal) {
        for (int i = 0; i < sps->MbHeightC; i++) {
            temp[i][0] = cof[i][0];
            for (int j = 1; j < sps->MbWidthC; j++)
                temp[i][j] = temp[i][j-1] + cof[i][j];
        }
        for (int i = 0; i < sps->MbHeightC; i++) {
            for (int j = 0; j < sps->MbWidthC; j++)
                mb_rres[i][j] = temp[i][j];
        }
    } else {
        for (int j = 0; j < sps->MbHeightC; j++) {
            for (int i = 0; i < sps->MbWidthC; i++)
                mb_rres[j][i] = cof[j][i];
        }
    }

    int max_pel_value_comp = (1 << sps->BitDepthC) - 1;
    for (int j = 0; j < sps->MbHeightC; ++j) {
        for (int i = 0; i < sps->MbWidthC; ++i)
            mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_pred[j][i] + mb_rres[j][i]);
    }
}


static void recon8x8(int **m7, imgpel **mb_rec, imgpel **mpr, int max_imgpel_value, int ioff)
{
    for (int j = 0; j < 8; j++) {
        int    *m_tr  = (*m7++) + ioff;
        imgpel *m_rec = (*mb_rec++) + ioff;
        imgpel *m_prd = (*mpr++) + ioff;

        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec++ = (imgpel) clip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
        *m_rec   = (imgpel) clip1(max_imgpel_value, (*m_prd  ) + rshift_rnd_sf(*m_tr  , DQ_BITS_8)); 
    }
}

static void recon8x8_lossless(int **m7, imgpel **mb_rec, imgpel **mpr, int max_imgpel_value, int ioff)
{
    for (int j = 0; j < 8; j++) {
        for (int i = ioff; i < ioff + 8; i++)
            (*mb_rec)[i] = (imgpel) clip1<int>(max_imgpel_value, ((*m7)[i] + (long)(*mpr)[i])); 
        mb_rec++;
        m7++;
        mpr++;
    }
}

static void sample_reconstruct(imgpel **curImg, imgpel **mpr, int **mb_rres, int mb_x, int opix_x, int width, int height, int max_imgpel_value, int dq_bits)
{
    for (int j = 0; j < height; j++) {
        imgpel *imgOrg = &curImg[j][opix_x];
        imgpel *imgPred = &mpr[j][mb_x];
        int *m7 = &mb_rres[j][mb_x]; 
        for (int i = 0; i < width; i++)
            *imgOrg++ = (imgpel) clip1( max_imgpel_value, rshift_rnd_sf(*m7++, dq_bits) + *imgPred++);
    }
}

static void icopy4x4(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    imgpel **mb_rec = &currSlice->mb_rec[pl][joff];
    imgpel **mpr    = &currSlice->mb_pred[pl][joff];

    for (int j = 0; j < 4; j++)
        memcpy((*mb_rec++) + ioff, (*mpr++) + ioff, 4 * sizeof(imgpel));
}

static void icopy8x8(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    imgpel **mb_rec = &currSlice->mb_rec[pl][joff];
    imgpel **mpr    = &currSlice->mb_pred[pl][joff];

    for (int j = 0; j < 8; j++)
        memcpy((*mb_rec++) + ioff, (*mpr++) + ioff, 8 * sizeof(imgpel));
}

void itrans4x4(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int  **mb_rres = currSlice->mb_rres[pl];
    int    max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    inverse4x4(currSlice->cof[pl], mb_rres, joff, ioff);

    sample_reconstruct(&currSlice->mb_rec[pl][joff], &currSlice->mb_pred[pl][joff], &mb_rres[joff], ioff, ioff, BLOCK_SIZE, BLOCK_SIZE, max_pel_value_comp, DQ_BITS);
}

static void itrans4x4_ls(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int    **cof     = currSlice->cof    [pl];
    imgpel **mb_pred = currSlice->mb_pred[pl];
    imgpel **mb_rec  = currSlice->mb_rec [pl];

    for (int j = joff; j < joff + BLOCK_SIZE; ++j) {
        for (int i = ioff; i < ioff + BLOCK_SIZE; ++i)
            mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_pred[j][i] + cof[j][i]);
    }
}

void itrans8x8(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int    **cof     = currSlice->cof    [pl];
    int    **mb_rres = currSlice->mb_rres[pl];
    imgpel **mb_pred = currSlice->mb_pred[pl];
    imgpel **mb_rec  = currSlice->mb_rec [pl];

    if (currMB->TransformBypassModeFlag == TRUE)
        recon8x8_lossless(&cof[joff], &mb_rec[joff], &mb_pred[joff], max_pel_value_comp, ioff);
    else {
        inverse8x8(&cof[joff], &mb_rres[joff], ioff);
        recon8x8  (&mb_rres[joff], &mb_rec[joff], &mb_pred[joff], max_pel_value_comp, ioff);
    }
}

void itrans16x16(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int **cof = currSlice->cof[pl];
    int **mb_rres = currSlice->mb_rres[pl];

    for (int jj = 0; jj < MB_BLOCK_SIZE; jj += BLOCK_SIZE) {
        inverse4x4(cof, mb_rres, jj, 0);
        inverse4x4(cof, mb_rres, jj, 4);
        inverse4x4(cof, mb_rres, jj, 8);
        inverse4x4(cof, mb_rres, jj, 12);
    }

    sample_reconstruct(currSlice->mb_rec[pl], currSlice->mb_pred[pl], mb_rres, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE, max_pel_value_comp, DQ_BITS);
}


void itrans_2(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int transform_pl = sps->separate_colour_plane_flag ? PLANE_Y : pl;
    int **cof = currSlice->cof[transform_pl];
    int qp_scaled = currMB->qp_scaled[transform_pl];
    int qp_per = qp_scaled / 6;
    int qp_rem = qp_scaled % 6;

    int invLevelScale = currSlice->InvLevelScale4x4_Intra[pl][qp_rem][0][0];
    int **M4;
    get_mem2Dint(&M4, BLOCK_SIZE, BLOCK_SIZE);

    // horizontal
    for (int j = 0; j < 4; ++j) {
        M4[j][0] = cof[j << 2][ 0];
        M4[j][1] = cof[j << 2][ 4];
        M4[j][2] = cof[j << 2][ 8];
        M4[j][3] = cof[j << 2][12];
    }

    ihadamard4x4(M4, M4);

    for (int j = 0; j < MB_BLOCK_SIZE; j += BLOCK_SIZE) {
        for (int i = 0; i < MB_BLOCK_SIZE; i += BLOCK_SIZE)
            cof[j][i] = rshift_rnd(((M4[j / 4][i / 4] * invLevelScale) << qp_per), 6);
    }

    free_mem2Dint(M4);
}

void itrans_420(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;

    int **cof = currSlice->cof[pl];
    int qp_scaled = currMB->qp_scaled[pl];
    int qp_per = qp_scaled / 6;
    int qp_rem = qp_scaled % 6;

    int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
        currSlice->InvLevelScale4x4_Intra[pl][qp_rem] :
        currSlice->InvLevelScale4x4_Inter[pl][qp_rem];
    int invLevelScale = InvLevelScale4x4[0][0];

    int M4[4];
    M4[0] = cof[0][0];
    M4[1] = cof[0][4];
    M4[2] = cof[4][0];
    M4[3] = cof[4][4];

    ihadamard2x2(M4, M4);

    cof[0][0] = ((M4[0] * invLevelScale) << qp_per) >> 5;
    cof[0][4] = ((M4[1] * invLevelScale) << qp_per) >> 5;
    cof[4][0] = ((M4[2] * invLevelScale) << qp_per) >> 5;
    cof[4][4] = ((M4[3] * invLevelScale) << qp_per) >> 5;
}

void itrans_422(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int **cof = currSlice->cof[pl];
    int qp_scaled = currMB->qpc[pl - 1] + 3 + sps->QpBdOffsetC;
    int qp_per = qp_scaled / 6;
    int qp_rem = qp_scaled % 6;

    int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
        currSlice->InvLevelScale4x4_Intra[pl][qp_rem] :
        currSlice->InvLevelScale4x4_Inter[pl][qp_rem];
    int invLevelScale = InvLevelScale4x4[0][0];
    int **M4;
    get_mem2Dint(&M4, BLOCK_SIZE, 2);

    for (int j = 0; j < 4; j++) {
        M4[j][0] = cof[j << 2][0];
        M4[j][1] = cof[j << 2][4];
    }

    ihadamard2x4(M4, M4);

    for (int j = 0; j < sps->MbHeightC; j += BLOCK_SIZE) {
        for (int i = 0; i < sps->MbWidthC; i += BLOCK_SIZE)
            cof[j][i] = rshift_rnd_sf((M4[j / 4][i / 4] * invLevelScale) << qp_per, 6);
    }

    free_mem2Dint(M4);
}


static void itrans_sp(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int qp = (currSlice->slice_type == SI_SLICE) ? currSlice->QsY : currSlice->SliceQpY;
    int qp_per = qp / 6;
    int qp_rem = qp % 6;

    int qp_per_sp = currSlice->QsY / 6;
    int qp_rem_sp = currSlice->QsY % 6;
    int q_bits_sp = Q_BITS + qp_per_sp;

    int    **cof     = currSlice->cof    [pl];
    int    **mb_rres = currSlice->mb_rres[pl];
    imgpel **mb_pred = currSlice->mb_pred[pl];
    imgpel **mb_rec  = currSlice->mb_rec [pl];
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    const int (*InvLevelScale4x4)  [4] = dequant_coef[qp_rem];
    const int (*InvLevelScale4x4SP)[4] = dequant_coef[qp_rem_sp];  
    int **PBlock;  

    get_mem2Dint(&PBlock, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    for (int j = 0; j < BLOCK_SIZE; ++j) {
        PBlock[j][0] = mb_pred[joff + j][ioff    ];
        PBlock[j][1] = mb_pred[joff + j][ioff + 1];
        PBlock[j][2] = mb_pred[joff + j][ioff + 2];
        PBlock[j][3] = mb_pred[joff + j][ioff + 3];
    }

    forward4x4(PBlock, PBlock, 0, 0);

    for (int j = 0; j < BLOCK_SIZE; ++j) {
        for (int i = 0; i < BLOCK_SIZE; ++i) {
            // recovering coefficient since they are already dequantized earlier
            int icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
            int ilev;
            if (currSlice->sp_for_switch_flag || currSlice->slice_type == SI_SLICE) {
                ilev = rshift_rnd_sf(abs(PBlock[j][i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
                ilev = isignab(ilev, PBlock[j][i]) + icof;
            } else {
                ilev = PBlock[j][i] + ((icof * InvLevelScale4x4[j][i] * A[j][i] <<  qp_per) >> 6);
                ilev = isign(ilev) * rshift_rnd_sf(abs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
            }
            cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
        }
    }

    inverse4x4(cof, mb_rres, joff, ioff);

    for (int j = joff; j < joff + BLOCK_SIZE; ++j) {
        mb_rec[j][ioff    ] = (imgpel)clip1(max_pel_value_comp, rshift_rnd_sf(mb_rres[j][ioff    ], DQ_BITS));
        mb_rec[j][ioff + 1] = (imgpel)clip1(max_pel_value_comp, rshift_rnd_sf(mb_rres[j][ioff + 1], DQ_BITS));
        mb_rec[j][ioff + 2] = (imgpel)clip1(max_pel_value_comp, rshift_rnd_sf(mb_rres[j][ioff + 2], DQ_BITS));
        mb_rec[j][ioff + 3] = (imgpel)clip1(max_pel_value_comp, rshift_rnd_sf(mb_rres[j][ioff + 3], DQ_BITS));
    }

    free_mem2Dint(PBlock);
}


static void itrans_sp_cr(mb_t *currMB, int uv)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int i,j,ilev, icof, n2,n1;
    int mp1[BLOCK_SIZE];
    imgpel **mb_pred = currSlice->mb_pred[uv + 1];
    int    **cof = currSlice->cof[uv + 1];
    int **PBlock = new_mem2Dint(MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    int qp_per    = (currSlice->SliceQpY < 0 ? currSlice->SliceQpY : QP_SCALE_CR[currSlice->SliceQpY]) / 6;
    int qp_rem    = (currSlice->SliceQpY < 0 ? currSlice->SliceQpY : QP_SCALE_CR[currSlice->SliceQpY]) % 6;

    int qp_per_sp = (currSlice->QsY < 0 ? currSlice->QsY : QP_SCALE_CR[currSlice->QsY]) / 6;
    int qp_rem_sp = (currSlice->QsY < 0 ? currSlice->QsY : QP_SCALE_CR[currSlice->QsY]) % 6;
    int q_bits_sp = Q_BITS + qp_per_sp;  

    if (currSlice->slice_type == SI_SLICE) {
        qp_per = qp_per_sp;
        qp_rem = qp_rem_sp;
    }

    for (j = 0; j < mb_cr_size_y; ++j) {
        for (i = 0; i < mb_cr_size_x; ++i) {
            PBlock[j][i] = mb_pred[j][i];
            mb_pred[j][i] = 0;
        }
    }

    for (n2 = 0; n2 < mb_cr_size_y; n2 += BLOCK_SIZE) {
        for (n1 = 0; n1 < mb_cr_size_x; n1 += BLOCK_SIZE)
            forward4x4(PBlock, PBlock, n2, n1);
    }

    //     2X2 transform of DC coeffs.
    mp1[0] = (PBlock[0][0] + PBlock[4][0] + PBlock[0][4] + PBlock[4][4]);
    mp1[1] = (PBlock[0][0] - PBlock[4][0] + PBlock[0][4] - PBlock[4][4]);
    mp1[2] = (PBlock[0][0] + PBlock[4][0] - PBlock[0][4] - PBlock[4][4]);
    mp1[3] = (PBlock[0][0] - PBlock[4][0] - PBlock[0][4] + PBlock[4][4]);

    if (currSlice->sp_for_switch_flag || currSlice->slice_type == SI_SLICE) {
        for (n2 = 0; n2 < 2; ++n2) {
            for (n1 = 0; n1 < 2; ++n1) {
                //quantization fo predicted block
                ilev = rshift_rnd_sf(abs (mp1[n1+n2*2]) * quant_coef[qp_rem_sp][0][0], q_bits_sp + 1);
                //addition
                ilev = isignab(ilev, mp1[n1+n2*2]) + cof[n2<<2][n1<<2];
                //dequantization
                mp1[n1+n2*2] =ilev * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
            }
        }

        for (n2 = 0; n2 < mb_cr_size_y; n2 += BLOCK_SIZE) {
            for (n1 = 0; n1 < mb_cr_size_x; n1 += BLOCK_SIZE) {
                for (j = 0; j < BLOCK_SIZE; ++j) {
                    for (i = 0; i < BLOCK_SIZE; ++i) {
                        // recovering coefficient since they are already dequantized earlier
                        cof[n2 + j][n1 + i] = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];
                        //quantization of the predicted block
                        ilev = rshift_rnd_sf(abs(PBlock[n2 + j][n1 + i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
                        //addition of the residual
                        ilev = isignab(ilev,PBlock[n2 + j][n1 + i]) + cof[n2 + j][n1 + i];
                        // Inverse quantization
                        cof[n2 + j][n1 + i] = ilev * dequant_coef[qp_rem_sp][j][i] << qp_per_sp;
                    }
                }
            }
        }
    } else {
        for (n2 = 0; n2 < 2; ++n2) {
            for (n1 = 0; n1 < 2; ++n1) {
                ilev = mp1[n1+n2*2] + (((cof[n2<<2][n1<<2] * dequant_coef[qp_rem][0][0] * A[0][0]) << qp_per) >> 5);
                ilev = isign(ilev) * rshift_rnd_sf(abs(ilev) * quant_coef[qp_rem_sp][0][0], q_bits_sp + 1);
                mp1[n1+n2*2] = ilev * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
            }
        }

        for (n2 = 0; n2 < mb_cr_size_y; n2 += BLOCK_SIZE) {
            for (n1 = 0; n1 < mb_cr_size_x; n1 += BLOCK_SIZE) {
                for (j = 0; j < BLOCK_SIZE; ++j) {
                    for (i = 0; i < BLOCK_SIZE; ++i) {
                        // recovering coefficient since they are already dequantized earlier
                        icof = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];
                        //dequantization and addition of the predicted block      
                        ilev = PBlock[n2 + j][n1 + i] + ((icof * dequant_coef[qp_rem][j][i] * A[j][i] << qp_per) >> 6);
                        //quantization and dequantization
                        ilev = isign(ilev) * rshift_rnd_sf(abs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
                        cof[n2 + j][n1 + i] = ilev * dequant_coef[qp_rem_sp][j][i] << qp_per_sp;
                    }
                }
            }
        }
    }

    cof[0][0] = (mp1[0] + mp1[1] + mp1[2] + mp1[3]) >> 1;
    cof[0][4] = (mp1[0] + mp1[1] - mp1[2] - mp1[3]) >> 1;
    cof[4][0] = (mp1[0] - mp1[1] + mp1[2] - mp1[3]) >> 1;
    cof[4][4] = (mp1[0] - mp1[1] - mp1[2] + mp1[3]) >> 1;

    free_mem2Dint(PBlock);
}


static void iMBtrans4x4(mb_t *currMB, ColorPlane pl, int smb)
{
    void (*itrans_4x4)(struct macroblock_t *currMB, ColorPlane pl, int ioff, int joff);
    itrans_4x4 = smb ? itrans_sp : Inv_Residual_trans_4x4;

    for (int y = 0; y < 16; y += 4) {
        for (int x = 0; x < 16; x += 4) {
            int block8x8 = (y / 8) * 2 + (x / 8);
            if (smb || currMB->TransformBypassModeFlag)
                itrans_4x4(currMB, pl, x, y);
            else if (currMB->cbp & (1 << block8x8))
                itrans4x4(currMB, pl, x, y);
            else
                icopy4x4(currMB, pl, x, y);
        }
    }
}

static void iMBtrans8x8(mb_t *currMB, ColorPlane pl)
{
    for (int y = 0; y < 16; y += 8) {
        for (int x = 0; x < 16; x += 8) {
            int block8x8 = (y / 8) * 2 + (x / 8);
            if (currMB->cbp & (1 << block8x8))
                itrans8x8(currMB, pl, x, y);
            else
                icopy8x8(currMB, pl, x, y);
        }
    }
}

void iTransform(mb_t *currMB, ColorPlane pl, int smb)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    StorablePicture *dec_picture = currSlice->dec_picture;
    imgpel **curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;

    if (smb || (currMB->cbp & 15)) {
        if (!currMB->transform_size_8x8_flag) // 4x4 inverse transform
            iMBtrans4x4(currMB, pl, smb); 
        else // 8x8 inverse transform
            iMBtrans8x8(currMB, pl);
        copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_rec[pl], currMB->pix_x, 0);
    } else {
        copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_pred[pl], currMB->pix_x, 0);
    }

    if (smb || (currMB->cbp & 15))
        currSlice->is_reset_coeff = FALSE;

    if (sps->chroma_format_idc == YUV400 || sps->chroma_format_idc == YUV444)
        return;

    void (*itrans_4x4)(struct macroblock_t *currMB, ColorPlane pl, int ioff, int joff);
    itrans_4x4 = (currMB->TransformBypassModeFlag == FALSE) ? itrans4x4 : itrans4x4_ls;

    for (int uv = 0; uv < 2; ++uv) {
        imgpel **curUV = &dec_picture->imgUV[uv][currMB->pix_c_y]; 
        imgpel **mb_rec  = currSlice->mb_rec [uv + 1];
        imgpel **mb_pred = currSlice->mb_pred[uv + 1];

        if (smb)
            itrans_sp_cr(currMB, uv);

        if (smb || (currMB->cbp >> 4)) {
            for (int joff = 0; joff < sps->MbHeightC; joff += BLOCK_SIZE) {
                for (int ioff = 0; ioff < sps->MbWidthC ;ioff += BLOCK_SIZE)
                    itrans_4x4(currMB, (ColorPlane)(uv + 1), ioff, joff);
            }
            copy_image_data(curUV, mb_rec, currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
        } else
            copy_image_data(curUV, mb_pred, currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
    }

    if (smb || (currMB->cbp >> 4))
        currSlice->is_reset_coeff_cr = FALSE;
}


void transform_t::inverse_transform_4x4(mb_t* mb)
{
}

void transform_t::inverse_transform_8x8(mb_t* mb)
{
}

void transform_t::inverse_transform_16x16(mb_t* mb)
{
}
