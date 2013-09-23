/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : transform.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

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


#include <functional>


using namespace vio::h264;


namespace vio  {
namespace h264 {


#define Q_BITS 15

static inline int isign(int x)
{
    return ( (x > 0) - (x < 0));
}

static inline int isignab(int a, int b)
{
    return ((b) < 0) ? -abs(a) : abs(a);
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

static const uint8_t QP_SCALE_CR[52] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35,
    35, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 39
};

static const int dequant_coef[6][4][4] = {
    {{ 10, 13, 10, 13 },
     { 13, 16, 13, 16 },
     { 10, 13, 10, 13 },
     { 13, 16, 13, 16 }},
    {{ 11, 14, 11, 14 },
     { 14, 18, 14, 18 },
     { 11, 14, 11, 14 },
     { 14, 18, 14, 18 }},
    {{ 13, 16, 13, 16 },
     { 16, 20, 16, 20 },
     { 13, 16, 13, 16 },
     { 16, 20, 16, 20 }},
    {{ 14, 18, 14, 18 },
     { 18, 23, 18, 23 },
     { 14, 18, 14, 18 },
     { 18, 23, 18, 23 }},
    {{ 16, 20, 16, 20 },
     { 20, 25, 20, 25 },
     { 16, 20, 16, 20 },
     { 20, 25, 20, 25 }},
    {{ 18, 23, 18, 23 },
     { 23, 29, 23, 29 },
     { 18, 23, 18, 23 },
     { 23, 29, 23, 29 }}
};

static const int quant_coef[6][4][4] = {
    {{ 13107,  8066, 13107,  8066 },
     {  8066,  5243,  8066,  5243 },
     { 13107,  8066, 13107,  8066 },
     {  8066,  5243,  8066,  5243 }},
    {{ 11916,  7490, 11916,  7490 },
     {  7490,  4660,  7490,  4660 },
     { 11916,  7490, 11916,  7490 },
     {  7490,  4660,  7490,  4660 }},
    {{ 10082,  6554, 10082,  6554 },
     {  6554,  4194,  6554,  4194 },
     { 10082,  6554, 10082,  6554 },
     {  6554,  4194,  6554,  4194 }},
    {{  9362,  5825,  9362,  5825 },
     {  5825,  3647,  5825,  3647 },
     {  9362,  5825,  9362,  5825 },
     {  5825,  3647,  5825,  3647 }},
    {{  8192,  5243,  8192,  5243 },
     {  5243,  3355,  5243,  3355 },
     {  8192,  5243,  8192,  5243 },
     {  5243,  3355,  5243,  3355 }},
    {{  7282,  4559,  7282,  4559 },
     {  4559,  2893,  4559,  2893 },
     {  7282,  4559,  7282,  4559 },
     {  4559,  2893,  4559,  2893 }}
};


static void copy_image_data_4x4(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2)
{
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1   + off1), (*imgBuf2   + off2), BLOCK_SIZE * sizeof (imgpel));
}

static void copy_image_data_8x8(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2)
{  
    for (int j = 0; j < BLOCK_SIZE_8x8; j += 4) {
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
        memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    }
}

static void copy_image_data_16x16(imgpel **imgBuf1, imgpel **imgBuf2, int off1, int off2)
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

void transform_t::inverse_4x4(int **d, int **r, int pos_y, int pos_x)
{
    int f[4][4];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 4; ++i) {
        int di0 = d[pos_y + i][pos_x + 0];
        int di1 = d[pos_y + i][pos_x + 1];
        int di2 = d[pos_y + i][pos_x + 2];
        int di3 = d[pos_y + i][pos_x + 3];

        int ei0 = di0 + di2;
        int ei1 = di0 - di2;
        int ei2 = (di1 >> 1) - di3;
        int ei3 = di1 + (di3 >> 1);

        f[i][0] = ei0 + ei3;
        f[i][1] = ei1 + ei2;
        f[i][2] = ei1 - ei2;
        f[i][3] = ei0 - ei3;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 4; ++j) {
        int f0j = f[0][j];
        int f1j = f[1][j];
        int f2j = f[2][j];
        int f3j = f[3][j];

        int g0j = f0j + f2j;
        int g1j = f0j - f2j;
        int g2j = (f1j >> 1) - f3j;
        int g3j = f1j + (f3j >> 1);

        int h0j = g0j + g3j;
        int h1j = g1j + g2j;
        int h2j = g1j - g2j;
        int h3j = g0j - g3j;

        r[pos_y + 0][pos_x + j] = (h0j + (1 << 5)) >> 6;
        r[pos_y + 1][pos_x + j] = (h1j + (1 << 5)) >> 6;
        r[pos_y + 2][pos_x + j] = (h2j + (1 << 5)) >> 6;
        r[pos_y + 3][pos_x + j] = (h3j + (1 << 5)) >> 6;
    }
}

void transform_t::inverse_8x8(int **d, int **r, int pos_y, int pos_x)
{
    int g[8][8];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 8; ++i) {
        int di0 = d[pos_y + i][pos_x + 0];
        int di1 = d[pos_y + i][pos_x + 1];
        int di2 = d[pos_y + i][pos_x + 2];
        int di3 = d[pos_y + i][pos_x + 3];
        int di4 = d[pos_y + i][pos_x + 4];
        int di5 = d[pos_y + i][pos_x + 5];
        int di6 = d[pos_y + i][pos_x + 6];
        int di7 = d[pos_y + i][pos_x + 7];

        int ei0 = di0 + di4;
        int ei1 = -di3 + di5 - di7 - (di7 >> 1);
        int ei2 = di0 - di4;
        int ei3 = di1 + di7 - di3 - (di3 >> 1);
        int ei4 = (di2 >> 1) - di6;
        int ei5 = -di1 + di7 + di5 + (di5 >> 1);
        int ei6 = di2 + (di6 >> 1);
        int ei7 = di3 + di5 + di1 + (di1 >> 1);

        int fi0 = ei0 + ei6;
        int fi1 = ei1 + (ei7 >> 2);
        int fi2 = ei2 + ei4;
        int fi3 = ei3 + (ei5 >> 2);
        int fi4 = ei2 - ei4;
        int fi5 = (ei3 >> 2) - ei5;
        int fi6 = ei0 - ei6;
        int fi7 = ei7 - (ei1 >> 2);

        g[i][0] = fi0 + fi7;
        g[i][1] = fi2 + fi5;
        g[i][2] = fi4 + fi3;
        g[i][3] = fi6 + fi1;
        g[i][4] = fi6 - fi1;
        g[i][5] = fi4 - fi3;
        g[i][6] = fi2 - fi5;
        g[i][7] = fi0 - fi7;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 8; ++j) {
        int g0j = g[0][j];
        int g1j = g[1][j];
        int g2j = g[2][j];
        int g3j = g[3][j];
        int g4j = g[4][j];
        int g5j = g[5][j];
        int g6j = g[6][j];
        int g7j = g[7][j];

        int h0j = g0j + g4j;
        int h1j = -g3j + g5j - g7j - (g7j >> 1);
        int h2j = g0j - g4j;
        int h3j = g1j + g7j - g3j - (g3j >> 1);
        int h4j = (g2j >> 1) - g6j;
        int h5j = -g1j + g7j + g5j + (g5j >> 1);
        int h6j = g2j + (g6j >> 1);
        int h7j = g3j + g5j + g1j + (g1j >> 1);

        int k0j = h0j + h6j;
        int k1j = h1j + (h7j >> 2);
        int k2j = h2j + h4j;
        int k3j = h3j + (h5j >> 2);
        int k4j = h2j - h4j;
        int k5j = (h3j >> 2) - h5j;
        int k6j = h0j - h6j;
        int k7j = h7j - (h1j >> 2);

        int m0j = k0j + k7j;
        int m1j = k2j + k5j;
        int m2j = k4j + k3j;
        int m3j = k6j + k1j;
        int m4j = k6j - k1j;
        int m5j = k4j - k3j;
        int m6j = k2j - k5j;
        int m7j = k0j - k7j;

        r[pos_y + 0][pos_x + j] = (m0j + (1 << 5)) >> 6;
        r[pos_y + 1][pos_x + j] = (m1j + (1 << 5)) >> 6;
        r[pos_y + 2][pos_x + j] = (m2j + (1 << 5)) >> 6;
        r[pos_y + 3][pos_x + j] = (m3j + (1 << 5)) >> 6;
        r[pos_y + 4][pos_x + j] = (m4j + (1 << 5)) >> 6;
        r[pos_y + 5][pos_x + j] = (m5j + (1 << 5)) >> 6;
        r[pos_y + 6][pos_x + j] = (m6j + (1 << 5)) >> 6;
        r[pos_y + 7][pos_x + j] = (m7j + (1 << 5)) >> 6;
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


void transform_t::bypass_4x4(int** r, int** f, int ioff, int joff, uint8_t pred_mode)
{
    if (pred_mode == intra_prediction_t::Intra_4x4_Vertical) {
        for (int i = 0; i < 4; ++i) {
            f[joff + 0][ioff + i] = r[joff + 0][ioff + i];
            for (int j = 1; j < 4; ++j)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j - 1][ioff + i];
        }
    } else if (pred_mode == intra_prediction_t::Intra_4x4_Horizontal) {
        for (int j = 0; j < 4; ++j) {
            f[joff + j][ioff + 0] = r[joff + j][ioff + 0];
            for (int i = 1; i < 4; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j][ioff + i - 1];
        }
    } else {
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i];
        }
    }
}

void transform_t::bypass_8x8(int** r, int** f, int ioff, int joff, uint8_t pred_mode)
{
    if (pred_mode == intra_prediction_t::Intra_8x8_Vertical) {
        for (int i = 0; i < 8; ++i) {
            f[joff + 0][ioff + i] = r[joff + 0][ioff + i];
            for (int j = 1; j < 8; ++j)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j - 1][ioff + i];
        }
    } else if (pred_mode == intra_prediction_t::Intra_8x8_Horizontal) {
        for (int j = 0; j < 8; ++j) {
            f[joff + j][ioff + 0] = r[joff + j][ioff + 0];
            for (int i = 1; i < 8; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j][ioff + i - 1];
        }
    } else {
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 8; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i];
        }
    }
}

void transform_t::bypass_16x16(int** r, int** f, int ioff, int joff, uint8_t pred_mode)
{
    if (pred_mode == intra_prediction_t::Intra_16x16_Vertical) {
        for (int i = 0; i < 16; ++i) {
            f[0][i] = r[0][i];
            for (int j = 1; j < 16; ++j)
                f[j][i] = r[j][i] + f[j - 1][i];
        }
    } else if (pred_mode == intra_prediction_t::Intra_16x16_Horizontal) {
        for (int j = 0; j < 16; ++j) {
            f[j][0] = r[j][0];
            for (int i = 1; i < 16; ++i)
                f[j][i] = r[j][i] + f[j][i - 1];
        }
    } else {
        for (int j = 0; j < 16; ++j) {
            for (int i = 0; i < 16; ++i)
                f[j][i] = r[j][i];
        }
    }
}

void transform_t::bypass_chroma(int** r, int** f, int nW, int nH, uint8_t pred_mode)
{
    if (pred_mode == intra_prediction_t::Intra_Chroma_Vertical) {
        for (int i = 0; i < nW; ++i) {
            f[0][i] = r[0][i];
            for (int j = 1; j < nH; ++j)
                f[j][i] = r[j][i] + f[j - 1][i];
        }
    } else if (pred_mode == intra_prediction_t::Intra_Chroma_Horizontal) {
        for (int j = 0; j < nH; ++j) {
            f[j][0] = r[j][0];
            for (int i = 1; i < nW; ++i)
                f[j][i] = r[j][i] + f[j][i - 1];
        }
    } else {
        for (int j = 0; j < nH; ++j) {
            for (int i = 0; i < nW; ++i)
                f[j][i] = r[j][i];
        }
    }
}


void transform_t::inverse_luma_dc(mb_t* mb, ColorPlane pl)
{
    if (!mb->TransformBypassModeFlag) {
        slice_t* slice = mb->p_Slice;
        sps_t* sps = slice->active_sps;

        int transform_pl = sps->separate_colour_plane_flag ? PLANE_Y : pl;
        int **cof = slice->cof[transform_pl];

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

        slice->quantization.inverse_itrans_2(mb, pl, M4);

        free_mem2Dint(M4);
    }
}

void transform_t::inverse_chroma_dc(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    bool smb = (slice->slice_type == SP_slice && !mb->is_intra_block) ||
               (slice->slice_type == SI_slice && mb->mb_type == SI4MB);
    int **cof = slice->cof[pl];

    if (!mb->TransformBypassModeFlag) {
        if (sps->ChromaArrayType == 1 && !smb) {
            int M4[4];
            M4[0] = cof[0][0];
            M4[1] = cof[0][4];
            M4[2] = cof[4][0];
            M4[3] = cof[4][4];

            ihadamard2x2(M4, M4);

            slice->quantization.inverse_itrans_420(mb, pl, M4);
        }
        if (sps->ChromaArrayType == 2) {
            int **M4;
            get_mem2Dint(&M4, BLOCK_SIZE, 2);

            for (int j = 0; j < 4; j++) {
                M4[j][0] = cof[j << 2][0];
                M4[j][1] = cof[j << 2][4];
            }

            ihadamard2x4(M4, M4);

            slice->quantization.inverse_itrans_422(mb, pl, M4);

            free_mem2Dint(M4);
        }
    }
}

void transform_t::inverse_transform_4x4(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    int block8x8 = (joff / 8) * 2 + (ioff / 8);

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        int i4x4 = ((joff / 4) / 2) * 8 + ((joff / 4) % 2) * 2 +
                   ((ioff / 4) / 2) * 4 + ((ioff / 4) % 2);
        uint8_t pred_mode = mb->Intra4x4PredMode[i4x4];
        if (mb->TransformBypassModeFlag)
            this->bypass_4x4(slice->cof[pl], slice->mb_rres[pl], ioff, joff, pred_mode);
        else
            this->inverse_4x4(slice->cof[pl], slice->mb_rres[pl], joff, ioff);
    }


    sps_t *sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int    **mb_rres = slice->mb_rres[pl];
    imgpel **mb_pred = slice->mb_pred[pl];
    imgpel **mb_rec  = slice->mb_rec [pl];

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        for (int j = joff; j < joff + 4; j++) {
            for (int i = ioff; i < ioff + 4; i++)
                mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_rres[j][i] + mb_pred[j][i]);
        }
    } else {
        for (int j = 0; j < 4; j++)
            memcpy(&mb_rec[joff + j][ioff], &mb_pred[joff + j][ioff], 4 * sizeof(imgpel));
    }

    copy_image_data_4x4(&curr_img[mb->mb.y * 16 + joff], &mb_rec[joff], mb->mb.x * 16 + ioff, ioff);
}

void transform_t::inverse_transform_8x8(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    int block8x8 = (joff / 8) * 2 + (ioff / 8);

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        uint8_t pred_mode = mb->Intra8x8PredMode[joff/8 * 2 + ioff/8];
        if (mb->TransformBypassModeFlag)
            this->bypass_8x8(slice->cof[pl], slice->mb_rres[pl], ioff, joff, pred_mode);
        else
            this->inverse_8x8(slice->cof[pl], slice->mb_rres[pl], joff, ioff);
    }


    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int    **mb_rres = slice->mb_rres[pl];
    imgpel **mb_pred = slice->mb_pred[pl];
    imgpel **mb_rec  = slice->mb_rec [pl];

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        for (int j = joff; j < joff + 8; j++) {
            for (int i = ioff; i < ioff + 8; i++)
                mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_pred[j][i] + mb_rres[j][i]);
        }
    } else {
        for (int j = 0; j < 8; j++)
            memcpy(&mb_rec[joff + j][ioff], &mb_pred[joff + j][ioff], 8 * sizeof(imgpel));
    }

    copy_image_data_8x8(&curr_img[mb->mb.y * 16 + joff], &mb_rec[joff], mb->mb.x * 16 + ioff, ioff);
}

void transform_t::inverse_transform_16x16(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;

    uint8_t pred_mode = mb->Intra16x16PredMode;
    if (mb->TransformBypassModeFlag)
        this->bypass_16x16(slice->cof[pl], slice->mb_rres[pl], ioff, joff, pred_mode);
    else {
        for (int j = 0; j < 16; j += 4) {
            for (int i = 0; i < 16; i += 4)
                this->inverse_4x4(slice->cof[pl], slice->mb_rres[pl], j, i);
        }
    }

    sps_t *sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int    **mb_rres = slice->mb_rres[pl];
    imgpel **mb_pred = slice->mb_pred[pl];
    imgpel **mb_rec  = slice->mb_rec [pl];

    for (int j = 0; j < 16; j++) {
        for (int i = 0; i < 16; i++)
            mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_rres[j][i] + mb_pred[j][i]);
    }

    copy_image_data_16x16(&curr_img[mb->mb.y * 16 + joff], &mb_rec[joff], mb->mb.x * 16 + ioff, ioff);
}

void transform_t::inverse_transform_chroma(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    uint8_t pred_mode = mb->intra_chroma_pred_mode;
    if (mb->TransformBypassModeFlag)
        this->bypass_chroma(slice->cof[pl], slice->mb_rres[pl], sps->MbWidthC, sps->MbHeightC, pred_mode);
    else {
        for (int joff = 0; joff < sps->MbHeightC; joff += 4) {
            for (int ioff = 0; ioff < sps->MbWidthC; ioff += 4)
                this->inverse_4x4(slice->cof[pl], slice->mb_rres[pl], joff, ioff);
        }
    }

    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int    **mb_rres = slice->mb_rres[pl];
    imgpel **mb_pred = slice->mb_pred[pl];
    imgpel **mb_rec  = slice->mb_rec [pl];

    for (int j = 0; j < sps->MbHeightC; j++) {
        for (int i = 0; i < sps->MbWidthC; i++)
            mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_rres[j][i] + mb_pred[j][i]);
    }

    for (int joff = 0; joff < sps->MbHeightC; joff += 4) {
        for (int ioff = 0; ioff < sps->MbWidthC; ioff += 4)
            copy_image_data_4x4(&curr_img[mb->mb.y * sps->MbHeightC + joff], &mb_rec[joff], mb->mb.x * sps->MbWidthC + ioff, ioff);
    }
}

void transform_t::inverse_transform_inter(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;

    if (mb->CodedBlockPatternLuma) {
        if (!mb->transform_size_8x8_flag) {
            for (int y = 0; y < 16; y += 4) {
                for (int x = 0; x < 16; x += 4)
                    this->inverse_transform_4x4(mb, pl, x, y);
            }
        } else {
            for (int y = 0; y < 16; y += 8) {
                for (int x = 0; x < 16; x += 8)
                    this->inverse_transform_8x8(mb, pl, x, y);
            }
        }
    } else
        copy_image_data_16x16(&curr_img[mb->mb.y * 16], slice->mb_pred[pl], mb->mb.x * 16, 0);

    if (mb->CodedBlockPatternLuma)
        slice->is_reset_coeff = false;

    if (sps->chroma_format_idc == YUV400 || sps->chroma_format_idc == YUV444)
        return;

    for (int uv = 0; uv < 2; ++uv) {
        imgpel **curUV = &dec_picture->imgUV[uv][mb->mb.y * sps->MbHeightC]; 
        imgpel **mb_pred = slice->mb_pred[uv + 1];

        if (mb->CodedBlockPatternChroma)
            this->inverse_transform_chroma(mb, (ColorPlane)(uv + 1));
        else
            copy_image_data(curUV, mb_pred, mb->mb.x * sps->MbWidthC, 0, sps->MbWidthC, sps->MbHeightC);
    }

    if (mb->CodedBlockPatternChroma)
        slice->is_reset_coeff_cr = false;
}


void transform_t::itrans_sp(mb_t *currMB, ColorPlane pl, int ioff, int joff)
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

    this->inverse_4x4(cof, mb_rres, joff, ioff);

    for (int j = joff; j < joff + BLOCK_SIZE; ++j) {
        mb_rec[j][ioff    ] = (imgpel)clip1(max_pel_value_comp, mb_rres[j][ioff    ]);
        mb_rec[j][ioff + 1] = (imgpel)clip1(max_pel_value_comp, mb_rres[j][ioff + 1]);
        mb_rec[j][ioff + 2] = (imgpel)clip1(max_pel_value_comp, mb_rres[j][ioff + 2]);
        mb_rec[j][ioff + 3] = (imgpel)clip1(max_pel_value_comp, mb_rres[j][ioff + 3]);
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

    int qp_per    = (currSlice->SliceQpY < 0 ? currSlice->SliceQpY : QP_SCALE_CR[currSlice->SliceQpY]) / 6;
    int qp_rem    = (currSlice->SliceQpY < 0 ? currSlice->SliceQpY : QP_SCALE_CR[currSlice->SliceQpY]) % 6;

    int qp_per_sp = (currSlice->QsY < 0 ? currSlice->QsY : QP_SCALE_CR[currSlice->QsY]) / 6;
    int qp_rem_sp = (currSlice->QsY < 0 ? currSlice->QsY : QP_SCALE_CR[currSlice->QsY]) % 6;
    int q_bits_sp = Q_BITS + qp_per_sp;  

    if (currSlice->slice_type == SI_SLICE) {
        qp_per = qp_per_sp;
        qp_rem = qp_rem_sp;
    }

    for (j = 0; j < sps->MbHeightC; ++j) {
        for (i = 0; i < sps->MbWidthC; ++i) {
            PBlock[j][i] = mb_pred[j][i];
            mb_pred[j][i] = 0;
        }
    }

    for (n2 = 0; n2 < sps->MbHeightC; n2 += BLOCK_SIZE) {
        for (n1 = 0; n1 < sps->MbWidthC; n1 += BLOCK_SIZE)
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

        for (n2 = 0; n2 < sps->MbHeightC; n2 += BLOCK_SIZE) {
            for (n1 = 0; n1 < sps->MbWidthC; n1 += BLOCK_SIZE) {
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

        for (n2 = 0; n2 < sps->MbHeightC; n2 += BLOCK_SIZE) {
            for (n1 = 0; n1 < sps->MbWidthC; n1 += BLOCK_SIZE) {
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

void transform_t::inverse_transform_sp(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;

    if (!mb->transform_size_8x8_flag) {
        for (int y = 0; y < 16; y += 4) {
            for (int x = 0; x < 16; x += 4)
                this->itrans_sp(mb, pl, x, y);
        }
    } else {
        for (int y = 0; y < 16; y += 8) {
            for (int x = 0; x < 16; x += 8)
                this->inverse_transform_8x8(mb, pl, x, y);
        }
    }
    copy_image_data_16x16(&curr_img[mb->mb.y * 16], slice->mb_rec[pl], mb->mb.x * 16, 0);

    slice->is_reset_coeff = false;

    if (sps->chroma_format_idc == YUV400 || sps->chroma_format_idc == YUV444)
        return;

    for (int uv = 0; uv < 2; ++uv) {
        itrans_sp_cr(mb, uv);
        this->inverse_transform_chroma(mb, (ColorPlane)(uv + 1));
    }

    slice->is_reset_coeff_cr = false;
}


}
}
