/*!
 ***********************************************************************
 * \file read_comp_cavlc.c
 *
 * \brief
 *     Read Coefficient Components (CAVLC version)
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

#include "global.h"
#include "slice.h"
#include "bitstream_elements.h"
#include "bitstream.h"
#include "macroblock.h"
#include "mb_read.h"
#include "transform.h"
#include "neighbour.h"

#include "mb_read_syntax.h"


#define IS_I16MB(MB)    ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)
#define IS_DIRECT(MB)   ((MB)->mb_type == 0 && (currSlice->slice_type == B_SLICE))

#define TOTRUN_NUM       15
#define RUNBEFORE_NUM_M1  6


static int readSyntaxElement_Level_VLC0(Bitstream *currStream)
{
    int level, sign;
    int len = 1;

    while (currStream->read_bits(1) == 0)
        len++;

    if (len < 15) {
        sign  = (len - 1) & 1;
        level = ((len - 1) >> 1) + 1;
    } else if (len == 15) {
        // escape code
        int code = 16 + currStream->read_bits(4);
        sign  = (code & 0x01);
        level = ((code >> 1) & 0x07) + 8;
        len  += 4;
    } else {
        // escape code
        int code = currStream->read_bits(len - 4);
        sign  = (code & 0x01);
        level = (code >> 1) + (2048 << (len - 16)) - 2032;
        len  += (len - 4);
    }

    return sign ? -level : level;
}

static int readSyntaxElement_Level_VLCN(Bitstream *currStream, int vlc)
{
    int level, sign;
    int len = 1;

    while (currStream->read_bits(1) == 0)
        len++;

    if (len < 16) {
        level = ((len - 1) << (vlc - 1)) + 1;
        // read (vlc-1) bits -> suffix
        if (vlc - 1 > 0)
            level += currStream->read_bits(vlc - 1);
        sign   = currStream->read_bits(1);
        len   += vlc;
    } else { // escape
        level  = (1 << (len - 5)) + (15 << (vlc - 1)) - 2047;
        level += currStream->read_bits(len - 5);
        sign   = currStream->read_bits(1);
        len   += (len - 4);
    }

    return sign ? -level : level;
}

static int readSyntaxElement_NumCoeffTrailingOnes(Bitstream *currStream, int *coeff, int *ones, int tab)
{
    static const byte lentab[3][4][17] = {
        {{ 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
         { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
         { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
         { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16}},
        {{ 2, 6, 6, 7, 8, 8, 9,11,11,12,12,12,13,13,13,14,14},
         { 0, 2, 5, 6, 6, 7, 8, 9,11,11,12,12,13,13,14,14,14},
         { 0, 0, 3, 6, 6, 7, 8, 9,11,11,12,12,13,13,13,14,14},
         { 0, 0, 0, 4, 4, 5, 6, 6, 7, 9,11,11,12,13,13,13,14}},
        {{ 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9,10,10,10,10},
         { 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10},
         { 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,10},
         { 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9,10,10,10}}
    };

    static const byte codtab[3][4][17] = {
        {{ 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7,4},
         { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10,6},
         { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9,5},
         { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12,8}},
        {{ 3,11, 7, 7, 7, 4, 7,15,11,15,11, 8,15,11, 7, 9,7},
         { 0, 2, 7,10, 6, 6, 6, 6,14,10,14,10,14,10,11, 8,6},
         { 0, 0, 3, 9, 5, 5, 5, 5,13, 9,13, 9,13, 9, 6,10,5},
         { 0, 0, 0, 5, 4, 6, 8, 4, 4, 4,12, 8,12,12, 8, 1,4}},
        {{15,15,11, 8,15,11, 9, 8,15,11,15,11, 8,13, 9, 5,1},
         { 0,14,15,12,10, 8,14,10,14,14,10,14,10, 7,12, 8,4},
         { 0, 0,13,14,11, 9,13, 9,13,10,13, 9,13, 9,11, 7,3},
         { 0, 0, 0,12,11,10, 9, 8,13,12,12,12, 8,12,10, 6,2}}
    };

    // tab is the index of Table used to code coeff_token
    // tab==3 means (8<=nC) which uses 6bit FLC
    if (tab == 3) {
        int code = currStream->read_bits(6);
        *coeff = (code >> 2);
        *ones  = (code & 3);
        if (*coeff == 0 && *ones == 3)
            // #c = 0, #t1 = 3 =>  #c = 0
            *ones = 0;
        else
            (*coeff)++;
        return 0;
    }

    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 17; i++) {
            int len = lentab[tab][j][i];
            int cod = codtab[tab][j][i];
            if (len > 0 && currStream->next_bits(len) == cod) {
                *coeff = i;
                *ones  = j;
                currStream->read_bits(len);
                return 0;
            }
        }
    }

    printf("ERROR: failed to find NumCoeff/TrailingOnes\n");
    exit(-1);
    return -1;
}

static int readSyntaxElement_NumCoeffTrailingOnesChromaDC(Bitstream *currStream, int *coeff, int *ones, int yuv)
{
    static const byte lentab[3][4][17] = {
        //YUV420
        {{ 2, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 1, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 3, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
        //YUV422
        {{ 1, 7, 7, 9, 9,10,11,12,13, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 2, 7, 7, 9,10,11,12,12, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 3, 7, 7, 9,10,11,12, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 0, 5, 6, 7, 7,10,11, 0, 0, 0, 0, 0, 0, 0, 0}},
        //YUV444
        {{ 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
         { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
         { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
         { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16}}
    };

    static const byte codtab[3][4][17] = {
        //YUV420
        {{ 1, 7, 4, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 1, 6, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
        //YUV422
        {{ 1,15,14, 7, 6, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 1,13,12, 5, 6, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 1,11,10, 4, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0},
         { 0, 0, 0, 1, 1, 9, 8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0}},
        //YUV444
        {{ 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7, 4},
         { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10, 6},
         { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9, 5},
         { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12, 8}}
    };

    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 17; i++) {
            int len = lentab[yuv][j][i];
            int cod = codtab[yuv][j][i];
            if (len > 0 && currStream->next_bits(len) == cod) {
                *coeff = i;
                *ones  = j;
                currStream->read_bits(len);
                return 0;
            }
        }
    }

    printf("ERROR: failed to find NumCoeff/TrailingOnes ChromaDC\n");
    exit(-1);
    return -1;
}

static int readSyntaxElement_TotalZeros(Bitstream *currStream, int *coeff, int tab)
{
    static const byte lentab[TOTRUN_NUM][16] = {
        { 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},
        { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},
        { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},
        { 5,3,4,4,3,3,3,4,3,4,5,5,5},
        { 4,4,4,3,3,3,3,3,4,5,4,5},
        { 6,5,3,3,3,3,3,3,4,3,6},
        { 6,5,3,3,3,2,3,4,3,6},
        { 6,4,5,3,2,2,3,3,6},
        { 6,6,4,2,2,3,2,5},
        { 5,5,3,2,2,2,4},
        { 4,4,3,3,1,3},
        { 4,4,2,1,3},
        { 3,3,1,2},
        { 2,2,1},
        { 1,1},
    };

    static const byte codtab[TOTRUN_NUM][16] = {
        {1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
        {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
        {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
        {3,7,5,4,6,5,4,3,3,2,2,1,0},
        {5,4,3,7,6,5,4,3,2,1,1,0},
        {1,1,7,6,5,4,3,2,1,1,0},
        {1,1,5,4,3,3,2,1,1,0},
        {1,1,1,3,3,2,2,1,0},
        {1,0,1,3,2,1,1,1,},
        {1,0,1,3,2,1,1,},
        {0,1,1,2,1,3},
        {0,1,1,1,1},
        {0,1,1,1},
        {0,1,1},
        {0,1},
    };

    for (int i = 0; i < 16; i++) {
        int len = lentab[tab][i];
        int cod = codtab[tab][i];
        if (len > 0 && currStream->next_bits(len) == cod) {
            *coeff = i;
            currStream->read_bits(len);
            return 0;
        }
    }

    printf("ERROR: failed to find Total Zeros !cdc\n");
    exit(-1);
    return -1;
}

static int readSyntaxElement_TotalZerosChromaDC(Bitstream *currStream, int *coeff, int yuv, int tab)
{
    static const byte lentab[3][TOTRUN_NUM][16] = {
        //YUV420
        {{ 1,2,3,3},
         { 1,2,2},
         { 1,1}},
        //YUV422
        {{ 1,3,3,4,4,4,5,5},
         { 3,2,3,3,3,3,3},
         { 3,3,2,2,3,3},
         { 3,2,2,2,3},
         { 2,2,2,2},
         { 2,2,1},
         { 1,1}},
        //YUV444
        {{ 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},
         { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},
         { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},
         { 5,3,4,4,3,3,3,4,3,4,5,5,5},
         { 4,4,4,3,3,3,3,3,4,5,4,5},
         { 6,5,3,3,3,3,3,3,4,3,6},
         { 6,5,3,3,3,2,3,4,3,6},
         { 6,4,5,3,2,2,3,3,6},
         { 6,6,4,2,2,3,2,5},
         { 5,5,3,2,2,2,4},
         { 4,4,3,3,1,3},
         { 4,4,2,1,3},
         { 3,3,1,2},
         { 2,2,1},
         { 1,1}}
    };

    static const byte codtab[3][TOTRUN_NUM][16] = {
        //YUV420
        {{ 1,1,1,0},
         { 1,1,0},
         { 1,0}},
        //YUV422
        {{ 1,2,3,2,3,1,1,0},
         { 0,1,1,4,5,6,7},
         { 0,1,1,2,6,7},
         { 6,0,1,2,7},
         { 0,1,2,3},
         { 0,1,1},
         { 0,1}},
        //YUV444
        {{1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
         {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
         {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
         {3,7,5,4,6,5,4,3,3,2,2,1,0},
         {5,4,3,7,6,5,4,3,2,1,1,0},
         {1,1,7,6,5,4,3,2,1,1,0},
         {1,1,5,4,3,3,2,1,1,0},
         {1,1,1,3,3,2,2,1,0},
         {1,0,1,3,2,1,1,1,},
         {1,0,1,3,2,1,1,},
         {0,1,1,2,1,3},
         {0,1,1,1,1},
         {0,1,1,1},
         {0,1,1},
         {0,1}}
    };

    for (int i = 0; i < 16; i++) {
        int len = lentab[yuv][tab][i];
        int cod = codtab[yuv][tab][i];
        if (len > 0 && currStream->next_bits(len) == cod) {
            *coeff = i;
            currStream->read_bits(len);
            return 0;
        }
    }

    printf("ERROR: failed to find Total Zeros\n");
    exit(-1);
    return -1;
}

static int readSyntaxElement_Run(Bitstream *currStream, int *coeff, int tab)
{
    static const byte lentab[TOTRUN_NUM][16] = {
        {1,1},
        {1,2,2},
        {2,2,2,2},
        {2,2,2,3,3},
        {2,2,3,3,3,3},
        {2,3,3,3,3,3,3},
        {3,3,3,3,3,3,3,4,5,6,7,8,9,10,11},
    };

    static const byte codtab[TOTRUN_NUM][16] = {
        {1,0},
        {1,1,0},
        {3,2,1,0},
        {3,2,1,1,0},
        {3,2,3,2,1,0},
        {3,0,1,3,2,5,4},
        {7,6,5,4,3,2,1,1,1,1,1,1,1,1,1},
    };

    for (int i = 0; i < 16; i++) {
        int len = lentab[tab][i];
        int cod = codtab[tab][i];
        if (len > 0 && currStream->next_bits(len) == cod) {
            *coeff = i;
            currStream->read_bits(len);
            return 0;
        }
    }

    printf("ERROR: failed to find Run\n");
    exit(-1);
    return -1;
}


static void read_coeff_4x4_CAVLC(mb_t *currMB, int block_type, int i, int j,
                                 int levarr[16], int runarr[16], int *number_coefficients)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;
    int mb_nr = currMB->mbAddrX;

    int k, code, vlcnum;
    int numcoeff = 0, numtrailingones;
    int level_two_or_higher;
    int numones, totzeros, abslevel, cdc=0, cac=0;
    int zerosleft, ntr, dptype = 0;
    int max_coeff_num = 0, nnz;
    static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768};    // maximum vlc = 6

    int num_cdc_coeff;
    if (sps->chroma_format_idc != YUV400)
        num_cdc_coeff = (((1 << sps->chroma_format_idc) & (~0x1)) << 1);
    else
        num_cdc_coeff = 0;

    int pl;
    switch (block_type) {
    case LUMA:
    case CB:
    case CR:
        max_coeff_num = 16;
        dptype = currMB->is_intra_block ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
        pl = block_type == LUMA ? 0 : block_type == CB ? 1 : 2;
        p_Vid->nz_coeff[mb_nr][pl][j][i] = 0;
        break;
    case LUMA_INTRA16x16DC:
    case CB_INTRA16x16DC:
    case CR_INTRA16x16DC:
        max_coeff_num = 16;
        dptype = SE_LUM_DC_INTRA;
        pl = block_type == LUMA_INTRA16x16DC ? 0 : block_type == CB_INTRA16x16DC ? 1 : 2;
        p_Vid->nz_coeff[mb_nr][pl][j][i] = 0; 
        break;
    case LUMA_INTRA16x16AC:
    case CB_INTRA16x16AC:
    case CR_INTRA16x16AC:
        max_coeff_num = 15;
        dptype = SE_LUM_AC_INTRA;
        pl = block_type == LUMA_INTRA16x16AC ? 0 : block_type == CB_INTRA16x16AC ? 1 : 2;
        p_Vid->nz_coeff[mb_nr][pl][j][i] = 0; 
        break;
    case CHROMA_DC:
        max_coeff_num = num_cdc_coeff;
        cdc = 1;
        dptype = currMB->is_intra_block ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
        p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
        break;
    case CHROMA_AC:
        max_coeff_num = 15;
        cac = 1;
        dptype = currMB->is_intra_block ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
        p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
        break;
    default:
        error ("read_coeff_4x4_CAVLC: invalid block type", 600);
        p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
        break;
    }

    DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][dptype]];
    Bitstream *currStream = dP->bitstream;

    if (!cdc) {
        // luma or chroma AC    
        if (block_type == LUMA ||
            block_type == LUMA_INTRA16x16DC || block_type == LUMA_INTRA16x16AC ||
            block_type == CHROMA_AC) {
            nnz = (!cac) ? predict_nnz(currMB, LUMA, i<<2, j<<2) :
                           predict_nnz_chroma(currMB, i, ((j-4)<<2));
        } else if (block_type == CB || block_type == CB_INTRA16x16DC || block_type == CB_INTRA16x16AC)
            nnz = predict_nnz(currMB, CB, i<<2, j<<2);
        else
            nnz = predict_nnz(currMB, CR, i<<2, j<<2);

        int tab = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));
        readSyntaxElement_NumCoeffTrailingOnes(currStream, &numcoeff, &numtrailingones, tab);

        if (block_type == LUMA ||
            block_type == LUMA_INTRA16x16DC || block_type == LUMA_INTRA16x16AC ||
            block_type == CHROMA_AC)
            p_Vid->nz_coeff[mb_nr][0][j][i] = (byte) numcoeff;
        else if (block_type == CB || block_type == CB_INTRA16x16DC || block_type == CB_INTRA16x16AC)
            p_Vid->nz_coeff[mb_nr][1][j][i] = (byte) numcoeff;
        else
            p_Vid->nz_coeff[mb_nr][2][j][i] = (byte) numcoeff;
    } else {
        // chroma DC
        int yuv = p_Vid->active_sps->chroma_format_idc - 1;
        readSyntaxElement_NumCoeffTrailingOnesChromaDC(currStream, &numcoeff, &numtrailingones, yuv);
    }

    memset(levarr, 0, max_coeff_num * sizeof(int));
    memset(runarr, 0, max_coeff_num * sizeof(int));

    numones = numtrailingones;
    *number_coefficients = numcoeff;

    if (numcoeff) {
        if (numtrailingones) {
            code = currStream->f(numtrailingones);
            ntr = numtrailingones;
            for (k = numcoeff - 1; k > numcoeff - 1 - numtrailingones; k--) {
                ntr--;
                levarr[k] = (code>>ntr)&1 ? -1 : 1;
            }
        }

        // decode levels
        level_two_or_higher = (numcoeff > 3 && numtrailingones == 3 ? 0 : 1);
        vlcnum              = (numcoeff > 10 && numtrailingones < 3 ? 1 : 0);

        for (k = numcoeff - 1 - numtrailingones; k >= 0; k--) {
            if (vlcnum == 0)
                levarr[k] = readSyntaxElement_Level_VLC0(currStream);
            else
                levarr[k] = readSyntaxElement_Level_VLCN(currStream, vlcnum);

            if (level_two_or_higher) {
                levarr[k] += (levarr[k] > 0 ? 1 : -1);
                level_two_or_higher = 0;
            }

            abslevel = iabs(levarr[k]);
            if (abslevel == 1)
                ++numones;

            // update VLC table
            if (abslevel > incVlc[vlcnum])
                ++vlcnum;

            if (k == numcoeff - 1 - numtrailingones && abslevel > 3)
                vlcnum = 2;      
        }

        if (numcoeff < max_coeff_num) {
            // decode total run
            vlcnum = numcoeff - 1;
            int yuv = p_Vid->active_sps->chroma_format_idc - 1;
            if (cdc)
                readSyntaxElement_TotalZerosChromaDC(currStream, &totzeros, yuv, vlcnum);
            else
                readSyntaxElement_TotalZeros(currStream, &totzeros, vlcnum);
        } else
            totzeros = 0;

        zerosleft = totzeros;
        i = numcoeff - 1;

        while (zerosleft > 0 && i > 0) {
            vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);
            readSyntaxElement_Run(currStream, &runarr[i], vlcnum);
            zerosleft -= runarr[i];
            i--;
        }
        runarr[i] = zerosleft;    
    }
}


static void read_tc_luma(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    const byte (*pos_scan8x8)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN8x8 : FIELD_SCAN8x8;

    if (IS_I16MB(currMB) && !currMB->dpl_flag) {
        int block_type = pl == PLANE_Y ? LUMA_INTRA16x16DC :
                         pl == PLANE_U ? CB_INTRA16x16DC : CR_INTRA16x16DC;
        int levarr[16], runarr[16], numcoeff;
        read_coeff_4x4_CAVLC(currMB, block_type, 0, 0, levarr, runarr, &numcoeff);

        int coef_ctr = -1;
        for (int k = 0; k < numcoeff; ++k) {
            if (levarr[k] != 0) {
                coef_ctr += runarr[k] + 1;
                int i0 = pos_scan4x4[coef_ctr][0];
                int j0 = pos_scan4x4[coef_ctr][1];
                currSlice->cof[pl][j0 << 2][i0 << 2] = levarr[k];
            }
        }

        if (!currMB->TransformBypassModeFlag) {
            int transform_pl = sps->separate_colour_plane_flag ? currSlice->colour_plane_id : pl;
            itrans_2(currMB, (ColorPlane)transform_pl);
        }
    }

    int qp_per = currMB->qp_scaled[pl] / 6;
    int qp_rem = currMB->qp_scaled[pl] % 6;
    int transform_pl = sps->separate_colour_plane_flag ? currSlice->colour_plane_id : pl;
    int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
        currSlice->InvLevelScale4x4_Intra[transform_pl][qp_rem] :
        currSlice->InvLevelScale4x4_Inter[transform_pl][qp_rem];
    int (*InvLevelScale8x8)[8] = currMB->is_intra_block ?
        currSlice->InvLevelScale8x8_Intra[transform_pl][qp_rem] :
        currSlice->InvLevelScale8x8_Inter[transform_pl][qp_rem];

    int start_scan = IS_I16MB(currMB) ? 1 : 0;
    byte **nzcoeff = currSlice->p_Vid->nz_coeff[currMB->mbAddrX][pl];

    int cur_context; 
    if (IS_I16MB(currMB)) {
        if (pl == PLANE_Y)
            cur_context = LUMA_INTRA16x16AC;
        else if (pl == PLANE_U)
            cur_context = CB_INTRA16x16AC;
        else
            cur_context = CR_INTRA16x16AC;
    } else {
        if (pl == PLANE_Y)
            cur_context = LUMA;
        else if (pl == PLANE_U)
            cur_context = CB;
        else
            cur_context = CR;
    }

    for (int block_y = 0; block_y < 4; block_y += 2) {
        for (int block_x = 0; block_x < 4; block_x += 2) {
            int b8 = 2*(block_y>>1) + (block_x>>1);

            if (currMB->cbp & (1 << b8)) {
                for (int j = block_y; j < block_y+2; ++j) {
                    for (int i = block_x; i < block_x+2; ++i) {
                        int levarr[16] = {0}, runarr[16] = {0}, numcoeff;
                        read_coeff_4x4_CAVLC(currMB, cur_context, i, j, levarr, runarr, &numcoeff);

                        int coef_ctr = start_scan - 1;

                        for (int k = 0; k < numcoeff; ++k) {
                            if (levarr[k] != 0) {
                                coef_ctr += runarr[k]+1;

                                if (!currMB->transform_size_8x8_flag) {
                                    currMB->s_cbp[pl].blk |= i64_power2((j<<2) + i);
                                    int i0 = pos_scan4x4[coef_ctr][0];
                                    int j0 = pos_scan4x4[coef_ctr][1];

                                    if (!currMB->TransformBypassModeFlag)
                                        currSlice->cof[pl][(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0]) << qp_per, 4);
                                    else
                                        currSlice->cof[pl][(j<<2) + j0][(i<<2) + i0] = levarr[k];
                                } else {
                                    // do same as CABAC for deblocking: any coeff in the 8x8 marks all the 4x4s
                                    //as containing coefficients
                                    currMB->s_cbp[pl].blk |= 51 << ((block_y << 2) + block_x);
                                    int b4 = 4 * coef_ctr + 2 * (j - block_y) + (i - block_x);
                                    int i0 = pos_scan8x8[b4][0];
                                    int j0 = pos_scan8x8[b4][1];

                                    if (!currMB->TransformBypassModeFlag)
                                        currSlice->cof[pl][block_y*4 +j0][block_x*4 +i0] = rshift_rnd_sf((levarr[k] * InvLevelScale8x8[j0][i0]) << qp_per, 6);
                                    else
                                        currSlice->cof[pl][block_y*4 +j0][block_x*4 +i0] = levarr[k];
                                }
                            }
                        }
                    }
                }
            } else {
                nzcoeff[block_y    ][block_x    ] = 0;
                nzcoeff[block_y    ][block_x + 1] = 0;
                nzcoeff[block_y + 1][block_x    ] = 0;
                nzcoeff[block_y + 1][block_x + 1] = 0;
            }
        }
    }
}

static void read_tc_chroma_420(mb_t *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int mb_nr = currMB->mbAddrX;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;

    int qp_per_uv[2];
    int qp_rem_uv[2];
    for (int i = 0; i < 2; ++i) {
        qp_per_uv[i] = currMB->qp_scaled[i + 1] / 6;
        qp_rem_uv[i] = currMB->qp_scaled[i + 1] % 6;
    }

    int num_blk8x8_uv = 0;
    if (sps->chroma_format_idc != YUV400)
        num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    int num_uv_blocks = num_blk8x8_uv >> 1;

    if (currMB->cbp > 15) {
        for (int ll = 0; ll < 3; ll += 2) {
            int uv = ll >> 1;

            int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
                currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv[uv]] :
                currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv[uv]];

            int levarr[16], runarr[16], numcoeff;
            read_coeff_4x4_CAVLC(currMB, CHROMA_DC, 0, 0, levarr, runarr, &numcoeff);

            int cofu[4] = { 0 };

            int coef_ctr = -1;
            for (int k = 0; k < numcoeff; ++k) {
                if (levarr[k] != 0) {
                    currMB->s_cbp[0].blk |= 0xf0000 << (ll << 1);
                    coef_ctr += runarr[k] + 1;
                    cofu[coef_ctr] = levarr[k];
                }
            }

            int smb = (currSlice->slice_type == SP_SLICE && !currMB->is_intra_block) ||
                      (currSlice->slice_type == SI_SLICE && currMB->mb_type == SI4MB);
            if (smb || currMB->TransformBypassModeFlag) {
                currSlice->cof[PLANE_U + uv][0][0] = cofu[0];
                currSlice->cof[PLANE_U + uv][0][4] = cofu[1];
                currSlice->cof[PLANE_U + uv][4][0] = cofu[2];
                currSlice->cof[PLANE_U + uv][4][4] = cofu[3];
            } else {
                int temp[4];
                ihadamard2x2(cofu, temp);

                currSlice->cof[PLANE_U + uv][0][0] = ((temp[0] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
                currSlice->cof[PLANE_U + uv][0][4] = ((temp[1] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
                currSlice->cof[PLANE_U + uv][4][0] = ((temp[2] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
                currSlice->cof[PLANE_U + uv][4][4] = ((temp[3] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
            }
        }
    }

    if (currMB->cbp > 31) {
        for (int b8 = 0; b8 < num_blk8x8_uv; ++b8) {
            int uv = b8 > (num_uv_blocks - 1);
            currMB->is_v_block = uv;

            int (*InvLevelScale4x4)[4] = NULL;
            if (!currMB->TransformBypassModeFlag)
                InvLevelScale4x4 = currMB->is_intra_block ?
                    currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv[uv]] :
                    currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv[uv]];

            for (int b4 = 0; b4 < 4; ++b4) {
                int i = cofuv_blk_x[0][b8][b4];
                int j = cofuv_blk_y[0][b8][b4];

                int levarr[16], runarr[16], numcoeff;
                read_coeff_4x4_CAVLC(currMB, CHROMA_AC, i + 2*uv, j + 4, levarr, runarr, &numcoeff);

                int coef_ctr = 0;
                for (int k = 0; k < numcoeff; ++k) {
                    if (levarr[k] != 0) {
                        currMB->s_cbp[0].blk |= i64_power2(cbp_blk_chroma[b8][b4]);
                        coef_ctr += runarr[k] + 1;

                        int i0 = pos_scan4x4[coef_ctr][0];
                        int j0 = pos_scan4x4[coef_ctr][1];

                        if (!currMB->TransformBypassModeFlag)
                            currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0])<<qp_per_uv[uv], 4);
                        else
                            currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = levarr[k];
                    }
                }
            }
        }
    } else
        memset(p_Vid->nz_coeff[mb_nr][1][0], 0, 2 * BLOCK_PIXELS * sizeof(byte));
}

static void read_tc_chroma_422(mb_t *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int mb_nr = currMB->mbAddrX;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;

    int qp_per_uv[2];
    int qp_rem_uv[2];
    for (int i = 0; i < 2; ++i) {
        qp_per_uv[i] = currMB->qp_scaled[i + 1] / 6;
        qp_rem_uv[i] = currMB->qp_scaled[i + 1] % 6;
    }

    int num_blk8x8_uv = 0;
    if (sps->chroma_format_idc != YUV400)
        num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    int num_uv_blocks = num_blk8x8_uv >> 1;

    if (currMB->cbp > 15) {
        for (int ll = 0; ll < 3; ll += 2) {
            int uv = ll >> 1;

            int qp_per_uv_dc = (currMB->qpc[uv] + 3 + sps->QpBdOffsetC) / 6;       //for YUV422 only
            int qp_rem_uv_dc = (currMB->qpc[uv] + 3 + sps->QpBdOffsetC) % 6;       //for YUV422 only
            int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
                currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv_dc] :
                currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv_dc];

            int levarr[16], runarr[16], numcoeff;
            read_coeff_4x4_CAVLC(currMB, CHROMA_DC, 0, 0, levarr, runarr, &numcoeff);

            int m3[2][4] = { { 0 }, { 0 } };

            int coef_ctr = -1;
            for (int k = 0; k < numcoeff; ++k) {
                if (levarr[k] != 0) {
                    currMB->s_cbp[0].blk |= ((int64)0xff0000) << (ll<<2);
                    coef_ctr += runarr[k] + 1;
                    int i0 = SCAN_YUV422[coef_ctr][0];
                    int j0 = SCAN_YUV422[coef_ctr][1];
                    m3[i0][j0] = levarr[k];
                }
            }

            if (!currMB->TransformBypassModeFlag) {
                int m4[2][4];

                m4[0][0] = m3[0][0] + m3[1][0];
                m4[0][1] = m3[0][1] + m3[1][1];
                m4[0][2] = m3[0][2] + m3[1][2];
                m4[0][3] = m3[0][3] + m3[1][3];

                m4[1][0] = m3[0][0] - m3[1][0];
                m4[1][1] = m3[0][1] - m3[1][1];
                m4[1][2] = m3[0][2] - m3[1][2];
                m4[1][3] = m3[0][3] - m3[1][3];

                int temp[2][4];
                for (int i = 0; i < 2; ++i) {
                    int m6[4];

                    m6[0] = m4[i][0] + m4[i][2];
                    m6[1] = m4[i][0] - m4[i][2];
                    m6[2] = m4[i][1] - m4[i][3];
                    m6[3] = m4[i][1] + m4[i][3];

                    temp[i][0] = m6[0] + m6[3];
                    temp[i][1] = m6[1] + m6[2];
                    temp[i][2] = m6[1] - m6[2];
                    temp[i][3] = m6[0] - m6[3];
                }

                for (int j = 0; j < sps->MbHeightC; j += BLOCK_SIZE) {
                    for (int i = 0; i < sps->MbWidthC; i += BLOCK_SIZE)
                        currSlice->cof[PLANE_U + uv][j][i] = rshift_rnd_sf((temp[i / 4][j / 4] * InvLevelScale4x4[0][0]) << qp_per_uv_dc, 6);
                }
            } else {
                for (int j = 0; j < sps->MbHeightC; j += BLOCK_SIZE) {
                    for (int i = 0; i < sps->MbWidthC; i += BLOCK_SIZE)
                        currSlice->cof[PLANE_U + uv][j][i] = m3[i / 4][j / 4];
                }
            }
        }
    }

    if (currMB->cbp > 31) {
        for (int b8 = 0; b8 < num_blk8x8_uv; ++b8) {
            int uv = b8 > (num_uv_blocks - 1);
            currMB->is_v_block = uv;

            int (*InvLevelScale4x4)[4] = NULL;
            if (!currMB->TransformBypassModeFlag)
                InvLevelScale4x4 = currMB->is_intra_block ?
                    currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv[uv]] :
                    currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv[uv]];

            for (int b4 = 0; b4 < 4; ++b4) {
                int i = cofuv_blk_x[1][b8][b4];
                int j = cofuv_blk_y[1][b8][b4];

                int levarr[16], runarr[16], numcoeff;
                read_coeff_4x4_CAVLC(currMB, CHROMA_AC, i + 2*uv, j + 4, levarr, runarr, &numcoeff);

                int coef_ctr = 0;
                for (int k = 0; k < numcoeff; ++k) {
                    if (levarr[k] != 0) {
                        currMB->s_cbp[0].blk |= i64_power2(cbp_blk_chroma[b8][b4]);
                        coef_ctr += runarr[k] + 1;
                        int i0 = pos_scan4x4[coef_ctr][0];
                        int j0 = pos_scan4x4[coef_ctr][1];

                        if (!currMB->TransformBypassModeFlag)
                            currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0])<<qp_per_uv[uv], 4);
                        else
                            currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = levarr[k];
                    }
                }
            }
        }
    } else
        memset(p_Vid->nz_coeff[mb_nr][1][0], 0, 2 * BLOCK_PIXELS * sizeof(byte));
}


void macroblock_t::read_CBP_and_coeffs_from_NAL_CAVLC()
{
    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;

    read_tc_luma(this, PLANE_Y);

    if (sps->chroma_format_idc == YUV420)
        read_tc_chroma_420(this);
    if (sps->chroma_format_idc == YUV422)
        read_tc_chroma_422(this);
    if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
        read_tc_luma(this, PLANE_U);
        read_tc_luma(this, PLANE_V);
    }
}
