/*!
 ***********************************************************************
 * \file read_comp_cabac.c
 *
 * \brief
 *     Read Coefficient Components
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

#include "global.h"
#include "slice.h"
#include "bitstream_elements.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "macroblock.h"
#include "mb_read.h"
#include "neighbour.h"
#include "transform.h"

#include "mb_read_syntax.h"


#define IS_I16MB(MB)    ((MB)->mb_type==I16MB  || (MB)->mb_type==IPCM)
#define IS_DIRECT(MB)   ((MB)->mb_type==0     && (currSlice->slice_type == B_SLICE ))


//! for the linfo_levrun_inter routine
static const byte NTAB1[4][8][2] = {
  {{1,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
  {{1,1},{1,2},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
  {{2,0},{1,3},{1,4},{1,5},{0,0},{0,0},{0,0},{0,0}},
  {{3,0},{2,1},{2,2},{1,6},{1,7},{1,8},{1,9},{4,0}},
};

static const byte LEVRUN1[16] = {
  4,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,
};

//! for the linfo_levrun__c2x2 routine
static const byte LEVRUN3[4] = { 2, 1, 0, 0 };

static const byte NTAB3[2][2][2] = {
  {{1,0},{0,0}},
  {{2,0},{1,1}},
};

/*!
 ************************************************************************
 * \par Input:
 *    length and info
 * \par Output:
 *    level, run
 ************************************************************************
 */
static void linfo_levrun_inter(int len, int info, int *level, int *irun)
{
  //assert (((len >> 1) - 5) < 32);

  if (len <= 9)
  {
    int l2     = imax(0,(len >> 1)-1);
    int inf    = info >> 1;

    *level = NTAB1[l2][inf][0];
    *irun  = NTAB1[l2][inf][1];
    if ((info & 0x01) == 1)
      *level = -*level;                   // make sign
  }
  else                                  // if len > 9, skip using the array
  {
    *irun  = (info & 0x1e) >> 1;
    *level = LEVRUN1[*irun] + (info >> 5) + ( 1 << ((len >> 1) - 5));
    if ((info & 0x01) == 1)
      *level = -*level;
  }
  
  if (len == 1) // EOB
    *level = 0;
}

/*!
 ************************************************************************
 * \par Input:
 *    length and info
 * \par Output:
 *    level, run
 ************************************************************************
 */
static void linfo_levrun_c2x2(int len, int info, int *level, int *irun)
{
  if (len<=5)
  {
    int l2     = imax(0, (len >> 1) - 1);
    int inf    = info >> 1;
    *level = NTAB3[l2][inf][0];
    *irun  = NTAB3[l2][inf][1];
    if ((info & 0x01) == 1)
      *level = -*level;                 // make sign
  }
  else                                  // if len > 5, skip using the array
  {
    *irun  = (info & 0x06) >> 1;
    *level = LEVRUN3[*irun] + (info >> 3) + (1 << ((len >> 1) - 3));
    if ((info & 0x01) == 1)
      *level = -*level;
  }
  
  if (len == 1) // EOB
    *level = 0;
}


static void read_comp_coeff_4x4_smb_CABAC(mb_t *currMB, SyntaxElement *currSE, ColorPlane pl, int block_y, int block_x, int start_scan, int64 *cbp_blk)
{
    slice_t *currSlice = currMB->p_Slice;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    DataPartition *dP;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    const byte *pos_scan_4x4 = pos_scan4x4[0];
    int **cof = currSlice->cof[pl];

    for (int j = block_y; j < block_y + BLOCK_SIZE_8x8; j += 4) {
        currMB->subblock_y = j; // position for coeff_count ctx

        for (int i = block_x; i < block_x + BLOCK_SIZE_8x8; i += 4) {
            currMB->subblock_x = i; // position for coeff_count ctx
            pos_scan_4x4 = pos_scan4x4[start_scan];
            int level = 1;

            if (start_scan == 0) {
                /*
                * make distinction between INTRA and INTER coded
                * luminance coefficients
                */
                currSE->type = currMB->is_intra_block ? SE_LUM_DC_INTRA : SE_LUM_DC_INTER;
                dP = &(currSlice->partArr[partMap[currSE->type]]);
                if (dP->bitstream->ei_flag)  
                    currSE->mapping = linfo_levrun_inter;
                else                                                     
                    currSE->reading = readRunLevel_CABAC;

                dP->readSyntaxElement(currMB, currSE, dP);
                level = currSE->value1;

                if (level != 0) {   /* leave if level == 0 */
                    pos_scan_4x4 += 2 * currSE->value2;

                    int i0 = *pos_scan_4x4++;
                    int j0 = *pos_scan_4x4++;

                    *cbp_blk |= i64_power2(j + (i >> 2)) ;
                    cof[j + j0][i + i0]= level;
                }
            }

            if (level != 0) {
                // make distinction between INTRA and INTER coded luminance coefficients
                currSE->type = currMB->is_intra_block ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
                dP = &(currSlice->partArr[partMap[currSE->type]]);
                if (dP->bitstream->ei_flag)  
                    currSE->mapping = linfo_levrun_inter;
                else                                                     
                    currSE->reading = readRunLevel_CABAC;

                for (int k = 1; k < 17 && level != 0; ++k) {
                    dP->readSyntaxElement(currMB, currSE, dP);
                    level = currSE->value1;

                    if (level != 0) {   /* leave if level == 0 */
                        pos_scan_4x4 += 2 * currSE->value2;

                        int i0 = *pos_scan_4x4++;
                        int j0 = *pos_scan_4x4++;

                        cof[j + j0][i + i0] = level;
                    }
                }
            }
        }
    }
}

static void read_comp_coeff_4x4_CABAC(mb_t *currMB, SyntaxElement *currSE, ColorPlane pl, int (*InvLevelScale4x4)[4], int qp_per, int cbp)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int start_scan = IS_I16MB (currMB)? 1 : 0; 
    int64 *cbp_blk = &currMB->s_cbp[pl].blk;

    if (pl == PLANE_Y || sps->separate_colour_plane_flag)
        currSE->context = IS_I16MB(currMB) ? LUMA_16AC : LUMA_4x4;
    else if (pl == PLANE_U)
        currSE->context = IS_I16MB(currMB) ? CB_16AC : CB_4x4;
    else
        currSE->context = IS_I16MB(currMB) ? CR_16AC : CR_4x4;

    for (int block_y = 0; block_y < MB_BLOCK_SIZE; block_y += BLOCK_SIZE_8x8) { /* all modes */
        int **cof = &currSlice->cof[pl][block_y];
        for (int block_x = 0; block_x < MB_BLOCK_SIZE; block_x += BLOCK_SIZE_8x8) {
            if (cbp & (1 << ((block_y >> 2) + (block_x >> 3)))) { // are there any coeff in current block at all
                read_comp_coeff_4x4_smb_CABAC(currMB, currSE, pl, block_y, block_x, start_scan, cbp_blk);

                if (!currMB->TransformBypassModeFlag) {
                    if (start_scan == 0) {
                        for (int j = 0; j < BLOCK_SIZE_8x8; ++j) {
                            int *coef = &cof[j][block_x];
                            int jj = j & 0x03;
                            for (int i = 0; i < BLOCK_SIZE_8x8; i += 4) {
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][0]) << qp_per, 4);
                                coef++;
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][1]) << qp_per, 4);
                                coef++;
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][2]) << qp_per, 4);
                                coef++;
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][3]) << qp_per, 4);
                                coef++;
                            }
                        }
                    } else {
                        for (int j = 0; j < BLOCK_SIZE_8x8; ++j) {
                            int *coef = &cof[j][block_x];
                            int jj = j & 0x03;
                            for (int i = 0; i < BLOCK_SIZE_8x8; i += 4) {
                                if ((jj != 0) && *coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][0]) << qp_per, 4);
                                coef++;
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][1]) << qp_per, 4);
                                coef++;
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][2]) << qp_per, 4);
                                coef++;
                                if (*coef)
                                    *coef = rshift_rnd_sf((*coef * InvLevelScale4x4[jj][3]) << qp_per, 4);
                                coef++;
                            }
                        }
                    }
                }
            }
        }
    }
}

static void read_comp_coeff_8x8_CABAC(mb_t *currMB, SyntaxElement *currSE, ColorPlane pl)
{
    for (int b8 = 0; b8 < 4; b8++) {
        if (currMB->cbp & (1 << b8)) { // are there any coefficients in the current block
            slice_t *currSlice = currMB->p_Slice;
            sps_t *sps = currSlice->active_sps;
            const byte *partMap = assignSE2partition[currSlice->dp_mode];
            DataPartition *dP;

            int64 cbp_mask = (int64) 51 << (4 * b8 - 2 * (b8 & 0x01)); // corresponds to 110011, as if all four 4x4 blocks contain coeff, shifted to block position            
            int64 *cur_cbp = &currMB->s_cbp[pl].blk;

            // select scan type
            const byte (*pos_scan8x8) = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN8x8[0] : FIELD_SCAN8x8[0];

            int qp_per = currMB->qp_scaled[pl] / 6;
            int qp_rem = currMB->qp_scaled[pl] % 6;
            int transform_pl = sps->separate_colour_plane_flag ? currSlice->colour_plane_id : pl;
            int (*InvLevelScale8x8)[8] = currMB->is_intra_block ?
                currSlice->InvLevelScale8x8_Intra[transform_pl][qp_rem] :
                currSlice->InvLevelScale8x8_Inter[transform_pl][qp_rem];

            // === set offset in current macroblock ===
            int boff_x = (b8&0x01) << 3;
            int boff_y = (b8 >> 1) << 3;
            int **tcoeffs = &currSlice->cof[pl][boff_y];

            currMB->subblock_x = boff_x; // position for coeff_count ctx
            currMB->subblock_y = boff_y; // position for coeff_count ctx

            if (pl == PLANE_Y || sps->separate_colour_plane_flag)
                currSE->context = LUMA_8x8;
            else if (pl == PLANE_U)
                currSE->context = CB_8x8;
            else
                currSE->context = CR_8x8;  
            currSE->reading = readRunLevel_CABAC;

            int level = 1;
            for (int k = 0; k < 65 && level != 0; ++k) {
                //============ read =============
                /*
                * make distinction between INTRA and INTER coded
                * luminance coefficients
                */
                currSE->type = currMB->is_intra_block ?
                    (k == 0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) :
                    (k == 0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER);
                currSE->reading = readRunLevel_CABAC;
                dP = &(currSlice->partArr[partMap[currSE->type]]);

                dP->readSyntaxElement(currMB, currSE, dP);
                level = currSE->value1;

                //============ decode =============
                if (level != 0) {   /* leave if level == 0 */
                    pos_scan8x8 += 2 * (currSE->value2);

                    int i = *pos_scan8x8++;
                    int j = *pos_scan8x8++;

                    *cur_cbp |= cbp_mask;

                    if (!currMB->TransformBypassModeFlag)
                        tcoeffs[j][boff_x + i] = rshift_rnd_sf((level * InvLevelScale8x8[j][i]) << qp_per, 6); // dequantization
                    else
                        tcoeffs[j][boff_x + i] = level;
                }
            }
        }
    }
}


static void read_tc_luma(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    DataPartition *dP = NULL;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    SyntaxElement currSE;

    // select scan type
    const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
    const byte *pos_scan_4x4 = pos_scan4x4[0];

    if (IS_I16MB(currMB)) {
        if (!currMB->dpl_flag) {
            int **cof = currSlice->cof[0];
            pos_scan_4x4 = pos_scan4x4[0];

            currSE.context = LUMA_16DC;
            currSE.type    = SE_LUM_DC_INTRA;
            dP = &(currSlice->partArr[partMap[currSE.type]]);
            if (dP->bitstream->ei_flag)
                currSE.mapping = linfo_levrun_inter;
            else
                currSE.reading = readRunLevel_CABAC;

            int level = 1; // just to get inside the loop

            for (int k = 0; k < 17 && level != 0; ++k) {
                dP->readSyntaxElement(currMB, &currSE, dP);
                level = currSE.value1;

                if (level != 0) { /* leave if level == 0 */
                    pos_scan_4x4 += 2 * currSE.value2;

                    int i0 = (*pos_scan_4x4++) << 2;
                    int j0 = (*pos_scan_4x4++) << 2;

                    cof[j0][i0] = level; // add new intra DC coeff
                }
            }

            if (!currMB->TransformBypassModeFlag)
                itrans_2(currMB, (ColorPlane)currSlice->colour_plane_id); // transform new intra DC
        }
    }

    int qp_per = currMB->qp_scaled[currSlice->colour_plane_id] / 6;
    int qp_rem = currMB->qp_scaled[currSlice->colour_plane_id] % 6;

    // luma coefficients
    //======= Other Modes & CABAC ========
    //------------------------------------          
    if (currMB->cbp) {
        if (currMB->transform_size_8x8_flag)
            //======= 8x8 transform size & CABAC ========
            read_comp_coeff_8x8_CABAC(currMB, &currSE, PLANE_Y); 
        else {
            int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
                currSlice->InvLevelScale4x4_Intra[currSlice->colour_plane_id][qp_rem] :
                currSlice->InvLevelScale4x4_Inter[currSlice->colour_plane_id][qp_rem];
            read_comp_coeff_4x4_CABAC(currMB, &currSE, PLANE_Y, InvLevelScale4x4, qp_per, currMB->cbp);        
        }
    }
}

static void read_tc_chroma_420(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    DataPartition *dP = NULL;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    SyntaxElement currSE;

    // select scan type
    const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
    const byte *pos_scan_4x4 = pos_scan4x4[0];

    int qp_per_uv[2];
    int qp_rem_uv[2];
    //init quant parameters for chroma 
    for (int i = 0; i < 2; ++i) {
        qp_per_uv[i] = currMB->qp_scaled[i + 1] / 6;
        qp_rem_uv[i] = currMB->qp_scaled[i + 1] % 6;
    }

    int num_blk8x8_uv;
    if (sps->chroma_format_idc != YUV400)
        num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    else
        num_blk8x8_uv = 0;
    int num_uv_blocks;
    if (sps->chroma_format_idc != YUV400)
        num_uv_blocks = (((1 << sps->chroma_format_idc) & (~0x1)) >> 1);
    else
        num_uv_blocks = 0;
    int num_cdc_coeff;
    if (sps->chroma_format_idc != YUV400)
        num_cdc_coeff = (((1 << sps->chroma_format_idc) & (~0x1)) << 1);
    else
        num_cdc_coeff = 0;

    //========================== CHROMA DC ============================
    //-----------------------------------------------------------------
    // chroma DC coeff
    if (currMB->cbp > 15) {
        CBPStructure *s_cbp = &currMB->s_cbp[0];

        for (int ll = 0; ll < 3; ll += 2) {
            int uv = ll >> 1;
            currMB->is_v_block = ll;

            int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
                currSlice->InvLevelScale4x4_Intra[uv + 1][qp_rem_uv[uv]] :
                currSlice->InvLevelScale4x4_Inter[uv + 1][qp_rem_uv[uv]];
            //===================== CHROMA DC YUV420 ======================
            int cofu[16];
            memset(cofu, 0, 4 * sizeof(int));

            currSE.context     = CHROMA_DC;
            currSE.type        = currMB->is_intra_block ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
            dP = &(currSlice->partArr[partMap[currSE.type]]);
            if (dP->bitstream->ei_flag)
                currSE.mapping = linfo_levrun_c2x2;
            else
                currSE.reading = readRunLevel_CABAC;

            int coef_ctr = -1;
            int level = 1;
            for (int k = 0; k < num_cdc_coeff + 1 && level != 0; ++k) {
                dP->readSyntaxElement(currMB, &currSE, dP);
                level = currSE.value1;
                if (level != 0) {
                    s_cbp->blk |= 0xf0000 << (ll << 1);
                    coef_ctr += currSE.value2 + 1;

                    // Bug: currSlice->cofu has only 4 entries, hence coef_ctr MUST be <4 (which is
                    // caught by the assert().  If it is bigger than 4, it starts patching the
                    // p_Vid->predmode pointer, which leads to bugs later on.
                    //
                    // This assert() should be left in the code, because it captures a very likely
                    // bug early when testing in error prone environments (or when testing NAL
                    // functionality).

                    assert(coef_ctr < num_cdc_coeff);
                    cofu[coef_ctr] = level;
                }
            }

            VideoParameters *p_Vid = currMB->p_Vid;
            int smb = (p_Vid->type == SP_SLICE && currMB->is_intra_block == FALSE) ||
                      (p_Vid->type == SI_SLICE && currMB->mb_type == SI4MB);
            if (smb || currMB->TransformBypassModeFlag) { // check to see if MB type is SPred or SIntra4x4
                currSlice->cof[uv + 1][0][0] = cofu[0];
                currSlice->cof[uv + 1][0][4] = cofu[1];
                currSlice->cof[uv + 1][4][0] = cofu[2];
                currSlice->cof[uv + 1][4][4] = cofu[3];
            } else {
                int temp[4];
                int scale_dc = InvLevelScale4x4[0][0];
                int **cof = currSlice->cof[uv + 1];

                ihadamard2x2(cofu, temp);

                cof[0][0] = (((temp[0] * scale_dc) << qp_per_uv[uv]) >> 5);
                cof[0][4] = (((temp[1] * scale_dc) << qp_per_uv[uv]) >> 5);
                cof[4][0] = (((temp[2] * scale_dc) << qp_per_uv[uv]) >> 5);
                cof[4][4] = (((temp[3] * scale_dc) << qp_per_uv[uv]) >> 5);
            }
        }
    }

    //========================== CHROMA AC ============================
    //-----------------------------------------------------------------
    // chroma AC coeff, all zero fram start_scan
    if (currMB->cbp > 31) {
        currSE.context = CHROMA_AC;
        currSE.type    = currMB->is_intra_block ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
        dP = &(currSlice->partArr[partMap[currSE.type]]);
        if (dP->bitstream->ei_flag)
            currSE.mapping = linfo_levrun_inter;
        else
            currSE.reading = readRunLevel_CABAC;

        CBPStructure *s_cbp = &currMB->s_cbp[0];

        int yuv = sps->chroma_format_idc - 1;

        for (int b8 = 0; b8 < num_blk8x8_uv; ++b8) {
            int uv = b8 > (num_uv_blocks - 1);
            currMB->is_v_block = uv;
            int **cof = currSlice->cof[uv + 1];
            int (*InvLevelScale4x4)[4] = NULL;
            if (!currMB->TransformBypassModeFlag)
                InvLevelScale4x4 = currMB->is_intra_block ?
                    currSlice->InvLevelScale4x4_Intra[uv + 1][qp_rem_uv[uv]] :
                    currSlice->InvLevelScale4x4_Inter[uv + 1][qp_rem_uv[uv]];

            for (int b4 = 0; b4 < 4; ++b4) {
                int i = cofuv_blk_x[yuv][b8][b4];
                int j = cofuv_blk_y[yuv][b8][b4];

                pos_scan_4x4 = pos_scan4x4[1];
                int level = 1;

                currMB->subblock_y = subblk_offset_y[yuv][b8][b4];
                currMB->subblock_x = subblk_offset_x[yuv][b8][b4];

                for (int k = 0; k < 16 && level != 0; ++k) {
                    dP->readSyntaxElement(currMB, &currSE, dP);
                    level = currSE.value1;

                    if (level != 0) {
                        s_cbp->blk |= i64_power2(cbp_blk_chroma[b8][b4]);
                        pos_scan_4x4 += currSE.value2 << 1;

                        int i0 = *pos_scan_4x4++;
                        int j0 = *pos_scan_4x4++;

                        if (!currMB->TransformBypassModeFlag)
                            cof[(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((level * InvLevelScale4x4[j0][i0])<<qp_per_uv[uv], 4);
                        else
                            cof[(j<<2) + j0][(i<<2) + i0] = level;
                    }
                }
            }
        }
    }
}

static void read_tc_chroma_422(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    DataPartition *dP = NULL;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    SyntaxElement currSE;

    // select scan type
    const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
    const byte *pos_scan_4x4 = pos_scan4x4[0];

    int qp_per_uv[2];
    int qp_rem_uv[2];
    //init quant parameters for chroma 
    for (int i = 0; i < 2; ++i) {
        qp_per_uv[i] = currMB->qp_scaled[i + 1] / 6;
        qp_rem_uv[i] = currMB->qp_scaled[i + 1] % 6;
    }

    int num_blk8x8_uv;
    if (sps->chroma_format_idc != YUV400)
        num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    else
        num_blk8x8_uv = 0;
    int num_uv_blocks;
    if (sps->chroma_format_idc != YUV400)
        num_uv_blocks = (((1 << sps->chroma_format_idc) & (~0x1)) >> 1);
    else
        num_uv_blocks = 0;
    int num_cdc_coeff;
    if (sps->chroma_format_idc != YUV400)
        num_cdc_coeff = (((1 << sps->chroma_format_idc) & (~0x1)) << 1);
    else
        num_cdc_coeff = 0;

    //========================== CHROMA DC ============================
    //-----------------------------------------------------------------
    // chroma DC coeff
    if (currMB->cbp > 15) {      
        for (int ll = 0; ll < 3; ll += 2) {
            int uv = ll >> 1;
            currMB->is_v_block = ll;

            int **imgcof = currSlice->cof[uv + 1];
            int m3[2][4] = {{0,0,0,0},{0,0,0,0}};
            int m4[2][4] = {{0,0,0,0},{0,0,0,0}};
            int m6[4];
            int qp_per_uv_dc = (currMB->qpc[uv] + 3 + sps->QpBdOffsetC) / 6;       //for YUV422 only
            int qp_rem_uv_dc = (currMB->qpc[uv] + 3 + sps->QpBdOffsetC) % 6;       //for YUV422 only
            int (*InvLevelScale4x4)[4] = NULL;
            if (currMB->is_intra_block)
                InvLevelScale4x4 = currSlice->InvLevelScale4x4_Intra[uv + 1][qp_rem_uv_dc];
            else 
                InvLevelScale4x4 = currSlice->InvLevelScale4x4_Inter[uv + 1][qp_rem_uv_dc];

            //===================== CHROMA DC YUV422 ======================
            CBPStructure *s_cbp = &currMB->s_cbp[0];
            int coef_ctr = -1;
            int level = 1;
            for (int k = 0; k < 9 && level != 0; ++k) {
                currSE.context = CHROMA_DC_2x4;
                currSE.type    = currMB->is_intra_block ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
                dP = &(currSlice->partArr[partMap[currSE.type]]);
                if (dP->bitstream->ei_flag)
                    currSE.mapping = linfo_levrun_c2x2;
                else
                    currSE.reading = readRunLevel_CABAC;
                dP->readSyntaxElement(currMB, &currSE, dP);
                level = currSE.value1;

                if (level != 0) {
                    s_cbp->blk |= ((int64)0xff0000) << (ll<<2);
                    coef_ctr += currSE.value2 + 1;
                    assert(coef_ctr < num_cdc_coeff);
                    int i0 = SCAN_YUV422[coef_ctr][0];
                    int j0 = SCAN_YUV422[coef_ctr][1];

                    m3[i0][j0] = level;
                }
            }

            // inverse CHROMA DC YUV422 transform
            // horizontal
            if (!currMB->TransformBypassModeFlag) {
                m4[0][0] = m3[0][0] + m3[1][0];
                m4[0][1] = m3[0][1] + m3[1][1];
                m4[0][2] = m3[0][2] + m3[1][2];
                m4[0][3] = m3[0][3] + m3[1][3];

                m4[1][0] = m3[0][0] - m3[1][0];
                m4[1][1] = m3[0][1] - m3[1][1];
                m4[1][2] = m3[0][2] - m3[1][2];
                m4[1][3] = m3[0][3] - m3[1][3];

                for (int i = 0; i < 2; ++i) {
                    m6[0] = m4[i][0] + m4[i][2];
                    m6[1] = m4[i][0] - m4[i][2];
                    m6[2] = m4[i][1] - m4[i][3];
                    m6[3] = m4[i][1] + m4[i][3];

                    imgcof[ 0][i << 2] = m6[0] + m6[3];
                    imgcof[ 4][i << 2] = m6[1] + m6[2];
                    imgcof[ 8][i << 2] = m6[1] - m6[2];
                    imgcof[12][i << 2] = m6[0] - m6[3];
                }

                for (int j = 0; j < sps->MbHeightC; j += BLOCK_SIZE) {
                    for (int i = 0; i < sps->MbWidthC; i+= BLOCK_SIZE)
                        imgcof[j][i] = rshift_rnd_sf((imgcof[j][i] * InvLevelScale4x4[0][0]) << qp_per_uv_dc, 6);
                }
            } else {
                for (int j = 0; j < 4; ++j) {
                    for (int i = 0; i < 2; ++i)
                        currSlice->cof[uv + 1][j<<2][i<<2] = m3[i][j];
                }
            }
        }
    }

    //========================== CHROMA AC ============================
    //-----------------------------------------------------------------
    // chroma AC coeff, all zero fram start_scan
    if (currMB->cbp > 31) {
        currSE.context = CHROMA_AC;
        currSE.type    = currMB->is_intra_block ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
        dP = &(currSlice->partArr[partMap[currSE.type]]);
        if (dP->bitstream->ei_flag)
            currSE.mapping = linfo_levrun_inter;
        else
            currSE.reading = readRunLevel_CABAC;

        int yuv = sps->chroma_format_idc - 1;

        CBPStructure *s_cbp = &currMB->s_cbp[0];
        for (int b8 = 0; b8 < num_blk8x8_uv; ++b8) {
            int uv = (b8 > ((num_uv_blocks) - 1 ));
            currMB->is_v_block = uv;

            int (*InvLevelScale4x4)[4] = NULL;
            if (!currMB->TransformBypassModeFlag)
                InvLevelScale4x4 = currMB->is_intra_block ?
                    currSlice->InvLevelScale4x4_Intra[uv + 1][qp_rem_uv[uv]] :
                    currSlice->InvLevelScale4x4_Inter[uv + 1][qp_rem_uv[uv]];

            for (int b4 = 0; b4 < 4; ++b4) {
                int i = cofuv_blk_x[yuv][b8][b4];
                int j = cofuv_blk_y[yuv][b8][b4];

                currMB->subblock_y = subblk_offset_y[yuv][b8][b4];
                currMB->subblock_x = subblk_offset_x[yuv][b8][b4];

                pos_scan_4x4 = pos_scan4x4[1];
                int level = 1;
                for (int k = 0; k < 16 && level != 0; ++k) {
                    dP->readSyntaxElement(currMB, &currSE, dP);
                    level = currSE.value1;

                    if (level != 0) {
                        s_cbp->blk |= i64_power2(cbp_blk_chroma[b8][b4]);
                        pos_scan_4x4 += (currSE.value2 << 1);

                        int i0 = *pos_scan_4x4++;
                        int j0 = *pos_scan_4x4++;

                        if (!currMB->TransformBypassModeFlag)
                            currSlice->cof[uv + 1][(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((level * InvLevelScale4x4[j0][i0]) << qp_per_uv[uv], 4);
                        else
                            currSlice->cof[uv + 1][(j<<2) + j0][(i<<2) + i0] = level;
                    }
                }
            }
        }
    }
}

static void read_tc_chroma_444(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    DataPartition *dP = NULL;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    SyntaxElement currSE;

    // select scan type
    const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;

    int qp_per_uv[2];
    int qp_rem_uv[2];

    for (int uv = 0; uv < 2; ++uv) {
    /*----------------------16x16DC Luma_Add----------------------*/
        if (IS_I16MB (currMB)) { // read DC coeffs for new intra modes       
            currSE.type = SE_LUM_DC_INTRA;
            if (!sps->separate_colour_plane_flag)
                currSE.context = LUMA_16DC; 
            else
                currSE.context = (uv==0) ? CB_16DC : CR_16DC;
            dP = &(currSlice->partArr[partMap[currSE.type]]);
            if (dP->bitstream->ei_flag)
                currSE.mapping = linfo_levrun_inter;
            else
                currSE.reading = readRunLevel_CABAC;

            int coef_ctr = -1;
            int level = 1;                            // just to get inside the loop
            for (int k = 0; k < 17 && level != 0; ++k) {
                dP->readSyntaxElement(currMB, &currSE, dP);
                level = currSE.value1;

                if (level != 0) {                    // leave if level == 0
                    coef_ctr += currSE.value2 + 1;

                    int i0 = pos_scan4x4[coef_ctr][0];
                    int j0 = pos_scan4x4[coef_ctr][1];
                    currSlice->cof[uv + 1][j0<<2][i0<<2] = level;
                }
            }

            if (!currMB->TransformBypassModeFlag)
                itrans_2(currMB, (ColorPlane) (uv + 1)); // transform new intra DC
        }

        //init constants for every chroma qp offset
        qp_per_uv[uv] = (currMB->qpc[uv] + sps->QpBdOffsetC) / 6;
        qp_rem_uv[uv] = (currMB->qpc[uv] + sps->QpBdOffsetC) % 6;

        int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
            currSlice->InvLevelScale4x4_Intra[uv + 1][qp_rem_uv[uv]] :
            currSlice->InvLevelScale4x4_Inter[uv + 1][qp_rem_uv[uv]];

        if (currMB->cbp) {
            if (currMB->transform_size_8x8_flag) 
                //======= 8x8 transform size & CABAC ========
                read_comp_coeff_8x8_CABAC(currMB, &currSE, (ColorPlane) (PLANE_U + uv)); 
            else //4x4
                read_comp_coeff_4x4_CABAC(currMB, &currSE, (ColorPlane) (PLANE_U + uv), InvLevelScale4x4, qp_per_uv[uv], currMB->cbp);
        }
    } 

}


void macroblock_t::read_CBP_and_coeffs_from_NAL_CABAC()
{
    mb_t *currMB = this;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    pps_t *pps = currSlice->active_pps;
    DataPartition *dP = NULL;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    SyntaxElement currSE;

    if (!IS_I16MB(currMB)) {
        //=====   C B P   =====
        //---------------------
        currSE.type = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)
                      ? SE_CBP_INTRA : SE_CBP_INTER;
        dP = &(currSlice->partArr[partMap[currSE.type]]);
        if (dP->bitstream->ei_flag) {
            currSE.mapping =
                (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB) ?
                (sps->chroma_format_idc == 0 || sps->chroma_format_idc == 3 ?
                 linfo_cbp_intra_other : linfo_cbp_intra_normal) :
                (sps->chroma_format_idc == 0 || sps->chroma_format_idc == 3 ?
                 linfo_cbp_inter_other : linfo_cbp_inter_normal);
        } else
            currSE.reading = read_CBP_CABAC;

        dP->readSyntaxElement(currMB, &currSE, dP);
        currMB->cbp = currSE.value1;

        //============= Transform size flag for INTER MBs =============
        //-------------------------------------------------------------
        int need_transform_size_flag =
            ((currMB->mb_type >= 1 && currMB->mb_type <= 3) ||
             (IS_DIRECT(currMB) && sps->direct_8x8_inference_flag) ||
             currMB->NoMbPartLessThan8x8Flag) &&
            (currMB->mb_type != I8MB) && (currMB->mb_type != I4MB) &&
            (currMB->cbp & 15) && pps->transform_8x8_mode_flag;

        if (need_transform_size_flag)
            currMB->transform_size_8x8_flag = parse_transform_size_8x8_flag(this);
    }

    //=====   DQUANT   =====
    //----------------------
    // Delta quant only if nonzero coeffs
    if (IS_I16MB(currMB) || currMB->cbp != 0) {
        currMB->parse_mb_qp();

        if (currSlice->dp_mode) {
            if (!currMB->is_intra_block && currSlice->dpC_NotPresent)
                currMB->dpl_flag = 1;
            if (currMB->is_intra_block && currSlice->dpB_NotPresent) {
                currMB->ei_flag  = 1;
                currMB->dpl_flag = 1;
            }
            check_dp_neighbors(currMB);
            if (currMB->dpl_flag)
                currMB->cbp = 0;
        }
    }

    currMB->update_qp(currSlice->SliceQpY);

    read_tc_luma(currMB);

    if (sps->chroma_format_idc == YUV420)
        read_tc_chroma_420(currMB);
    if (sps->chroma_format_idc == YUV422)
        read_tc_chroma_422(currMB);
    if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag)
        read_tc_chroma_444(currMB);
}
