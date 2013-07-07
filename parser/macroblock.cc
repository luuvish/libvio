/*!
 ***********************************************************************
 * \file macroblock.c
 *
 * \brief
 *     Decode a Macroblock
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Detlev Marpe
 *    - Gabi Blaettermann
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Lowell Winger                   <lwinger@lsil.com>
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

#include <math.h>

#include "global.h"
#include "slice.h"
#include "dpb.h"
#include "dpb_mvc.h"
#include "bitstream_elements.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "macroblock.h"
#include "mb_read.h"
#include "fmo.h"
#include "image.h"
#include "neighbour.h"
#include "biaridecod.h"
#include "transform.h"
#include "mv_prediction.h"
#include "inter_prediction.h"
#include "intra_prediction.h"

#include "mb.h"

#define MB_PIXELS            256 // MB_BLOCK_SIZE * MB_BLOCK_SIZE
#define MB_BLOCK_SHIFT         4
#define MB_BLOCK_PARTITIONS   16 // (BLOCK_MULTIPLE * BLOCK_MULTIPLE)


/*!
************************************************************************
* \brief
*    updates chroma QP according to luma QP and bit depth
************************************************************************
*/
void update_qp(Macroblock *currMB, int qp)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    sps_t *sps = p_Vid->active_sps;
    StorablePicture *dec_picture = currMB->p_Slice->dec_picture;
    currMB->qp = qp;
    currMB->qp_scaled[0] = qp + sps->QpBdOffsetY;

    for (int i = 0; i < 2; i++) {
        currMB->qpc[i] = iClip3 (-(sps->QpBdOffsetC), 51, currMB->qp + dec_picture->chroma_qp_offset[i]);
        currMB->qpc[i] = currMB->qpc[i] < 0 ? currMB->qpc[i] : QP_SCALE_CR[currMB->qpc[i]];
        currMB->qp_scaled[i + 1] = currMB->qpc[i] + sps->QpBdOffsetC;
    }

    currMB->is_lossless = (currMB->qp_scaled[0] == 0 && sps->qpprime_y_zero_transform_bypass_flag);
    set_read_comp_coeff_cavlc(currMB);
    set_read_comp_coeff_cabac(currMB);
}

/*!
 ************************************************************************
 * \brief
 *    initializes the current macroblock
 ************************************************************************
 */
void start_macroblock(Macroblock *currMB)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currSlice->p_Vid;
    int mb_nr = currMB->mbAddrX;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_cr_size = mb_cr_size_x * mb_cr_size_y;

    /* Update coordinates of the current macroblock */
    if (currSlice->MbaffFrameFlag) {
        currMB->mb.x = (short) (   (mb_nr) % (2 * sps->PicWidthInMbs));
        currMB->mb.y = (short) (2*((mb_nr) / (2 * sps->PicWidthInMbs)));

        currMB->mb.y += (currMB->mb.x & 0x01);
        currMB->mb.x >>= 1;
    } else
        currMB->mb = p_Vid->PicPos[mb_nr];

    /* Define pixel/block positions */
    int mb_x = currMB->mb.x;
    int mb_y = currMB->mb.y;
    currMB->block_x     = mb_x << BLOCK_SHIFT;                 /* horizontal block position */
    currMB->block_y     = mb_y << BLOCK_SHIFT;                 /* vertical block position */
    currMB->block_y_aff = mb_y << BLOCK_SHIFT;                 /* interlace relative vertical position */
    currMB->pix_x       = mb_x << MB_BLOCK_SHIFT;              /* horizontal luma pixel position */
    currMB->pix_y       = mb_y << MB_BLOCK_SHIFT;              /* vertical luma pixel position */
    currMB->pix_c_x     = mb_x * mb_cr_size_x;  /* horizontal chroma pixel position */
    currMB->pix_c_y     = mb_y * mb_cr_size_y;  /* vertical chroma pixel position */

    // reset intra mode
    currMB->is_intra_block = FALSE;
    // reset mode info
    currMB->mb_type         = 0;
    currMB->delta_quant     = 0;
    currMB->cbp             = 0;    
    currMB->c_ipred_mode    = DC_PRED_8;

    // Save the slice number of this macroblock. When the macroblock below
    // is coded it will use this to decide if prediction for above is possible
    currMB->slice_nr = (short) currSlice->current_slice_nr;

    CheckAvailabilityOfNeighbors(currMB);

    set_read_and_store_CBP(currMB, currSlice->active_sps->chroma_format_idc);

    // Reset syntax element entries in MB struct
    if (currSlice->slice_type != I_SLICE) {
        if (currSlice->slice_type != B_SLICE)
            memset(currMB->mvd[0][0][0], 0, MB_BLOCK_PARTITIONS * 2 * sizeof(short));
        else
            memset(currMB->mvd[0][0][0], 0, 2 * MB_BLOCK_PARTITIONS * 2 * sizeof(short));
    }

    memset(currMB->s_cbp, 0, 3 * sizeof(CBPStructure));

    // initialize currSlice->mb_rres
    if (currSlice->is_reset_coeff == FALSE) {
        memset( currSlice->mb_rres[0][0], 0, MB_PIXELS * sizeof(int));
        memset( currSlice->mb_rres[1][0], 0, mb_cr_size * sizeof(int));
        memset( currSlice->mb_rres[2][0], 0, mb_cr_size * sizeof(int));
        if (currSlice->is_reset_coeff_cr == FALSE) {
            memset( currSlice->cof[0][0], 0, 3 * MB_PIXELS * sizeof(int));
            currSlice->is_reset_coeff_cr = TRUE;
        } else
            memset( currSlice->cof[0][0], 0, MB_PIXELS * sizeof(int));

        currSlice->is_reset_coeff = TRUE;
    }

    // store filtering parameters for this MB
    currMB->DFDisableIdc      = currSlice->disable_deblocking_filter_idc;
    currMB->DFAlphaC0Offset   = currSlice->FilterOffsetA;
    currMB->DFBetaOffset      = currSlice->FilterOffsetB;
    currMB->list_offset       = 0;
    currMB->mixedModeEdgeFlag = 0;
}

/*!
 ************************************************************************
 * \brief
 *    set coordinates of the next macroblock
 *    check end_of_slice condition
 ************************************************************************
 */
bool exit_macroblock(Slice *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    int eos_bit = (!currSlice->MbaffFrameFlag || currSlice->current_mb_nr % 2);
    int startcode_follows;
    sps_t *sps = currSlice->active_sps;
    int PicSizeInMbs = sps->PicWidthInMbs * (sps->FrameHeightInMbs / (1 + currSlice->field_pic_flag));

    //! The if() statement below resembles the original code, which tested
    //! p_Vid->current_mb_nr == p_Vid->PicSizeInMbs.  Both is, of course, nonsense
    //! In an error prone environment, one can only be sure to have a new
    //! picture by checking the tr of the next slice header!

    ++(currSlice->num_dec_mb);

    if (currSlice->current_mb_nr == PicSizeInMbs - 1)
        return 1;

    currSlice->current_mb_nr = FmoGetNextMBNr(p_Vid, currSlice->current_mb_nr);

    if (currSlice->active_pps->entropy_coding_mode_flag)
        startcode_follows = cabac_startcode_follows(currSlice, eos_bit);
    else
        startcode_follows = uvlc_startcode_follows(currSlice, eos_bit);

    if (currSlice->current_mb_nr == -1) { // End of Slice group, MUST be end of slice
        assert(startcode_follows);
        return 1;
    }

    if (!startcode_follows)
        return 0;

    if (currSlice->slice_type == I_SLICE || currSlice->slice_type == SI_SLICE)
        return 1;
    if (p_Vid->active_pps->entropy_coding_mode_flag)
        return 1;
    if (currSlice->cod_counter <= 0)
        return 1;

    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for P-Frames
 ************************************************************************
 */
static void interpret_mb_mode_P(Macroblock *currMB)
{
    static const short ICBPTAB[6] = {0,16,32,15,31,47};
    short  mbmode = currMB->mb_type;

    if (mbmode < 4) {
        currMB->mb_type = mbmode;
        memset(currMB->b8mode, mbmode, 4 * sizeof(char));
        memset(currMB->b8pdir, 0, 4 * sizeof(char));
    } else if ((mbmode == 4 || mbmode == 5)) {
        currMB->mb_type = P8x8;
        currMB->p_Slice->allrefzero = (mbmode == 5);
    } else if (mbmode == 6) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I4MB;
        memset(currMB->b8mode, IBLOCK, 4 * sizeof(char));
        memset(currMB->b8pdir,     -1, 4 * sizeof(char));
    } else if (mbmode == 31) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = IPCM;
        currMB->cbp = -1;
        currMB->i16mode = 0;

        memset(currMB->b8mode, 0, 4 * sizeof(char));
        memset(currMB->b8pdir,-1, 4 * sizeof(char));
    } else {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I16MB;
        currMB->cbp = ICBPTAB[((mbmode-7))>>2];
        currMB->i16mode = ((mbmode-7)) & 0x03;
        memset(currMB->b8mode, 0, 4 * sizeof(char));
        memset(currMB->b8pdir,-1, 4 * sizeof(char));
    }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for I-Frames
 ************************************************************************
 */
static void interpret_mb_mode_I(Macroblock *currMB)
{
    static const short ICBPTAB[6] = {0,16,32,15,31,47};
    short mbmode = currMB->mb_type;

    if (mbmode == 0) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I4MB;
        memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else if (mbmode == 25) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type=IPCM;
        currMB->cbp= -1;
        currMB->i16mode = 0;

        memset(currMB->b8mode, 0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I16MB;
        currMB->cbp= ICBPTAB[(mbmode-1)>>2];
        currMB->i16mode = (mbmode-1) & 0x03;
        memset(currMB->b8mode, 0, 4 * sizeof(char));
        memset(currMB->b8pdir,-1, 4 * sizeof(char));
    }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for B-Frames
 ************************************************************************
 */
static void interpret_mb_mode_B(Macroblock *currMB)
{
    static const char offset2pdir16x16[12]   = {0, 0, 1, 2, 0,0,0,0,0,0,0,0};
    static const char offset2pdir16x8[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},{1,0},
                                                {0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2},{0,0}};
    static const char offset2pdir8x16[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},
                                                {1,0},{0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2}};

    static const char ICBPTAB[6] = {0,16,32,15,31,47};

    short i, mbmode;
    short mbtype = currMB->mb_type;

    //--- set mbtype, b8type, and b8pdir ---
    if (mbtype == 0) { // direct
        mbmode = 0;
        memset(currMB->b8mode, 0, 4 * sizeof(char));
        memset(currMB->b8pdir, 2, 4 * sizeof(char));
    } else if (mbtype == 23) { // intra4x4
        currMB->is_intra_block = TRUE;
        mbmode = I4MB;
        memset(currMB->b8mode, IBLOCK,4 * sizeof(char));
        memset(currMB->b8pdir, -1,4 * sizeof(char));
    } else if ((mbtype > 23) && (mbtype < 48) ) { // intra16x16
        currMB->is_intra_block = TRUE;
        mbmode = I16MB;
        memset(currMB->b8mode,  0, 4 * sizeof(char));
        memset(currMB->b8pdir, -1, 4 * sizeof(char));

        currMB->cbp     = (int) ICBPTAB[(mbtype-24)>>2];
        currMB->i16mode = (mbtype-24) & 0x03;
    } else if (mbtype == 22) { // 8x8(+split)
        mbmode = P8x8;       // b8mode and pdir is transmitted in additional codewords
    } else if (mbtype < 4) { // 16x16
        mbmode = 1;
        memset(currMB->b8mode, 1,4 * sizeof(char));
        memset(currMB->b8pdir, offset2pdir16x16[mbtype], 4 * sizeof(char));
    } else if(mbtype == 48) {
        currMB->is_intra_block = TRUE;
        mbmode=IPCM;
        memset(currMB->b8mode, 0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));

        currMB->cbp= -1;
        currMB->i16mode = 0;
    } else if ((mbtype & 0x01) == 0) { // 16x8
        mbmode = 2;
        memset(currMB->b8mode, 2,4 * sizeof(char));
        for(i=0;i<4;++i)
            currMB->b8pdir[i] = offset2pdir16x8 [mbtype][i>>1];
    } else {
        mbmode=3;
        memset(currMB->b8mode, 3,4 * sizeof(char));
        for(i=0;i<4; ++i)
            currMB->b8pdir[i] = offset2pdir8x16 [mbtype][i&0x01];
    }
    currMB->mb_type = mbmode;
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for SI-Frames
 ************************************************************************
 */
static void interpret_mb_mode_SI(Macroblock *currMB)
{
    //VideoParameters *p_Vid = currMB->p_Vid;
    const int ICBPTAB[6] = {0,16,32,15,31,47};
    short mbmode = currMB->mb_type;

    if (mbmode == 0) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = SI4MB;
        memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
        currMB->p_Slice->siblock[currMB->mb.y][currMB->mb.x]=1;
    } else if (mbmode == 1) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I4MB;
        memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else if (mbmode == 26) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type=IPCM;
        currMB->cbp= -1;
        currMB->i16mode = 0;
        memset(currMB->b8mode,0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I16MB;
        currMB->cbp= ICBPTAB[(mbmode-2)>>2];
        currMB->i16mode = (mbmode-2) & 0x03;
        memset(currMB->b8mode,0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    }
}

void interpret_mb_mode(Macroblock *currMB)
{
    Slice *currSlice = currMB->p_Slice;
    switch (currSlice->slice_type) {
    case P_SLICE: 
    case SP_SLICE:
        interpret_mb_mode_P(currMB);
        break;
    case B_SLICE:
        interpret_mb_mode_B(currMB);
        break;
    case I_SLICE: 
        interpret_mb_mode_I(currMB);
        break;
    case SI_SLICE: 
        interpret_mb_mode_SI(currMB);
        break;
    default:
        printf("Unsupported slice type\n");
        break;
    }
}
