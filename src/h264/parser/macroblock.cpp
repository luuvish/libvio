#include <math.h>

#include "global.h"
#include "slice.h"
#include "dpb.h"
#include "bitstream_cabac.h"
#include "data_partition.h"
#include "macroblock.h"
#include "mb_read.h"
#include "fmo.h"
#include "image.h"
#include "neighbour.h"
#include "transform.h"
#include "mv_prediction.h"
#include "inter_prediction.h"
#include "intra_prediction.h"

#define MB_PIXELS            256 // MB_BLOCK_SIZE * MB_BLOCK_SIZE
#define MB_BLOCK_SHIFT         4
#define MB_BLOCK_PARTITIONS   16 // (BLOCK_MULTIPLE * BLOCK_MULTIPLE)



void macroblock_t::init(slice_t *slice)
{
    this->p_Slice = slice;
    this->p_Vid   = slice->p_Vid;
    this->mbAddrX = slice->current_mb_nr;

    sps_t *sps = slice->active_sps;
    int mb_cr_size = sps->MbWidthC * sps->MbHeightC;

    /* Update coordinates of the current macroblock */
    if (slice->MbaffFrameFlag) {
        this->mb.x = (this->mbAddrX / 2) % sps->PicWidthInMbs;
        this->mb.y = (this->mbAddrX / 2) / sps->PicWidthInMbs * 2 + (this->mbAddrX % 2);
    } else {
        this->mb.x = this->mbAddrX % sps->PicWidthInMbs;
        this->mb.y = this->mbAddrX / sps->PicWidthInMbs;
    }

    /* Define pixel/block positions */
    int mb_x = this->mb.x;
    int mb_y = this->mb.y;
    this->block_x     = mb_x << BLOCK_SHIFT;    /* horizontal block position */
    this->block_y     = mb_y << BLOCK_SHIFT;    /* vertical block position */
    this->pix_x       = mb_x * 16;              /* horizontal luma pixel position */
    this->pix_y       = mb_y * 16;              /* vertical luma pixel position */
    this->pix_c_x     = mb_x * sps->MbWidthC;   /* horizontal chroma pixel position */
    this->pix_c_y     = mb_y * sps->MbHeightC;  /* vertical chroma pixel position */

    this->mb_qp_delta = 0;
    // reset intra mode
    this->is_intra_block = 0;
    // reset mode info
    this->mb_skip_flag   = 0;
    this->mb_type        = 0;
    this->cbp            = 0;    
    this->intra_chroma_pred_mode = Intra_Chroma_DC;

    // Save the slice number of this macroblock. When the macroblock below
    // is coded it will use this to decide if prediction for above is possible
    this->slice_nr = (short) slice->current_slice_nr;

    CheckAvailabilityOfNeighbors(this);

    // Reset syntax element entries in MB struct
    if (slice->slice_type != I_SLICE) {
        if (slice->slice_type != B_SLICE)
            memset(this->mvd[0][0][0], 0, MB_BLOCK_PARTITIONS * 2 * sizeof(short));
        else
            memset(this->mvd[0][0][0], 0, 2 * MB_BLOCK_PARTITIONS * 2 * sizeof(short));
    }

    memset(this->s_cbp, 0, 3 * sizeof(CBPStructure));

    // initialize slice->mb_rres
    if (!slice->is_reset_coeff) {
        memset(slice->mb_rres[0][0], 0, MB_PIXELS * sizeof(int));
        memset(slice->mb_rres[1][0], 0, mb_cr_size * sizeof(int));
        memset(slice->mb_rres[2][0], 0, mb_cr_size * sizeof(int));
        if (!slice->is_reset_coeff_cr) {
            memset(slice->cof[0][0], 0, 3 * MB_PIXELS * sizeof(int));
            slice->is_reset_coeff_cr = 1;
        } else
            memset(slice->cof[0][0], 0, MB_PIXELS * sizeof(int));

        slice->is_reset_coeff = 1;
    }

    // store filtering parameters for this MB
    this->DFDisableIdc      = slice->disable_deblocking_filter_idc;
    this->DFAlphaC0Offset   = slice->FilterOffsetA;
    this->DFBetaOffset      = slice->FilterOffsetB;
    this->mixedModeEdgeFlag = 0;

    pps_t *pps = slice->active_pps;
    int CurrMbAddr = this->mbAddrX;

    this->update_qp(slice->SliceQpY);

    this->mb_field_decoding_flag = 0;
    if (slice->MbaffFrameFlag) {
        bool prevMbSkipped = (CurrMbAddr % 2 == 1) ?
            slice->mb_data[CurrMbAddr - 1].mb_skip_flag : 0;
        if (CurrMbAddr % 2 == 0 || prevMbSkipped) {
            if (this->mbAvailA)
                this->mb_field_decoding_flag = slice->mb_data[this->mbAddrA].mb_field_decoding_flag;
            else if (this->mbAvailB)
                this->mb_field_decoding_flag = slice->mb_data[this->mbAddrB].mb_field_decoding_flag;
        } else
            this->mb_field_decoding_flag = slice->mb_data[CurrMbAddr - 1].mb_field_decoding_flag;
    }

    if (pps->entropy_coding_mode_flag)
        CheckAvailabilityOfNeighborsCABAC(this);
}

bool macroblock_t::close(slice_t *slice)
{
    slice_t *currSlice = slice;
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
        startcode_follows = eos_bit && currSlice->partArr[0].de_cabac.decode_terminate();
    else
        startcode_follows = !currSlice->partArr[0].more_rbsp_data();

    if (currSlice->current_mb_nr == -1) { // End of slice_t group, MUST be end of slice
        assert(startcode_follows);
        return 1;
    }

    if (!startcode_follows)
        return 0;

    if (currSlice->slice_type == I_SLICE || currSlice->slice_type == SI_SLICE)
        return 1;
    if (p_Vid->active_pps->entropy_coding_mode_flag)
        return 1;
    if (currSlice->mb_skip_run <= 0)
        return 1;

    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for P-Frames
 ************************************************************************
 */
static void interpret_mb_mode_P(mb_t *currMB)
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
        currMB->Intra16x16PredMode = 0;

        memset(currMB->b8mode, 0, 4 * sizeof(char));
        memset(currMB->b8pdir,-1, 4 * sizeof(char));
    } else {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I16MB;
        currMB->cbp = ICBPTAB[((mbmode-7))>>2];
        currMB->Intra16x16PredMode = ((mbmode-7)) & 0x03;
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
static void interpret_mb_mode_I(mb_t *currMB)
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
        currMB->Intra16x16PredMode = 0;

        memset(currMB->b8mode, 0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I16MB;
        currMB->cbp= ICBPTAB[(mbmode-1)>>2];
        currMB->Intra16x16PredMode = (mbmode-1) & 0x03;
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
static void interpret_mb_mode_B(mb_t *currMB)
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
        currMB->Intra16x16PredMode = (mbtype-24) & 0x03;
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
        currMB->Intra16x16PredMode = 0;
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
static void interpret_mb_mode_SI(mb_t *currMB)
{
    //VideoParameters *p_Vid = currMB->p_Vid;
    const int ICBPTAB[6] = {0,16,32,15,31,47};
    short mbmode = currMB->mb_type;

    if (mbmode == 0) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = SI4MB;
        memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else if (mbmode == 1) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I4MB;
        memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else if (mbmode == 26) {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = IPCM;
        currMB->cbp = -1;
        currMB->Intra16x16PredMode = 0;
        memset(currMB->b8mode,0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    } else {
        currMB->is_intra_block = TRUE;
        currMB->mb_type = I16MB;
        currMB->cbp = ICBPTAB[(mbmode - 2) >> 2];
        currMB->Intra16x16PredMode = (mbmode - 2) & 0x03;
        memset(currMB->b8mode,0,4 * sizeof(char));
        memset(currMB->b8pdir,-1,4 * sizeof(char));
    }
}

void macroblock_t::interpret_mb_mode()
{
    slice_t *slice = this->p_Slice;
    switch (slice->slice_type) {
    case P_SLICE: 
    case SP_SLICE:
        interpret_mb_mode_P(this);
        break;
    case B_SLICE:
        interpret_mb_mode_B(this);
        break;
    case I_SLICE: 
        interpret_mb_mode_I(this);
        break;
    case SI_SLICE: 
        interpret_mb_mode_SI(this);
        break;
    default:
        printf("Unsupported slice type\n");
        break;
    }
}
