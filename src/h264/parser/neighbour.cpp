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
 *  File      : neighbour.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "data_partition.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "neighbour.h"


using namespace vio::h264;


neighbour_t neighbour;


void neighbour_t::get_mb2pos(slice_t* slice, int mbAddr, int& xI, int& yI)
{
    sps_t* sps = slice->active_sps;

    if (slice->MbaffFrameFlag == 0) {
        xI = mbAddr % sps->PicWidthInMbs * 16;
        yI = mbAddr / sps->PicWidthInMbs * 16;
    } else {
        xI = (mbAddr / 2) % sps->PicWidthInMbs * 16;
        yI = (mbAddr / 2) / sps->PicWidthInMbs * 32;
        mb_t* mb = &slice->mb_data[mbAddr];
        if (mb->mb_field_decoding_flag == 0)
            yI += mbAddr % 2 * 16;
        else
            yI += mbAddr % 2;
    }
}

mb_t* neighbour_t::get_pos2mb(slice_t* slice, int xP, int yP, int& mbAddr)
{
    sps_t* sps = slice->active_sps;

    if (xP < 0 || xP >= sps->PicWidthInMbs * 16)
        return nullptr;
    if (yP < 0 || yP >= slice->PicHeightInMbs * 16)
        return nullptr;

    mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((yP / 16) * sps->PicWidthInMbs + (xP / 16)) :
        ((yP / 32) * sps->PicWidthInMbs + (xP / 16)) * 2;
    if (mbAddr < 0 || mbAddr >= sps->PicWidthInMbs * slice->PicHeightInMbs)
        return nullptr;

    mb_t* mb = &slice->mb_data[mbAddr];
    if (slice->MbaffFrameFlag) {
        if ((mb->mb_field_decoding_flag == 0) ? (yP & 16) : (yP & 1)) {
            ++mbAddr;
            ++mb;
        }
    }
    return mb;
}

int neighbour_t::predict_nnz(mb_t* mb, int pl, int i, int j)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    PixelPos pixA, pixB;
    int mb_size[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };
    int plc = sps->separate_colour_plane_flag || pl == 0 ? 0 : 1;
    get4x4Neighbour(mb, i - 1, j, mb_size[plc], &pixA);
    get4x4Neighbour(mb, i, j - 1, mb_size[plc], &pixB);

    bool availableFlagA = pixA.available;
    bool availableFlagB = pixB.available;
    if (pps->constrained_intra_pred_flag && slice->dp_mode == PAR_DP_3) {
        if (mb->is_intra_block) {
            availableFlagA &= slice->mb_data[pixA.mb_addr].is_intra_block;
            availableFlagB &= slice->mb_data[pixB.mb_addr].is_intra_block;
        }
    }

    uint8_t nA = 0;
    if (availableFlagA) {
        mb_t *mb = &slice->mb_data[pixA.mb_addr];
        //if (mb->mb_type == PSKIP || mb->mb_type == BSKIP_DIRECT)
        //    nA = 0;
        //else if (mb->mb_type != IPCM && (mb->cbp & 15) == 0)
        //    nA = 0;
        //else if (mb->mb_type == IPCM)
        //    nA = 16;
        //else
            nA = mb->nz_coeff[pl][pixA.y][pixA.x];
    }

    uint8_t nB = 0;
    if (availableFlagB) {
        mb_t *mb = &slice->mb_data[pixB.mb_addr];
        //if (mb->mb_type == PSKIP || mb->mb_type == BSKIP_DIRECT)
        //    nB = 0;
        //else if (mb->mb_type != IPCM && (mb->cbp & 15) == 0)
        //    nB = 0;
        //else if (mb->mb_type == IPCM)
        //    nB = 16;
        //else
            nB = mb->nz_coeff[pl][pixB.y][pixB.x];
    }

    uint8_t nC = nA + nB;
    if (availableFlagA && availableFlagB)
        nC = (nC + 1) >> 1;
    return nC;
}



static void getNonAffNeighbour(mb_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
{
    int maxW = mb_size[0], maxH = mb_size[1];

    if (xN < 0) {
        if (yN < 0) {
            pix->mb_addr   = currMB->mbAddrD;
            pix->available = currMB->mbAvailD;
        } else if (yN < maxH) {
            pix->mb_addr   = currMB->mbAddrA;
            pix->available = currMB->mbAvailA;
        } else
            pix->available = false;
    } else if (xN < maxW) {
        if (yN < 0) {
            pix->mb_addr   = currMB->mbAddrB;
            pix->available = currMB->mbAvailB;
        } else if (yN < maxH) {
            pix->mb_addr   = currMB->mbAddrX;
            pix->available = true;
        } else
            pix->available = false;
    } else if (xN >= maxW && yN < 0) {
        pix->mb_addr   = currMB->mbAddrC;
        pix->available = currMB->mbAvailC;
    } else
        pix->available = false;

    if (pix->available) {
        BlockPos *CurPos = &(currMB->p_Vid->PicPos[pix->mb_addr]);
        pix->x     = (short)(xN & (maxW - 1));
        pix->y     = (short)(yN & (maxH - 1));    
        pix->pos_x = (short)(pix->x + CurPos->x * maxW);
        pix->pos_y = (short)(pix->y + CurPos->y * maxH);
    }
}

static void getAffNeighbour(mb_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    int maxW, maxH;
    int yM = -1;

    maxW = mb_size[0];
    maxH = mb_size[1];

    // initialize to "not available"
    pix->available = false;

    if (yN > maxH - 1)
        return;
    if (xN > maxW - 1 && yN >= 0 && yN < maxH)
        return;

    if (xN < 0) {
        if (yN < 0) {
            if (!currMB->mb_field_decoding_flag) { // frame
                if ((currMB->mbAddrX & 0x01) == 0) { // top
                    pix->mb_addr   = currMB->mbAddrD  + 1;
                    pix->available = currMB->mbAvailD;
                    yM = yN;
                } else { // bottom
                    pix->mb_addr   = currMB->mbAddrA;
                    pix->available = currMB->mbAvailA;
                    if (currMB->mbAvailA) {
                        if (!p_Vid->mb_data[currMB->mbAddrA].mb_field_decoding_flag)
                            yM = yN;
                        else {
                            (pix->mb_addr)++;
                            yM = (yN + maxH) >> 1;
                        }
                    }
                }
            } else { // field
                if ((currMB->mbAddrX & 0x01) == 0) { // top
                    pix->mb_addr   = currMB->mbAddrD;
                    pix->available = currMB->mbAvailD;
                    if (currMB->mbAvailD) {
                        if (!p_Vid->mb_data[currMB->mbAddrD].mb_field_decoding_flag) {
                            (pix->mb_addr)++;
                            yM = 2 * yN;
                        } else
                            yM = yN;
                    }
                } else { // bottom
                    pix->mb_addr   = currMB->mbAddrD + 1;
                    pix->available = currMB->mbAvailD;
                    yM = yN;
                }
            }
        } else { // xN < 0 && yN >= 0
            if (yN >= 0 && yN <maxH) {
                if (!currMB->mb_field_decoding_flag) { // frame
                    if ((currMB->mbAddrX & 0x01) == 0) { // top
                        pix->mb_addr   = currMB->mbAddrA;
                        pix->available = currMB->mbAvailA;
                        if (currMB->mbAvailA) {
                            if (!p_Vid->mb_data[currMB->mbAddrA].mb_field_decoding_flag)
                                yM = yN;
                            else {
                                (pix->mb_addr)+= ((yN & 0x01) != 0);
                                yM = yN >> 1;
                            }
                        }
                    } else { // bottom
                        pix->mb_addr   = currMB->mbAddrA;
                        pix->available = currMB->mbAvailA;
                        if (currMB->mbAvailA) {
                            if (!p_Vid->mb_data[currMB->mbAddrA].mb_field_decoding_flag) {
                                (pix->mb_addr)++;
                                yM = yN;
                            } else {
                                (pix->mb_addr)+= ((yN & 0x01) != 0);
                                yM = (yN + maxH) >> 1;
                            }
                        }
                    }
                } else { // field
                    if ((currMB->mbAddrX & 0x01) == 0) { // top
                        pix->mb_addr  = currMB->mbAddrA;
                        pix->available = currMB->mbAvailA;
                        if (currMB->mbAvailA) {
                            if (!p_Vid->mb_data[currMB->mbAddrA].mb_field_decoding_flag) {
                                if (yN < (maxH >> 1))
                                    yM = yN << 1;
                                else {
                                    (pix->mb_addr)++;
                                    yM = (yN << 1 ) - maxH;
                                }
                            } else
                                yM = yN;
                        }
                    } else { // bottom
                        pix->mb_addr  = currMB->mbAddrA;
                        pix->available = currMB->mbAvailA;
                        if (currMB->mbAvailA) {
                            if (!p_Vid->mb_data[currMB->mbAddrA].mb_field_decoding_flag) {
                                if (yN < (maxH >> 1))
                                    yM = (yN << 1) + 1;
                                else {
                                    (pix->mb_addr)++;
                                    yM = (yN << 1 ) + 1 - maxH;
                                }
                            } else {
                                (pix->mb_addr)++;
                                yM = yN;
                            }
                        }
                    }
                }
            }
        }
    } else { // xN >= 0
        if (xN >= 0 && xN < maxW) {
            if (yN < 0) {
                if (!currMB->mb_field_decoding_flag) { //frame
                    if ((currMB->mbAddrX & 0x01) == 0) { //top
                        pix->mb_addr  = currMB->mbAddrB;
                        // for the deblocker if the current MB is a frame and the one above is a field
                        // then the neighbor is the top MB of the pair
                        if (currMB->mbAvailB)
                            pix->mb_addr  += 1;

                        pix->available = currMB->mbAvailB;
                        yM = yN;
                    } else { // bottom
                        pix->mb_addr   = currMB->mbAddrX - 1;
                        pix->available = true;
                        yM = yN;
                    }
                } else { // field
                    if ((currMB->mbAddrX & 0x01) == 0) { // top
                        pix->mb_addr   = currMB->mbAddrB;
                        pix->available = currMB->mbAvailB;
                        if (currMB->mbAvailB) {
                            if (!p_Vid->mb_data[currMB->mbAddrB].mb_field_decoding_flag) {
                                (pix->mb_addr)++;
                                yM = 2* yN;
                            } else
                                yM = yN;
                        }
                    } else { // bottom
                        pix->mb_addr   = currMB->mbAddrB + 1;
                        pix->available = currMB->mbAvailB;
                        yM = yN;
                    }
                }
            } else { // yN >=0
                // for the deblocker if this is the extra edge then do this special stuff
                if (yN >= 0 && yN < maxH) {
                    pix->mb_addr   = currMB->mbAddrX;
                    pix->available = true;
                    yM = yN;
                }
            }
        } else { // xN >= maxW
            if (yN < 0) {
                if (!currMB->mb_field_decoding_flag) { // frame
                    if ((currMB->mbAddrX & 0x01) == 0) { // top
                        pix->mb_addr  = currMB->mbAddrC + 1;
                        pix->available = currMB->mbAvailC;
                        yM = yN;
                    } else // bottom
                        pix->available = false;
                } else { // field
                    if ((currMB->mbAddrX & 0x01) == 0) { // top
                        pix->mb_addr   = currMB->mbAddrC;
                        pix->available = currMB->mbAvailC;
                        if (currMB->mbAvailC) {
                            if (!p_Vid->mb_data[currMB->mbAddrC].mb_field_decoding_flag) {
                                (pix->mb_addr)++;
                                yM = 2 * yN;
                            } else
                                yM = yN;
                        }
                    } else { // bottom
                        pix->mb_addr   = currMB->mbAddrC + 1;
                        pix->available = currMB->mbAvailC;
                        yM = yN;
                    }
                }
            }
        }
    }

    if (pix->available) {
        pix->x = (short)(xN & (maxW - 1));
        pix->y = (short)(yM & (maxH - 1));

        if (p_Vid->dec_picture->mb_aff_frame_flag) {
            BlockPos *pPos = &p_Vid->PicPos[pix->mb_addr >> 1];
            pix->pos_x = (short)pPos->x;
            pix->pos_y = (short)((pPos->y << 1) + (pix->mb_addr & 0x01));
        } else {
            BlockPos *pPos = &p_Vid->PicPos[pix->mb_addr];
            pix->pos_x = (short)pPos->x;
            pix->pos_y = (short)pPos->y;
        }

        pix->pos_x = pix->x + (short)(pix->pos_x * mb_size[0]);
        pix->pos_y = pix->y + (short)(pix->pos_y * mb_size[1]);
    }
}

void getNeighbour(mb_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
{
    if (currMB->p_Slice->MbaffFrameFlag)
        getAffNeighbour(currMB, xN, yN, mb_size, pix);
    else
        getNonAffNeighbour(currMB, xN, yN, mb_size, pix);
}

static bool mb_is_available(int mbAddr, mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    if (mbAddr < 0 || mbAddr > slice->PicSizeInMbs - 1)
        return false;

    // the following line checks both: slice number and if the mb has been decoded
    if (slice->mb_data[mbAddr].slice_nr != mb->slice_nr)
        return false;

    return true;
}

void CheckAvailabilityOfNeighbors(mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    BlockPos* PicPos = mb->p_Vid->PicPos;

    if (slice->MbaffFrameFlag) {
        const int mb_pair = mb->mbAddrX >> 1;
        mb->mbAddrA = 2 * (mb_pair - 1);
        mb->mbAddrB = 2 * (mb_pair - sps->PicWidthInMbs);
        mb->mbAddrC = 2 * (mb_pair - sps->PicWidthInMbs + 1);
        mb->mbAddrD = 2 * (mb_pair - sps->PicWidthInMbs - 1);

        mb->mbAvailA = mb_is_available(mb->mbAddrA, mb) && (PicPos[mb_pair    ].x != 0);
        mb->mbAvailB = mb_is_available(mb->mbAddrB, mb);
        mb->mbAvailC = mb_is_available(mb->mbAddrC, mb) && (PicPos[mb_pair + 1].x != 0);
        mb->mbAvailD = mb_is_available(mb->mbAddrD, mb) && (PicPos[mb_pair    ].x != 0);
    } else {
        const int mb_nr = mb->mbAddrX;
        mb->mbAddrA = mb_nr - 1;
        mb->mbAddrB = mb_nr - sps->PicWidthInMbs;
        mb->mbAddrC = mb_nr - sps->PicWidthInMbs + 1;
        mb->mbAddrD = mb_nr - sps->PicWidthInMbs - 1;

        mb->mbAvailA = mb_is_available(mb->mbAddrA, mb) && (PicPos[mb_nr    ].x != 0);
        mb->mbAvailB = mb_is_available(mb->mbAddrB, mb);
        mb->mbAvailC = mb_is_available(mb->mbAddrC, mb) && (PicPos[mb_nr + 1].x != 0);
        mb->mbAvailD = mb_is_available(mb->mbAddrD, mb) && (PicPos[mb_nr    ].x != 0);
    }

    mb->mb_left = mb->mbAvailA ? &slice->mb_data[mb->mbAddrA] : NULL;
    mb->mb_up   = mb->mbAvailB ? &slice->mb_data[mb->mbAddrB] : NULL;
}

void CheckAvailabilityOfNeighborsCABAC(mb_t *currMB)
{
    PixelPos up, left;
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    getNeighbour(currMB, -1,  0, mb_size, &left);
    getNeighbour(currMB,  0, -1, mb_size, &up);

    if (up.available)
        currMB->mb_up = &currMB->p_Slice->mb_data[up.mb_addr];
    else
        currMB->mb_up = NULL;

    if (left.available)
        currMB->mb_left = &currMB->p_Slice->mb_data[left.mb_addr];
    else
        currMB->mb_left = NULL;
}



void get4x4Neighbour(mb_t *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix)
{
    getNeighbour(currMB, block_x, block_y, mb_size, pix);

    if (pix->available) {
        pix->x >>= 2;
        pix->y >>= 2;
        pix->pos_x >>= 2;
        pix->pos_y >>= 2;
    }
}

void get4x4NeighbourBase(mb_t *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix)
{
    getNeighbour(currMB, block_x, block_y, mb_size, pix);

    if (pix->available) {
        pix->x >>= 2;
        pix->y >>= 2;
    }
}


void get_neighbors(mb_t *currMB, PixelPos *block, int mb_x, int mb_y, int blockshape_x)
{
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
  
    get4x4Neighbour(currMB, mb_x - 1,            mb_y    , mb_size, block    );
    get4x4Neighbour(currMB, mb_x,                mb_y - 1, mb_size, block + 1);
    get4x4Neighbour(currMB, mb_x + blockshape_x, mb_y - 1, mb_size, block + 2);  

    if (mb_y > 0) {
        if (mb_x < 8) { // first column of 8x8 blocks
            if (mb_y == 8) {
                if (blockshape_x == MB_BLOCK_SIZE)      
                    block[2].available  = 0;
            } else if (mb_x + blockshape_x == 8)
                block[2].available = 0;
        } else if (mb_x + blockshape_x == MB_BLOCK_SIZE)
            block[2].available = 0;
    }

    if (!block[2].available) {
        get4x4Neighbour(currMB, mb_x - 1, mb_y - 1, mb_size, block + 3);
        block[2] = block[3];
    }
}


void check_dp_neighbors(mb_t *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;

    PixelPos pixA, pixB;
    int mb_size[2] = { sps->MbWidthC, sps->MbHeightC };

    getNeighbour(currMB, -1,  0, mb_size, &pixA);
    getNeighbour(currMB,  0, -1, mb_size, &pixB);

    if (!(pps->constrained_intra_pred_flag && currMB->is_intra_block)) {
        if (pixA.available)
            currMB->dpl_flag |= p_Vid->mb_data[pixA.mb_addr].dpl_flag;
        if (pixB.available)
            currMB->dpl_flag |= p_Vid->mb_data[pixB.mb_addr].dpl_flag;
    }
}

