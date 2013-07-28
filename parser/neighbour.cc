
/*!
 *************************************************************************************
 * \file mb_access.c
 *
 * \brief
 *    Functions for macroblock neighborhoods
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten Suehring
 *************************************************************************************
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "neighbour.h"

/*!
 ************************************************************************
 * \brief
 *    returns 1 if the macroblock at the given address is available
 ************************************************************************
 */
bool mb_is_available(int mbAddr, mb_t *currMB)
{
  slice_t *currSlice = currMB->p_Slice;
  if ((mbAddr < 0) || (mbAddr > ((int)currSlice->dec_picture->PicSizeInMbs - 1)))
    return FALSE;

  // the following line checks both: slice number and if the mb has been decoded
  if (!currMB->DeblockCall)
  {
    if (currSlice->mb_data[mbAddr].slice_nr != currMB->slice_nr)
      return FALSE;
  }

  return TRUE;
}


/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 ************************************************************************
 */
void CheckAvailabilityOfNeighbors(mb_t *currMB)
{
  slice_t *currSlice = currMB->p_Slice;
  StorablePicture *dec_picture = currSlice->dec_picture; //p_Vid->dec_picture;
  const int mb_nr = currMB->mbAddrX;
  BlockPos *PicPos = currMB->p_Vid->PicPos;

  if (dec_picture->mb_aff_frame_flag)
  {
    int cur_mb_pair = mb_nr >> 1;
    currMB->mbAddrA = 2 * (cur_mb_pair - 1);
    currMB->mbAddrB = 2 * (cur_mb_pair - dec_picture->PicWidthInMbs);
    currMB->mbAddrC = 2 * (cur_mb_pair - dec_picture->PicWidthInMbs + 1);
    currMB->mbAddrD = 2 * (cur_mb_pair - dec_picture->PicWidthInMbs - 1);

    currMB->mbAvailA = (Boolean) (mb_is_available(currMB->mbAddrA, currMB) && ((PicPos[cur_mb_pair    ].x)!=0));
    currMB->mbAvailB = (Boolean) (mb_is_available(currMB->mbAddrB, currMB));
    currMB->mbAvailC = (Boolean) (mb_is_available(currMB->mbAddrC, currMB) && ((PicPos[cur_mb_pair + 1].x)!=0));
    currMB->mbAvailD = (Boolean) (mb_is_available(currMB->mbAddrD, currMB) && ((PicPos[cur_mb_pair    ].x)!=0));
  }
  else
  {
    BlockPos *p_pic_pos = &PicPos[mb_nr    ];
    currMB->mbAddrA = mb_nr - 1;
    currMB->mbAddrD = currMB->mbAddrA - dec_picture->PicWidthInMbs;
    currMB->mbAddrB = currMB->mbAddrD + 1;
    currMB->mbAddrC = currMB->mbAddrB + 1;


    currMB->mbAvailA = (Boolean) (mb_is_available(currMB->mbAddrA, currMB) && ((p_pic_pos->x)!=0));
    currMB->mbAvailD = (Boolean) (mb_is_available(currMB->mbAddrD, currMB) && ((p_pic_pos->x)!=0));
    currMB->mbAvailC = (Boolean) (mb_is_available(currMB->mbAddrC, currMB) && (((p_pic_pos + 1)->x)!=0));
    currMB->mbAvailB = (Boolean) (mb_is_available(currMB->mbAddrB, currMB));        
  }

  currMB->mb_left = (currMB->mbAvailA) ? &(currSlice->mb_data[currMB->mbAddrA]) : NULL;
  currMB->mb_up   = (currMB->mbAvailB) ? &(currSlice->mb_data[currMB->mbAddrB]) : NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 ************************************************************************
 */
void CheckAvailabilityOfNeighborsNormal(mb_t *currMB)
{
  slice_t *currSlice = currMB->p_Slice;
  StorablePicture *dec_picture = currSlice->dec_picture; //p_Vid->dec_picture;
  const int mb_nr = currMB->mbAddrX;
  BlockPos *PicPos = currMB->p_Vid->PicPos;

  BlockPos *p_pic_pos = &PicPos[mb_nr    ];
  currMB->mbAddrA = mb_nr - 1;
  currMB->mbAddrD = currMB->mbAddrA - dec_picture->PicWidthInMbs;
  currMB->mbAddrB = currMB->mbAddrD + 1;
  currMB->mbAddrC = currMB->mbAddrB + 1;


  currMB->mbAvailA = (Boolean) (mb_is_available(currMB->mbAddrA, currMB) && ((p_pic_pos->x)!=0));
  currMB->mbAvailD = (Boolean) (mb_is_available(currMB->mbAddrD, currMB) && ((p_pic_pos->x)!=0));
  currMB->mbAvailC = (Boolean) (mb_is_available(currMB->mbAddrC, currMB) && (((p_pic_pos + 1)->x)!=0));
  currMB->mbAvailB = (Boolean) (mb_is_available(currMB->mbAddrB, currMB));        


  currMB->mb_left = (currMB->mbAvailA) ? &(currSlice->mb_data[currMB->mbAddrA]) : NULL;
  currMB->mb_up   = (currMB->mbAvailB) ? &(currSlice->mb_data[currMB->mbAddrB]) : NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 ************************************************************************
 */
void CheckAvailabilityOfNeighborsMBAFF(mb_t *currMB)
{
  slice_t *currSlice = currMB->p_Slice;
  StorablePicture *dec_picture = currSlice->dec_picture;
  const int mb_nr = currMB->mbAddrX;
  BlockPos *PicPos = currMB->p_Vid->PicPos;

  int cur_mb_pair = mb_nr >> 1;
  currMB->mbAddrA = 2 * (cur_mb_pair - 1);
  currMB->mbAddrB = 2 * (cur_mb_pair - dec_picture->PicWidthInMbs);
  currMB->mbAddrC = 2 * (cur_mb_pair - dec_picture->PicWidthInMbs + 1);
  currMB->mbAddrD = 2 * (cur_mb_pair - dec_picture->PicWidthInMbs - 1);

  currMB->mbAvailA = (Boolean) (mb_is_available(currMB->mbAddrA, currMB) && ((PicPos[cur_mb_pair    ].x)!=0));
  currMB->mbAvailB = (Boolean) (mb_is_available(currMB->mbAddrB, currMB));
  currMB->mbAvailC = (Boolean) (mb_is_available(currMB->mbAddrC, currMB) && ((PicPos[cur_mb_pair + 1].x)!=0));
  currMB->mbAvailD = (Boolean) (mb_is_available(currMB->mbAddrD, currMB) && ((PicPos[cur_mb_pair    ].x)!=0));

  currMB->mb_left = (currMB->mbAvailA) ? &(currSlice->mb_data[currMB->mbAddrA]) : NULL;
  currMB->mb_up   = (currMB->mbAvailB) ? &(currSlice->mb_data[currMB->mbAddrB]) : NULL;
}

void CheckAvailabilityOfNeighborsCABAC(mb_t *currMB)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  PixelPos up, left;
  int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

  p_Vid->getNeighbour(currMB, -1,  0, mb_size, &left);
  p_Vid->getNeighbour(currMB,  0, -1, mb_size, &up);

  if (up.available)
    currMB->mb_up = &currMB->p_Slice->mb_data[up.mb_addr];
  else
    currMB->mb_up = NULL;

  if (left.available)
    currMB->mb_left = &currMB->p_Slice->mb_data[left.mb_addr];
  else
    currMB->mb_left = NULL;
}


/*!
 ************************************************************************
 * \brief
 *    returns the x and y macroblock coordinates for a given MbAddress
 ************************************************************************
 */
void get_mb_block_pos_normal (BlockPos *PicPos, int mb_addr, short *x, short *y)
{
  BlockPos *pPos = &PicPos[ mb_addr ];
  *x = (short) pPos->x;
  *y = (short) pPos->y;
}

/*!
 ************************************************************************
 * \brief
 *    returns the x and y macroblock coordinates for a given MbAddress
 *    for mbaff type slices
 ************************************************************************
 */
void get_mb_block_pos_mbaff (BlockPos *PicPos, int mb_addr, short *x, short *y)
{
  BlockPos *pPos = &PicPos[ mb_addr >> 1 ];
  *x = (short)  pPos->x;
  *y = (short) ((pPos->y << 1) + (mb_addr & 0x01));
}

/*!
 ************************************************************************
 * \brief
 *    returns the x and y sample coordinates for a given MbAddress
 ************************************************************************
 */
void get_mb_pos (VideoParameters *p_Vid, int mb_addr, int mb_size[2], short *x, short *y)
{
  p_Vid->get_mb_block_pos(p_Vid->PicPos, mb_addr, x, y);

  (*x) = (short) ((*x) * mb_size[0]);
  (*y) = (short) ((*y) * mb_size[1]);
}


/*!
 ************************************************************************
 * \brief
 *    get neighbouring positions for non-aff coding
 * \param currMB
 *   current macroblock
 * \param xN
 *    input x position
 * \param yN
 *    input y position
 * \param mb_size
 *    mb_t size in pixel (according to luma or chroma MB access)
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void getNonAffNeighbour(mb_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
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
            pix->available = FALSE;
    } else if (xN < maxW) {
        if (yN < 0) {
            pix->mb_addr   = currMB->mbAddrB;
            pix->available = currMB->mbAvailB;
        } else if (yN < maxH) {
            pix->mb_addr   = currMB->mbAddrX;
            pix->available = TRUE;
        } else
            pix->available = FALSE;
    } else if (xN >= maxW && yN < 0) {
        pix->mb_addr   = currMB->mbAddrC;
        pix->available = currMB->mbAvailC;
    } else
        pix->available = FALSE;

    if (pix->available || currMB->DeblockCall) {
        BlockPos *CurPos = &(currMB->p_Vid->PicPos[pix->mb_addr]);
        pix->x     = (short)(xN & (maxW - 1));
        pix->y     = (short)(yN & (maxH - 1));    
        pix->pos_x = (short)(pix->x + CurPos->x * maxW);
        pix->pos_y = (short)(pix->y + CurPos->y * maxH);
    }
}

/*!
 ************************************************************************
 * \brief
 *    get neighboring positions for aff coding
 * \param currMB
 *   current macroblock
 * \param xN
 *    input x position
 * \param yN
 *    input y position
 * \param mb_size
 *    mb_t size in pixel (according to luma or chroma MB access)
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void getAffNeighbour(mb_t *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    int maxW, maxH;
    int yM = -1;

    maxW = mb_size[0];
    maxH = mb_size[1];

    // initialize to "not available"
    pix->available = FALSE;

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
                        if (currMB->mbAvailB) {
                            if (!(currMB->DeblockCall == 1 && (p_Vid->mb_data[currMB->mbAddrB]).mb_field_decoding_flag))
                                pix->mb_addr  += 1;
                        }

                        pix->available = currMB->mbAvailB;
                        yM = yN;
                    } else { // bottom
                        pix->mb_addr   = currMB->mbAddrX - 1;
                        pix->available = TRUE;
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
                if (yN == 0 && currMB->DeblockCall == 2) {
                    pix->mb_addr  = currMB->mbAddrB + 1;
                    pix->available = TRUE;
                    yM = yN - 1;
                } else if (yN >= 0 && yN < maxH) {
                    pix->mb_addr   = currMB->mbAddrX;
                    pix->available = TRUE;
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
                        pix->available = FALSE;
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

    if (pix->available || currMB->DeblockCall) {
        pix->x = (short)(xN & (maxW - 1));
        pix->y = (short)(yM & (maxH - 1));
        get_mb_pos(p_Vid, pix->mb_addr, mb_size, &(pix->pos_x), &(pix->pos_y));
        pix->pos_x = pix->pos_x + pix->x;
        pix->pos_y = pix->pos_y + pix->y;
    }
}


/*!
 ************************************************************************
 * \brief
 *    get neighboring 4x4 block
 * \param currMB
 *   current macroblock
 * \param block_x
 *    input x block position
 * \param block_y
 *    input y block position
 * \param mb_size
 *    mb_t size in pixel (according to luma or chroma MB access)
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void get4x4Neighbour (mb_t *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix)
{
  currMB->p_Vid->getNeighbour(currMB, block_x, block_y, mb_size, pix);

  if (pix->available)
  {
    pix->x >>= 2;
    pix->y >>= 2;
    pix->pos_x >>= 2;
    pix->pos_y >>= 2;
  }
}

/*!
 ************************************************************************
 * \brief
 *    get neighboring 4x4 block
 * \param currMB
 *   current macroblock
 * \param block_x
 *    input x block position
 * \param block_y
 *    input y block position
 * \param mb_size
 *    mb_t size in pixel (according to luma or chroma MB access)
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void get4x4NeighbourBase (mb_t *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix)
{
  currMB->p_Vid->getNeighbour(currMB, block_x, block_y, mb_size, pix);

  if (pix->available)
  {
    pix->x >>= 2;
    pix->y >>= 2;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Get current block spatial neighbors
 ************************************************************************
 */
void get_neighbors(mb_t *currMB,       // <--  current mb_t
                   PixelPos   *block,     // <--> neighbor blocks
                   int         mb_x,         // <--  block x position
                   int         mb_y,         // <--  block y position
                   int         blockshape_x  // <--  block width
                   )
{
  int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
  
  get4x4Neighbour(currMB, mb_x - 1,            mb_y    , mb_size, block    );
  get4x4Neighbour(currMB, mb_x,                mb_y - 1, mb_size, block + 1);
  get4x4Neighbour(currMB, mb_x + blockshape_x, mb_y - 1, mb_size, block + 2);  

  if (mb_y > 0)
  {
    if (mb_x < 8)  // first column of 8x8 blocks
    {
      if (mb_y == 8 )
      {
        if (blockshape_x == MB_BLOCK_SIZE)      
          block[2].available  = 0;
      }
      else if (mb_x + blockshape_x == 8)
      {
        block[2].available = 0;
      }
    }
    else if (mb_x + blockshape_x == MB_BLOCK_SIZE)
    {
      block[2].available = 0;
    }
  }

  if (!block[2].available)
  {
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

    p_Vid->getNeighbour(currMB, -1,  0, mb_size, &pixA);
    p_Vid->getNeighbour(currMB,  0, -1, mb_size, &pixB);

    if (!(pps->constrained_intra_pred_flag && currMB->is_intra_block)) {
        if (pixA.available)
            currMB->dpl_flag |= p_Vid->mb_data[pixA.mb_addr].dpl_flag;
        if (pixB.available)
            currMB->dpl_flag |= p_Vid->mb_data[pixB.mb_addr].dpl_flag;
    }
}



int predict_nnz(mb_t *currMB, int pl, int i, int j)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    pps_t *pps = currSlice->active_pps;

    PixelPos pixA, pixB;
    int mb_size[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };
    int plc = sps->separate_colour_plane_flag || pl == 0 ? 0 : 1;
    get4x4Neighbour(currMB, i - 1, j, mb_size[plc], &pixA);
    get4x4Neighbour(currMB, i, j - 1, mb_size[plc], &pixB);

    bool availableFlagA = pixA.available;
    bool availableFlagB = pixB.available;
    if (pps->constrained_intra_pred_flag && currSlice->dp_mode == PAR_DP_3) {
        if (currMB->is_intra_block) {
            availableFlagA &= p_Vid->mb_data[pixA.mb_addr].is_intra_block;
            availableFlagB &= p_Vid->mb_data[pixB.mb_addr].is_intra_block;
        }
    }

    uint8_t nA = 0;
    if (availableFlagA) {
        mb_t *mb = &p_Vid->mb_data[pixA.mb_addr];
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
        mb_t *mb = &p_Vid->mb_data[pixB.mb_addr];
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
