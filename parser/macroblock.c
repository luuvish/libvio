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

#include "block.h"
#include "global.h"
#include "mbuffer.h"
#include "mbuffer_mvc.h"
#include "bitstream_elements.h"
#include "bitstream_cabac.h"
#include "bitstream_vlc.h"
#include "macroblock.h"
#include "mb_read.h"
#include "fmo.h"
#include "image.h"
#include "neighbour.h"
#include "biaridecod.h"
#include "transform.h"
#include "mc_prediction.h"
#include "quant.h"
#include "mv_prediction.h"
#include "mb_prediction.h"
#include "fast_memory.h"
#include "filehandle.h"

#include "mb.h"


extern void update_direct_types                (Slice *currSlice);
extern void set_intra_prediction_modes         (Slice *currSlice);

void set_chroma_qp(Macroblock* currMB)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  StorablePicture *dec_picture = currMB->p_Slice->dec_picture;
  int i;

  for (i=0; i<2; ++i)
  {
    currMB->qpc[i] = iClip3 ( -p_Vid->bitdepth_chroma_qp_scale, 51, currMB->qp + dec_picture->chroma_qp_offset[i] );
    currMB->qpc[i] = currMB->qpc[i] < 0 ? currMB->qpc[i] : QP_SCALE_CR[currMB->qpc[i]];
    currMB->qp_scaled[i + 1] = currMB->qpc[i] + p_Vid->bitdepth_chroma_qp_scale;
  }
}

/*!
************************************************************************
* \brief
*    updates chroma QP according to luma QP and bit depth
************************************************************************
*/
void update_qp(Macroblock *currMB, int qp)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  currMB->qp = qp;
  currMB->qp_scaled[0] = qp + p_Vid->bitdepth_luma_qp_scale;
  set_chroma_qp(currMB);
  currMB->is_lossless = (Boolean) ((currMB->qp_scaled[0] == 0) && (p_Vid->lossless_qpprime_flag == 1));
  set_read_comp_coeff_cavlc(currMB);
  set_read_comp_coeff_cabac(currMB);
}


void invScaleCoeff(Macroblock *currMB, int level, int run, int qp_per, int i, int j, int i0, int j0, int coef_ctr, const byte (*pos_scan4x4)[2], int (*InvLevelScale4x4)[4])
{
  if (level != 0)    /* leave if level == 0 */
  {
    coef_ctr += run + 1;

    i0 = pos_scan4x4[coef_ctr][0];
    j0 = pos_scan4x4[coef_ctr][1];

    currMB->s_cbp[0].blk |= i64_power2((j << 2) + i) ;
    currMB->p_Slice->cof[0][(j<<2) + j0][(i<<2) + i0]= rshift_rnd_sf((level * InvLevelScale4x4[j0][i0]) << qp_per, 4);
  }
}

static inline void setup_mb_pos_info(Macroblock *currMB)
{
  int mb_x = currMB->mb.x;
  int mb_y = currMB->mb.y;
  currMB->block_x     = mb_x << BLOCK_SHIFT;           /* horizontal block position */
  currMB->block_y     = mb_y << BLOCK_SHIFT;           /* vertical block position */
  currMB->block_y_aff = currMB->block_y;                       /* interlace relative vertical position */
  currMB->pix_x       = mb_x << MB_BLOCK_SHIFT;        /* horizontal luma pixel position */
  currMB->pix_y       = mb_y << MB_BLOCK_SHIFT;        /* vertical luma pixel position */
  currMB->pix_c_x     = mb_x * currMB->p_Vid->mb_cr_size_x;    /* horizontal chroma pixel position */
  currMB->pix_c_y     = mb_y * currMB->p_Vid->mb_cr_size_y;    /* vertical chroma pixel position */
}

/*!
 ************************************************************************
 * \brief
 *    initializes the current macroblock
 ************************************************************************
 */
void start_macroblock(Slice *currSlice, Macroblock **currMB)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  int mb_nr = currSlice->current_mb_nr;
  
  *currMB = &currSlice->mb_data[mb_nr]; 

  (*currMB)->p_Slice = currSlice;
  (*currMB)->p_Vid   = p_Vid;  
  (*currMB)->mbAddrX = mb_nr;

  //assert (mb_nr < (int) p_Vid->PicSizeInMbs);

  /* Update coordinates of the current macroblock */
  if (currSlice->mb_aff_frame_flag)
  {
    (*currMB)->mb.x = (short) (   (mb_nr) % ((2*p_Vid->width) / MB_BLOCK_SIZE));
    (*currMB)->mb.y = (short) (2*((mb_nr) / ((2*p_Vid->width) / MB_BLOCK_SIZE)));

    (*currMB)->mb.y += ((*currMB)->mb.x & 0x01);
    (*currMB)->mb.x >>= 1;
  }
  else
  {
    (*currMB)->mb = p_Vid->PicPos[mb_nr];
  }

  /* Define pixel/block positions */
  setup_mb_pos_info(*currMB);

  // reset intra mode
  (*currMB)->is_intra_block = FALSE;
  // reset mode info
  (*currMB)->mb_type         = 0;
  (*currMB)->delta_quant     = 0;
  (*currMB)->cbp             = 0;    
  (*currMB)->c_ipred_mode    = DC_PRED_8;

  // Save the slice number of this macroblock. When the macroblock below
  // is coded it will use this to decide if prediction for above is possible
  (*currMB)->slice_nr = (short) currSlice->current_slice_nr;

  CheckAvailabilityOfNeighbors(*currMB);

  // Select appropriate MV predictor function
  init_motion_vector_prediction(*currMB, currSlice->mb_aff_frame_flag);

  set_read_and_store_CBP(currMB, currSlice->active_sps->chroma_format_idc);

  // Reset syntax element entries in MB struct

  if (currSlice->slice_type != I_SLICE)
  {
    if (currSlice->slice_type != B_SLICE)
      fast_memset((*currMB)->mvd[0][0][0], 0, MB_BLOCK_PARTITIONS * 2 * sizeof(short));
    else
      fast_memset((*currMB)->mvd[0][0][0], 0, 2 * MB_BLOCK_PARTITIONS * 2 * sizeof(short));
  }
  
  fast_memset((*currMB)->s_cbp, 0, 3 * sizeof(CBPStructure));

  // initialize currSlice->mb_rres
  if (currSlice->is_reset_coeff == FALSE)
  {
    fast_memset_zero( currSlice->mb_rres[0][0], MB_PIXELS * sizeof(int));
    fast_memset_zero( currSlice->mb_rres[1][0], p_Vid->mb_cr_size * sizeof(int));
    fast_memset_zero( currSlice->mb_rres[2][0], p_Vid->mb_cr_size * sizeof(int));
    if (currSlice->is_reset_coeff_cr == FALSE)
    {
      fast_memset_zero( currSlice->cof[0][0], 3 * MB_PIXELS * sizeof(int));
      currSlice->is_reset_coeff_cr = TRUE;
    }
    else
    {
      fast_memset_zero( currSlice->cof[0][0], MB_PIXELS * sizeof(int));
    }
    //fast_memset_zero( currSlice->cof[0][0], MB_PIXELS * sizeof(int));
    //fast_memset_zero( currSlice->cof[1][0], p_Vid->mb_cr_size * sizeof(int));
    //fast_memset_zero( currSlice->cof[2][0], p_Vid->mb_cr_size * sizeof(int));

    //fast_memset(currSlice->fcf[0][0], 0, MB_PIXELS * sizeof(int)); // reset luma coeffs   
    //fast_memset(currSlice->fcf[1][0], 0, MB_PIXELS * sizeof(int));
    //fast_memset(currSlice->fcf[2][0], 0, MB_PIXELS * sizeof(int));

    currSlice->is_reset_coeff = TRUE;
  }

  // store filtering parameters for this MB
  (*currMB)->DFDisableIdc    = currSlice->DFDisableIdc;
  (*currMB)->DFAlphaC0Offset = currSlice->DFAlphaC0Offset;
  (*currMB)->DFBetaOffset    = currSlice->DFBetaOffset;
  (*currMB)->list_offset     = 0;
  (*currMB)->mixedModeEdgeFlag = 0;
}

/*!
 ************************************************************************
 * \brief
 *    set coordinates of the next macroblock
 *    check end_of_slice condition
 ************************************************************************
 */
Boolean exit_macroblock(Slice *currSlice, int eos_bit)
{
  VideoParameters *p_Vid = currSlice->p_Vid;

 //! The if() statement below resembles the original code, which tested
  //! p_Vid->current_mb_nr == p_Vid->PicSizeInMbs.  Both is, of course, nonsense
  //! In an error prone environment, one can only be sure to have a new
  //! picture by checking the tr of the next slice header!

// printf ("exit_macroblock: FmoGetLastMBOfPicture %d, p_Vid->current_mb_nr %d\n", FmoGetLastMBOfPicture(), p_Vid->current_mb_nr);
  ++(currSlice->num_dec_mb);

  if(currSlice->current_mb_nr == p_Vid->PicSizeInMbs - 1) //if (p_Vid->num_dec_mb == p_Vid->PicSizeInMbs)
  {
    return TRUE;
  }
  // ask for last mb in the slice  CAVLC
  else
  {

    currSlice->current_mb_nr = FmoGetNextMBNr (p_Vid, currSlice->current_mb_nr);

    if (currSlice->current_mb_nr == -1)     // End of Slice group, MUST be end of slice
    {
      assert (currSlice->nal_startcode_follows (currSlice, eos_bit) == TRUE);
      return TRUE;
    }

    if(currSlice->nal_startcode_follows(currSlice, eos_bit) == FALSE)
      return FALSE;

    if(currSlice->slice_type == I_SLICE  || currSlice->slice_type == SI_SLICE || p_Vid->active_pps->entropy_coding_mode_flag == (Boolean) CABAC)
      return TRUE;
    if(currSlice->cod_counter <= 0)
      return TRUE;
    return FALSE;
  }
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

  if(mbmode < 4)
  {
    currMB->mb_type = mbmode;
    memset(currMB->b8mode, mbmode, 4 * sizeof(char));
    memset(currMB->b8pdir, 0, 4 * sizeof(char));
  }
  else if((mbmode == 4 || mbmode == 5))
  {
    currMB->mb_type = P8x8;
    currMB->p_Slice->allrefzero = (mbmode == 5);
  }
  else if(mbmode == 6)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type = I4MB;
    memset(currMB->b8mode, IBLOCK, 4 * sizeof(char));
    memset(currMB->b8pdir,     -1, 4 * sizeof(char));
  }
  else if(mbmode == 31)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type = IPCM;
    currMB->cbp = -1;
    currMB->i16mode = 0;

    memset(currMB->b8mode, 0, 4 * sizeof(char));
    memset(currMB->b8pdir,-1, 4 * sizeof(char));
  }
  else
  {
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
  short mbmode   = currMB->mb_type;

  if (mbmode == 0)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type = I4MB;
    memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));
  }
  else if(mbmode == 25)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type=IPCM;
    currMB->cbp= -1;
    currMB->i16mode = 0;

    memset(currMB->b8mode, 0,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));
  }
  else
  {
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
  short mbtype  = currMB->mb_type;

  //--- set mbtype, b8type, and b8pdir ---
  if (mbtype == 0)       // direct
  {
    mbmode=0;
    memset(currMB->b8mode, 0, 4 * sizeof(char));
    memset(currMB->b8pdir, 2, 4 * sizeof(char));
  }
  else if (mbtype == 23) // intra4x4
  {
    currMB->is_intra_block = TRUE;
    mbmode=I4MB;
    memset(currMB->b8mode, IBLOCK,4 * sizeof(char));
    memset(currMB->b8pdir, -1,4 * sizeof(char));
  }
  else if ((mbtype > 23) && (mbtype < 48) ) // intra16x16
  {
    currMB->is_intra_block = TRUE;
    mbmode=I16MB;
    memset(currMB->b8mode,  0, 4 * sizeof(char));
    memset(currMB->b8pdir, -1, 4 * sizeof(char));

    currMB->cbp     = (int) ICBPTAB[(mbtype-24)>>2];
    currMB->i16mode = (mbtype-24) & 0x03;
  }
  else if (mbtype == 22) // 8x8(+split)
  {
    mbmode=P8x8;       // b8mode and pdir is transmitted in additional codewords
  }
  else if (mbtype < 4)   // 16x16
  {
    mbmode = 1;
    memset(currMB->b8mode, 1,4 * sizeof(char));
    memset(currMB->b8pdir, offset2pdir16x16[mbtype], 4 * sizeof(char));
  }
  else if(mbtype == 48)
  {
    currMB->is_intra_block = TRUE;
    mbmode=IPCM;
    memset(currMB->b8mode, 0,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));

    currMB->cbp= -1;
    currMB->i16mode = 0;
  }
  else if ((mbtype & 0x01) == 0) // 16x8
  {
    mbmode = 2;
    memset(currMB->b8mode, 2,4 * sizeof(char));
    for(i=0;i<4;++i)
    {
      currMB->b8pdir[i] = offset2pdir16x8 [mbtype][i>>1];
    }
  }
  else
  {
    mbmode=3;
    memset(currMB->b8mode, 3,4 * sizeof(char));
    for(i=0;i<4; ++i)
    {
      currMB->b8pdir[i] = offset2pdir8x16 [mbtype][i&0x01];
    }
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
  short mbmode   = currMB->mb_type;

  if (mbmode == 0)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type = SI4MB;
    memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));
    currMB->p_Slice->siblock[currMB->mb.y][currMB->mb.x]=1;
  }
  else if (mbmode == 1)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type = I4MB;
    memset(currMB->b8mode,IBLOCK,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));
  }
  else if(mbmode == 26)
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type=IPCM;
    currMB->cbp= -1;
    currMB->i16mode = 0;
    memset(currMB->b8mode,0,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));
  }

  else
  {
    currMB->is_intra_block = TRUE;
    currMB->mb_type = I16MB;
    currMB->cbp= ICBPTAB[(mbmode-2)>>2];
    currMB->i16mode = (mbmode-2) & 0x03;
    memset(currMB->b8mode,0,4 * sizeof(char));
    memset(currMB->b8pdir,-1,4 * sizeof(char));
  }
}


/*!
 ************************************************************************
 * \brief
 *    Set mode interpretation based on slice type
 ************************************************************************
 */
void setup_slice_methods(Slice *currSlice)
{
  setup_read_macroblock (currSlice);

  switch (currSlice->slice_type)
  {
  case P_SLICE: 
    currSlice->interpret_mb_mode         = interpret_mb_mode_P;
    currSlice->decode_one_component      = decode_one_component_p_slice;
    currSlice->update_direct_mv_info     = NULL;
#if (MVC_EXTENSION_ENABLE)
    currSlice->init_lists                = currSlice->view_id ? init_lists_p_slice_mvc : init_lists_p_slice;
#else
    currSlice->init_lists                = init_lists_p_slice;
#endif
    break;
  case SP_SLICE:
    currSlice->interpret_mb_mode         = interpret_mb_mode_P;
    currSlice->decode_one_component      = decode_one_component_sp_slice;
    currSlice->update_direct_mv_info     = NULL;
#if (MVC_EXTENSION_ENABLE)
    currSlice->init_lists                = currSlice->view_id ? init_lists_p_slice_mvc : init_lists_p_slice;
#else
    currSlice->init_lists                = init_lists_p_slice;
#endif
    break;
  case B_SLICE:
    currSlice->interpret_mb_mode         = interpret_mb_mode_B;
    currSlice->decode_one_component      = decode_one_component_b_slice;
    update_direct_types(currSlice);
#if (MVC_EXTENSION_ENABLE)
    currSlice->init_lists                = currSlice->view_id ? init_lists_b_slice_mvc : init_lists_b_slice;
#else
    currSlice->init_lists                = init_lists_b_slice;
#endif
    break;
  case I_SLICE: 
    currSlice->interpret_mb_mode         = interpret_mb_mode_I;
    currSlice->decode_one_component      = decode_one_component_i_slice;
    currSlice->update_direct_mv_info     = NULL;
#if (MVC_EXTENSION_ENABLE)
    currSlice->init_lists                = currSlice->view_id ? init_lists_i_slice_mvc : init_lists_i_slice;
#else
    currSlice->init_lists                = init_lists_i_slice;
#endif
    break;
  case SI_SLICE: 
    currSlice->interpret_mb_mode         = interpret_mb_mode_SI;
    currSlice->decode_one_component      = decode_one_component_i_slice;
    currSlice->update_direct_mv_info     = NULL;
#if (MVC_EXTENSION_ENABLE)
    currSlice->init_lists                = currSlice->view_id ? init_lists_i_slice_mvc : init_lists_i_slice;
#else
    currSlice->init_lists                = init_lists_i_slice;
#endif
    break;
  default:
    printf("Unsupported slice type\n");
    break;
  }

  set_intra_prediction_modes(currSlice);
}

