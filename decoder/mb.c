
#include <math.h>

#include "block.h"
#include "global.h"
#include "slice.h"
#include "mbuffer.h"
#include "mbuffer_mvc.h"
#include "bitstream_elements.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "macroblock.h"
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

/*!
 ************************************************************************
 * \brief
 *    decode one color component in an I slice
 ************************************************************************
 */

int decode_one_component_i_slice(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  //For residual DPCM
  currMB->ipmode_DPCM = NO_INTRA_PMODE; 
  if(currMB->mb_type == IPCM)
    mb_pred_ipcm(currMB);
  else if (currMB->mb_type==I16MB)
    mb_pred_intra16x16(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == I4MB)
    mb_pred_intra4x4(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == I8MB) 
    mb_pred_intra8x8(currMB, curr_plane, currImg, dec_picture);

  return 1;
}

/*!
 ************************************************************************
 * \brief
 *    decode one color component for a p slice
 ************************************************************************
 */
int decode_one_component_p_slice(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  //For residual DPCM
  currMB->ipmode_DPCM = NO_INTRA_PMODE; 
  if(currMB->mb_type == IPCM)
    mb_pred_ipcm(currMB);
  else if (currMB->mb_type==I16MB)
    mb_pred_intra16x16(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == I4MB)
    mb_pred_intra4x4(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == I8MB) 
    mb_pred_intra8x8(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == PSKIP)
    mb_pred_skip(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == P16x16)
    mb_pred_p_inter16x16(currMB, curr_plane, dec_picture);  
  else if (currMB->mb_type == P16x8)
    mb_pred_p_inter16x8(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == P8x16)
    mb_pred_p_inter8x16(currMB, curr_plane, dec_picture);
  else
    mb_pred_p_inter8x8(currMB, curr_plane, dec_picture);

  return 1;
}


/*!
 ************************************************************************
 * \brief
 *    decode one color component for a sp slice
 ************************************************************************
 */
int decode_one_component_sp_slice(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{   
  //For residual DPCM
  currMB->ipmode_DPCM = NO_INTRA_PMODE; 

  if (currMB->mb_type == IPCM)
    mb_pred_ipcm(currMB);
  else if (currMB->mb_type==I16MB)
    mb_pred_intra16x16(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == I4MB)
    mb_pred_intra4x4(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == I8MB) 
    mb_pred_intra8x8(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == PSKIP)
    mb_pred_sp_skip(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == P16x16)
    mb_pred_p_inter16x16(currMB, curr_plane, dec_picture);  
  else if (currMB->mb_type == P16x8)
    mb_pred_p_inter16x8(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == P8x16)
    mb_pred_p_inter8x16(currMB, curr_plane, dec_picture);
  else
    mb_pred_p_inter8x8(currMB, curr_plane, dec_picture);

  return 1;
}



/*!
 ************************************************************************
 * \brief
 *    decode one color component for a b slice
 ************************************************************************
 */

int decode_one_component_b_slice(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{  
  //For residual DPCM
  currMB->ipmode_DPCM = NO_INTRA_PMODE; 

  if(currMB->mb_type == IPCM)
    mb_pred_ipcm(currMB);
  else if (currMB->mb_type==I16MB)
    mb_pred_intra16x16(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == I4MB)
    mb_pred_intra4x4(currMB, curr_plane, currImg, dec_picture);
  else if (currMB->mb_type == I8MB) 
    mb_pred_intra8x8(currMB, curr_plane, currImg, dec_picture);  
  else if (currMB->mb_type == P16x16)
    mb_pred_p_inter16x16(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == P16x8)
    mb_pred_p_inter16x8(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == P8x16)
    mb_pred_p_inter8x16(currMB, curr_plane, dec_picture);
  else if (currMB->mb_type == BSKIP_DIRECT)
  {
    Slice *currSlice = currMB->p_Slice;
    if (currSlice->direct_spatial_mv_pred_flag == 0)
    {
      if (currSlice->active_sps->direct_8x8_inference_flag)
        mb_pred_b_d8x8temporal (currMB, curr_plane, currImg, dec_picture);
      else
        mb_pred_b_d4x4temporal (currMB, curr_plane, currImg, dec_picture);
    }
    else
    {
      if (currSlice->active_sps->direct_8x8_inference_flag)
        mb_pred_b_d8x8spatial (currMB, curr_plane, currImg, dec_picture);
      else
        mb_pred_b_d4x4spatial (currMB, curr_plane, currImg, dec_picture);
    }
  }
  else
    mb_pred_b_inter8x8 (currMB, curr_plane, dec_picture);

 return 1;
}

// probably a better way (or place) to do this, but I'm not sure what (where) it is [CJV]
// this is intended to make get_block_luma faster, but I'm still performing
// this at the MB level, and it really should be done at the slice level
static void init_cur_imgy(VideoParameters *p_Vid,Slice *currSlice,int pl)
{
  int i,j;
  if (p_Vid->separate_colour_plane_flag == 0)
  {
    StorablePicture *vidref = p_Vid->no_reference_picture;
    int noref = (currSlice->framepoc < p_Vid->recovery_poc);    
    if (pl==PLANE_Y) 
    {
      for (j = 0; j < 6; j++)    // for (j = 0; j < (currSlice->slice_type==B_SLICE?2:1); j++) 
      {
        for (i = 0; i < currSlice->listXsize[j] ; i++) 
        {
          StorablePicture *curr_ref = currSlice->listX[j][i];
          if (curr_ref) 
          {
            curr_ref->no_ref = noref && (curr_ref == vidref);
            curr_ref->cur_imgY = curr_ref->imgY;
          }
        }
      }
    }
    else 
    {
      for (j = 0; j < 6; j++)  //for (j = 0; j < (currSlice->slice_type==B_SLICE?2:1); j++)
      {
        for (i = 0; i < currSlice->listXsize[j]; i++) 
        {
          StorablePicture *curr_ref = currSlice->listX[j][i];
          if (curr_ref) 
          {
            curr_ref->no_ref = noref && (curr_ref == vidref);
            curr_ref->cur_imgY = curr_ref->imgUV[pl-1]; 
          }
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    decode one macroblock
 ************************************************************************
 */

int decode_one_macroblock(Macroblock *currMB, StorablePicture *dec_picture)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;  

  // macroblock decoding **************************************************
  if (currSlice->chroma444_not_separate)  
  {
    if (!currMB->is_intra_block)
    {
      init_cur_imgy(p_Vid, currSlice, PLANE_Y);
      currSlice->decode_one_component(currMB, PLANE_Y, dec_picture->imgY, dec_picture);
      init_cur_imgy(p_Vid, currSlice, PLANE_U);
      currSlice->decode_one_component(currMB, PLANE_U, dec_picture->imgUV[0], dec_picture);
      init_cur_imgy(p_Vid, currSlice, PLANE_V);
      currSlice->decode_one_component(currMB, PLANE_V, dec_picture->imgUV[1], dec_picture);
    }
    else
    {
      currSlice->decode_one_component(currMB, PLANE_Y, dec_picture->imgY, dec_picture);
      currSlice->decode_one_component(currMB, PLANE_U, dec_picture->imgUV[0], dec_picture);
      currSlice->decode_one_component(currMB, PLANE_V, dec_picture->imgUV[1], dec_picture);      
    }
    currSlice->is_reset_coeff = FALSE;
    currSlice->is_reset_coeff_cr = FALSE;
  }
  else
  {
    currSlice->decode_one_component(currMB, PLANE_Y, dec_picture->imgY, dec_picture);
  }

  return 0;
}
