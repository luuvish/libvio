
/*!
 *************************************************************************************
 * \file mc_direct.c
 *
 * \brief
 *    Direct Prediction functions
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *      - Yuwen He                 <yhe@dolby.com>
 *
 *************************************************************************************
 */
#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "transform.h"
#include "inter_prediction.h"
#include "intra_prediction.h"
#include "neighbour.h"
#include "mv_prediction.h"


#define INVALIDINDEX  (-135792468)


static inline int RSD(int x)
{
 return ((x&2)?(x|1):(x&(~1)));
}

//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

static void update_direct_mv_info_temporal(Macroblock *currMB)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;
  int j,k;
  int partmode        = ((currMB->mb_type == P8x8) ? 4 : currMB->mb_type);
  int step_h0         = BLOCK_STEP [partmode][0];
  int step_v0         = BLOCK_STEP [partmode][1];

  int i0, j0, j6;

  int j4, i4;
  StorablePicture *dec_picture = currSlice->dec_picture;

  int list_offset = currMB->list_offset;
  StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
  StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

  Boolean has_direct = (currMB->b8mode[0] == 0) | (currMB->b8mode[1] == 0) | (currMB->b8mode[2] == 0) | (currMB->b8mode[3] == 0);

  if (has_direct)
  {
    int mv_scale = 0;

    for (k = 0; k < 4; ++k) // Scan all blocks
    {
      if (currMB->b8mode[k] == 0)
      {
        currMB->b8pdir[k] = 2;
        for(j0 = 2 * (k >> 1); j0 < 2 * (k >> 1) + 2; j0 += step_v0)
        {
          for(i0 = currMB->block_x + 2*(k & 0x01); i0 < currMB->block_x + 2 * (k & 0x01)+2; i0 += step_h0)
          {
            PicMotionParams *colocated;
            int refList;
            int ref_idx;
            int mapped_idx = -1, iref;

            colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
              &list1[0]->mv_info[RSD(currMB->block_y_aff + j0)][RSD(i0)] : &list1[0]->mv_info[currMB->block_y_aff + j0][i0];
            if(currSlice->MbaffFrameFlag)
            {
              assert(p_Vid->active_sps->direct_8x8_inference_flag);
              if(!currMB->mb_field_decoding_flag && ((currSlice->listX[LIST_1][0]->iCodingType==FRAME_MB_PAIR_CODING && currSlice->listX[LIST_1][0]->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
                (currSlice->listX[LIST_1][0]->iCodingType==FIELD_CODING)))
              {
                if (iabs(dec_picture->poc - currSlice->listX[LIST_1+4][0]->poc)> iabs(dec_picture->poc -currSlice->listX[LIST_1+2][0]->poc) )
                {
                  colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                    &currSlice->listX[LIST_1+2][0]->mv_info[RSD(currMB->block_y_aff + j0)>>1][RSD(i0)] : &currSlice->listX[LIST_1+2][0]->mv_info[(currMB->block_y_aff + j0)>>1][i0];
                }
                else
                {
                  colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                    &currSlice->listX[LIST_1+4][0]->mv_info[RSD(currMB->block_y_aff + j0)>>1][RSD(i0)] : &currSlice->listX[LIST_1+4][0]->mv_info[(currMB->block_y_aff + j0)>>1][i0];
                }
              }
            }
            else if(!p_Vid->active_sps->frame_mbs_only_flag && !currSlice->field_pic_flag && currSlice->listX[LIST_1][0]->iCodingType != FRAME_CODING)
            {
              if (iabs(dec_picture->poc - list1[0]->bottom_field->poc)> iabs(dec_picture->poc -list1[0]->top_field->poc) )
              {
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                  &list1[0]->top_field->mv_info[RSD(currMB->block_y_aff + j0)>>1][RSD(i0)] : &list1[0]->top_field->mv_info[(currMB->block_y_aff + j0)>>1][i0];
              }
              else
              {
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                  &list1[0]->bottom_field->mv_info[RSD(currMB->block_y_aff + j0)>>1][RSD(i0)] : &list1[0]->bottom_field->mv_info[(currMB->block_y_aff + j0)>>1][i0];
              }
            }
            else if(!p_Vid->active_sps->frame_mbs_only_flag && currSlice->field_pic_flag && currSlice->structure!=list1[0]->structure && list1[0]->coded_frame)
            {
              if (!currSlice->bottom_field_flag)
              {
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                  &list1[0]->frame->top_field->mv_info[RSD(currMB->block_y_aff + j0)][RSD(i0)] : &list1[0]->frame->top_field->mv_info[currMB->block_y_aff + j0][i0];
              }
              else
              {
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                  &list1[0]->frame->bottom_field->mv_info[RSD(currMB->block_y_aff + j0)][RSD(i0)] : &list1[0]->frame->bottom_field->mv_info[currMB->block_y_aff + j0][i0];
              }
            }

            refList = colocated->ref_idx[LIST_0 ]== -1 ? LIST_1 : LIST_0;
            ref_idx = colocated->ref_idx[refList];

            if (ref_idx == -1)
            {
              for (j4 = currMB->block_y + j0; j4 < currMB->block_y + j0 + step_v0; ++j4)
              {
                for (i4 = i0; i4 < i0 + step_h0; ++i4)
                {
                  PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];
                  mv_info->ref_pic[LIST_0] = list0[0];
                  mv_info->ref_pic[LIST_1] = list1[0];
                  mv_info->mv [LIST_0] = zero_mv;
                  mv_info->mv [LIST_1] = zero_mv;
                  mv_info->ref_idx [LIST_0] = 0;
                  mv_info->ref_idx [LIST_1] = 0;
                }
              }
            }
            else
            {
              if( (currSlice->MbaffFrameFlag && ( (currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure==FRAME) || 
                (!currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure!=FRAME))) ||
                (!currSlice->MbaffFrameFlag && ((currSlice->field_pic_flag==0 && colocated->ref_pic[refList]->structure!=FRAME)||
                (currSlice->field_pic_flag==1 && colocated->ref_pic[refList]->structure==FRAME))) )
              {
                //! Frame with field co-located
                for (iref = 0; iref < imin(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]); ++iref)
                {
                  if (currSlice->listX[LIST_0 + list_offset][iref]->top_field == colocated->ref_pic[refList] ||
                    currSlice->listX[LIST_0 + list_offset][iref]->bottom_field == colocated->ref_pic[refList] ||
                    currSlice->listX[LIST_0 + list_offset][iref]->frame == colocated->ref_pic[refList] ) 
                  {
                    if (currSlice->field_pic_flag && (currSlice->listX[LIST_0 + list_offset][iref]->structure != currSlice->structure))
                    {
                      mapped_idx=INVALIDINDEX;
                    }
                    else
                    {
                      mapped_idx = iref;            
                      break;
                    }
                  }
                  else //! invalid index. Default to zero even though this case should not happen
                    mapped_idx=INVALIDINDEX;
                }
              }
              else
              {
                for (iref = 0; iref < imin(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]); ++iref)
                {
                  if (currSlice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList])
                  {
                    mapped_idx=iref;
                    break;
                  }
                  else //! invalid index. Default to zero even though this case should not happen
                    mapped_idx=INVALIDINDEX;
                }
              }

              if (mapped_idx != INVALIDINDEX)
              {
                for (j = j0; j < j0 + step_v0; ++j)
                {
                  j4 = currMB->block_y + j;
                  j6 = currMB->block_y_aff + j;

                  for (i4 = i0; i4 < i0 + step_h0; ++i4)
                  {
                    PicMotionParams *colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                      &list1[0]->mv_info[RSD(j6)][RSD(i4)] : &list1[0]->mv_info[j6][i4];
                    PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];
                    int mv_y;
                    if(currSlice->MbaffFrameFlag)
                    {
                      if(!currMB->mb_field_decoding_flag && ((currSlice->listX[LIST_1][0]->iCodingType==FRAME_MB_PAIR_CODING && currSlice->listX[LIST_1][0]->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
                        (currSlice->listX[LIST_1][0]->iCodingType==FIELD_CODING)))
                      {
                        if (iabs(dec_picture->poc - currSlice->listX[LIST_1+4][0]->poc)> iabs(dec_picture->poc -currSlice->listX[LIST_1+2][0]->poc) )
                        {
                          colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                            &currSlice->listX[LIST_1+2][0]->mv_info[RSD(j6)>>1][RSD(i4)] : &currSlice->listX[LIST_1+2][0]->mv_info[j6>>1][i4];
                        }
                        else
                        {
                          colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                            &currSlice->listX[LIST_1+4][0]->mv_info[RSD(j6)>>1][RSD(i4)] : &currSlice->listX[LIST_1+4][0]->mv_info[j6>>1][i4];
                        }
                      }
                    }
                    else if(!p_Vid->active_sps->frame_mbs_only_flag && !currSlice->field_pic_flag && currSlice->listX[LIST_1][0]->iCodingType!=FRAME_CODING)
                    {
                      if (iabs(dec_picture->poc - list1[0]->bottom_field->poc)> iabs(dec_picture->poc -list1[0]->top_field->poc) )
                      {
                        colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                          &list1[0]->top_field->mv_info[RSD(j6)>>1][RSD(i4)] : &list1[0]->top_field->mv_info[(j6)>>1][i4];
                      }
                      else
                      {
                        colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                          &list1[0]->bottom_field->mv_info[RSD(j6)>>1][RSD(i4)] : &list1[0]->bottom_field->mv_info[(j6)>>1][i4];
                      }
                    }
                    else if(!p_Vid->active_sps->frame_mbs_only_flag && currSlice->field_pic_flag && currSlice->structure!=list1[0]->structure && list1[0]->coded_frame)
                    {
                      if (!currSlice->bottom_field_flag)
                      {
                        colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                          &list1[0]->frame->top_field->mv_info[RSD(j6)][RSD(i4)] : &list1[0]->frame->top_field->mv_info[j6][i4];
                      }
                      else
                      {
                        colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                          &list1[0]->frame->bottom_field->mv_info[RSD(j6)][RSD(i4)] : &list1[0]->frame->bottom_field->mv_info[j6][i4];
                      }
                    }

                    mv_y = colocated->mv[refList].mv_y; 
                    if((currSlice->MbaffFrameFlag && !currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure!=FRAME) ||
                      (!currSlice->MbaffFrameFlag && currSlice->field_pic_flag==0 && colocated->ref_pic[refList]->structure!=FRAME))
                      mv_y *= 2;
                    else if((currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure==FRAME) ||
                      (!currSlice->MbaffFrameFlag && currSlice->field_pic_flag==1 && colocated->ref_pic[refList]->structure==FRAME))
                      mv_y /= 2;

                    mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];

                    mv_info->ref_idx [LIST_0] = (char) mapped_idx;
                    mv_info->ref_idx [LIST_1] = 0;

                    mv_info->ref_pic[LIST_0] = list0[mapped_idx];
                    mv_info->ref_pic[LIST_1] = list1[0];

                    if (mv_scale == 9999 || currSlice->listX[LIST_0+list_offset][mapped_idx]->is_long_term)
                    {
                      mv_info->mv[LIST_0].mv_x = colocated->mv[refList].mv_x;
                      mv_info->mv[LIST_0].mv_y = (short) mv_y;
                      mv_info->mv[LIST_1] = zero_mv;
                    }
                    else
                    {
                      mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * colocated->mv[refList].mv_x + 128 ) >> 8);
                      mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * mv_y/*colocated->mv[refList].mv_y*/ + 128 ) >> 8);
                      mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - colocated->mv[refList].mv_x);
                      mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - mv_y/*colocated->mv[refList].mv_y*/);
                    }
                  }
                }
              }
              else if (INVALIDINDEX == mapped_idx)
              {
                error("temporal direct error: colocated block has ref that is unavailable",-1111);
              }
            }
          }
        }
      }
    }    
  }
}


static inline void update_neighbor_mvs(PicMotionParams **motion, const PicMotionParams *mv_info, int i4)
{
  (*motion++)[i4 + 1] = *mv_info;
  (*motion  )[i4    ] = *mv_info;
  (*motion  )[i4 + 1] = *mv_info;
}

/*!
*************************************************************************************
* \brief
*    Colocated info <= direct_inference is disabled. 
*************************************************************************************
*/
static int get_colocated_info_4x4(Macroblock *currMB, StorablePicture *list1, int i, int j)
{
  if (list1->is_long_term)
    return 1;
  else
  {
    PicMotionParams *fs = &list1->mv_info[j][i];

    int moving = !(((
          (fs->ref_idx[LIST_0] == 0)
      &&  (iabs(fs->mv[LIST_0].mv_x)>>1 == 0)
      &&  (iabs(fs->mv[LIST_0].mv_y)>>1 == 0)))
      || ((fs->ref_idx[LIST_0] == -1)
      &&  (fs->ref_idx[LIST_1] == 0)
      &&  (iabs(fs->mv[LIST_1].mv_x)>>1 == 0)
      &&  (iabs(fs->mv[LIST_1].mv_y)>>1 == 0)));

    return moving;  
  }
}

/*!
*************************************************************************************
* \brief
*    Temporary function for colocated info when direct_inference is enabled. Will be replaced with 
*    function that will access directly motion information
*************************************************************************************
*/
static int get_colocated_info_8x8(Macroblock *currMB, StorablePicture *list1, int i, int j)
{
  if (list1->is_long_term)
    return 1;
  else
  {
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    if( (currSlice->MbaffFrameFlag) ||
      (!p_Vid->active_sps->frame_mbs_only_flag && ((!currSlice->field_pic_flag && list1->iCodingType == FIELD_CODING)||(currSlice->structure!=list1->structure && list1->coded_frame))))
    {
      int jj = RSD(j);
      int ii = RSD(i);
      int jdiv = (jj>>1);
      int moving;
      PicMotionParams *fs = &list1->mv_info[jj][ii];
      
      if(currSlice->field_pic_flag && currSlice->structure!=list1->structure && list1->coded_frame)
      {
         if(!currSlice->bottom_field_flag)
           fs = list1->top_field->mv_info[jj] + ii;
         else
           fs = list1->bottom_field->mv_info[jj] + ii;
      }
      else
      {
        if( (currSlice->MbaffFrameFlag && ((!currMB->mb_field_decoding_flag && list1->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
          (!currMB->mb_field_decoding_flag && list1->iCodingType == FIELD_CODING))) 
          || (!currSlice->MbaffFrameFlag))
        {
          if (iabs(currSlice->dec_picture->poc - list1->bottom_field->poc)> iabs(currSlice->dec_picture->poc -list1->top_field->poc) )
          {
            fs = list1->top_field->mv_info[jdiv] + ii;
          }
          else
          {
            fs = list1->bottom_field->mv_info[jdiv] + ii;
          }
        }
      }
      moving = !((((fs->ref_idx[LIST_0] == 0)
        &&  (iabs(fs->mv[LIST_0].mv_x)>>1 == 0)
        &&  (iabs(fs->mv[LIST_0].mv_y)>>1 == 0)))
        || ((fs->ref_idx[LIST_0] == -1)
        &&  (fs->ref_idx[LIST_1] == 0)
        &&  (iabs(fs->mv[LIST_1].mv_x)>>1 == 0)
        &&  (iabs(fs->mv[LIST_1].mv_y)>>1 == 0)));
      return moving;
    }
    else
    {
      PicMotionParams *fs = &list1->mv_info[RSD(j)][RSD(i)];
      int moving;
      if(currMB->p_Vid->active_sps->separate_colour_plane_flag && currMB->p_Vid->yuv_format==YUV444)
        fs = &list1->JVmv_info[currMB->p_Slice->colour_plane_id][RSD(j)][RSD(i)];
      moving= !((((fs->ref_idx[LIST_0] == 0)
        &&  (iabs(fs->mv[LIST_0].mv_x)>>1 == 0)
        &&  (iabs(fs->mv[LIST_0].mv_y)>>1 == 0)))
        || ((fs->ref_idx[LIST_0] == -1)
        &&  (fs->ref_idx[LIST_1] == 0)
        &&  (iabs(fs->mv[LIST_1].mv_x)>>1 == 0)
        &&  (iabs(fs->mv[LIST_1].mv_y)>>1 == 0)));
      return moving;  
    }
  }
}


static void update_direct_mv_info_spatial_8x8(Macroblock *currMB)
{
  Boolean has_direct = (currMB->b8mode[0] == 0) | (currMB->b8mode[1] == 0) | (currMB->b8mode[2] == 0) | (currMB->b8mode[3] == 0);

  if (has_direct)
  {
    //VideoParameters *p_Vid = currMB->p_Vid;
    Slice *currSlice = currMB->p_Slice;
    int i,j,k;

    int j4, i4;
    StorablePicture *dec_picture = currSlice->dec_picture;

    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    char  l0_rFrame, l1_rFrame;
    MotionVector pmvl0, pmvl1;
    int is_not_moving;
    PicMotionParams *mv_info = NULL;

    prepare_direct_params(currMB, dec_picture, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    for (k = 0; k < 4; ++k)
    {
      if (currMB->b8mode[k] == 0)
      {
        i = 2 * (k & 0x01);
        j = 2 * (k >> 1);

        //j6 = currMB->block_y_aff + j;
        j4 = currMB->block_y     + j;
        i4 = currMB->block_x     + i;

        mv_info = &dec_picture->mv_info[j4][i4];

        is_not_moving = (get_colocated_info_8x8(currMB, list1[0], i4, currMB->block_y_aff + j) == 0);

        if (is_not_moving && (l0_rFrame == 0 || l1_rFrame == 0))
        {            
          if (l1_rFrame == -1)
          {
            if  (l0_rFrame == 0)
            {
              mv_info->ref_pic[LIST_0] = list0[0];
              mv_info->ref_pic[LIST_1] = list1[0];
              mv_info->mv[LIST_0] = zero_mv;
              mv_info->mv[LIST_1] = zero_mv;
              mv_info->ref_idx[LIST_0] = 0;
              mv_info->ref_idx[LIST_1] = -1;
            }
            else
            {
              mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
              mv_info->ref_pic[LIST_1] = NULL;
              mv_info->mv[LIST_0] = pmvl0;
              mv_info->mv[LIST_1] = zero_mv;                    
              mv_info->ref_idx[LIST_0] = l0_rFrame;
              mv_info->ref_idx[LIST_1] = -1;
            }
          }
          else if (l0_rFrame == -1)
          {
            if  (l1_rFrame == 0)
            {
              mv_info->ref_pic[LIST_0] = NULL;
              mv_info->ref_pic[LIST_1] = list1[0];
              mv_info->mv[LIST_0] = zero_mv;
              mv_info->mv[LIST_1] = zero_mv;                    
              mv_info->ref_idx[LIST_0] = -1;
              mv_info->ref_idx[LIST_1] = 0;
            }
            else
            {
              mv_info->ref_pic[LIST_0] = NULL;
              mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
              mv_info->mv[LIST_0] = zero_mv;
              mv_info->mv[LIST_1] = pmvl1;                    
              mv_info->ref_idx[LIST_0] = -1;
              mv_info->ref_idx[LIST_1] = l1_rFrame;
            }
          }
          else
          {
            if  (l0_rFrame == 0)
            {
              mv_info->ref_pic[LIST_0] = list0[0];
              mv_info->mv[LIST_0] = zero_mv;
              mv_info->ref_idx[LIST_0] = 0;
            }
            else
            {
              mv_info->ref_pic[LIST_1] = list1[(short) l0_rFrame];
              mv_info->mv[LIST_0] = pmvl0;
              mv_info->ref_idx[LIST_0] = l0_rFrame;
            }

            if  (l1_rFrame == 0)
            {
              mv_info->ref_pic[LIST_1] = list1[0];
              mv_info->mv[LIST_1] = zero_mv;
              mv_info->ref_idx[LIST_1] = 0;
            }
            else
            {                    
              mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
              mv_info->mv[LIST_1] = pmvl1;
              mv_info->ref_idx[LIST_1] = l1_rFrame;
            }
          }
        }
        else
        {
          if (l0_rFrame < 0 && l1_rFrame < 0)
          {
            mv_info->ref_pic[LIST_0] = list0[0];
            mv_info->ref_pic[LIST_1] = list1[0];
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;
            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = 0;
          }
          else if (l0_rFrame < 0)
          {
            mv_info->ref_pic[LIST_0] = NULL;
            mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = pmvl1;
            mv_info->ref_idx[LIST_0] = -1;
            mv_info->ref_idx[LIST_1] = l1_rFrame;
          }
          else  if (l1_rFrame < 0)
          {
            mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
            mv_info->ref_pic[LIST_1] = NULL;

            mv_info->mv[LIST_0] = pmvl0;
            mv_info->mv[LIST_1] = zero_mv;
            mv_info->ref_idx[LIST_0] = l0_rFrame;
            mv_info->ref_idx[LIST_1] = -1;
          }
          else
          {
            mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
            mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
            mv_info->mv[LIST_0] = pmvl0;
            mv_info->mv[LIST_1] = pmvl1;
            mv_info->ref_idx[LIST_0] = l0_rFrame;
            mv_info->ref_idx[LIST_1] = l1_rFrame;
          }
        }
        update_neighbor_mvs(&dec_picture->mv_info[j4], mv_info, i4);              
      }
    }
  }
}

static void update_direct_mv_info_spatial_4x4(Macroblock *currMB)
{
  Boolean has_direct = (currMB->b8mode[0] == 0) | (currMB->b8mode[1] == 0) | (currMB->b8mode[2] == 0) | (currMB->b8mode[3] == 0);

  if (has_direct)
  {   
    VideoParameters *p_Vid = currMB->p_Vid;
    Slice *currSlice = currMB->p_Slice;
    int i,j,k;

    int j4, i4;
    StorablePicture *dec_picture = p_Vid->dec_picture;

    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    char  l0_rFrame, l1_rFrame;
    MotionVector pmvl0, pmvl1;

    prepare_direct_params(currMB, dec_picture, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);
    for (k = 0; k < 4; ++k)
    {
      if (currMB->b8mode[k] == 0)
      {

        i = 2 * (k & 0x01);
        for(j = 2 * (k >> 1); j < 2 * (k >> 1)+2;++j)
        {
          j4 = currMB->block_y     + j;

          for(i4 = currMB->block_x + i; i4 < currMB->block_x + i + 2; ++i4)
          {
            PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];
            //===== DIRECT PREDICTION =====      
            if (l0_rFrame == 0 || l1_rFrame == 0)
            {
              int is_not_moving = (get_colocated_info_4x4(currMB, list1[0], i4, currMB->block_y_aff + j) == 0);

              if (l1_rFrame == -1)
              {
                if (is_not_moving)
                {
                  mv_info->ref_pic[LIST_0] = list0[0];
                  mv_info->ref_pic[LIST_1] = NULL;
                  mv_info->mv[LIST_0] = zero_mv;
                  mv_info->mv[LIST_1] = zero_mv;
                  mv_info->ref_idx[LIST_0] = 0;
                  mv_info->ref_idx[LIST_1] = -1;
                }
                else
                {
                  mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                  mv_info->ref_pic[LIST_1] = NULL;
                  mv_info->mv[LIST_0] = pmvl0;
                  mv_info->mv[LIST_1] = zero_mv;
                  mv_info->ref_idx[LIST_0] = l0_rFrame;
                  mv_info->ref_idx[LIST_1] = -1;
                }
              }
              else if (l0_rFrame == -1) 
              {
                if  (is_not_moving)
                {
                  mv_info->ref_pic[LIST_0] = NULL;
                  mv_info->ref_pic[LIST_1] = list1[0];
                  mv_info->mv[LIST_0] = zero_mv;
                  mv_info->mv[LIST_1] = zero_mv;
                  mv_info->ref_idx[LIST_0] = -1;
                  mv_info->ref_idx[LIST_1] = 0;
                }
                else
                {
                  mv_info->ref_pic[LIST_0] = NULL;
                  mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                  mv_info->mv[LIST_0] = zero_mv;            
                  mv_info->mv[LIST_1] = pmvl1;
                  mv_info->ref_idx[LIST_0] = -1;
                  mv_info->ref_idx[LIST_1] = l1_rFrame;
                }
              }
              else
              {
                if (l0_rFrame == 0 && ((is_not_moving)))
                {
                  mv_info->ref_pic[LIST_0] = list0[0];
                  mv_info->mv[LIST_0] = zero_mv;
                  mv_info->ref_idx[LIST_0] = 0;
                }
                else
                {
                  mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                  mv_info->mv[LIST_0] = pmvl0;
                  mv_info->ref_idx[LIST_0] = l0_rFrame;
                }

                if  (l1_rFrame == 0 && ((is_not_moving)))
                {
                  mv_info->ref_pic[LIST_1] = list1[0];
                  mv_info->mv[LIST_1] = zero_mv;
                  mv_info->ref_idx[LIST_1]    = 0;
                }
                else
                {
                  mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                  mv_info->mv[LIST_1] = pmvl1;
                  mv_info->ref_idx[LIST_1] = l1_rFrame;              
                }            
              }
            }
            else 
            {       
              mv_info = &dec_picture->mv_info[j4][i4];

              if (l0_rFrame < 0 && l1_rFrame < 0)
              {
                mv_info->ref_pic[LIST_0] = list0[0];
                mv_info->ref_pic[LIST_1] = list1[0];
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = zero_mv;
                mv_info->ref_idx[LIST_0] = 0;
                mv_info->ref_idx[LIST_1] = 0;
              }
              else if (l1_rFrame == -1)
              {
                mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                mv_info->ref_pic[LIST_1] = NULL;
                mv_info->mv[LIST_0] = pmvl0;
                mv_info->mv[LIST_1] = zero_mv;
                mv_info->ref_idx[LIST_0] = l0_rFrame;
                mv_info->ref_idx[LIST_1] = -1;
              }
              else if (l0_rFrame == -1) 
              {
                mv_info->ref_pic[LIST_0] = NULL;
                mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = pmvl1;
                mv_info->ref_idx[LIST_0] = -1;
                mv_info->ref_idx[LIST_1] = l1_rFrame;
              }
              else
              {
                mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                mv_info->mv[LIST_0] = pmvl0;
                mv_info->mv[LIST_1] = pmvl1;
                mv_info->ref_idx[LIST_0] = l0_rFrame;
                mv_info->ref_idx[LIST_1] = l1_rFrame;            
              }
            }
          }
        }
      }
    }        
  }
}


void update_direct_mv_info(Macroblock *currMB)
{
    Slice *currSlice = currMB->p_Slice;
    if (currSlice->active_sps->direct_8x8_inference_flag) {
        if (currSlice->direct_spatial_mv_pred_flag)
            update_direct_mv_info_spatial_8x8(currMB);
        else
            update_direct_mv_info_temporal(currMB);
    } else {
        if (currSlice->direct_spatial_mv_pred_flag)
            update_direct_mv_info_spatial_4x4(currMB);
        else
            update_direct_mv_info_temporal(currMB);
    }
}


int get_direct8x8temporal(Macroblock *currMB, StorablePicture *dec_picture, int block8x8)
{
    short ref_idx;
    int refList;

    int k, i, j, i4, j4, j6;
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    PicMotionParams *mv_info = NULL, *colocated = NULL;
  
    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    int pred_dir = currMB->b8pdir[block8x8];

    int k_start = (block8x8 << 2);

    for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++) {
        i         =  (decode_block_scan[k] & 3);
        j         = ((decode_block_scan[k] >> 2) & 3);
        i4        = currMB->block_x + i;
        j4        = currMB->block_y + j;
        j6        = currMB->block_y_aff + j;
        mv_info   = &dec_picture->mv_info[j4][i4];
        colocated = &list1[0]->mv_info[RSD(j6)][RSD(i4)];

        if (currMB->p_Vid->active_sps->separate_colour_plane_flag && currMB->p_Vid->yuv_format==YUV444)
            colocated = &list1[0]->JVmv_info[currMB->p_Slice->colour_plane_id][RSD(j6)][RSD(i4)];
        if (currSlice->MbaffFrameFlag) {
            assert(p_Vid->active_sps->direct_8x8_inference_flag);
            if (!currMB->mb_field_decoding_flag && ((currSlice->listX[LIST_1][0]->iCodingType==FRAME_MB_PAIR_CODING && currSlice->listX[LIST_1][0]->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
                (currSlice->listX[LIST_1][0]->iCodingType==FIELD_CODING))) {
                if (iabs(dec_picture->poc - currSlice->listX[LIST_1+4][0]->poc)> iabs(dec_picture->poc -currSlice->listX[LIST_1+2][0]->poc) )
                    colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                                &currSlice->listX[LIST_1+2][0]->mv_info[RSD(j6)>>1][RSD(i4)] : &currSlice->listX[LIST_1+2][0]->mv_info[j6>>1][i4];
                else
                    colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                                &currSlice->listX[LIST_1+4][0]->mv_info[RSD(j6)>>1][RSD(i4)] : &currSlice->listX[LIST_1+4][0]->mv_info[j6>>1][i4];
            }
        } else if (!p_Vid->active_sps->frame_mbs_only_flag && 
                   (!currSlice->field_pic_flag && currSlice->listX[LIST_1][0]->iCodingType!=FRAME_CODING)) {
            if (iabs(dec_picture->poc - list1[0]->bottom_field->poc)> iabs(dec_picture->poc -list1[0]->top_field->poc) )
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                            &list1[0]->top_field->mv_info[RSD(j6)>>1][RSD(i4)] : &list1[0]->top_field->mv_info[j6>>1][i4];
            else
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                            &list1[0]->bottom_field->mv_info[RSD(j6)>>1][RSD(i4)] : &list1[0]->bottom_field->mv_info[j6>>1][i4];
        } else if (!p_Vid->active_sps->frame_mbs_only_flag && currSlice->field_pic_flag && currSlice->structure!=list1[0]->structure && list1[0]->coded_frame) {
            if (!currSlice->bottom_field_flag)
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                            &list1[0]->frame->top_field->mv_info[RSD(j6)][RSD(i4)] : &list1[0]->frame->top_field->mv_info[j6][i4];
            else
                colocated = p_Vid->active_sps->direct_8x8_inference_flag ? 
                            &list1[0]->frame->bottom_field->mv_info[RSD(j6)][RSD(i4)] : &list1[0]->frame->bottom_field->mv_info[j6][i4];
        }

        assert (pred_dir<=2);

        refList = (colocated->ref_idx[LIST_0]== -1 ? LIST_1 : LIST_0);
        ref_idx =  colocated->ref_idx[refList];

        if (ref_idx==-1) { // co-located is intra mode
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;

            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = 0;
        } else { // co-located skip or inter mode
            int mapped_idx=0;
            int iref;
            if ((currSlice->MbaffFrameFlag && ( (currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure==FRAME) || 
                (!currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure!=FRAME))) ||
                (!currSlice->MbaffFrameFlag && ((currSlice->field_pic_flag==0 && colocated->ref_pic[refList]->structure!=FRAME)
                ||(currSlice->field_pic_flag==1 && colocated->ref_pic[refList]->structure==FRAME))) ) {
                for (iref = 0; iref < imin(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]);iref++) {
                    if (currSlice->listX[LIST_0 + list_offset][iref]->top_field == colocated->ref_pic[refList] || 
                        currSlice->listX[LIST_0 + list_offset][iref]->bottom_field == colocated->ref_pic[refList] ||
                        currSlice->listX[LIST_0 + list_offset][iref]->frame == colocated->ref_pic[refList]) {
                        if ((currSlice->field_pic_flag==1) && (currSlice->listX[LIST_0 + list_offset][iref]->structure != currSlice->structure))
                            mapped_idx=INVALIDINDEX;
                        else {
                            mapped_idx = iref;            
                            break;
                        }
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx=INVALIDINDEX;
                }
            } else {
                for (iref = 0; iref < imin(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]);iref++) {
                    if (currSlice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                        mapped_idx = iref;            
                        break;
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx=INVALIDINDEX;
                }
            }

            if (INVALIDINDEX != mapped_idx) {
                int mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];
                int mv_y = colocated->mv[refList].mv_y; 
                if ((currSlice->MbaffFrameFlag && !currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure!=FRAME) ||
                   (!currSlice->MbaffFrameFlag && currSlice->field_pic_flag==0 && colocated->ref_pic[refList]->structure!=FRAME) )
                    mv_y *= 2;
                else if ((currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure==FRAME) ||
                        (!currSlice->MbaffFrameFlag && currSlice->field_pic_flag==1 && colocated->ref_pic[refList]->structure==FRAME) )
                    mv_y /= 2;

                //! In such case, an array is needed for each different reference.
                if (mv_scale == 9999 || currSlice->listX[LIST_0 + list_offset][mapped_idx]->is_long_term) {
                    mv_info->mv[LIST_0].mv_x = colocated->mv[refList].mv_x;
                    mv_info->mv[LIST_0].mv_y = (short) mv_y;
                    mv_info->mv[LIST_1] = zero_mv;
                } else {
                    mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * colocated->mv[refList].mv_x + 128 ) >> 8);
                    mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * mv_y/*colocated->mv[refList].mv_y*/ + 128 ) >> 8);

                    mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - colocated->mv[refList].mv_x);
                    mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - mv_y/*colocated->mv[refList].mv_y*/);
                }

                mv_info->ref_idx[LIST_0] = (char) mapped_idx; //colocated->ref_idx[refList];
                mv_info->ref_idx[LIST_1] = 0;
            } else if (INVALIDINDEX == mapped_idx)
                error("temporal direct error: colocated block has ref that is unavailable",-1111);
        }
        // store reference picture ID determined by direct mode
        mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
        mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
    }

    return pred_dir;
}

int get_direct4x4temporal(Macroblock *currMB, StorablePicture *dec_picture, int block8x8)
{
    short ref_idx;
    int refList;

    int k;
    Slice *currSlice = currMB->p_Slice;
  
    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    int pred_dir = currMB->b8pdir[block8x8];

    int k_start = (block8x8 << 2);
    int k_end = k_start + BLOCK_MULTIPLE;

    for (k = k_start; k < k_end; k ++) {
        int i =  (decode_block_scan[k] & 3);
        int j = ((decode_block_scan[k] >> 2) & 3);
        int i4   = currMB->block_x + i;
        int j4   = currMB->block_y + j;
        int j6   = currMB->block_y_aff + j;
        PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];
        PicMotionParams *colocated = &list1[0]->mv_info[j6][i4];
        if (currMB->p_Vid->active_sps->separate_colour_plane_flag && currMB->p_Vid->yuv_format==YUV444)
            colocated = &list1[0]->JVmv_info[currMB->p_Slice->colour_plane_id][RSD(j6)][RSD(i4)];
        assert (pred_dir<=2);

        refList = (colocated->ref_idx[LIST_0]== -1 ? LIST_1 : LIST_0);
        ref_idx =  colocated->ref_idx[refList];

        if (ref_idx==-1) { // co-located is intra mode
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;

            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = 0;
        } else { // co-located skip or inter mode
            int mapped_idx=0;
            int iref;

            for (iref=0;iref<imin(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]);iref++) {
                if (currSlice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                    mapped_idx=iref;
                    break;
                } else //! invalid index. Default to zero even though this case should not happen
                    mapped_idx=INVALIDINDEX;
            }
            if (INVALIDINDEX == mapped_idx)
                error("temporal direct error: colocated block has ref that is unavailable",-1111);
            else {
                int mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];

                //! In such case, an array is needed for each different reference.
                if (mv_scale == 9999 || currSlice->listX[LIST_0+list_offset][mapped_idx]->is_long_term) {
                    mv_info->mv[LIST_0] = colocated->mv[refList];
                    mv_info->mv[LIST_1] = zero_mv;
                } else {
                    mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * colocated->mv[refList].mv_x + 128 ) >> 8);
                    mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * colocated->mv[refList].mv_y + 128 ) >> 8);

                    mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - colocated->mv[refList].mv_x);
                    mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - colocated->mv[refList].mv_y);
                }

                mv_info->ref_idx[LIST_0] = (char) mapped_idx; //colocated->ref_idx[refList];
                mv_info->ref_idx[LIST_1] = 0;
            }
        }
        // store reference picture ID determined by direct mode
        mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
        mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
    }

    return pred_dir;
}

int get_direct8x8spatial_eq(Macroblock *currMB, StorablePicture *dec_picture, int block8x8, MotionVector *pmvl0, MotionVector *pmvl1, char l0_rFrame, char l1_rFrame)
{
    int i4, j4;
    Slice *currSlice = currMB->p_Slice;

    PicMotionParams *mv_info;
    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    int pred_dir = 0;

    int is_not_moving;

    int k_start = (block8x8 << 2);

    int i  =  (decode_block_scan[k_start] & 3);
    int j  = ((decode_block_scan[k_start] >> 2) & 3);
    i4  = currMB->block_x + i;
    j4  = currMB->block_y + j;

    is_not_moving = (get_colocated_info_8x8(currMB, list1[0], i4, currMB->block_y_aff + j) == 0);

    mv_info = &dec_picture->mv_info[j4][i4];

    //===== DIRECT PREDICTION =====
    if (l1_rFrame == -1) {
        if (is_not_moving) {
            mv_info->ref_pic[LIST_0] = list0[0];
            mv_info->ref_pic[LIST_1] = NULL;
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;
            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = -1;
        } else {
            mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
            mv_info->ref_pic[LIST_1] = NULL;
            mv_info->mv[LIST_0] = *pmvl0;
            mv_info->mv[LIST_1] = zero_mv;
            mv_info->ref_idx[LIST_0] = l0_rFrame;
            mv_info->ref_idx[LIST_1] = -1;
        }
        pred_dir = 0;
    } else if (l0_rFrame == -1) {
        if (is_not_moving) {
            mv_info->ref_pic[LIST_0] = NULL;
            mv_info->ref_pic[LIST_1] = list1[0];
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;
            mv_info->ref_idx[LIST_0] = -1;
            mv_info->ref_idx[LIST_1] = 0;
        } else {
            mv_info->ref_pic[LIST_0] = NULL;
            mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
            mv_info->mv[LIST_0] = zero_mv;            
            mv_info->mv[LIST_1] = *pmvl1;
            mv_info->ref_idx[LIST_0] = -1;
            mv_info->ref_idx[LIST_1] = l1_rFrame;
        }
        pred_dir = 1;
    } else {
        if (l0_rFrame == 0 && ((is_not_moving))) {
            mv_info->ref_pic[LIST_0] = list0[0];
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->ref_idx[LIST_0] = 0;
        } else {
            mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
            mv_info->mv[LIST_0] = *pmvl0;
            mv_info->ref_idx[LIST_0] = l0_rFrame;
        }

        if (l1_rFrame == 0 && ((is_not_moving))) {
            mv_info->ref_pic[LIST_1] = list1[0];
            mv_info->mv[LIST_1] = zero_mv;
            mv_info->ref_idx[LIST_1]    = 0;
        } else {
            mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
            mv_info->mv[LIST_1] = *pmvl1;
            mv_info->ref_idx[LIST_1] = l1_rFrame;              
        }
        pred_dir = 2;
    }

    update_neighbor_mvs(&dec_picture->mv_info[j4], mv_info, i4);

    return pred_dir;
}

int get_direct8x8spatial_ne(Macroblock *currMB, StorablePicture *dec_picture, int block8x8, MotionVector *pmvl0, MotionVector *pmvl1, char l0_rFrame, char l1_rFrame)
{
    int i4, j4;
    Slice *currSlice = currMB->p_Slice;

    PicMotionParams *mv_info;
    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    int pred_dir = 0;

    //===== DIRECT PREDICTION =====
    if (l0_rFrame < 0 && l1_rFrame < 0) {
        pred_dir = 2;
        for (j4 = currMB->block_y; j4 < currMB->block_y + BLOCK_MULTIPLE; j4 += 2) {
            for (i4 = currMB->block_x; i4 < currMB->block_x + BLOCK_MULTIPLE; i4 += 2) {
                mv_info = &dec_picture->mv_info[j4][i4];

                mv_info->ref_pic[LIST_0] = list0[0];
                mv_info->ref_pic[LIST_1] = list1[0];
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = zero_mv;
                mv_info->ref_idx[LIST_0] = 0;
                mv_info->ref_idx[LIST_1] = 0;            

                update_neighbor_mvs(&dec_picture->mv_info[j4], mv_info, i4);
            }
        }
    } else if (l1_rFrame == -1) {
        pred_dir = 0;

        for (j4 = currMB->block_y; j4 < currMB->block_y + BLOCK_MULTIPLE; j4 += 2) {
            for (i4 = currMB->block_x; i4 < currMB->block_x + BLOCK_MULTIPLE; i4 += 2) {
                mv_info = &dec_picture->mv_info[j4][i4];

                mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                mv_info->ref_pic[LIST_1] = NULL;
                mv_info->mv[LIST_0] = *pmvl0;
                mv_info->mv[LIST_1] = zero_mv;
                mv_info->ref_idx[LIST_0] = l0_rFrame;
                mv_info->ref_idx[LIST_1] = -1;

                update_neighbor_mvs(&dec_picture->mv_info[j4], mv_info, i4);
            }
        }
    } else if (l0_rFrame == -1) {
        pred_dir = 1;

        for (j4 = currMB->block_y; j4 < currMB->block_y + BLOCK_MULTIPLE; j4 += 2) {
            for (i4 = currMB->block_x; i4 < currMB->block_x + BLOCK_MULTIPLE; i4 += 2) {
                mv_info = &dec_picture->mv_info[j4][i4];

                mv_info->ref_pic[LIST_0] = NULL;
                mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = *pmvl1;
                mv_info->ref_idx[LIST_0] = -1;
                mv_info->ref_idx[LIST_1] = l1_rFrame;

                update_neighbor_mvs(&dec_picture->mv_info[j4], mv_info, i4);
            }
        }
    } else {
        pred_dir = 2;

        for (j4 = currMB->block_y; j4 < currMB->block_y + BLOCK_MULTIPLE; j4 += 2) {
            for (i4 = currMB->block_x; i4 < currMB->block_x + BLOCK_MULTIPLE; i4 += 2) {
                mv_info = &dec_picture->mv_info[j4][i4];

                mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                mv_info->mv[LIST_0] = *pmvl0;
                mv_info->mv[LIST_1] = *pmvl1;
                mv_info->ref_idx[LIST_0] = l0_rFrame;
                mv_info->ref_idx[LIST_1] = l1_rFrame;            

                update_neighbor_mvs(&dec_picture->mv_info[j4], mv_info, i4);
            }
        }
    }

    return pred_dir;
}

int get_direct4x4spatial(Macroblock *currMB, StorablePicture *dec_picture, int block8x8, MotionVector *pmvl0, MotionVector *pmvl1, char l0_rFrame, char l1_rFrame)
{
    int k;
    Slice *currSlice = currMB->p_Slice;

    PicMotionParams *mv_info;
    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    int pred_dir = 0;

    int k_start = (block8x8 << 2);
    int k_end = k_start + BLOCK_MULTIPLE;

    for (k = k_start; k < k_end; k ++) {
        int i  =  (decode_block_scan[k] & 3);
        int j  = ((decode_block_scan[k] >> 2) & 3);
        int i4  = currMB->block_x + i;
        int j4  = currMB->block_y + j;
  
        mv_info = &dec_picture->mv_info[j4][i4];
        //===== DIRECT PREDICTION =====      
        if (l0_rFrame == 0 || l1_rFrame == 0) {
            int is_not_moving = (get_colocated_info_4x4(currMB, list1[0], i4, currMB->block_y_aff + j) == 0);

            if (l1_rFrame == -1) {
                if (is_not_moving) {
                    mv_info->ref_pic[LIST_0] = list0[0];
                    mv_info->ref_pic[LIST_1] = NULL;
                    mv_info->mv[LIST_0] = zero_mv;
                    mv_info->mv[LIST_1] = zero_mv;
                    mv_info->ref_idx[LIST_0] = 0;
                    mv_info->ref_idx[LIST_1] = -1;
                } else {
                    mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                    mv_info->ref_pic[LIST_1] = NULL;
                    mv_info->mv[LIST_0] = *pmvl0;
                    mv_info->mv[LIST_1] = zero_mv;
                    mv_info->ref_idx[LIST_0] = l0_rFrame;
                    mv_info->ref_idx[LIST_1] = -1;
                }
                pred_dir = 0;
            } else if (l0_rFrame == -1) {
                if (is_not_moving) {
                    mv_info->ref_pic[LIST_0] = NULL;
                    mv_info->ref_pic[LIST_1] = list1[0];
                    mv_info->mv[LIST_0] = zero_mv;
                    mv_info->mv[LIST_1] = zero_mv;
                    mv_info->ref_idx[LIST_0] = -1;
                    mv_info->ref_idx[LIST_1] = 0;
                } else {
                    mv_info->ref_pic[LIST_0] = NULL;
                    mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                    mv_info->mv[LIST_0] = zero_mv;            
                    mv_info->mv[LIST_1] = *pmvl1;
                    mv_info->ref_idx[LIST_0] = -1;
                    mv_info->ref_idx[LIST_1] = l1_rFrame;
                }
                pred_dir = 1;
            } else {
                if (l0_rFrame == 0 && ((is_not_moving))) {
                    mv_info->ref_pic[LIST_0] = list0[0];
                    mv_info->mv[LIST_0] = zero_mv;
                    mv_info->ref_idx[LIST_0] = 0;
                } else {
                    mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                    mv_info->mv[LIST_0] = *pmvl0;
                    mv_info->ref_idx[LIST_0] = l0_rFrame;
                }

                if (l1_rFrame == 0 && ((is_not_moving))) {
                    mv_info->ref_pic[LIST_1] = list1[0];
                    mv_info->mv[LIST_1] = zero_mv;
                    mv_info->ref_idx[LIST_1] = 0;
                } else {
                    mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                    mv_info->mv[LIST_1] = *pmvl1;
                    mv_info->ref_idx[LIST_1] = l1_rFrame;              
                }
                pred_dir = 2;
            }
        } else {
            mv_info = &dec_picture->mv_info[j4][i4];

            if (l0_rFrame < 0 && l1_rFrame < 0) {
                pred_dir = 2;
                mv_info->ref_pic[LIST_0] = list0[0];
                mv_info->ref_pic[LIST_1] = list1[0];
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = zero_mv;
                mv_info->ref_idx[LIST_0] = 0;
                mv_info->ref_idx[LIST_1] = 0;
            } else if (l1_rFrame == -1) {
                pred_dir = 0;
                mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                mv_info->ref_pic[LIST_1] = NULL;
                mv_info->mv[LIST_0] = *pmvl0;
                mv_info->mv[LIST_1] = zero_mv;
                mv_info->ref_idx[LIST_0] = l0_rFrame;
                mv_info->ref_idx[LIST_1] = -1;
            } else if (l0_rFrame == -1) {
                pred_dir = 1;
                mv_info->ref_pic[LIST_0] = NULL;
                mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = *pmvl1;
                mv_info->ref_idx[LIST_0] = -1;
                mv_info->ref_idx[LIST_1] = l1_rFrame;
            } else {
                pred_dir = 2;
                mv_info->ref_pic[LIST_0] = list0[(short) l0_rFrame];
                mv_info->ref_pic[LIST_1] = list1[(short) l1_rFrame];
                mv_info->mv[LIST_0] = *pmvl0;
                mv_info->mv[LIST_1] = *pmvl1;
                mv_info->ref_idx[LIST_0] = l0_rFrame;
                mv_info->ref_idx[LIST_1] = l1_rFrame;            
            }
        }
    }

    return pred_dir;
}

int get_inter8x8(Macroblock *currMB, StorablePicture *dec_picture, int block8x8)
{
    int block_size_x, block_size_y;
    int k;
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    int list_offset = currMB->list_offset;
    StorablePicture **list0 = currSlice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = currSlice->listX[LIST_1 + list_offset];

    int mv_mode  = currMB->b8mode[block8x8];
    int pred_dir = currMB->b8pdir[block8x8];
    int k_start, k_end, k_inc;

    if ( mv_mode != 0 ) {
        k_start = (block8x8 << 2);
        k_inc = (mv_mode == SMB8x4) ? 2 : 1;
        k_end = (mv_mode == SMB8x8) ? k_start + 1 : ((mv_mode == SMB4x4) ? k_start + 4 : k_start + k_inc + 1);

        block_size_x = ( mv_mode == SMB8x4 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
        block_size_y = ( mv_mode == SMB4x8 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
    } else {
        k_start = (block8x8 << 2);
        k_end = k_start;
        k_inc = 1;

        if (p_Vid->active_sps->direct_8x8_inference_flag) {
            block_size_x = SMB_BLOCK_SIZE;
            block_size_y = SMB_BLOCK_SIZE;
            k_end ++;
        } else {
            block_size_x = BLOCK_SIZE;
            block_size_y = BLOCK_SIZE;
            k_end += BLOCK_MULTIPLE;
        }

        // Prepare mvs (needed for deblocking and mv prediction
        if (currSlice->direct_spatial_mv_pred_flag) {
            for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++) {
                int i  =  (decode_block_scan[k] & 3);
                int j  = ((decode_block_scan[k] >> 2) & 3);
                int i4 = currMB->block_x + i;
                int j4 = currMB->block_y + j;
                PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];

                assert (pred_dir<=2);

                //===== DIRECT PREDICTION =====
                // motion information should be already set 
                if (mv_info->ref_idx[LIST_1] == -1)
                    pred_dir = 0;
                else if (mv_info->ref_idx[LIST_0] == -1)
                    pred_dir = 1;
                else
                    pred_dir = 2;
            }
        } else {
            for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++) {
                int i  =  (decode_block_scan[k] & 3);
                int j  = ((decode_block_scan[k] >> 2) & 3);
                int i4 = currMB->block_x + i;
                int j4 = currMB->block_y + j;
                PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];

                assert (pred_dir<=2);

                // store reference picture ID determined by direct mode
                mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
                mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
            }
        }
    }

    return pred_dir;
}


static inline void set_direct_references(const PixelPos *mb, char *l0_rFrame, char *l1_rFrame, PicMotionParams **mv_info)
{
  if (mb->available)
  {
    char *ref_idx = mv_info[mb->pos_y][mb->pos_x].ref_idx;
    *l0_rFrame  = ref_idx[LIST_0];
    *l1_rFrame  = ref_idx[LIST_1];
  }
  else
  {
    *l0_rFrame  = -1;
    *l1_rFrame  = -1;
  }
}


static void set_direct_references_mb_field(const PixelPos *mb, char *l0_rFrame, char *l1_rFrame, PicMotionParams **mv_info, Macroblock *mb_data)
{
  if (mb->available)
  {
    char *ref_idx = mv_info[mb->pos_y][mb->pos_x].ref_idx;
    if (mb_data[mb->mb_addr].mb_field_decoding_flag)
    {
      *l0_rFrame  = ref_idx[LIST_0];
      *l1_rFrame  = ref_idx[LIST_1];
    }
    else
    {
      *l0_rFrame  = (ref_idx[LIST_0] < 0) ? ref_idx[LIST_0] : ref_idx[LIST_0] * 2;
      *l1_rFrame  = (ref_idx[LIST_1] < 0) ? ref_idx[LIST_1] : ref_idx[LIST_1] * 2;
    }
  }
  else
  {
    *l0_rFrame  = -1;
    *l1_rFrame  = -1;
  }
}

static void set_direct_references_mb_frame(const PixelPos *mb, char *l0_rFrame, char *l1_rFrame, PicMotionParams **mv_info, Macroblock *mb_data)
{
  if (mb->available)
  {
    char *ref_idx = mv_info[mb->pos_y][mb->pos_x].ref_idx;
    if (mb_data[mb->mb_addr].mb_field_decoding_flag)
    {
      *l0_rFrame  = (ref_idx[LIST_0] >> 1);
      *l1_rFrame  = (ref_idx[LIST_1] >> 1);
    }
    else
    {
      *l0_rFrame  = ref_idx[LIST_0];
      *l1_rFrame  = ref_idx[LIST_1];
    }
  }
  else
  {
    *l0_rFrame  = -1;
    *l1_rFrame  = -1;
  }
}

void prepare_direct_params(Macroblock *currMB, StorablePicture *dec_picture, MotionVector *pmvl0, MotionVector *pmvl1, char *l0_rFrame, char *l1_rFrame)
{
  Slice *currSlice = currMB->p_Slice;
  char l0_refA, l0_refB, l0_refC;
  char l1_refA, l1_refB, l1_refC;
  PicMotionParams **mv_info = dec_picture->mv_info;
  
  PixelPos mb[4];

  get_neighbors(currMB, mb, 0, 0, 16);

  if (!currSlice->MbaffFrameFlag)
  {
    set_direct_references(&mb[0], &l0_refA, &l1_refA, mv_info);
    set_direct_references(&mb[1], &l0_refB, &l1_refB, mv_info);
    set_direct_references(&mb[2], &l0_refC, &l1_refC, mv_info);
  }
  else
  {
    VideoParameters *p_Vid = currMB->p_Vid;
    if (currMB->mb_field_decoding_flag)
    {
      set_direct_references_mb_field(&mb[0], &l0_refA, &l1_refA, mv_info, p_Vid->mb_data);
      set_direct_references_mb_field(&mb[1], &l0_refB, &l1_refB, mv_info, p_Vid->mb_data);
      set_direct_references_mb_field(&mb[2], &l0_refC, &l1_refC, mv_info, p_Vid->mb_data);
    }
    else
    {
      set_direct_references_mb_frame(&mb[0], &l0_refA, &l1_refA, mv_info, p_Vid->mb_data);
      set_direct_references_mb_frame(&mb[1], &l0_refB, &l1_refB, mv_info, p_Vid->mb_data);
      set_direct_references_mb_frame(&mb[2], &l0_refC, &l1_refC, mv_info, p_Vid->mb_data);
    }
  }

  *l0_rFrame = (char) imin(imin((unsigned char) l0_refA, (unsigned char) l0_refB), (unsigned char) l0_refC);
  *l1_rFrame = (char) imin(imin((unsigned char) l1_refA, (unsigned char) l1_refB), (unsigned char) l1_refC);

  if (*l0_rFrame >=0)
    GetMVPredictor (currMB, mb, pmvl0, *l0_rFrame, mv_info, LIST_0, 0, 0, 16, 16);

  if (*l1_rFrame >=0)
    GetMVPredictor (currMB, mb, pmvl1, *l1_rFrame, mv_info, LIST_1, 0, 0, 16, 16);
}
