#include "global.h"
#include "memalloc.h"

using vio::h264::mb_t;
using namespace vio::h264;

#include "erc_api.h"
#include "slice.h"
#include "macroblock.h"

/*!
 ************************************************************************
 * \brief
 *    Initinize the error concealment module
 ************************************************************************
 */
void ercInit(VideoParameters *p_Vid, int pic_sizex, int pic_sizey, int flag)
{
  ercClose(p_Vid, p_Vid->erc_errorVar);
  p_Vid->erc_object_list = (objectBuffer_t *) calloc((pic_sizex * pic_sizey) >> 6, sizeof(objectBuffer_t));
  if (p_Vid->erc_object_list == NULL) no_mem_exit("ercInit: erc_object_list");

  // the error concealment instance is allocated
  p_Vid->erc_errorVar = ercOpen();

  // set error concealment ON
  ercSetErrorConcealment(p_Vid->erc_errorVar, flag);
}

/*!
 ************************************************************************
 * \brief
 *      Allocates data structures used in error concealment.
 *\return
 *      The allocated ercVariables_t is returned.
 ************************************************************************
 */
ercVariables_t *ercOpen( void )
{
  ercVariables_t *errorVar = NULL;

  errorVar = (ercVariables_t *)malloc( sizeof(ercVariables_t));
  if ( errorVar == NULL ) no_mem_exit("ercOpen: errorVar");

  errorVar->nOfMBs = 0;
  errorVar->segments = NULL;
  errorVar->currSegment = 0;
  errorVar->yCondition = NULL;
  errorVar->uCondition = NULL;
  errorVar->vCondition = NULL;
  errorVar->prevFrameYCondition = NULL;

  errorVar->concealment = 1;

  return errorVar;
}

/*!
 ************************************************************************
 * \brief
 *      Resets the variables used in error detection.
 *      Should be called always when starting to decode a new frame.
 * \param errorVar
 *      Variables for error concealment
 * \param nOfMBs
 *      Number of macroblocks in a frame
 * \param numOfSegments
 *    Estimated number of segments (memory reserved)
 * \param picSizeX
 *      Width of the frame in pixels.
 ************************************************************************
 */
void ercReset( ercVariables_t *errorVar, int nOfMBs, int numOfSegments, int picSizeX )
{
  char *tmp = NULL;
  int i = 0;

  if ( errorVar && errorVar->concealment )
  {
    ercSegment_t *segments = NULL;
    // If frame size has been changed
    if ( nOfMBs != errorVar->nOfMBs && errorVar->yCondition != NULL )
    {
      free( errorVar->yCondition );
      errorVar->yCondition = NULL;
      free( errorVar->prevFrameYCondition );
      errorVar->prevFrameYCondition = NULL;
      free( errorVar->uCondition );
      errorVar->uCondition = NULL;
      free( errorVar->vCondition );
      errorVar->vCondition = NULL;
      free( errorVar->segments );
      errorVar->segments = NULL;
    }

    // If the structures are uninitialized (first frame, or frame size is changed)
    if ( errorVar->yCondition == NULL )
    {
      errorVar->segments = (ercSegment_t *)malloc( numOfSegments*sizeof(ercSegment_t) );
      if ( errorVar->segments == NULL ) no_mem_exit("ercReset: errorVar->segments");
      memset( errorVar->segments, 0, numOfSegments*sizeof(ercSegment_t));
      errorVar->nOfSegments = numOfSegments;

      errorVar->yCondition = (char *)malloc( 4*nOfMBs*sizeof(char) );
      if ( errorVar->yCondition == NULL ) no_mem_exit("ercReset: errorVar->yCondition");
      errorVar->prevFrameYCondition = (char *)malloc( 4*nOfMBs*sizeof(char) );
      if ( errorVar->prevFrameYCondition == NULL ) no_mem_exit("ercReset: errorVar->prevFrameYCondition");
      errorVar->uCondition = (char *)malloc( nOfMBs*sizeof(char) );
      if ( errorVar->uCondition == NULL ) no_mem_exit("ercReset: errorVar->uCondition");
      errorVar->vCondition = (char *)malloc( nOfMBs*sizeof(char) );
      if ( errorVar->vCondition == NULL ) no_mem_exit("ercReset: errorVar->vCondition");
      errorVar->nOfMBs = nOfMBs;
    }
    else
    {
      // Store the yCondition struct of the previous frame
      tmp = errorVar->prevFrameYCondition;
      errorVar->prevFrameYCondition = errorVar->yCondition;
      errorVar->yCondition = tmp;
    }

    // Reset tables and parameters
    memset( errorVar->yCondition, 0, 4*nOfMBs*sizeof(*errorVar->yCondition));
    memset( errorVar->uCondition, 0,   nOfMBs*sizeof(*errorVar->uCondition));
    memset( errorVar->vCondition, 0,   nOfMBs*sizeof(*errorVar->vCondition));

    if (errorVar->nOfSegments != numOfSegments)
    {
      free( errorVar->segments );
      errorVar->segments = NULL;
      errorVar->segments = (ercSegment_t *)malloc( numOfSegments*sizeof(ercSegment_t) );
      if ( errorVar->segments == NULL ) no_mem_exit("ercReset: errorVar->segments");
      errorVar->nOfSegments = numOfSegments;
    }

    //memset( errorVar->segments, 0, errorVar->nOfSegments * sizeof(ercSegment_t));

    segments = errorVar->segments;
    for ( i = 0; i < errorVar->nOfSegments; i++ )
    {
      segments->startMBPos = 0;
      segments->endMBPos = (short) (nOfMBs - 1);
      (segments++)->fCorrupted = 1; //! mark segments as corrupted
    }

    errorVar->currSegment = 0;
    errorVar->nOfCorruptedSegments = 0;
  }
}

/*!
 ************************************************************************
 * \brief
 *      Resets the variables used in error detection.
 *      Should be called always when starting to decode a new frame.
 * \param p_Vid
 *      VideoParameters variable
 * \param errorVar
 *      Variables for error concealment
 ************************************************************************
 */
void ercClose(VideoParameters *p_Vid,  ercVariables_t *errorVar )
{
  if ( errorVar != NULL )
  {
    if (errorVar->yCondition != NULL)
    {
      free( errorVar->segments );
      free( errorVar->yCondition );
      free( errorVar->uCondition );
      free( errorVar->vCondition );
      free( errorVar->prevFrameYCondition );
    }
    free( errorVar );
    errorVar = NULL;
  }

  if (p_Vid->erc_object_list)
  {
    free(p_Vid->erc_object_list);
    p_Vid->erc_object_list=NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *      Sets error concealment ON/OFF. Can be invoked only between frames, not during a frame
 * \param errorVar
 *      Variables for error concealment
 * \param value
 *      New value
 ************************************************************************
 */
void ercSetErrorConcealment( ercVariables_t *errorVar, int value )
{
  if ( errorVar != NULL )
    errorVar->concealment = value;
}

/*!
 ************************************************************************
 * \brief
 *      Creates a new segment in the segment-list, and marks the start MB and bit position.
 *      If the end of the previous segment was not explicitly marked by "ercStopSegment",
 *      also marks the end of the previous segment.
 *      If needed, it reallocates the segment-list for a larger storage place.
 * \param currMBNum
 *      The MB number where the new slice/segment starts
 * \param segment
 *      Segment/slice_t No. counted by the caller
 * \param bitPos
 *      data_partition_t pointer: number of bits read from the buffer.
 * \param errorVar
 *      Variables for error detector
 ************************************************************************
 */
void ercStartSegment( int currMBNum, int segment, unsigned int bitPos, ercVariables_t *errorVar )
{
  if ( errorVar && errorVar->concealment )
  {
    errorVar->currSegmentCorrupted = 0;

    errorVar->segments[ segment ].fCorrupted = 0;
    errorVar->segments[ segment ].startMBPos = (short) currMBNum;

  }
}

/*!
 ************************************************************************
 * \brief
 *      Marks the end position of a segment.
 * \param currMBNum
 *      The last MB number of the previous segment
 * \param segment
 *      Segment/slice_t No. counted by the caller
 *      If (segment<0) the internal segment counter is used.
 * \param bitPos
 *      data_partition_t pointer: number of bits read from the buffer.
 * \param errorVar
 *      Variables for error detector
 ************************************************************************
 */
void ercStopSegment( int currMBNum, int segment, unsigned int bitPos, ercVariables_t *errorVar )
{
  if ( errorVar && errorVar->concealment )
  {
    errorVar->segments[ segment ].endMBPos = (short) currMBNum;
    errorVar->currSegment++;
  }
}

/*!
 ************************************************************************
 * \brief
 *      Marks the current segment (the one which has the "currMBNum" MB in it)
 *      as lost: all the blocks of the MBs in the segment as corrupted.
 * \param picSizeX
 *      Width of the frame in pixels.
 * \param errorVar
 *      Variables for error detector
 ************************************************************************
 */
void ercMarkCurrSegmentLost(int picSizeX, ercVariables_t *errorVar )
{
  int j = 0;
  int current_segment;

  current_segment = errorVar->currSegment-1;
  if ( errorVar && errorVar->concealment )
  {
    if (errorVar->currSegmentCorrupted == 0)
    {
      errorVar->nOfCorruptedSegments++;
      errorVar->currSegmentCorrupted = 1;
    }

    for ( j = errorVar->segments[current_segment].startMBPos; j <= errorVar->segments[current_segment].endMBPos; j++ )
    {
      errorVar->yCondition[MBNum2YBlock (j, 0, picSizeX)] = ERC_BLOCK_CORRUPTED;
      errorVar->yCondition[MBNum2YBlock (j, 1, picSizeX)] = ERC_BLOCK_CORRUPTED;
      errorVar->yCondition[MBNum2YBlock (j, 2, picSizeX)] = ERC_BLOCK_CORRUPTED;
      errorVar->yCondition[MBNum2YBlock (j, 3, picSizeX)] = ERC_BLOCK_CORRUPTED;
      errorVar->uCondition[j] = ERC_BLOCK_CORRUPTED;
      errorVar->vCondition[j] = ERC_BLOCK_CORRUPTED;
    }
    errorVar->segments[current_segment].fCorrupted = 1;
  }
}

/*!
 ************************************************************************
 * \brief
 *      Marks the current segment (the one which has the "currMBNum" MB in it)
 *      as OK: all the blocks of the MBs in the segment as OK.
 * \param picSizeX
 *      Width of the frame in pixels.
 * \param errorVar
 *      Variables for error detector
 ************************************************************************
 */
void ercMarkCurrSegmentOK(int picSizeX, ercVariables_t *errorVar )
{
  int j = 0;
  int current_segment;

  current_segment = errorVar->currSegment-1;
  if ( errorVar && errorVar->concealment )
  {
    // mark all the Blocks belonging to the segment as OK */
    for ( j = errorVar->segments[current_segment].startMBPos; j <= errorVar->segments[current_segment].endMBPos; j++ )
    {
      errorVar->yCondition[MBNum2YBlock (j, 0, picSizeX)] = ERC_BLOCK_OK;
      errorVar->yCondition[MBNum2YBlock (j, 1, picSizeX)] = ERC_BLOCK_OK;
      errorVar->yCondition[MBNum2YBlock (j, 2, picSizeX)] = ERC_BLOCK_OK;
      errorVar->yCondition[MBNum2YBlock (j, 3, picSizeX)] = ERC_BLOCK_OK;
      errorVar->uCondition[j] = ERC_BLOCK_OK;
      errorVar->vCondition[j] = ERC_BLOCK_OK;
    }
    errorVar->segments[current_segment].fCorrupted = 0;
  }
}

/*!
 ************************************************************************
 * \brief
 *      Marks the Blocks of the given component (YUV) of the current MB as concealed.
 * \param currMBNum
 *      Selects the segment where this MB number is in.
 * \param comp
 *      Component to mark (0:Y, 1:U, 2:V, <0:All)
 * \param picSizeX
 *      Width of the frame in pixels.
 * \param errorVar
 *      Variables for error detector
 ************************************************************************
 */
void ercMarkCurrMBConcealed( int currMBNum, int comp, int picSizeX, ercVariables_t *errorVar )
{
  int setAll = 0;

  if ( errorVar && errorVar->concealment )
  {
    if (comp < 0)
    {
      setAll = 1;
      comp = 0;
    }

    switch (comp)
    {
    case 0:
      errorVar->yCondition[MBNum2YBlock (currMBNum, 0, picSizeX)] = ERC_BLOCK_CONCEALED;
      errorVar->yCondition[MBNum2YBlock (currMBNum, 1, picSizeX)] = ERC_BLOCK_CONCEALED;
      errorVar->yCondition[MBNum2YBlock (currMBNum, 2, picSizeX)] = ERC_BLOCK_CONCEALED;
      errorVar->yCondition[MBNum2YBlock (currMBNum, 3, picSizeX)] = ERC_BLOCK_CONCEALED;
      if (!setAll)
        break;
    case 1:
      errorVar->uCondition[currMBNum] = ERC_BLOCK_CONCEALED;
      if (!setAll)
        break;
    case 2:
      errorVar->vCondition[currMBNum] = ERC_BLOCK_CONCEALED;
    }
  }
}

void erc_picture(VideoParameters *p_Vid, storable_picture **dec_picture)
{
    frame recfr;
    recfr.p_Vid = p_Vid;
    recfr.yptr = &(*dec_picture)->imgY[0][0];
    if ((*dec_picture)->chroma_format_idc != YUV400) {
        recfr.uptr = &(*dec_picture)->imgUV[0][0][0];
        recfr.vptr = &(*dec_picture)->imgUV[1][0][0];
    }

    //! this is always true at the beginning of a picture
    int ercStartMB = 0;
    int ercSegment = 0;

    //! mark the start of the first segment
    if (!(*dec_picture)->mb_aff_frame_flag) {
        int i;
        ercStartSegment(0, ercSegment, 0 , p_Vid->erc_errorVar);
        //! generate the segments according to the macroblock map
        for (i = 1; i < (int) (*dec_picture)->PicSizeInMbs; ++i) {
            if (p_Vid->mb_data[i].ei_flag != p_Vid->mb_data[i-1].ei_flag) {
                ercStopSegment(i-1, ercSegment, 0, p_Vid->erc_errorVar); //! stop current segment

                //! mark current segment as lost or OK
                if (p_Vid->mb_data[i-1].ei_flag)
                    ercMarkCurrSegmentLost((*dec_picture)->size_x, p_Vid->erc_errorVar);
                else
                    ercMarkCurrSegmentOK((*dec_picture)->size_x, p_Vid->erc_errorVar);

                ++ercSegment;  //! next segment
                ercStartSegment(i, ercSegment, 0 , p_Vid->erc_errorVar); //! start new segment
                ercStartMB = i;//! save start MB for this segment
            }
        }
        //! mark end of the last segment
        ercStopSegment((*dec_picture)->PicSizeInMbs-1, ercSegment, 0, p_Vid->erc_errorVar);
        if (p_Vid->mb_data[i-1].ei_flag)
            ercMarkCurrSegmentLost((*dec_picture)->size_x, p_Vid->erc_errorVar);
        else
            ercMarkCurrSegmentOK((*dec_picture)->size_x, p_Vid->erc_errorVar);

        //! call the right error concealment function depending on the frame type.
        p_Vid->erc_mvperMB /= (*dec_picture)->PicSizeInMbs;

        p_Vid->erc_img = p_Vid;

        if ((*dec_picture)->slice_type == I_slice || (*dec_picture)->slice_type == SI_slice) // I-frame
            ercConcealIntraFrame(p_Vid, &recfr, (*dec_picture)->size_x, (*dec_picture)->size_y, p_Vid->erc_errorVar);
        else
            ercConcealInterFrame(&recfr, p_Vid->erc_object_list, (*dec_picture)->size_x, (*dec_picture)->size_y, p_Vid->erc_errorVar, (*dec_picture)->chroma_format_idc);
    }
}

/*!
 ************************************************************************
 * \brief
 *    write the encoding mode and motion vectors of current
 *    MB to the buffer of the error concealment module.
 ************************************************************************
 */
void ercWriteMBMODEandMV(mb_t *currMB)
{
  VideoParameters *p_Vid = currMB->p_Slice->p_Vid;
  int i, ii, jj, currMBNum = currMB->mbAddrX;
  storable_picture *dec_picture = p_Vid->dec_picture;
  int mbx = xPosMB(currMBNum, dec_picture->size_x), mby = yPosMB(currMBNum, dec_picture->size_x);
  objectBuffer_t *currRegion, *pRegion;

  currRegion = p_Vid->erc_object_list + (currMBNum<<2);

  if(p_Vid->type != B_slice) //non-B frame
  {
    for (i=0; i<4; ++i)
    {
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type == I_16x16 ? REGMODE_INTRA :
        currMB->SubMbType[i]==I_4x4 ? REGMODE_INTRA_8x8  :
        currMB->SubMbType[i]==0    ? REGMODE_INTER_COPY :
        currMB->SubMbType[i]==1    ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);
      if (currMB->SubMbType[i]==0 || currMB->SubMbType[i]==I_4x4)  // INTRA OR COPY
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        ii              = 4*mbx + (i & 0x01)*2;
        jj              = 4*mby + (i >> 1  )*2;
        if (currMB->SubMbType[i]>=5 && currMB->SubMbType[i]<=7)  // SMALL BLOCKS
        {
          pRegion->mv[0]  = (dec_picture->mv_info[jj][ii].mv[LIST_0].mv_x + dec_picture->mv_info[jj][ii + 1].mv[LIST_0].mv_x + dec_picture->mv_info[jj + 1][ii].mv[LIST_0].mv_x + dec_picture->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_x + 2)/4;
          pRegion->mv[1]  = (dec_picture->mv_info[jj][ii].mv[LIST_0].mv_y + dec_picture->mv_info[jj][ii + 1].mv[LIST_0].mv_y + dec_picture->mv_info[jj + 1][ii].mv[LIST_0].mv_y + dec_picture->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_y + 2)/4;
        }
        else // 16x16, 16x8, 8x16, 8x8
        {
          pRegion->mv[0]  = dec_picture->mv_info[jj][ii].mv[LIST_0].mv_x;
          pRegion->mv[1]  = dec_picture->mv_info[jj][ii].mv[LIST_0].mv_y;
        }
        currMB->p_Slice->erc_mvperMB      += abs(pRegion->mv[0]) + abs(pRegion->mv[1]);
        pRegion->mv[2]    = dec_picture->mv_info[jj][ii].ref_idx[LIST_0];
      }
    }
  }
  else  //B-frame
  {
    for (i=0; i<4; ++i)
    {
      ii                  = 4*mbx + (i%2)*2;
      jj                  = 4*mby + (i/2)*2;
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I_16x16  ? REGMODE_INTRA      :
        currMB->SubMbType[i]==I_4x4 ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
      if (currMB->mb_type==I_16x16 || currMB->SubMbType[i]==I_4x4)  // INTRA
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        int idx = (dec_picture->mv_info[jj][ii].ref_idx[0] < 0) ? 1 : 0;
        pRegion->mv[0]    = (dec_picture->mv_info[jj][ii].mv[idx].mv_x + 
          dec_picture->mv_info[jj][ii+1].mv[idx].mv_x + 
          dec_picture->mv_info[jj+1][ii].mv[idx].mv_x + 
          dec_picture->mv_info[jj+1][ii+1].mv[idx].mv_x + 2)/4;
        pRegion->mv[1]    = (dec_picture->mv_info[jj][ii].mv[idx].mv_y + 
          dec_picture->mv_info[jj][ii+1].mv[idx].mv_y + 
          dec_picture->mv_info[jj+1][ii].mv[idx].mv_y + 
          dec_picture->mv_info[jj+1][ii+1].mv[idx].mv_y + 2)/4;
        currMB->p_Slice->erc_mvperMB      += abs(pRegion->mv[0]) + abs(pRegion->mv[1]);

        pRegion->mv[2]  = (dec_picture->mv_info[jj][ii].ref_idx[idx]);
      }
    }
  }
}
