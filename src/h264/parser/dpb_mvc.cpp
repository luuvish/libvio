#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "dpb.h"
#include "memalloc.h"
#include "output.h"

#if (MVC_EXTENSION_ENABLE)

static int GetViewIdx(VideoParameters *p_Vid, int iVOIdx)
{
  int iViewIdx = -1;
  int *piViewIdMap;

  if( p_Vid->active_subset_sps )
  {
    assert( p_Vid->active_subset_sps->num_views_minus1 >= iVOIdx && iVOIdx >= 0 );
    piViewIdMap = p_Vid->active_subset_sps->view_id;
    iViewIdx = piViewIdMap[iVOIdx];    
  }

  return iViewIdx;
}

static int get_maxViewIdx (VideoParameters *p_Vid, int view_id, int anchor_pic_flag, int listidx)
{
  int VOIdx;
  int maxViewIdx = 0;

  VOIdx = view_id; 
  if(VOIdx >= 0)
  {
    if(anchor_pic_flag)
      maxViewIdx = listidx? p_Vid->active_subset_sps->num_anchor_refs_l1[VOIdx] : p_Vid->active_subset_sps->num_anchor_refs_l0[VOIdx];
    else
      maxViewIdx = listidx? p_Vid->active_subset_sps->num_non_anchor_refs_l1[VOIdx] : p_Vid->active_subset_sps->num_non_anchor_refs_l0[VOIdx];
  }

  return maxViewIdx;
}

/*!
 ************************************************************************
 * \brief
 *    Returns inter-view prediction pic with given targetViewID
 *
 ************************************************************************
 */
static storable_picture*  get_inter_view_pic(VideoParameters *p_Vid, slice_t *currSlice, int targetViewID, int currPOC, int listidx)
{
  unsigned i;
  unsigned int listinterview_size;
  frame_store **fs_listinterview;

  if (listidx == 0)
  {
    fs_listinterview = currSlice->fs_listinterview0;
    listinterview_size = currSlice->listinterviewidx0; 
  }
  else
  {
    fs_listinterview = currSlice->fs_listinterview1;
    listinterview_size = currSlice->listinterviewidx1; 
  }

  for(i=0; i<listinterview_size; i++)
  {
    if (fs_listinterview[i]->layer_id == GetVOIdx( p_Vid, targetViewID ))
    {
      if(!currSlice->field_pic_flag && fs_listinterview[i]->frame->poc == currPOC)
      {
        return fs_listinterview[i]->frame;
      }
      else if(currSlice->field_pic_flag && !currSlice->bottom_field_flag && fs_listinterview[i]->top_field->poc == currPOC)
      {
        return fs_listinterview[i]->top_field;
      }
      else if(currSlice->field_pic_flag && currSlice->bottom_field_flag && fs_listinterview[i]->bottom_field->poc == currPOC)
      {
        return fs_listinterview[i]->bottom_field;
      }
    }
  }

  return NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Generates a alternating field list from a given frame_store inter-view list
 *
 ************************************************************************
 */
static void gen_pic_list_from_frame_interview_list(bool bottom_field_flag, frame_store **fs_list, int list_idx, storable_picture **list, char *list_size)
{
  int i;

  if (!bottom_field_flag)
  {
    for (i=0; i<list_idx; i++)
    {
      list[(int)(*list_size)] = fs_list[i]->top_field;
      (*list_size)++;
    }
  }
  else
  {
    for (i=0; i<list_idx; i++)
    {
      list[(int)(*list_size)] = fs_list[i]->bottom_field;
      (*list_size)++;
    }
  }
}


/*!
************************************************************************
* \brief
*    Initialize reference lists depending on current slice type
*
************************************************************************
*/
void init_lists_i_slice_mvc(slice_t *currSlice)
{
  //VideoParameters *p_Vid = currSlice->p_Vid;

  currSlice->listinterviewidx0 = 0;
  currSlice->listinterviewidx1 = 0;

  currSlice->listXsize[0] = 0;
  currSlice->listXsize[1] = 0;
}

/*!
************************************************************************
* \brief
*    Initialize reference lists for a P slice_t
*
************************************************************************
*/
void init_lists_p_slice_mvc(slice_t *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  dpb_t *p_Dpb = currSlice->p_Dpb;

  unsigned int i;

  int list0idx = 0;
  int listltidx = 0;

  frame_store **fs_list0;
  frame_store **fs_listlt;

  int currPOC = currSlice->ThisPOC;
  int anchor_pic_flag = currSlice->anchor_pic_flag;

  currSlice->listinterviewidx0 = 0;
  currSlice->listinterviewidx1 = 0;

  if (!currSlice->field_pic_flag)
  {
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_used==3)
      {
        if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
        {
          currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
        }
      }
    }
    // order list 0 by PicNum
    qsort((void *)currSlice->listX[0], list0idx, sizeof(storable_picture*), compare_pic_by_pic_num_desc);
    currSlice->listXsize[0] = (char) list0idx;
    //printf("listX[0] (PicNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", currSlice->listX[0][i]->pic_num);} printf("\n");

    // long term handling
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ltref[i]->is_used==3)
      {
        if (p_Dpb->fs_ltref[i]->frame->is_long_term)
        {
          currSlice->listX[0][list0idx++]=p_Dpb->fs_ltref[i]->frame;
        }
      }
    }
    qsort((void *)&currSlice->listX[0][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(storable_picture*), compare_pic_by_lt_pic_num_asc);
    currSlice->listXsize[0] = (char) list0idx;
  }
  else
  {
    fs_list0 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL==fs_list0)
      no_mem_exit("init_lists: fs_list0");
    fs_listlt = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL==fs_listlt)
      no_mem_exit("init_lists: fs_listlt");

    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_reference)
      {
        fs_list0[list0idx++] = p_Dpb->fs_ref[i];
      }
    }

    qsort((void *)fs_list0, list0idx, sizeof(frame_store*), compare_fs_by_frame_num_desc);

    //printf("fs_list0 (FrameNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->FrameNumWrap);} printf("\n");

    currSlice->listXsize[0] = 0;
    gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);

    //printf("listX[0] (PicNum): "); for (i=0; i < currSlice->listXsize[0]; i++){printf ("%d  ", currSlice->listX[0][i]->pic_num);} printf("\n");

    // long term handling
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
    {
      fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
    }

    qsort((void *)fs_listlt, listltidx, sizeof(frame_store*), compare_fs_by_lt_pic_idx_asc);

    gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);

    free(fs_list0);
    free(fs_listlt);
  }

  currSlice->listXsize[1] = 0;    

#if !SIMULCAST_ENABLE
  if (currSlice->svc_extension_flag == 0)
  {        
    int curr_view_id = currSlice->layer_id;
    currSlice->fs_listinterview0 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL==currSlice->fs_listinterview0)
      no_mem_exit("init_lists: fs_listinterview0");
    list0idx = currSlice->listXsize[0];
    if (!currSlice->field_pic_flag)
    {
      append_interview_list(p_Vid->p_Dpb_layer[1], false, false, 0, currSlice->fs_listinterview0, &currSlice->listinterviewidx0, currPOC, curr_view_id, anchor_pic_flag);
      for (i=0; i<(unsigned int)currSlice->listinterviewidx0; i++)
      {
        currSlice->listX[0][list0idx++]=currSlice->fs_listinterview0[i]->frame;
      }
      currSlice->listXsize[0] = (char) list0idx;
    }
    else
    {
      append_interview_list(p_Vid->p_Dpb_layer[1], currSlice->field_pic_flag, currSlice->bottom_field_flag, 0, currSlice->fs_listinterview0, &currSlice->listinterviewidx0, currPOC, curr_view_id, anchor_pic_flag);
      gen_pic_list_from_frame_interview_list(currSlice->bottom_field_flag, currSlice->fs_listinterview0, currSlice->listinterviewidx0, currSlice->listX[0], &currSlice->listXsize[0]);
    }
  }
#endif
  // set max size
  currSlice->listXsize[0] = (char) min<int>(currSlice->listXsize[0], currSlice->num_ref_idx_l0_active_minus1 + 1);
  currSlice->listXsize[1] = (char) min<int>(currSlice->listXsize[1], currSlice->num_ref_idx_l1_active_minus1 + 1);

  // set the unused list entries to NULL
  for (i=currSlice->listXsize[0]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[0][i] = p_Vid->no_reference_picture;
  }
  for (i=currSlice->listXsize[1]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[1][i] = p_Vid->no_reference_picture;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initialize reference lists for a B slice_t
 *
 ************************************************************************
 */
void init_lists_b_slice_mvc(slice_t *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  dpb_t *p_Dpb = currSlice->p_Dpb;

  unsigned int i;
  int j;

  int list0idx = 0;
  int list0idx_1 = 0;
  int listltidx = 0;

  frame_store **fs_list0;
  frame_store **fs_list1;
  frame_store **fs_listlt;

  int currPOC = currSlice->ThisPOC;
  int anchor_pic_flag = currSlice->anchor_pic_flag;

  currSlice->listinterviewidx0 = 0;
  currSlice->listinterviewidx1 = 0;

  {
    // B-slice_t
    if (!currSlice->field_pic_flag)
    {
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (currSlice->framepoc >= p_Dpb->fs_ref[i]->frame->poc) //!KS use >= for error concealment
              //            if (p_Vid->framepoc > p_Dpb->fs_ref[i]->frame->poc)
            {
              currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)currSlice->listX[0], list0idx, sizeof(storable_picture*), compare_pic_by_poc_desc);
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (currSlice->framepoc < p_Dpb->fs_ref[i]->frame->poc)
            {
              currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)&currSlice->listX[0][list0idx_1], list0idx-list0idx_1, sizeof(storable_picture*), compare_pic_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        currSlice->listX[1][list0idx-list0idx_1+j]=currSlice->listX[0][j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        currSlice->listX[1][j-list0idx_1]=currSlice->listX[0][j];
      }

      currSlice->listXsize[0] = currSlice->listXsize[1] = (char) list0idx;

      //      printf("currSlice->listX[0] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<currSlice->listXsize[0]; i++){printf ("%d  ", currSlice->listX[0][i]->poc);} printf("\n");
      //      printf("currSlice->listX[1] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<currSlice->listXsize[1]; i++){printf ("%d  ", currSlice->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ltref[i]->is_used==3)
        {
          if (p_Dpb->fs_ltref[i]->frame->is_long_term)
          {
            currSlice->listX[0][list0idx]   = p_Dpb->fs_ltref[i]->frame;
            currSlice->listX[1][list0idx++] = p_Dpb->fs_ltref[i]->frame;
          }
        }
      }
      qsort((void *)&currSlice->listX[0][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(storable_picture*), compare_pic_by_lt_pic_num_asc);
      qsort((void *)&currSlice->listX[1][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(storable_picture*), compare_pic_by_lt_pic_num_asc);
      currSlice->listXsize[0] = currSlice->listXsize[1] = (char) list0idx;
    }
    else
    {
      fs_list0 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
      if (NULL==fs_list0)
        no_mem_exit("init_lists: fs_list0");
      fs_list1 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
      if (NULL==fs_list1)
        no_mem_exit("init_lists: fs_list1");
      fs_listlt = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
      if (NULL==fs_listlt)
        no_mem_exit("init_lists: fs_listlt");

      currSlice->listXsize[0] = 0;
      currSlice->listXsize[1] = 1;

      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (currSlice->ThisPOC >= p_Dpb->fs_ref[i]->poc)
          {
            fs_list0[list0idx++] = p_Dpb->fs_ref[i];
          }
        }
      }
      qsort((void *)fs_list0, list0idx, sizeof(frame_store*), compare_fs_by_poc_desc);
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (currSlice->ThisPOC < p_Dpb->fs_ref[i]->poc)
          {
            fs_list0[list0idx++] = p_Dpb->fs_ref[i];
          }
        }
      }
      qsort((void *)&fs_list0[list0idx_1], list0idx-list0idx_1, sizeof(frame_store*), compare_fs_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        fs_list1[list0idx-list0idx_1+j]=fs_list0[j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        fs_list1[j-list0idx_1]=fs_list0[j];
      }

      currSlice->listXsize[0] = 0;
      currSlice->listXsize[1] = 0;
      gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);
      gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list1, list0idx, currSlice->listX[1], &currSlice->listXsize[1], 0);

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
          fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
      }

      qsort((void *)fs_listlt, listltidx, sizeof(frame_store*), compare_fs_by_lt_pic_idx_asc);

      gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);
      gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[1], &currSlice->listXsize[1], 1);

      free(fs_list0);
      free(fs_list1);
      free(fs_listlt);
    }
  }

  if ((currSlice->listXsize[0] == currSlice->listXsize[1]) && (currSlice->listXsize[0] > 1))
  {
    // check if lists are identical, if yes swap first two elements of currSlice->listX[1]
    int diff=0;
    for (j = 0; j< currSlice->listXsize[0]; j++)
    {
      if (currSlice->listX[0][j] != currSlice->listX[1][j])
      {
        diff = 1;
        break;
      }
    }
    if (!diff)
    {
      storable_picture *tmp_s = currSlice->listX[1][0];
      currSlice->listX[1][0]=currSlice->listX[1][1];
      currSlice->listX[1][1]=tmp_s;
    }
  }

#if !SIMULCAST_ENABLE
  if (currSlice->svc_extension_flag == 0)
  {
    int curr_view_id = currSlice->view_id;
    // B-slice_t
    currSlice->fs_listinterview0 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL==currSlice->fs_listinterview0)
      no_mem_exit("init_lists: fs_listinterview0");
    currSlice->fs_listinterview1 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL==currSlice->fs_listinterview1)
      no_mem_exit("init_lists: fs_listinterview1");
    list0idx = currSlice->listXsize[0];

    if (!currSlice->field_pic_flag)
    {
      append_interview_list(
        p_Vid->p_Dpb_layer[1],
        false, false, 0, currSlice->fs_listinterview0, &currSlice->listinterviewidx0, currPOC, curr_view_id, anchor_pic_flag);
      append_interview_list(
        p_Vid->p_Dpb_layer[1],
        false, false, 1, currSlice->fs_listinterview1, &currSlice->listinterviewidx1, currPOC, curr_view_id, anchor_pic_flag);
      for (i=0; i<(unsigned int)currSlice->listinterviewidx0; i++)
      {
        currSlice->listX[0][list0idx++]=currSlice->fs_listinterview0[i]->frame;
      }
      currSlice->listXsize[0] = (char) list0idx;
      list0idx = currSlice->listXsize[1];
      for (i=0; i<(unsigned int)currSlice->listinterviewidx1; i++)
      {
        currSlice->listX[1][list0idx++] = currSlice->fs_listinterview1[i]->frame;
      }
      currSlice->listXsize[1] = (char) list0idx;
    }
    else
    {
      append_interview_list(
        p_Vid->p_Dpb_layer[1],
        currSlice->field_pic_flag, currSlice->bottom_field_flag, 0, currSlice->fs_listinterview0, &currSlice->listinterviewidx0, currPOC, curr_view_id, anchor_pic_flag);
      gen_pic_list_from_frame_interview_list(currSlice->bottom_field_flag, currSlice->fs_listinterview0, currSlice->listinterviewidx0, currSlice->listX[0], &currSlice->listXsize[0]);
      append_interview_list(
        p_Vid->p_Dpb_layer[1],
        currSlice->field_pic_flag, currSlice->bottom_field_flag, 1, currSlice->fs_listinterview1, &currSlice->listinterviewidx1, currPOC, curr_view_id, anchor_pic_flag);
      gen_pic_list_from_frame_interview_list(currSlice->bottom_field_flag, currSlice->fs_listinterview1, currSlice->listinterviewidx1, currSlice->listX[1], &currSlice->listXsize[1]);
    }    
  }
#endif
  // set max size
  currSlice->listXsize[0] = (char) min<int>(currSlice->listXsize[0], currSlice->num_ref_idx_l0_active_minus1 + 1);
  currSlice->listXsize[1] = (char) min<int>(currSlice->listXsize[1], currSlice->num_ref_idx_l1_active_minus1 + 1);

  // set the unused list entries to NULL
  for (i=currSlice->listXsize[0]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[0][i] = p_Vid->no_reference_picture;
  }
  for (i=currSlice->listXsize[1]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[1][i] = p_Vid->no_reference_picture;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Reordering process for inter-view reference pictures
 *
 ************************************************************************
 */
static void reorder_interview(VideoParameters *p_Vid, slice_t *currSlice, storable_picture **RefPicListX, int num_ref_idx_lX_active_minus1, int *refIdxLX, int targetViewID, int currPOC, int listidx)
{
  int cIdx, nIdx;
  storable_picture *picLX;

  picLX = get_inter_view_pic(p_Vid, currSlice, targetViewID, currPOC, listidx);

  if (picLX)
  {
    for( cIdx = num_ref_idx_lX_active_minus1+1; cIdx > *refIdxLX; cIdx-- )
      RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

    RefPicListX[ (*refIdxLX)++ ] = picLX;

    nIdx = *refIdxLX;

    for( cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1+1; cIdx++ )
    {
      if((GetViewIdx( p_Vid, RefPicListX[cIdx]->view_id ) != targetViewID) || (RefPicListX[cIdx]->poc != currPOC))
        RefPicListX[ nIdx++ ] = RefPicListX[ cIdx ];
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reordering process for MVC reference picture lists
 *
 ************************************************************************
 */
void reorder_ref_pic_list_mvc(slice_t *currSlice, int cur_list, int **anchor_ref, int **non_anchor_ref, int view_id, int anchor_pic_flag, int currPOC, int listidx)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  uint8_t  *modification_of_pic_nums_idc = currSlice->modification_of_pic_nums_idc[cur_list];
  uint32_t *abs_diff_pic_num_minus1      = currSlice->abs_diff_pic_num_minus1[cur_list];
  uint32_t *long_term_pic_num            = currSlice->long_term_pic_num[cur_list];
  int num_ref_idx_lX_active_minus1 =
    cur_list == 0 ? currSlice->num_ref_idx_l0_active_minus1
                  : currSlice->num_ref_idx_l1_active_minus1;
  uint32_t *abs_diff_view_idx_minus1 = currSlice->abs_diff_view_idx_minus1[cur_list];

  int i;

  int maxPicNum, currPicNum, picNumLXNoWrap, picNumLXPred, picNumLX;
  int picViewIdxLX, targetViewID;
  int refIdxLX = 0;
  int maxViewIdx =0;
  int curr_VOIdx = -1;
  int picViewIdxLXPred=-1;

  if (!currSlice->field_pic_flag)
  {
    maxPicNum  = p_Vid->active_sps->MaxFrameNum;
    currPicNum = currSlice->frame_num;
  }
  else
  {
    maxPicNum  = 2 * p_Vid->active_sps->MaxFrameNum;
    currPicNum = 2 * currSlice->frame_num + 1;
  }

  if(currSlice->svc_extension_flag==0)
  {
    curr_VOIdx = view_id;
    maxViewIdx = get_maxViewIdx(p_Vid, view_id, anchor_pic_flag, 0);
    picViewIdxLXPred=-1;
  }

  picNumLXPred = currPicNum;

  for (i=0; modification_of_pic_nums_idc[i]!=3; i++)
  {
    if (modification_of_pic_nums_idc[i] > 5)
      error ("Invalid modification_of_pic_nums_idc command", 500);

    if (modification_of_pic_nums_idc[i] < 2)
    {
      if (modification_of_pic_nums_idc[i] == 0)
      {
        if( picNumLXPred < abs_diff_pic_num_minus1[i] + 1 )
          picNumLXNoWrap = picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 ) + maxPicNum;
        else
          picNumLXNoWrap = picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 );
      }
      else // (modification_of_pic_nums_idc[i] == 1)
      {
        if( picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 )  >=  maxPicNum )
          picNumLXNoWrap = picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 ) - maxPicNum;
        else
          picNumLXNoWrap = picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 );
      }
      picNumLXPred = picNumLXNoWrap;

      if( picNumLXNoWrap > currPicNum )
        picNumLX = picNumLXNoWrap - maxPicNum;
      else
        picNumLX = picNumLXNoWrap;

      reorder_short_term(currSlice, cur_list, num_ref_idx_lX_active_minus1, picNumLX, &refIdxLX, view_id);
    }
    else if (modification_of_pic_nums_idc[i] == 2) //(modification_of_pic_nums_idc[i] == 2)
    {
      reorder_long_term(currSlice, currSlice->listX[cur_list], num_ref_idx_lX_active_minus1, long_term_pic_num[i], &refIdxLX, view_id);
    }
    else 
    {
      if(modification_of_pic_nums_idc[i] == 4) //(modification_of_pic_nums_idc[i] == 4)
      {
        picViewIdxLX = picViewIdxLXPred - (abs_diff_view_idx_minus1[i] + 1);
        if( picViewIdxLX <0)
          picViewIdxLX += maxViewIdx;
      }
      else //(modification_of_pic_nums_idc[i] == 5)
      {
        picViewIdxLX = picViewIdxLXPred + (abs_diff_view_idx_minus1[i] + 1);
        if( picViewIdxLX >= maxViewIdx)
          picViewIdxLX -= maxViewIdx;
      }
      picViewIdxLXPred = picViewIdxLX;

      if (anchor_pic_flag)
        targetViewID = anchor_ref[curr_VOIdx][picViewIdxLX];
      else
        targetViewID = non_anchor_ref[curr_VOIdx][picViewIdxLX];

      reorder_interview(p_Vid, currSlice, currSlice->listX[cur_list], num_ref_idx_lX_active_minus1, &refIdxLX, targetViewID, currPOC, listidx);
    }
  }
  // that's a definition
  currSlice->listXsize[cur_list] = (char) (num_ref_idx_lX_active_minus1 + 1);
}

void reorder_lists_mvc(slice_t * currSlice, int currPOC)
{
  VideoParameters *p_Vid = currSlice->p_Vid;

  if ((currSlice->slice_type != I_SLICE)&&(currSlice->slice_type != SI_SLICE))
  {
    if (currSlice->ref_pic_list_modification_flag_l0)
    {
      reorder_ref_pic_list_mvc(currSlice, LIST_0,
        p_Vid->active_subset_sps->anchor_ref_l0,
        p_Vid->active_subset_sps->non_anchor_ref_l0,
        currSlice->view_id, currSlice->anchor_pic_flag, currPOC, 0);
    }
    if (p_Vid->no_reference_picture == currSlice->listX[0][currSlice->num_ref_idx_l0_active_minus1])
    {
      if (p_Vid->non_conforming_stream)
        printf("RefPicList0[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_l0_active_minus1);
      else
        error("RefPicList0[ num_ref_idx_l0_active_minus1 ] in MVC layer is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    currSlice->listXsize[0] = (char)currSlice->num_ref_idx_l0_active_minus1 + 1;
  }
  if (currSlice->slice_type == B_SLICE)
  {
    if (currSlice->ref_pic_list_modification_flag_l1)
    {
      reorder_ref_pic_list_mvc(currSlice, LIST_1,
        p_Vid->active_subset_sps->anchor_ref_l1,
        p_Vid->active_subset_sps->non_anchor_ref_l1,
        currSlice->view_id, currSlice->anchor_pic_flag, currPOC, 1);
    }
    if (p_Vid->no_reference_picture == currSlice->listX[1][currSlice->num_ref_idx_l1_active_minus1])
    {
      if (p_Vid->non_conforming_stream)
        printf("RefPicList1[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_l1_active_minus1);
      else
        error("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    currSlice->listXsize[1] = (char)currSlice->num_ref_idx_l1_active_minus1 + 1;
  }
}

#endif
 
