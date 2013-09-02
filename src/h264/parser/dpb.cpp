#include <limits.h>

#include "global.h"
#include "erc_api.h"
#include "slice.h"
#include "image.h"
#include "dpb.h"
#include "memalloc.h"
#include "output.h"

static inline int RSD(int x)
{
    return ((x&2)?(x|1):(x&(~1)));
}


#define MAX_LIST_SIZE 33


/*!
 ************************************************************************
 * \brief
 *    Check if one of the frames/fields in frame store is used for reference
 ************************************************************************
 */
static int is_used_for_reference(FrameStore* fs)
{
  if (fs->is_reference)
  {
    return 1;
  }

  if (fs->is_used == 3) // frame
  {
    if (fs->frame->used_for_reference)
    {
      return 1;
    }
  }

  if (fs->is_used & 1) // top field
  {
    if (fs->top_field)
    {
      if (fs->top_field->used_for_reference)
      {
        return 1;
      }
    }
  }

  if (fs->is_used & 2) // bottom field
  {
    if (fs->bottom_field)
    {
      if (fs->bottom_field->used_for_reference)
      {
        return 1;
      }
    }
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    find smallest POC in the DPB.
 ************************************************************************
 */
void get_smallest_poc(dpb_t *p_Dpb, int *poc,int * pos)
{
  uint32 i;

  if (p_Dpb->used_size<1)
  {
    error("Cannot determine smallest POC, DPB empty.",150);
  }

  *pos=-1;
  *poc = INT_MAX;
  for (i = 0; i < p_Dpb->used_size; i++)
  {
    if ((*poc > p_Dpb->fs[i]->poc)&&(!p_Dpb->fs[i]->is_output))
    {
      *poc = p_Dpb->fs[i]->poc;
      *pos = i;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Remove a picture from DPB which is no longer needed.
 ************************************************************************
 */
int remove_unused_frame_from_dpb(dpb_t *p_Dpb)
{
  uint32 i;

  // check for frames that were already output and no longer used for reference
  for (i = 0; i < p_Dpb->used_size; i++)
  {
    if (p_Dpb->fs[i]->is_output && (!is_used_for_reference(p_Dpb->fs[i])))
    {
      remove_frame_from_dpb(p_Dpb, i);
      return 1;
    }
  }
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Returns the size of the dpb depending on level and picture size
 *
 *
 ************************************************************************
 */
int getDpbSize(VideoParameters *p_Vid, sps_t *active_sps)
{
  int pic_size = (active_sps->pic_width_in_mbs_minus1 + 1) * (active_sps->pic_height_in_map_units_minus1 + 1) * (active_sps->frame_mbs_only_flag?1:2) * 384;

  int size = 0;

  switch (active_sps->level_idc)
  {
  case 9:
    size = 152064;
    break;
  case 10:
    size = 152064;
    break;
  case 11:
    if (!is_FREXT_profile(active_sps->profile_idc) && (active_sps->constraint_set3_flag == 1))
      size = 152064;
    else
      size = 345600;
    break;
  case 12:
    size = 912384;
    break;
  case 13:
    size = 912384;
    break;
  case 20:
    size = 912384;
    break;
  case 21:
    size = 1824768;
    break;
  case 22:
    size = 3110400;
    break;
  case 30:
    size = 3110400;
    break;
  case 31:
    size = 6912000;
    break;
  case 32:
    size = 7864320;
    break;
  case 40:
    size = 12582912;
    break;
  case 41:
    size = 12582912;
    break;
  case 42:
    size = 13369344;
    break;
  case 50:
    size = 42393600;
    break;
  case 51:
    size = 70778880;
    break;
  case 52:
    size = 70778880;
    break;
  default:
    error("undefined level", 500);
    break;
  }

  size /= pic_size;
#if MVC_EXTENSION_ENABLE
  if(p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH)
  {
    int num_views = p_Vid->active_subset_sps->num_views_minus1+1;
    size = min(2 * size, max<int>(1, round(log2(num_views))) * 16) / num_views;
  }
  else
#endif
  {
    size = min(size, 16);
  }

  if (active_sps->vui_parameters_present_flag && active_sps->vui_parameters.bitstream_restriction_flag)
  {
    int size_vui;
    if ((int)active_sps->vui_parameters.max_dec_frame_buffering > size)
    {
      error("max_dec_frame_buffering larger than MaxDpbSize", 500);
    }
    size_vui = max<int>(1, active_sps->vui_parameters.max_dec_frame_buffering);
#ifdef _DEBUG
    if(size_vui < size)
    {
      printf("Warning: max_dec_frame_buffering(%d) is less than DPB size(%d) calculated from Profile/Level.\n", size_vui, size);
    }
#endif
    size = size_vui;    
  }

  return size;
}

/*!
 ************************************************************************
 * \brief
 *    Check then number of frames marked "used for reference" and break
 *    if maximum is exceeded
 *
 ************************************************************************
 */


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for decoded picture buffer and initialize with sane values.
 *
 ************************************************************************
 */
void init_dpb(VideoParameters *p_Vid, dpb_t *p_Dpb, int type)
{
  unsigned i; 
  sps_t *sps = p_Vid->active_sps;

  p_Dpb->p_Vid = p_Vid;
  if (p_Dpb->init_done)
  {
    free_dpb(p_Dpb);
  }

  p_Dpb->size = getDpbSize(p_Vid, sps) + p_Vid->p_Inp->dpb_plus[type==2? 1: 0];
  p_Dpb->num_ref_frames = sps->max_num_ref_frames; 

#if (MVC_EXTENSION_ENABLE)
  if ((unsigned int)sps->max_dec_frame_buffering < sps->max_num_ref_frames)
#else
  if (p_Dpb->size < sps->max_num_ref_frames)
#endif
  {
    error ("DPB size at specified level is smaller than the specified number of reference frames. This is not allowed.\n", 1000);
  }

  p_Dpb->used_size = 0;
  p_Dpb->last_picture = NULL;

  p_Dpb->ref_frames_in_buffer = 0;
  p_Dpb->ltref_frames_in_buffer = 0;

  p_Dpb->fs = (FrameStore **)calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs)
    no_mem_exit("init_dpb: p_Dpb->fs");

  p_Dpb->fs_ref = (FrameStore **)calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ref)
    no_mem_exit("init_dpb: p_Dpb->fs_ref");

  p_Dpb->fs_ltref = (FrameStore **)calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ltref)
    no_mem_exit("init_dpb: p_Dpb->fs_ltref");

#if (MVC_EXTENSION_ENABLE)
  p_Dpb->fs_ilref = (FrameStore **)calloc(1, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ilref)
    no_mem_exit("init_dpb: p_Dpb->fs_ilref");
#endif

  for (i = 0; i < p_Dpb->size; i++)
  {
    p_Dpb->fs[i]       = alloc_frame_store();
    p_Dpb->fs_ref[i]   = NULL;
    p_Dpb->fs_ltref[i] = NULL;
    p_Dpb->fs[i]->layer_id = MVC_INIT_VIEW_ID;
#if (MVC_EXTENSION_ENABLE)
    p_Dpb->fs[i]->view_id = MVC_INIT_VIEW_ID;
    p_Dpb->fs[i]->inter_view_flag[0] = p_Dpb->fs[i]->inter_view_flag[1] = 0;
    p_Dpb->fs[i]->anchor_pic_flag[0] = p_Dpb->fs[i]->anchor_pic_flag[1] = 0;
#endif
  }
#if (MVC_EXTENSION_ENABLE)
  if (type == 2)
  {
    p_Dpb->fs_ilref[0] = alloc_frame_store();
    // These may need some cleanups
    p_Dpb->fs_ilref[0]->view_id = MVC_INIT_VIEW_ID;
    p_Dpb->fs_ilref[0]->inter_view_flag[0] = p_Dpb->fs_ilref[0]->inter_view_flag[1] = 0;
    p_Dpb->fs_ilref[0]->anchor_pic_flag[0] = p_Dpb->fs_ilref[0]->anchor_pic_flag[1] = 0;
    // given that this is in a different buffer, do we even need proc_flag anymore?    
  }
  else
    p_Dpb->fs_ilref[0] = NULL;
#endif

  /* allocate a dummy storable picture */
  if(!p_Vid->no_reference_picture)
  {
    p_Vid->no_reference_picture = alloc_storable_picture (p_Vid, FRAME,
        sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
        sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);
    p_Vid->no_reference_picture->top_field    = p_Vid->no_reference_picture;
    p_Vid->no_reference_picture->bottom_field = p_Vid->no_reference_picture;
    p_Vid->no_reference_picture->frame        = p_Vid->no_reference_picture;
  }
  p_Dpb->last_output_poc = INT_MIN;

#if (MVC_EXTENSION_ENABLE)
  p_Dpb->last_output_view_id = -1;
#endif

  p_Vid->last_has_mmco_5 = 0;

  p_Dpb->init_done = 1;

  // picture error concealment
  if(p_Vid->conceal_mode !=0 && !p_Vid->last_out_fs)
    p_Vid->last_out_fs = alloc_frame_store();

}

void re_init_dpb(VideoParameters *p_Vid, dpb_t *p_Dpb, int type)
{
  int i; 
  sps_t *active_sps = p_Vid->active_sps;
  int iDpbSize;

  iDpbSize = getDpbSize(p_Vid, active_sps)+p_Vid->p_Inp->dpb_plus[type==2? 1: 0];
  p_Dpb->num_ref_frames = active_sps->max_num_ref_frames;
  if( iDpbSize > (int)p_Dpb->size)
  {
#if (MVC_EXTENSION_ENABLE)
    if ((unsigned int)active_sps->max_dec_frame_buffering < active_sps->max_num_ref_frames)
#else
    if (p_Dpb->size < active_sps->max_num_ref_frames)
#endif
    {
      error ("DPB size at specified level is smaller than the specified number of reference frames. This is not allowed.\n", 1000);
    }

    p_Dpb->fs = (FrameStore **)realloc(p_Dpb->fs, iDpbSize * sizeof (FrameStore*));
    if (NULL==p_Dpb->fs)
      no_mem_exit("re_init_dpb: p_Dpb->fs");

    p_Dpb->fs_ref = (FrameStore **)realloc(p_Dpb->fs_ref, iDpbSize * sizeof (FrameStore*));
    if (NULL==p_Dpb->fs_ref)
      no_mem_exit("re_init_dpb: p_Dpb->fs_ref");

    p_Dpb->fs_ltref = (FrameStore **)realloc(p_Dpb->fs_ltref, iDpbSize * sizeof (FrameStore*));
    if (NULL==p_Dpb->fs_ltref)
      no_mem_exit("re_init_dpb: p_Dpb->fs_ltref");

#if (MVC_EXTENSION_ENABLE)
    if(!p_Dpb->fs_ilref)
    {
      p_Dpb->fs_ilref = (FrameStore **)calloc(1, sizeof (FrameStore*));
      if (NULL==p_Dpb->fs_ilref)
        no_mem_exit("init_dpb: p_Dpb->fs_ilref");
    }
#endif

    for (i = p_Dpb->size; i < iDpbSize; i++)
    {
      p_Dpb->fs[i]       = alloc_frame_store();
      p_Dpb->fs_ref[i]   = NULL;
      p_Dpb->fs_ltref[i] = NULL;
#if (MVC_EXTENSION_ENABLE)
      p_Dpb->fs[i]->view_id = MVC_INIT_VIEW_ID;
      p_Dpb->fs[i]->inter_view_flag[0] = p_Dpb->fs[i]->inter_view_flag[1] = 0;
      p_Dpb->fs[i]->anchor_pic_flag[0] = p_Dpb->fs[i]->anchor_pic_flag[1] = 0;
#endif
    }

#if (MVC_EXTENSION_ENABLE)
  if (type == 2 && !p_Dpb->fs_ilref[0])
  {
    p_Dpb->fs_ilref[0] = alloc_frame_store();
    // These may need some cleanups
    p_Dpb->fs_ilref[0]->view_id = MVC_INIT_VIEW_ID;
    p_Dpb->fs_ilref[0]->inter_view_flag[0] = p_Dpb->fs_ilref[0]->inter_view_flag[1] = 0;
    p_Dpb->fs_ilref[0]->anchor_pic_flag[0] = p_Dpb->fs_ilref[0]->anchor_pic_flag[1] = 0;
    // given that this is in a different buffer, do we even need proc_flag anymore?    
  }
  else
    p_Dpb->fs_ilref[0] = NULL;
#endif

    p_Dpb->size = iDpbSize;
    p_Dpb->last_output_poc = INT_MIN;
#if (MVC_EXTENSION_ENABLE)
    p_Dpb->last_output_view_id = -1;
#endif
    p_Dpb->init_done = 1;

  }
}

/*!
 ************************************************************************
 * \brief
 *    Free memory for decoded picture buffer.
 ************************************************************************
 */
void free_dpb(dpb_t *p_Dpb)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  unsigned i;
  if (p_Dpb->fs)
  {
    for (i=0; i<p_Dpb->size; i++)
    {
      free_frame_store(p_Dpb->fs[i]);
    }
    free (p_Dpb->fs);
    p_Dpb->fs=NULL;
  }

  if (p_Dpb->fs_ref)
  {
    free (p_Dpb->fs_ref);
  }
  if (p_Dpb->fs_ltref)
  {
    free (p_Dpb->fs_ltref);
  }

#if (MVC_EXTENSION_ENABLE)
  if (p_Dpb->fs_ilref)
  {
    for (i=0; i<1; i++)
    {
      free_frame_store(p_Dpb->fs_ilref[i]);
    }
    free (p_Dpb->fs_ilref);
    p_Dpb->fs_ilref=NULL;
  }

  p_Dpb->last_output_view_id = -1;
#endif

  p_Dpb->last_output_poc = INT_MIN;

  p_Dpb->init_done = 0;

  // picture error concealment
  if(p_Vid->conceal_mode != 0 || p_Vid->last_out_fs)
      free_frame_store(p_Vid->last_out_fs);

  if(p_Vid->no_reference_picture)
  {
    free_storable_picture(p_Vid->no_reference_picture);
    p_Vid->no_reference_picture = NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for decoded picture buffer frame stores and initialize with sane values.
 *
 * \return
 *    the allocated FrameStore structure
 ************************************************************************
 */
FrameStore* alloc_frame_store(void)
{
  FrameStore *f;

  f = (FrameStore *)calloc (1, sizeof(FrameStore));
  if (NULL==f)
    no_mem_exit("alloc_frame_store: f");

  f->is_used      = 0;
  f->is_reference = 0;
  f->is_long_term = 0;
  f->is_orig_reference = 0;

  f->is_output = 0;

  f->frame        = NULL;;
  f->top_field    = NULL;
  f->bottom_field = NULL;

  return f;
}

void alloc_pic_motion(PicMotionParamsOld *motion, int size_y, int size_x)
{
  motion->mb_field_decoding_flag = (byte *)calloc (size_y * size_x, sizeof(byte));
  if (motion->mb_field_decoding_flag == NULL)
    no_mem_exit("alloc_storable_picture: motion->mb_field_decoding_flag");
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for a stored picture.
 *
 * \param p_Vid
 *    VideoParameters
 * \param structure
 *    picture structure
 * \param size_x
 *    horizontal luma size
 * \param size_y
 *    vertical luma size
 * \param size_x_cr
 *    horizontal chroma size
 * \param size_y_cr
 *    vertical chroma size
 *
 * \return
 *    the allocated StorablePicture structure
 ************************************************************************
 */
StorablePicture* alloc_storable_picture(VideoParameters *p_Vid, PictureStructure structure, int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output)
{
  sps_t *sps = p_Vid->active_sps;  

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (sps->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (sps->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

  StorablePicture *s;
  int   nplane;

  s = (StorablePicture *)calloc (1, sizeof(StorablePicture));
  if (NULL==s)
    no_mem_exit("alloc_storable_picture: s");

  if (structure!=FRAME)
  {
    size_y    /= 2;
    size_y_cr /= 2;
  }

  s->PicSizeInMbs = (size_x*size_y)/256;
  s->imgUV = NULL;

  get_mem2Dpel_pad (&(s->imgY), size_y, size_x, MCBUF_LUMA_PAD_Y, MCBUF_LUMA_PAD_X);
  s->iLumaStride = size_x+2*MCBUF_LUMA_PAD_X;
  s->iLumaExpandedHeight = size_y+2*MCBUF_LUMA_PAD_Y;

  if (sps->chroma_format_idc != YUV400)
  {
    get_mem3Dpel_pad(&(s->imgUV), 2, size_y_cr, size_x_cr, iChromaPadY, iChromaPadX);
  }

  s->chroma_format_idc = sps->chroma_format_idc;
  s->iChromaStride =size_x_cr + 2*iChromaPadX;
  s->iChromaExpandedHeight = size_y_cr + 2*iChromaPadY;
  s->iChromaPadY = iChromaPadY;
  s->iChromaPadX = iChromaPadX;

  s->separate_colour_plane_flag = sps->separate_colour_plane_flag;

  get_mem2Dmp     ( &s->mv_info, (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));
  alloc_pic_motion( &s->motion , (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));

  if( (sps->separate_colour_plane_flag != 0) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      get_mem2Dmp      (&s->JVmv_info[nplane], (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));
      alloc_pic_motion(&s->JVmotion[nplane] , (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));
    }
  }

  s->PicNum   = 0;
  s->frame_num = 0;
  s->LongTermFrameIdx = 0;
  s->LongTermPicNum   = 0;
  s->used_for_reference  = 0;
  s->is_long_term        = 0;
  s->non_existing        = 0;
  s->is_output           = 0;
#if (MVC_EXTENSION_ENABLE)
  s->view_id = -1;
#endif

  s->structure=structure;

  s->size_x = size_x;
  s->size_y = size_y;
  s->size_x_cr = size_x_cr;
  s->size_y_cr = size_y_cr;
  s->size_x_m1 = size_x - 1;
  s->size_y_m1 = size_y - 1;
  s->size_x_cr_m1 = size_x_cr - 1;
  s->size_y_cr_m1 = size_y_cr - 1;

  s->top_field    = p_Vid->no_reference_picture;
  s->bottom_field = p_Vid->no_reference_picture;
  s->frame        = p_Vid->no_reference_picture;

  s->dec_ref_pic_marking_buffer = NULL;

  s->coded_frame  = 0;
  s->mb_aff_frame_flag  = 0;

  s->top_poc = s->bottom_poc = s->poc = 0;
  s->seiHasTone_mapping = 0;

  if(!sps->frame_mbs_only_flag && structure != FRAME)
  {
    int i, j;
    for(j = 0; j < MAX_NUM_SLICES; j++)
    {
      for (i = 0; i < 2; i++)
      {
        s->listX[j][i] = (StorablePicture **)calloc(MAX_LIST_SIZE, sizeof (StorablePicture*)); // +1 for reordering
        if (NULL==s->listX[j][i])
        no_mem_exit("alloc_storable_picture: s->listX[i]");
      }
    }
  }

  return s;
}

/*!
 ************************************************************************
 * \brief
 *    Free frame store memory.
 *
 * \param p_Vid
 *    VideoParameters
 * \param f
 *    FrameStore to be freed
 *
 ************************************************************************
 */
void free_frame_store(FrameStore* f)
{
  if (f)
  {
    if (f->frame)
    {
      free_storable_picture(f->frame);
      f->frame=NULL;
    }
    if (f->top_field)
    {
      free_storable_picture(f->top_field);
      f->top_field=NULL;
    }
    if (f->bottom_field)
    {
      free_storable_picture(f->bottom_field);
      f->bottom_field=NULL;
    }
    free(f);
  }
}

void free_pic_motion(PicMotionParamsOld *motion)
{
  if (motion->mb_field_decoding_flag)
  {
    free(motion->mb_field_decoding_flag);
    motion->mb_field_decoding_flag = NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Free picture memory.
 *
 * \param p
 *    Picture to be freed
 *
 ************************************************************************
 */
void free_storable_picture(StorablePicture* p)
{
    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (p->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (p->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

  int nplane;
  if (p)
  {
    if (p->mv_info)
    {
      free_mem2Dmp(p->mv_info);
      p->mv_info = NULL;
    }
    free_pic_motion(&p->motion);

    if( (p->separate_colour_plane_flag != 0) )
    {
      for( nplane=0; nplane<MAX_PLANE; nplane++ )
      {
        if (p->JVmv_info[nplane])
        {
          free_mem2Dmp(p->JVmv_info[nplane]);
          p->JVmv_info[nplane] = NULL;
        }
        free_pic_motion(&p->JVmotion[nplane]);
      }
    }

    if (p->imgY)
    {
      free_mem2Dpel_pad(p->imgY, MCBUF_LUMA_PAD_Y, MCBUF_LUMA_PAD_X);
      p->imgY = NULL;
    }

    if (p->imgUV)
    {
      free_mem3Dpel_pad(p->imgUV, 2, p->iChromaPadY, p->iChromaPadX);
      p->imgUV=NULL;
    }


    if (p->seiHasTone_mapping)
      free(p->tone_mapping_lut);

    {
      int i, j;
      for(j = 0; j < MAX_NUM_SLICES; j++)
      {
        for(i=0; i<2; i++)
        {
          if(p->listX[j][i])
          {
            free(p->listX[j][i]);
            p->listX[j][i] = NULL;
          }
        }
      }
    }
    free(p);
    p = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    mark FrameStore unused for reference
 *
 ************************************************************************
 */
void unmark_for_reference(FrameStore* fs)
{
  if (fs->is_used & 1)
  {
    if (fs->top_field)
    {
      fs->top_field->used_for_reference = 0;
    }
  }
  if (fs->is_used & 2)
  {
    if (fs->bottom_field)
    {
      fs->bottom_field->used_for_reference = 0;
    }
  }
  if (fs->is_used == 3)
  {
    if (fs->top_field && fs->bottom_field)
    {
      fs->top_field->used_for_reference = 0;
      fs->bottom_field->used_for_reference = 0;
    }
    fs->frame->used_for_reference = 0;
  }

  fs->is_reference = 0;

  if(fs->frame)
  {
    free_pic_motion(&fs->frame->motion);
  }

  if (fs->top_field)
  {
    free_pic_motion(&fs->top_field->motion);
  }

  if (fs->bottom_field)
  {
    free_pic_motion(&fs->bottom_field->motion);
  }
}


/*!
 ************************************************************************
 * \brief
 *    mark FrameStore unused for reference and reset long term flags
 *
 ************************************************************************
 */
void unmark_for_long_term_reference(FrameStore* fs)
{
  if (fs->is_used & 1)
  {
    if (fs->top_field)
    {
      fs->top_field->used_for_reference = 0;
      fs->top_field->is_long_term = 0;
    }
  }
  if (fs->is_used & 2)
  {
    if (fs->bottom_field)
    {
      fs->bottom_field->used_for_reference = 0;
      fs->bottom_field->is_long_term = 0;
    }
  }
  if (fs->is_used == 3)
  {
    if (fs->top_field && fs->bottom_field)
    {
      fs->top_field->used_for_reference = 0;
      fs->top_field->is_long_term = 0;
      fs->bottom_field->used_for_reference = 0;
      fs->bottom_field->is_long_term = 0;
    }
    fs->frame->used_for_reference = 0;
    fs->frame->is_long_term = 0;
  }

  fs->is_reference = 0;
  fs->is_long_term = 0;
}


static void gen_field_ref_ids(VideoParameters *p_Vid, StorablePicture *p)
{
  int i,j;
   //! Generate Frame parameters from field information.

  //copy the list;
  for(j=0; j<p_Vid->iSliceNumOfCurrPic; j++)
  {
    if(p->listX[j][LIST_0])
    {
      p->listXsize[j][LIST_0] =  p_Vid->ppSliceList[j]->listXsize[LIST_0];
      for(i=0; i<p->listXsize[j][LIST_0]; i++)
        p->listX[j][LIST_0][i] = p_Vid->ppSliceList[j]->listX[LIST_0][i];
    }
    if(p->listX[j][LIST_1])
    {
      p->listXsize[j][LIST_1] =  p_Vid->ppSliceList[j]->listXsize[LIST_1];
      for(i=0; i<p->listXsize[j][LIST_1]; i++)
        p->listX[j][LIST_1][i] = p_Vid->ppSliceList[j]->listX[LIST_1][i];
    }
  }
}

void insert_picture_in_dpb(VideoParameters *p_Vid, FrameStore* fs, StorablePicture* p)
{
  assert (p!=NULL);
  assert (fs!=NULL);
  switch (p->structure)
  {
  case FRAME:
    fs->frame = p;
    fs->is_used = 3;
    if (p->used_for_reference)
    {
      fs->is_reference = 3;
      fs->is_orig_reference = 3;
      if (p->is_long_term)
      {
        fs->is_long_term = 3;
        fs->LongTermFrameIdx = p->LongTermFrameIdx;
      }
    }
    fs->layer_id = p->layer_id;
#if (MVC_EXTENSION_ENABLE)
    fs->view_id = p->view_id;
    fs->inter_view_flag[0] = fs->inter_view_flag[1] = p->inter_view_flag;
    fs->anchor_pic_flag[0] = fs->anchor_pic_flag[1] = p->anchor_pic_flag;
#endif
    // generate field views
    dpb_split_field(p_Vid, fs);
    break;
  case TOP_FIELD:
    fs->top_field = p;
    fs->is_used |= 1;
    fs->layer_id = p->layer_id;
#if (MVC_EXTENSION_ENABLE)
    fs->view_id = p->view_id;
    fs->inter_view_flag[0] = p->inter_view_flag;
    fs->anchor_pic_flag[0] = p->anchor_pic_flag;
#endif
    if (p->used_for_reference)
    {
      fs->is_reference |= 1;
      fs->is_orig_reference |= 1;
      if (p->is_long_term)
      {
        fs->is_long_term |= 1;
        fs->LongTermFrameIdx = p->LongTermFrameIdx;
      }
    }
    if (fs->is_used == 3)
    {
      // generate frame view
      dpb_combine_field(p_Vid, fs);
    }
    else
    {
      fs->poc = p->poc;
    }
    gen_field_ref_ids(p_Vid, p);
    break;
  case BOTTOM_FIELD:
    fs->bottom_field = p;
    fs->is_used |= 2;
    fs->layer_id = p->layer_id;
#if (MVC_EXTENSION_ENABLE)
    fs->view_id = p->view_id;
    fs->inter_view_flag[1] = p->inter_view_flag;
    fs->anchor_pic_flag[1] = p->anchor_pic_flag;
#endif
    if (p->used_for_reference)
    {
      fs->is_reference |= 2;
      fs->is_orig_reference |= 2;
      if (p->is_long_term)
      {
        fs->is_long_term |= 2;
        fs->LongTermFrameIdx = p->LongTermFrameIdx;
      }
    }
    if (fs->is_used == 3)
    {
      // generate frame view
      dpb_combine_field(p_Vid, fs);
    }
    else
    {
      fs->poc = p->poc;
    }
    gen_field_ref_ids(p_Vid, p);
    break;
  }
  fs->FrameNum = p->PicNum;
  fs->recovery_frame = p->recovery_frame;

  fs->is_output = p->is_output;

  if (fs->is_used==3)
  {
    calculate_frame_no(p_Vid, p);
  }
}


/*!
 ************************************************************************
 * \brief
 *    remove one frame from DPB
 ************************************************************************
 */
void remove_frame_from_dpb(dpb_t *p_Dpb, int pos)
{
  FrameStore* fs = p_Dpb->fs[pos];
  FrameStore* tmp;
  unsigned i;

  //printf ("remove frame with frame_num #%d\n", fs->frame_num);
  switch (fs->is_used)
  {
  case 3:
    free_storable_picture(fs->frame);
    free_storable_picture(fs->top_field);
    free_storable_picture(fs->bottom_field);
    fs->frame=NULL;
    fs->top_field=NULL;
    fs->bottom_field=NULL;
    break;
  case 2:
    free_storable_picture(fs->bottom_field);
    fs->bottom_field=NULL;
    break;
  case 1:
    free_storable_picture(fs->top_field);
    fs->top_field=NULL;
    break;
  case 0:
    break;
  default:
    error("invalid frame store type",500);
  }
  fs->is_used = 0;
  fs->is_long_term = 0;
  fs->is_reference = 0;
  fs->is_orig_reference = 0;

  // move empty framestore to end of buffer
  tmp = p_Dpb->fs[pos];

  for (i=pos; i<p_Dpb->used_size-1;i++)
  {
    p_Dpb->fs[i] = p_Dpb->fs[i+1];
  }
  p_Dpb->fs[p_Dpb->used_size-1] = tmp;
  p_Dpb->used_size--;
}




/*!
 ************************************************************************
 * \brief
 *    Output one picture stored in the DPB.
 ************************************************************************
 */
int output_one_frame_from_dpb(dpb_t *p_Dpb)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  int poc, pos;
  //diagnostics
  if (p_Dpb->used_size < 1)
  {
    error("Cannot output frame, DPB empty.",150);
  }

  // find smallest POC
  get_smallest_poc(p_Dpb, &poc, &pos);

  if(pos==-1)
  {
    return 0;
  }

  // call the output function
//  printf ("output frame with frame_num #%d, poc %d (dpb. p_Dpb->size=%d, p_Dpb->used_size=%d)\n", p_Dpb->fs[pos]->frame_num, p_Dpb->fs[pos]->frame->poc, p_Dpb->size, p_Dpb->used_size);

#if (DISABLE_ERC == 0)
  // picture error concealment
  if(p_Vid->conceal_mode != 0)
  {
    if(p_Dpb->last_output_poc == 0)
    {
      write_lost_ref_after_idr(p_Dpb, pos);
    }
#if (MVC_EXTENSION_ENABLE)
    write_lost_non_ref_pic(p_Dpb, poc, p_Vid->p_out_mvc[p_Dpb->layer_id]);
#else
    write_lost_non_ref_pic(p_Dpb, poc, p_Vid->p_out);
#endif
  }
#endif
// JVT-P072 ends

#if (MVC_EXTENSION_ENABLE)
  write_stored_frame(p_Vid, p_Dpb->fs[pos], p_Vid->p_out_mvc[p_Dpb->layer_id]);
#else
  write_stored_frame(p_Vid, p_Dpb->fs[pos], p_Vid->p_out);
#endif

  // picture error concealment
  if(p_Vid->conceal_mode == 0)
  {
    if (p_Dpb->last_output_poc >= poc)
    {
      error ("output POC must be in ascending order", 150);
    }
  }

  p_Dpb->last_output_poc = poc;

  // free frame store and move empty store to end of buffer
  if (!is_used_for_reference(p_Dpb->fs[pos]))
  {
    remove_frame_from_dpb(p_Dpb, pos);
  }
  return 1;
}



/*!
 ************************************************************************
 * \brief
 *    All stored picture are output. Should be called to empty the buffer
 ************************************************************************
 */
void flush_dpb(dpb_t *p_Dpb)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  uint32 i;

  if(!p_Dpb->init_done)
    return;
#if (DISABLE_ERC == 0)
  if (p_Vid->conceal_mode != 0)
    conceal_non_ref_pics(p_Dpb, 0);
#endif
  // mark all frames unused
  for (i=0; i<p_Dpb->used_size; i++)
  {
#if MVC_EXTENSION_ENABLE
    assert( p_Dpb->fs[i]->view_id == p_Dpb->layer_id);
#endif
    unmark_for_reference (p_Dpb->fs[i]);
  }

  while (remove_unused_frame_from_dpb(p_Dpb)) ;

  // output frames in POC order
  while (p_Dpb->used_size && output_one_frame_from_dpb(p_Dpb)) ;

  p_Dpb->last_output_poc = INT_MIN;
}

#if (MVC_EXTENSION_ENABLE)
void flush_dpbs(dpb_t **p_Dpb_layers, int nLayers)
{
  VideoParameters *p_Vid = p_Dpb_layers[0]->p_Vid;
  dpb_t *p_Dpb;
  int i, j, used_size;

#if (DISABLE_ERC == 0)
  if (p_Vid->conceal_mode != 0)
  {
    conceal_non_ref_pics(p_Dpb_layers[0], 0);
  }
#endif
  // mark all frames unused
  for(j=0; j<nLayers; j++)
  {
    p_Dpb = p_Dpb_layers[j];
    if(p_Dpb->init_done)
    {
      for (i=0; i<(int)p_Dpb->used_size; i++)
      {
        unmark_for_reference (p_Dpb->fs[i]);
      }
      while (remove_unused_frame_from_dpb(p_Dpb)) ;
    }
  }
  // output frames in POC order
  used_size = p_Dpb_layers[0]->used_size;
  for(j=1; j<nLayers; j++)
  {
    p_Dpb = p_Dpb_layers[j];
    if(p_Dpb->init_done)
    {
      if(p_Dpb->used_size &&  (p_Dpb->used_size != used_size))
      {
        assert(!"DPB used_size is not equal!");
      }
      if((int)p_Dpb->used_size > used_size)
      {
        used_size = (int)p_Dpb->used_size;
      }
    }
  }  
  while (used_size)
  {
    for(j=0; j<nLayers; j++)
    {
      p_Dpb = p_Dpb_layers[j];
      if(p_Dpb->used_size)
        output_one_frame_from_dpb(p_Dpb);
    }
    used_size--;
  }
  for(j=0; j<nLayers; j++)
  {
    p_Dpb = p_Dpb_layers[j];
    p_Dpb->last_output_poc = INT_MIN;
  }  
}
#endif


/*!
 ************************************************************************
 * \brief
 *    Extract top field from a frame
 ************************************************************************
 */
void dpb_split_field(VideoParameters *p_Vid, FrameStore *fs)
{
  int i, j, ii, jj, jj4;
  int idiv,jdiv;
  int currentmb;
  int twosz16 = 2 * (fs->frame->size_x >> 4);
  StorablePicture *fs_top = NULL, *fs_btm = NULL; 
  StorablePicture *frame = fs->frame;

  fs->poc = frame->poc;

  if (!frame->frame_mbs_only_flag)
  {
    fs_top = fs->top_field    = alloc_storable_picture(p_Vid, TOP_FIELD,    frame->size_x, frame->size_y, frame->size_x_cr, frame->size_y_cr, 1);
    fs_btm = fs->bottom_field = alloc_storable_picture(p_Vid, BOTTOM_FIELD, frame->size_x, frame->size_y, frame->size_x_cr, frame->size_y_cr, 1);

    for (i = 0; i < (frame->size_y >> 1); i++)
    {
      memcpy(fs_top->imgY[i], frame->imgY[i*2], frame->size_x*sizeof(imgpel));
    }

    for (i = 0; i< (frame->size_y_cr >> 1); i++)
    {
      memcpy(fs_top->imgUV[0][i], frame->imgUV[0][i*2], frame->size_x_cr*sizeof(imgpel));
      memcpy(fs_top->imgUV[1][i], frame->imgUV[1][i*2], frame->size_x_cr*sizeof(imgpel));
    }

    for (i = 0; i < (frame->size_y>>1); i++)
    {
      memcpy(fs_btm->imgY[i], frame->imgY[i*2 + 1], frame->size_x*sizeof(imgpel));
    }

    for (i = 0; i < (frame->size_y_cr>>1); i++)
    {
      memcpy(fs_btm->imgUV[0][i], frame->imgUV[0][i*2 + 1], frame->size_x_cr*sizeof(imgpel));
      memcpy(fs_btm->imgUV[1][i], frame->imgUV[1][i*2 + 1], frame->size_x_cr*sizeof(imgpel));
    }

    fs_top->poc = frame->top_poc;
    fs_btm->poc = frame->bottom_poc;

#if (MVC_EXTENSION_ENABLE)
    fs_top->view_id = frame->view_id;
    fs_btm->view_id = frame->view_id;
#endif

    fs_top->frame_poc =  frame->frame_poc;

    fs_top->bottom_poc = fs_btm->bottom_poc =  frame->bottom_poc;
    fs_top->top_poc    = fs_btm->top_poc    =  frame->top_poc;
    fs_btm->frame_poc  = frame->frame_poc;

    fs_top->used_for_reference = fs_btm->used_for_reference
                                      = frame->used_for_reference;
    fs_top->is_long_term = fs_btm->is_long_term
                                = frame->is_long_term;
    fs->LongTermFrameIdx = fs_top->LongTermFrameIdx
                            = fs_btm->LongTermFrameIdx
                            = frame->LongTermFrameIdx;

    fs_top->coded_frame = fs_btm->coded_frame = 1;
    fs_top->mb_aff_frame_flag = fs_btm->mb_aff_frame_flag
                        = frame->mb_aff_frame_flag;

    frame->top_field    = fs_top;
    frame->bottom_field = fs_btm;
    frame->frame         = frame;
    fs_top->bottom_field = fs_btm;
    fs_top->frame        = frame;
    fs_top->top_field = fs_top;
    fs_btm->top_field = fs_top;
    fs_btm->frame     = frame;
    fs_btm->bottom_field = fs_btm;

#if (MVC_EXTENSION_ENABLE)
    fs_top->view_id = fs_btm->view_id = fs->view_id;
    fs_top->inter_view_flag = fs->inter_view_flag[0];
    fs_btm->inter_view_flag = fs->inter_view_flag[1];
#endif

    fs_top->chroma_format_idc = fs_btm->chroma_format_idc = frame->chroma_format_idc;
    fs_top->iCodingType = fs_btm->iCodingType = frame->iCodingType;
    if(frame->used_for_reference)
    {
      pad_dec_picture(p_Vid, fs_top);
      pad_dec_picture(p_Vid, fs_btm);
    }
  }
  else
  {
    fs->top_field       = NULL;
    fs->bottom_field    = NULL;
    frame->top_field    = NULL;
    frame->bottom_field = NULL;
    frame->frame = frame;
  }

  if (!frame->frame_mbs_only_flag)
  {
    if (frame->mb_aff_frame_flag)
    {
      PicMotionParamsOld *frm_motion = &frame->motion;
      for (j=0 ; j< (frame->size_y >> 3); j++)
      {
        jj = (j >> 2)*8 + (j & 0x03);
        jj4 = jj + 4;
        jdiv = (j >> 1);
        for (i=0 ; i < (frame->size_x>>2); i++)
        {
          idiv = (i >> 2);

          currentmb = twosz16*(jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);
          // Assign field mvs attached to MB-Frame buffer to the proper buffer
          if (frm_motion->mb_field_decoding_flag[currentmb])
          {
            fs_btm->mv_info[j][i].mv[LIST_0] = frame->mv_info[jj4][i].mv[LIST_0];
            fs_btm->mv_info[j][i].mv[LIST_1] = frame->mv_info[jj4][i].mv[LIST_1];
            fs_btm->mv_info[j][i].ref_idx[LIST_0] = frame->mv_info[jj4][i].ref_idx[LIST_0];
            if(fs_btm->mv_info[j][i].ref_idx[LIST_0] >=0)
              fs_btm->mv_info[j][i].ref_pic[LIST_0] = p_Vid->ppSliceList[frame->mv_info[jj4][i].slice_no]->listX[4][(short) fs_btm->mv_info[j][i].ref_idx[LIST_0]];
            else
              fs_btm->mv_info[j][i].ref_pic[LIST_0] = NULL;
            fs_btm->mv_info[j][i].ref_idx[LIST_1] = frame->mv_info[jj4][i].ref_idx[LIST_1];
            if(fs_btm->mv_info[j][i].ref_idx[LIST_1] >=0)
              fs_btm->mv_info[j][i].ref_pic[LIST_1] = p_Vid->ppSliceList[frame->mv_info[jj4][i].slice_no]->listX[5][(short) fs_btm->mv_info[j][i].ref_idx[LIST_1]];
            else
              fs_btm->mv_info[j][i].ref_pic[LIST_1] = NULL;
          
            fs_top->mv_info[j][i].mv[LIST_0] = frame->mv_info[jj][i].mv[LIST_0];
            fs_top->mv_info[j][i].mv[LIST_1] = frame->mv_info[jj][i].mv[LIST_1];
            fs_top->mv_info[j][i].ref_idx[LIST_0] = frame->mv_info[jj][i].ref_idx[LIST_0];
            if(fs_top->mv_info[j][i].ref_idx[LIST_0] >=0)
              fs_top->mv_info[j][i].ref_pic[LIST_0] = p_Vid->ppSliceList[frame->mv_info[jj][i].slice_no]->listX[2][(short) fs_top->mv_info[j][i].ref_idx[LIST_0]];
            else
              fs_top->mv_info[j][i].ref_pic[LIST_0] = NULL;
            fs_top->mv_info[j][i].ref_idx[LIST_1] = frame->mv_info[jj][i].ref_idx[LIST_1];
            if(fs_top->mv_info[j][i].ref_idx[LIST_1] >=0)
              fs_top->mv_info[j][i].ref_pic[LIST_1] = p_Vid->ppSliceList[frame->mv_info[jj][i].slice_no]->listX[3][(short) fs_top->mv_info[j][i].ref_idx[LIST_1]];
            else
              fs_top->mv_info[j][i].ref_pic[LIST_1] = NULL;
          }
        }
      }
    }
  
      //! Generate field MVs from Frame MVs
    for (j=0 ; j < (frame->size_y >> 3) ; j++)
    {
      jj = 2* RSD(j);
      jdiv = (j >> 1);
      for (i=0 ; i < (frame->size_x >> 2) ; i++)
      {
        ii = RSD(i);
        idiv = (i >> 2);

        currentmb = twosz16 * (jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);

        if (!frame->mb_aff_frame_flag  || !frame->motion.mb_field_decoding_flag[currentmb])
        {
          fs_top->mv_info[j][i].mv[LIST_0] = fs_btm->mv_info[j][i].mv[LIST_0] = frame->mv_info[jj][ii].mv[LIST_0];
          fs_top->mv_info[j][i].mv[LIST_1] = fs_btm->mv_info[j][i].mv[LIST_1] = frame->mv_info[jj][ii].mv[LIST_1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)
          if (frame->mv_info[jj][ii].ref_idx[LIST_0] == -1)
          {
            fs_top->mv_info[j][i].ref_idx[LIST_0] = fs_btm->mv_info[j][i].ref_idx[LIST_0] = - 1;
            fs_top->mv_info[j][i].ref_pic[LIST_0] = fs_btm->mv_info[j][i].ref_pic[LIST_0] = NULL;
          }
          else
          {
            fs_top->mv_info[j][i].ref_idx[LIST_0] = fs_btm->mv_info[j][i].ref_idx[LIST_0] = frame->mv_info[jj][ii].ref_idx[LIST_0];
            fs_top->mv_info[j][i].ref_pic[LIST_0] = fs_btm->mv_info[j][i].ref_pic[LIST_0] = p_Vid->ppSliceList[frame->mv_info[jj][ii].slice_no]->listX[LIST_0][(short) frame->mv_info[jj][ii].ref_idx[LIST_0]];
          }

          if (frame->mv_info[jj][ii].ref_idx[LIST_1] == -1)
          {
            fs_top->mv_info[j][i].ref_idx[LIST_1] = fs_btm->mv_info[j][i].ref_idx[LIST_1] = - 1;
            fs_top->mv_info[j][i].ref_pic[LIST_1] = fs_btm->mv_info[j][i].ref_pic[LIST_1] = NULL;
          }
          else
          {
            fs_top->mv_info[j][i].ref_idx[LIST_1] = fs_btm->mv_info[j][i].ref_idx[LIST_1] = frame->mv_info[jj][ii].ref_idx[LIST_1];
            fs_top->mv_info[j][i].ref_pic[LIST_1] = fs_btm->mv_info[j][i].ref_pic[LIST_1] = p_Vid->ppSliceList[frame->mv_info[jj][ii].slice_no]->listX[LIST_1][(short) frame->mv_info[jj][ii].ref_idx[LIST_1]];
          }
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields,
 *    YUV components and display information only
 ************************************************************************
 */
void dpb_combine_field_yuv(VideoParameters *p_Vid, FrameStore *fs)
{
  int i, j;

  if (!fs->frame)
  {
    fs->frame = alloc_storable_picture(p_Vid, FRAME, fs->top_field->size_x, fs->top_field->size_y*2, fs->top_field->size_x_cr, fs->top_field->size_y_cr*2, 1);
  }

  for (i=0; i<fs->top_field->size_y; i++)
  {
    memcpy(fs->frame->imgY[i*2],     fs->top_field->imgY[i]   , fs->top_field->size_x * sizeof(imgpel));     // top field
    memcpy(fs->frame->imgY[i*2 + 1], fs->bottom_field->imgY[i], fs->bottom_field->size_x * sizeof(imgpel)); // bottom field
  }

  for (j = 0; j < 2; j++)
  {
    for (i=0; i<fs->top_field->size_y_cr; i++)
    {
      memcpy(fs->frame->imgUV[j][i*2],     fs->top_field->imgUV[j][i],    fs->top_field->size_x_cr*sizeof(imgpel));
      memcpy(fs->frame->imgUV[j][i*2 + 1], fs->bottom_field->imgUV[j][i], fs->bottom_field->size_x_cr*sizeof(imgpel));
    }
  }
  fs->poc=fs->frame->poc =fs->frame->frame_poc = min (fs->top_field->poc, fs->bottom_field->poc);

  fs->bottom_field->frame_poc=fs->top_field->frame_poc=fs->frame->poc;

  fs->bottom_field->top_poc=fs->frame->top_poc=fs->top_field->poc;
  fs->top_field->bottom_poc=fs->frame->bottom_poc=fs->bottom_field->poc;

  fs->frame->used_for_reference = (fs->top_field->used_for_reference && fs->bottom_field->used_for_reference );
  fs->frame->is_long_term = (fs->top_field->is_long_term && fs->bottom_field->is_long_term );

  if (fs->frame->is_long_term)
    fs->frame->LongTermFrameIdx = fs->LongTermFrameIdx;

  fs->frame->top_field    = fs->top_field;
  fs->frame->bottom_field = fs->bottom_field;
  fs->frame->frame = fs->frame;

  fs->frame->coded_frame = 0;

  fs->frame->chroma_format_idc = fs->top_field->chroma_format_idc;
  fs->frame->frame_cropping_flag = fs->top_field->frame_cropping_flag;
  if (fs->frame->frame_cropping_flag)
  {
    fs->frame->frame_crop_top_offset = fs->top_field->frame_crop_top_offset;
    fs->frame->frame_crop_bottom_offset = fs->top_field->frame_crop_bottom_offset;
    fs->frame->frame_crop_left_offset = fs->top_field->frame_crop_left_offset;
    fs->frame->frame_crop_right_offset = fs->top_field->frame_crop_right_offset;
  }

  fs->top_field->frame = fs->bottom_field->frame = fs->frame;
  fs->top_field->top_field = fs->top_field;
  fs->top_field->bottom_field = fs->bottom_field;
  fs->bottom_field->top_field = fs->top_field;
  fs->bottom_field->bottom_field = fs->bottom_field;
  if(fs->top_field->used_for_reference || fs->bottom_field->used_for_reference)
  {
    pad_dec_picture(p_Vid, fs->frame);
  }

}


/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields
 ************************************************************************
 */
void dpb_combine_field(VideoParameters *p_Vid, FrameStore *fs)
{
  int i,j, jj, jj4, k, l;

  dpb_combine_field_yuv(p_Vid, fs);

#if (MVC_EXTENSION_ENABLE)
  fs->frame->view_id = fs->view_id;
#endif
  fs->frame->iCodingType = fs->top_field->iCodingType; //FIELD_CODING;
   //! Use inference flag to remap mvs/references

  //! Generate Frame parameters from field information.

  for (j=0 ; j < (fs->top_field->size_y >> 2) ; j++)
  {
    jj = (j<<1);
    jj4 = jj + 1;
    for (i=0 ; i< (fs->top_field->size_x >> 2) ; i++)
    {
      fs->frame->mv_info[jj][i].mv[LIST_0] = fs->top_field->mv_info[j][i].mv[LIST_0];
      fs->frame->mv_info[jj][i].mv[LIST_1] = fs->top_field->mv_info[j][i].mv[LIST_1];

      fs->frame->mv_info[jj][i].ref_idx[LIST_0] = fs->top_field->mv_info[j][i].ref_idx[LIST_0];
      fs->frame->mv_info[jj][i].ref_idx[LIST_1] = fs->top_field->mv_info[j][i].ref_idx[LIST_1];

      /* bug: top field list doesnot exist.*/
      l = fs->top_field->mv_info[j][i].slice_no;
      k = fs->top_field->mv_info[j][i].ref_idx[LIST_0];
      fs->frame->mv_info[jj][i].ref_pic[LIST_0] = k>=0? fs->top_field->listX[l][LIST_0][k]: NULL;  
      k = fs->top_field->mv_info[j][i].ref_idx[LIST_1];
      fs->frame->mv_info[jj][i].ref_pic[LIST_1] = k>=0? fs->top_field->listX[l][LIST_1][k]: NULL;

      //! association with id already known for fields.
      fs->frame->mv_info[jj4][i].mv[LIST_0] = fs->bottom_field->mv_info[j][i].mv[LIST_0];
      fs->frame->mv_info[jj4][i].mv[LIST_1] = fs->bottom_field->mv_info[j][i].mv[LIST_1];

      fs->frame->mv_info[jj4][i].ref_idx[LIST_0]  = fs->bottom_field->mv_info[j][i].ref_idx[LIST_0];
      fs->frame->mv_info[jj4][i].ref_idx[LIST_1]  = fs->bottom_field->mv_info[j][i].ref_idx[LIST_1];
      l = fs->bottom_field->mv_info[j][i].slice_no;

      k = fs->bottom_field->mv_info[j][i].ref_idx[LIST_0];
      fs->frame->mv_info[jj4][i].ref_pic[LIST_0] = k>=0? fs->bottom_field->listX[l][LIST_0][k]: NULL;
      k = fs->bottom_field->mv_info[j][i].ref_idx[LIST_1];
      fs->frame->mv_info[jj4][i].ref_pic[LIST_1] = k>=0? fs->bottom_field->listX[l][LIST_1][k]: NULL;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *      Tian Dong
 *          June 13, 2002, Modified on July 30, 2003
 *
 *      If a gap in frame_num is found, try to fill the gap
 * \param p_Vid
 *    VideoParameters structure
 *
 ************************************************************************
 */
void fill_frame_num_gap(VideoParameters *p_Vid, slice_t *currSlice)
{
  sps_t *sps = p_Vid->active_sps;
  
  int CurrFrameNum;
  int UnusedShortTermFrameNum;
  StorablePicture *picture = NULL;
  int tmp1 = currSlice->delta_pic_order_cnt[0];
  int tmp2 = currSlice->delta_pic_order_cnt[1];
  currSlice->delta_pic_order_cnt[0] = currSlice->delta_pic_order_cnt[1] = 0;

  printf("A gap in frame number is found, try to fill it.\n");

  UnusedShortTermFrameNum = (p_Vid->pre_frame_num + 1) % sps->MaxFrameNum;
  CurrFrameNum = currSlice->frame_num;

  while (CurrFrameNum != UnusedShortTermFrameNum)
  {
    picture = alloc_storable_picture (p_Vid, FRAME,
        sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
        sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);
    picture->coded_frame = 1;
    picture->PicNum = UnusedShortTermFrameNum;
    picture->frame_num = UnusedShortTermFrameNum;
    picture->non_existing = 1;
    picture->is_output = 1;
    picture->used_for_reference = 1;
    picture->adaptive_ref_pic_buffering_flag = 0;
#if (MVC_EXTENSION_ENABLE)
    picture->view_id = currSlice->view_id;
#endif

    currSlice->frame_num = UnusedShortTermFrameNum;
    if (sps->pic_order_cnt_type!=0)
    {
      decode_poc(p_Vid, p_Vid->ppSliceList[0]);
    }
    picture->top_poc    = currSlice->TopFieldOrderCnt;
    picture->bottom_poc = currSlice->BottomFieldOrderCnt;
    picture->frame_poc  = currSlice->framepoc;
    picture->poc        = currSlice->framepoc;

    store_picture_in_dpb(currSlice->p_Dpb, picture);

    picture=NULL;
    p_Vid->pre_frame_num = UnusedShortTermFrameNum;
    UnusedShortTermFrameNum = (UnusedShortTermFrameNum + 1) % sps->MaxFrameNum;
  }
  currSlice->delta_pic_order_cnt[0] = tmp1;
  currSlice->delta_pic_order_cnt[1] = tmp2;
  currSlice->frame_num = CurrFrameNum;
}



#if (MVC_EXTENSION_ENABLE)
int GetMaxDecFrameBuffering(VideoParameters *p_Vid)
{
  int i, j, iMax, iMax_1 = 0, iMax_2 = 0;
  sub_sps_t *curr_subset_sps;
  sps_t *curr_sps;

  curr_subset_sps = p_Vid->SubsetSeqParSet;
  curr_sps = p_Vid->SeqParSet;
  for(i=0; i<MAXSPS; i++)
  {
    if(curr_subset_sps->Valid && curr_subset_sps->sps.seq_parameter_set_id < MAXSPS)
    {
      j = curr_subset_sps->sps.max_dec_frame_buffering;

      if (curr_subset_sps->sps.vui_parameters_present_flag && curr_subset_sps->sps.vui_parameters.bitstream_restriction_flag)
      {
        if ((int)curr_subset_sps->sps.vui_parameters.max_dec_frame_buffering > j)
        {
          error ("max_dec_frame_buffering larger than MaxDpbSize", 500);
        }
        j = max<int>(1, curr_subset_sps->sps.vui_parameters.max_dec_frame_buffering);
      }

      if(j > iMax_2)
        iMax_2 = j;
    }
    
    if(curr_sps->Valid)
    {
      j = curr_sps->max_dec_frame_buffering;

      if (curr_sps->vui_parameters_present_flag && curr_sps->vui_parameters.bitstream_restriction_flag)
      {
        if ((int)curr_sps->vui_parameters.max_dec_frame_buffering > j)
        {
          error ("max_dec_frame_buffering larger than MaxDpbSize", 500);
        }
        j = max<int>(1, curr_sps->vui_parameters.max_dec_frame_buffering);
      }

      if(j > iMax_1)
        iMax_1 = j;
    }
    curr_subset_sps++;
    curr_sps++;
  }  
      
  if (iMax_1 > 0 && iMax_2 > 0)
    iMax = iMax_1 + iMax_2;
  else
    iMax = (iMax_1 >0? iMax_1*2 : iMax_2*2);
  return iMax;
}

static int is_view_id_in_ref_view_list(int view_id, int *ref_view_id, int num_ref_views)
{
   int i;
   for(i=0; i<num_ref_views; i++)
   {
     if(view_id == ref_view_id[i])
       break;
   }

   return (num_ref_views && (i<num_ref_views));
}

void append_interview_list(dpb_t *p_Dpb, 
                           bool field_pic_flag,
                           bool bottom_field_flag,
                           int list_idx, 
                           FrameStore **list,
                           int *listXsize, 
                           int currPOC, 
                           int curr_view_id, 
                           int anchor_pic_flag)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  int iVOIdx = curr_view_id;
  int pic_avail;
  int poc = 0;
  int fld_idx;
  int num_ref_views, *ref_view_id;
  FrameStore *fs = p_Dpb->fs_ilref[0];


  if(iVOIdx <0)
    printf("Error: iVOIdx: %d is not less than 0\n", iVOIdx);

  if(anchor_pic_flag)
  {
    num_ref_views = list_idx? p_Vid->active_subset_sps->num_anchor_refs_l1[iVOIdx] : p_Vid->active_subset_sps->num_anchor_refs_l0[iVOIdx];
    ref_view_id   = list_idx? p_Vid->active_subset_sps->anchor_ref_l1[iVOIdx]:p_Vid->active_subset_sps->anchor_ref_l0[iVOIdx];
  }
  else
  {
    num_ref_views = list_idx? p_Vid->active_subset_sps->num_non_anchor_refs_l1[iVOIdx] : p_Vid->active_subset_sps->num_non_anchor_refs_l0[iVOIdx];
    ref_view_id = list_idx? p_Vid->active_subset_sps->non_anchor_ref_l1[iVOIdx]:p_Vid->active_subset_sps->non_anchor_ref_l0[iVOIdx];
  }

  if(bottom_field_flag)
    fld_idx = 1;
  else
    fld_idx = 0;

    if(!field_pic_flag)
    {
      pic_avail = (fs->is_used == 3);
      if (pic_avail)
        poc = fs->frame->poc;
    }
    else if(!bottom_field_flag)
    {
      pic_avail = fs->is_used & 1;
      if (pic_avail)
        poc = fs->top_field->poc;
    }
    else
    {
      pic_avail = fs->is_used & 2;
      if (pic_avail)
        poc = fs->bottom_field->poc;
    }

    if(pic_avail && fs->inter_view_flag[fld_idx])
    {
      if(poc == currPOC)
      {
        if(is_view_id_in_ref_view_list(fs->view_id, ref_view_id, num_ref_views))
        {
          //add one inter-view reference;
          list[*listXsize] = fs; 
          //next;
          (*listXsize)++;
        }
      }
    }
}

#endif

void process_picture_in_dpb_s(VideoParameters *p_Vid, StorablePicture *p_pic)
{
  ImageData *p_img_out = &p_Vid->tempData3;
  imgpel***  d_img;
  int i;

  if(p_Vid->tempData3.frm_data[0] == NULL)
    init_img_data( p_Vid, &(p_Vid->tempData3), p_Vid->active_sps);

  if (p_pic->structure == FRAME)
  {
    d_img = p_img_out->frm_data;
  }
  else //If reference picture is a field, then frm_data will actually contain field data and therefore top/bottom stride is set accordingly.
  {
    if (p_pic->structure == TOP_FIELD)
    {
      d_img = p_img_out->top_data;
    }
    else
    {
      d_img = p_img_out->bot_data;
    }
  }

  for(i=0; i<p_pic->size_y; i++)
    memcpy(d_img[0][i], p_pic->imgY[i], p_pic->size_x*sizeof(imgpel));
  if (p_Vid->active_sps->chroma_format_idc != YUV400)
  {
    for(i=0; i<p_pic->size_y_cr; i++)
      memcpy(d_img[1][i], p_pic->imgUV[0][i], p_pic->size_x_cr * sizeof(imgpel));
    for(i=0; i<p_pic->size_y_cr; i++)
      memcpy(d_img[2][i], p_pic->imgUV[1][i], p_pic->size_x_cr * sizeof(imgpel));
  }
}

int init_img_data(VideoParameters *p_Vid, ImageData *p_ImgData, sps_t *sps)
{
  int memory_size = 0;
  int nplane;
  
  // allocate memory for reference frame buffers: p_ImgData->frm_data
  p_ImgData->yuv_format    = (ColorFormat) sps->chroma_format_idc;
  p_ImgData->frm_stride[0] = sps->PicWidthInMbs * 16;
  p_ImgData->frm_stride[1] = p_ImgData->frm_stride[2] = sps->PicWidthInMbs * sps->MbWidthC;
  p_ImgData->top_stride[0] = p_ImgData->bot_stride[0] = p_ImgData->frm_stride[0] << 1;
  p_ImgData->top_stride[1] = p_ImgData->top_stride[2] = p_ImgData->bot_stride[1] = p_ImgData->bot_stride[2] = p_ImgData->frm_stride[1] << 1;

  if( sps->separate_colour_plane_flag )
  {
    for( nplane=0; nplane < MAX_PLANE; nplane++ )
    {
      memory_size += get_mem2Dpel(&(p_ImgData->frm_data[nplane]), sps->FrameHeightInMbs * 16, sps->PicWidthInMbs * 16);
    }
  }
  else
  {
    memory_size += get_mem2Dpel(&(p_ImgData->frm_data[0]), sps->FrameHeightInMbs * 16, sps->PicWidthInMbs * 16);

    if (sps->chroma_format_idc != YUV400)
    {
      int i, j, k;
      memory_size += get_mem2Dpel(&(p_ImgData->frm_data[1]), sps->FrameHeightInMbs * sps->MbHeightC, sps->PicWidthInMbs * sps->MbWidthC);
      memory_size += get_mem2Dpel(&(p_ImgData->frm_data[2]), sps->FrameHeightInMbs * sps->MbHeightC, sps->PicWidthInMbs * sps->MbWidthC);

      if (sizeof(imgpel) == sizeof(unsigned char))
      {
        for (k = 1; k < 3; k++)
          memset(p_ImgData->frm_data[k][0], 128, sps->FrameHeightInMbs * sps->MbHeightC * sps->PicWidthInMbs * sps->MbWidthC * sizeof(imgpel));
      }
      else
      {
        imgpel mean_val;

        for (k = 1; k < 3; k++)
        {
          mean_val = (imgpel) (sps->BitDepthC);

          for (j = 0; j < sps->FrameHeightInMbs * sps->MbHeightC; j++)
          {
            for (i = 0; i < sps->PicWidthInMbs * sps->MbWidthC; i++)
              p_ImgData->frm_data[k][j][i] = mean_val;
          }
        }
      }
    }
  }

  if (!sps->frame_mbs_only_flag)
  {
    // allocate memory for field reference frame buffers
    memory_size += init_top_bot_planes(p_ImgData->frm_data[0], sps->FrameHeightInMbs * 16, &(p_ImgData->top_data[0]), &(p_ImgData->bot_data[0]));

    if (sps->chroma_format_idc != YUV400)
    {
      memory_size += 4*(sizeof(imgpel**));

      memory_size += init_top_bot_planes(p_ImgData->frm_data[1], sps->FrameHeightInMbs * sps->MbHeightC, &(p_ImgData->top_data[1]), &(p_ImgData->bot_data[1]));
      memory_size += init_top_bot_planes(p_ImgData->frm_data[2], sps->FrameHeightInMbs * sps->MbHeightC, &(p_ImgData->top_data[2]), &(p_ImgData->bot_data[2]));
    }
  }

  return memory_size;
}

void free_img_data(VideoParameters *p_Vid, ImageData *p_ImgData)
{
  if ( p_Vid->active_sps->separate_colour_plane_flag )
  {
    int nplane;

    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      if (p_ImgData->frm_data[nplane])
      {
        free_mem2Dpel(p_ImgData->frm_data[nplane]);      // free ref frame buffers
        p_ImgData->frm_data[nplane] = NULL;
      }
    }
  }
  else
  {
    if (p_ImgData->frm_data[0])
    {
      free_mem2Dpel(p_ImgData->frm_data[0]);      // free ref frame buffers
      p_ImgData->frm_data[0] = NULL;
    }
    
    if (p_ImgData->yuv_format != YUV400)
    {
      if (p_ImgData->frm_data[1])
      {
        free_mem2Dpel(p_ImgData->frm_data[1]);
        p_ImgData->frm_data[1] = NULL;
      }
      if (p_ImgData->frm_data[2])
      {
        free_mem2Dpel(p_ImgData->frm_data[2]);
        p_ImgData->frm_data[2] = NULL;
      }
    }
  }
  
  if (!p_Vid->active_sps->frame_mbs_only_flag)
  {
    free_top_bot_planes(p_ImgData->top_data[0], p_ImgData->bot_data[0]);

    if (p_ImgData->yuv_format != YUV400)
    {
      free_top_bot_planes(p_ImgData->top_data[1], p_ImgData->bot_data[1]);
      free_top_bot_planes(p_ImgData->top_data[2], p_ImgData->bot_data[2]);
    }
  }
}

static inline void copy_img_data(imgpel *out_img, imgpel *in_img, int ostride, int istride, unsigned int size_y, unsigned int size_x)
{
  unsigned int i;
  for(i = 0; i < size_y; i++)
  {
    memcpy(out_img, in_img, size_x);
    out_img += ostride;
    in_img += istride;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Remove a picture from DPB which is no longer needed.
 ************************************************************************
 */
#if !MVC_EXTENSION_ENABLE
int remove_unused_proc_pic_from_dpb(dpb_t *p_Dpb)
{
  assert(!"The function is not available\n");
  return 0;
}
#endif

/*!
 ************************************************************************
 * \brief
 *    Store a processed picture in DPB. This includes cheking for space in DPB and
 *    flushing frames.
 *    If we received a frame, we need to check for a new store, if we
 *    got a field, check if it's the second field of an already allocated
 *    store.
 *
 * \param p_Vid
 *    VideoParameters
 * \param p
 *    Picture to be stored
 *
 ************************************************************************
 */
#if (MVC_EXTENSION_ENABLE)
void store_proc_picture_in_dpb(dpb_t *p_Dpb, StorablePicture* p)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  FrameStore *fs = p_Dpb->fs_ilref[0];
  if(p_Dpb->used_size_il>0 && fs->is_used==3)
  {
    //checking;
#ifdef _DEBUG
    if(p->structure==FRAME)
      assert(fs->frame->frame_poc != p->poc);
    else if(p->structure==TOP_FIELD)
      assert(fs->top_field->top_poc != p->poc);
    else if(p->structure==BOTTOM_FIELD)
      assert(fs->bottom_field->bottom_poc != p->poc);
#endif
    if(fs->frame)
    {
      free_storable_picture(fs->frame);
      fs->frame = NULL;
    }
    if(fs->top_field)
    {
      free_storable_picture(fs->top_field);
      fs->top_field = NULL;
    }
    if(fs->bottom_field)
    {
      free_storable_picture(fs->bottom_field);
      fs->bottom_field = NULL;
    }
    fs->is_used = 0;
    fs->is_reference = 0;
    p_Dpb->used_size_il--;   
  }
#ifdef _DEBUG  
  if(fs->is_used>0)
  {
    //checking;
    if(p->structure==FRAME)
      assert(fs->frame == NULL);
    else if(p->structure==TOP_FIELD)
      assert(fs->top_field == NULL);
    else if(p->structure==BOTTOM_FIELD)
      assert(fs->bottom_field == NULL);
  }
#endif

  insert_picture_in_dpb(p_Vid, fs, p);
  if((p->structure==FRAME && fs->is_used == 3) || (p->structure!=FRAME && fs->is_used && fs->is_used <3))
   p_Dpb->used_size_il++;  
}

/*!
 ************************************************************************
 * \brief
 *    Clone an encoded frame picture structure
 ************************************************************************
 */
StorablePicture * clone_storable_picture( VideoParameters *p_Vid, StorablePicture *p_pic )
{
  int i, j;
  int nplane;
  int *istride = NULL;
  int ostride[2];
  imgpel ***img_in = NULL;

  sps_t *sps = p_Vid->active_sps;
  StorablePicture *p_stored_pic = alloc_storable_picture (p_Vid, (PictureStructure) p_Vid->structure,
      sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
      sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 0);

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (sps->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (sps->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

  p_stored_pic->PicNum = p_pic->PicNum;
  p_stored_pic->frame_num = p_pic->frame_num;
  p_stored_pic->LongTermFrameIdx = p_pic->LongTermFrameIdx;
  p_stored_pic->LongTermPicNum = p_pic->LongTermPicNum;
  p_stored_pic->is_long_term = 0;
  p_stored_pic->non_existing = p_pic->non_existing;
  p_stored_pic->structure = p_pic->structure;
  p_stored_pic->size_x = p_pic->size_x;
  p_stored_pic->size_y = p_pic->size_y;
  p_stored_pic->size_x_cr = p_pic->size_x_cr;
  p_stored_pic->size_y_cr = p_pic->size_y_cr;
  p_stored_pic->size_x_m1 = p_pic->size_x - 1;
  p_stored_pic->size_y_m1 = p_pic->size_y - 1;
  p_stored_pic->size_x_cr_m1 = p_pic->size_x_cr - 1;
  p_stored_pic->size_y_cr_m1 = p_pic->size_y_cr - 1;
  
  p_stored_pic->mb_aff_frame_flag = p_pic->mb_aff_frame_flag;
  p_stored_pic->seiHasTone_mapping = p_pic->seiHasTone_mapping;
  p_stored_pic->poc         = p_pic->poc;
  p_stored_pic->top_poc     = p_pic->top_poc;
  p_stored_pic->bottom_poc  = p_pic->bottom_poc;
  p_stored_pic->frame_poc   = p_pic->frame_poc;
  p_stored_pic->PicNum     = p_pic->PicNum;
  p_stored_pic->frame_num   = p_pic->frame_num;
  p_stored_pic->coded_frame = 1;
  p_stored_pic->qp = p_pic->qp;
  p_stored_pic->slice_qp_delta = p_pic->slice_qp_delta;
  p_stored_pic->chroma_qp_offset[0] = p_pic->chroma_qp_offset[0];
  p_stored_pic->chroma_qp_offset[1] = p_pic->chroma_qp_offset[1];

  p_stored_pic->slice_type = p_pic->slice_type;
  p_stored_pic->idr_flag = p_pic->idr_flag;
  p_stored_pic->no_output_of_prior_pics_flag = p_pic->no_output_of_prior_pics_flag;
  p_stored_pic->long_term_reference_flag = 0;
  p_stored_pic->adaptive_ref_pic_buffering_flag = 0;
  p_stored_pic->dec_ref_pic_marking_buffer = NULL;
  p_stored_pic->PicWidthInMbs = p_pic->PicWidthInMbs;
  p_stored_pic->recovery_frame = p_pic->recovery_frame;
  p_stored_pic->chroma_format_idc = p_pic->chroma_format_idc;
  p_stored_pic->frame_mbs_only_flag = p_pic->frame_mbs_only_flag;
  p_stored_pic->frame_cropping_flag = p_pic->frame_cropping_flag;

  if (p_stored_pic->frame_cropping_flag)
  {
    p_stored_pic->frame_crop_left_offset   = p_pic->frame_crop_left_offset;
    p_stored_pic->frame_crop_right_offset  = p_pic->frame_crop_right_offset;
    p_stored_pic->frame_crop_top_offset    = p_pic->frame_crop_top_offset;
    p_stored_pic->frame_crop_bottom_offset = p_pic->frame_crop_bottom_offset;
  }
  
  // store BL reconstruction

  ostride[0] = p_stored_pic->iLumaStride;
  ostride[1] = p_stored_pic->iChromaStride;
  if (p_stored_pic->structure == FRAME)
  {
    istride = p_Vid->tempData3.frm_stride;
    img_in  = p_Vid->tempData3.frm_data;
  }
  else if (p_stored_pic->structure == TOP_FIELD)
  {
    istride = p_Vid->tempData3.top_stride;
    img_in  = p_Vid->tempData3.top_data;
  }
  else
  {
    istride = p_Vid->tempData3.bot_stride;
    img_in  = p_Vid->tempData3.bot_data;
  }

  copy_img_data(&p_stored_pic->imgY[0][0], &img_in[0][0][0], ostride[0], istride[0], p_pic->size_y, p_pic->size_x * sizeof(imgpel)); 

  pad_buf(*p_stored_pic->imgY, p_stored_pic->size_x, p_stored_pic->size_y, p_stored_pic->iLumaStride, MCBUF_LUMA_PAD_X, MCBUF_LUMA_PAD_Y);

  if (p_Vid->active_sps->chroma_format_idc != YUV400)
  {    
    copy_img_data(&p_stored_pic->imgUV[0][0][0], &img_in[1][0][0], ostride[1], istride[1], p_pic->size_y_cr, p_pic->size_x_cr*sizeof(imgpel));
    pad_buf(*p_stored_pic->imgUV[0], p_stored_pic->size_x_cr, p_stored_pic->size_y_cr, p_stored_pic->iChromaStride, iChromaPadX, iChromaPadY);
    copy_img_data(&p_stored_pic->imgUV[1][0][0], &img_in[2][0][0], ostride[1], istride[2], p_pic->size_y_cr, p_pic->size_x_cr*sizeof(imgpel));
    pad_buf(*p_stored_pic->imgUV[1], p_stored_pic->size_x_cr, p_stored_pic->size_y_cr, p_stored_pic->iChromaStride, iChromaPadX, iChromaPadY);
  }

  for (j = 0; j < (p_pic->size_y >> BLOCK_SHIFT); j++)
  {
    char *ref_idx = p_stored_pic->mv_info[j][0].ref_idx;
    for (i = 0; i < (p_pic->size_x >> BLOCK_SHIFT); i++)
    {          
      *((short *) ref_idx) = -1;
      ref_idx += sizeof(PicMotionParams);
    }
  }

  if( (p_Vid->active_sps->separate_colour_plane_flag != 0) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      for (j = 0; j < (p_pic->size_y >> BLOCK_SHIFT); j++)
      {
        for (i = 0; i < (p_pic->size_x >> BLOCK_SHIFT); i++)
        {
          p_stored_pic->JVmv_info[nplane][j][i].ref_idx[LIST_0] = -1;
          p_stored_pic->JVmv_info[nplane][j][i].ref_idx[LIST_1] = -1;
        }
      }
    }
  }

  // MVC-related parameters
  p_stored_pic->inter_view_flag = p_pic->inter_view_flag;
  p_stored_pic->anchor_pic_flag = 0;
  p_stored_pic->view_id = 0;
  p_stored_pic->proc_flag = 1;
  p_stored_pic->is_output = 1;
  p_stored_pic->used_for_reference = 1;

  return p_stored_pic;
}
#endif



