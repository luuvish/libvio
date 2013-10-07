#include <limits.h>

#include "global.h"
#include "input_parameters.h"

#include "slice.h"
#include "image.h"
#include "dpb.h"
#include "memalloc.h"
#include "output.h"

using vio::h264::mb_t;

#include "erc_api.h"

static inline int RSD(int x)
{
    return ((x&2)?(x|1):(x&(~1)));
}


#define MAX_LIST_SIZE 33


static inline int is_FREXT_profile(unsigned int profile_idc) 
{
    return ( profile_idc >= FREXT_HP || profile_idc == FREXT_CAVLC444 );
}

/*!
 ************************************************************************
 * \brief
 *    Check if one of the frames/fields in frame store is used for reference
 ************************************************************************
 */
static int is_used_for_reference(frame_store* fs)
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
  uint32_t i;

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
  uint32_t i;

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


int getDpbSize(VideoParameters* p_Vid, sps_t* active_sps)
{
    int pic_size = (active_sps->pic_width_in_mbs_minus1 + 1) * (active_sps->pic_height_in_map_units_minus1 + 1) * (active_sps->frame_mbs_only_flag ? 1 : 2) * 384;

    int size = 0;

    switch (active_sps->level_idc) {
    case 9:
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
    case 13:
    case 20:
        size = 912384;
        break;
    case 21:
        size = 1824768;
        break;
    case 22:
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
    case 52:
        size = 70778880;
        break;
    default:
        error("undefined level", 500);
        break;
    }

    size /= pic_size;
#if MVC_EXTENSION_ENABLE
    if (p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH) {
        int num_views = p_Vid->active_subset_sps->num_views_minus1+1;
        size = min(2 * size, max<int>(1, round(log2(num_views))) * 16) / num_views;
    } else
#endif
        size = min(size, 16);

    if (active_sps->vui_parameters_present_flag && active_sps->vui_parameters.bitstream_restriction_flag) {
        int size_vui;
        if ((int)active_sps->vui_parameters.max_dec_frame_buffering > size)
            error("max_dec_frame_buffering larger than MaxDpbSize", 500);
        size_vui = max<int>(1, active_sps->vui_parameters.max_dec_frame_buffering);
#ifdef _DEBUG
        if (size_vui < size)
            printf("Warning: max_dec_frame_buffering(%d) is less than DPB size(%d) calculated from Profile/Level.\n", size_vui, size);
#endif
        size = size_vui;    
    }

    return size;
}


void init_dpb(VideoParameters* p_Vid, dpb_t* p_Dpb, int type)
{
  unsigned i; 
    sps_t* sps = p_Vid->active_sps;

    p_Dpb->p_Vid = p_Vid;
    if (p_Dpb->init_done)
        free_dpb(p_Dpb);

    p_Dpb->size = getDpbSize(p_Vid, sps) + p_Vid->p_Inp->dpb_plus[type == 2 ? 1 : 0];
    p_Dpb->num_ref_frames = sps->max_num_ref_frames; 

#if (MVC_EXTENSION_ENABLE)
    if ((unsigned int)sps->max_dec_frame_buffering < sps->max_num_ref_frames)
#else
    if (p_Dpb->size < sps->max_num_ref_frames)
#endif
        error ("DPB size at specified level is smaller than the specified number of reference frames. This is not allowed.\n", 1000);

    p_Dpb->used_size = 0;
    p_Dpb->last_picture = NULL;

    p_Dpb->ref_frames_in_buffer = 0;
    p_Dpb->ltref_frames_in_buffer = 0;

    p_Dpb->fs = (frame_store**)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL == p_Dpb->fs)
        no_mem_exit("init_dpb: p_Dpb->fs");

    p_Dpb->fs_ref = (frame_store**)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL == p_Dpb->fs_ref)
        no_mem_exit("init_dpb: p_Dpb->fs_ref");

    p_Dpb->fs_ltref = (frame_store**)calloc(p_Dpb->size, sizeof (frame_store*));
    if (NULL == p_Dpb->fs_ltref)
        no_mem_exit("init_dpb: p_Dpb->fs_ltref");

#if (MVC_EXTENSION_ENABLE)
    p_Dpb->fs_ilref = (frame_store**)calloc(1, sizeof (frame_store*));
    if (NULL == p_Dpb->fs_ilref)
        no_mem_exit("init_dpb: p_Dpb->fs_ilref");
#endif

    for (i = 0; i < p_Dpb->size; i++) {
        p_Dpb->fs[i]       = new frame_store {};
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
    if (type == 2) {
        p_Dpb->fs_ilref[0] = new frame_store {};
        // These may need some cleanups
        p_Dpb->fs_ilref[0]->view_id = MVC_INIT_VIEW_ID;
        p_Dpb->fs_ilref[0]->inter_view_flag[0] = p_Dpb->fs_ilref[0]->inter_view_flag[1] = 0;
        p_Dpb->fs_ilref[0]->anchor_pic_flag[0] = p_Dpb->fs_ilref[0]->anchor_pic_flag[1] = 0;
        // given that this is in a different buffer, do we even need proc_flag anymore?    
    } else
        p_Dpb->fs_ilref[0] = NULL;
#endif

    /* allocate a dummy storable picture */
    if (!p_Vid->no_reference_picture) {
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
    if (p_Vid->conceal_mode != 0 && !p_Vid->last_out_fs)
        p_Vid->last_out_fs = new frame_store {};
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

    p_Dpb->fs = (frame_store **)realloc(p_Dpb->fs, iDpbSize * sizeof (frame_store*));
    if (NULL==p_Dpb->fs)
      no_mem_exit("re_init_dpb: p_Dpb->fs");

    p_Dpb->fs_ref = (frame_store **)realloc(p_Dpb->fs_ref, iDpbSize * sizeof (frame_store*));
    if (NULL==p_Dpb->fs_ref)
      no_mem_exit("re_init_dpb: p_Dpb->fs_ref");

    p_Dpb->fs_ltref = (frame_store **)realloc(p_Dpb->fs_ltref, iDpbSize * sizeof (frame_store*));
    if (NULL==p_Dpb->fs_ltref)
      no_mem_exit("re_init_dpb: p_Dpb->fs_ltref");

#if (MVC_EXTENSION_ENABLE)
    if(!p_Dpb->fs_ilref)
    {
      p_Dpb->fs_ilref = (frame_store **)calloc(1, sizeof (frame_store*));
      if (NULL==p_Dpb->fs_ilref)
        no_mem_exit("init_dpb: p_Dpb->fs_ilref");
    }
#endif

    for (i = p_Dpb->size; i < iDpbSize; i++)
    {
      p_Dpb->fs[i]       = new frame_store {};
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
    p_Dpb->fs_ilref[0] = new frame_store {};
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

    if (p_Dpb->fs) {
        for (int i = 0; i < p_Dpb->size; i++)
            delete p_Dpb->fs[i];
        free (p_Dpb->fs);
        p_Dpb->fs=NULL;
    }

    if (p_Dpb->fs_ref)
        free (p_Dpb->fs_ref);
    if (p_Dpb->fs_ltref)
        free (p_Dpb->fs_ltref);

#if (MVC_EXTENSION_ENABLE)
    if (p_Dpb->fs_ilref) {
        for (int i = 0; i < 1; i++)
            delete p_Dpb->fs_ilref[i];
        free (p_Dpb->fs_ilref);
        p_Dpb->fs_ilref=NULL;
    }

    p_Dpb->last_output_view_id = -1;
#endif

    p_Dpb->last_output_poc = INT_MIN;

    p_Dpb->init_done = 0;

    // picture error concealment
    if (p_Vid->conceal_mode != 0 || p_Vid->last_out_fs)
        delete p_Vid->last_out_fs;

    if (p_Vid->no_reference_picture) {
        free_storable_picture(p_Vid->no_reference_picture);
        p_Vid->no_reference_picture = NULL;
    }
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
 *    the allocated storable_picture structure
 ************************************************************************
 */
storable_picture* alloc_storable_picture(VideoParameters *p_Vid, PictureStructure structure, int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output)
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

    storable_picture* s;
    int nplane;

    s = (storable_picture *)calloc (1, sizeof(storable_picture));
    if (NULL == s)
        no_mem_exit("alloc_storable_picture: s");

    if (structure != FRAME) {
        size_y    /= 2;
        size_y_cr /= 2;
    }

    s->PicSizeInMbs = (size_x * size_y) / 256;
    s->imgUV = NULL;

    get_mem2Dpel_pad(&s->imgY, size_y, size_x, MCBUF_LUMA_PAD_Y, MCBUF_LUMA_PAD_X);
    s->iLumaStride = size_x + 2 * MCBUF_LUMA_PAD_X;
    s->iLumaExpandedHeight = size_y + 2 * MCBUF_LUMA_PAD_Y;

    if (sps->chroma_format_idc != YUV400)
        get_mem3Dpel_pad(&s->imgUV, 2, size_y_cr, size_x_cr, iChromaPadY, iChromaPadX);

    s->chroma_format_idc     = sps->chroma_format_idc;
    s->iChromaStride         = size_x_cr + 2 * iChromaPadX;
    s->iChromaExpandedHeight = size_y_cr + 2 * iChromaPadY;
    s->iChromaPadY           = iChromaPadY;
    s->iChromaPadX           = iChromaPadX;

    s->separate_colour_plane_flag = sps->separate_colour_plane_flag;

    get_mem2Dmp(&s->mv_info, (size_y / 4), (size_x / 4));
    s->motion.mb_field_decoding_flag = new bool[(size_x / 4) * (size_y / 4)];

    if (sps->separate_colour_plane_flag) {
        for (nplane = 0; nplane < 3; ++nplane) {
            get_mem2Dmp(&s->JVmv_info[nplane], (size_y / 4), (size_x / 4));
            s->JVmotion[nplane].mb_field_decoding_flag = new bool[(size_x / 4) * (size_y / 4)];
        }
    }

    s->PicNum             = 0;
    s->frame_num          = 0;
    s->LongTermFrameIdx   = 0;
    s->LongTermPicNum     = 0;
    s->used_for_reference = 0;
    s->is_long_term       = 0;
    s->non_existing       = 0;
    s->is_output          = 0;
#if (MVC_EXTENSION_ENABLE)
    s->view_id            = -1;
#endif

    s->structure    = structure;

    s->size_x       = size_x;
    s->size_y       = size_y;
    s->size_x_cr    = size_x_cr;
    s->size_y_cr    = size_y_cr;
    s->size_x_m1    = size_x - 1;
    s->size_y_m1    = size_y - 1;
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

    return s;
}


frame_store::~frame_store()
{
    if (this->frame)
        free_storable_picture(this->frame);
    if (this->top_field)
        free_storable_picture(this->top_field);
    if (this->bottom_field)
        free_storable_picture(this->bottom_field);
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
void free_storable_picture(storable_picture* p)
{
  int nplane;
  if (p)
  {
    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (p->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (p->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

    if (p->mv_info)
    {
      free_mem2Dmp(p->mv_info);
      p->mv_info = NULL;
    }
    delete []p->motion.mb_field_decoding_flag;

    if( (p->separate_colour_plane_flag != 0) )
    {
      for( nplane=0; nplane<3; nplane++ )
      {
        if (p->JVmv_info[nplane])
        {
          free_mem2Dmp(p->JVmv_info[nplane]);
          p->JVmv_info[nplane] = NULL;
        }
        delete []p->JVmotion[nplane].mb_field_decoding_flag;
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

    free(p);
    p = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    mark frame_store unused for reference
 *
 ************************************************************************
 */
void unmark_for_reference(frame_store* fs)
{
    if ((fs->is_used & 1) && fs->top_field) {
        fs->top_field->used_for_reference = 0;
    }
    if ((fs->is_used & 2) && fs->bottom_field) {
        fs->bottom_field->used_for_reference = 0;
    }
    if (fs->is_used == 3) {
        if (fs->top_field && fs->bottom_field) {
            fs->top_field->used_for_reference = 0;
            fs->bottom_field->used_for_reference = 0;
        }
        fs->frame->used_for_reference = 0;
    }
    fs->is_reference = 0;

    if (fs->frame) {
        delete []fs->frame->motion.mb_field_decoding_flag;
        fs->frame->motion.mb_field_decoding_flag = nullptr;
    }

    if (fs->top_field) {
        delete []fs->top_field->motion.mb_field_decoding_flag;
        fs->top_field->motion.mb_field_decoding_flag = nullptr;
    }

    if (fs->bottom_field) {
        delete []fs->bottom_field->motion.mb_field_decoding_flag;
        fs->bottom_field->motion.mb_field_decoding_flag = nullptr;
    }
}


/*!
 ************************************************************************
 * \brief
 *    mark frame_store unused for reference and reset long term flags
 *
 ************************************************************************
 */
void unmark_for_long_term_reference(frame_store* fs)
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


static void gen_field_ref_ids(VideoParameters *p_Vid, storable_picture *p)
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

/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields
 ************************************************************************
 */
static void dpb_combine_field(VideoParameters *p_Vid, frame_store *fs)
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

void insert_picture_in_dpb(VideoParameters *p_Vid, frame_store* fs, storable_picture* p)
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
    p_Vid->calculate_frame_no(p);
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
  frame_store* fs = p_Dpb->fs[pos];
  frame_store* tmp;
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
        error("Cannot output frame, DPB empty.",150);

    // find smallest POC
    get_smallest_poc(p_Dpb, &poc, &pos);

    if (pos == -1)
        return 0;

  // call the output function

#if (DISABLE_ERC == 0)
    // picture error concealment
    if (p_Vid->conceal_mode != 0) {
        if (p_Dpb->last_output_poc == 0)
            write_lost_ref_after_idr(p_Dpb, pos);
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
    if (p_Vid->conceal_mode == 0) {
        if (p_Dpb->last_output_poc >= poc)
            error ("output POC must be in ascending order", 150);
    }

    p_Dpb->last_output_poc = poc;

    // free frame store and move empty store to end of buffer
    if (!is_used_for_reference(p_Dpb->fs[pos]))
        remove_frame_from_dpb(p_Dpb, pos);
    return 1;
}


void flush_dpb(dpb_t *p_Dpb)
{
    VideoParameters *p_Vid = p_Dpb->p_Vid;

    if (!p_Dpb->init_done)
        return;
#if (DISABLE_ERC == 0)
    if (p_Vid->conceal_mode != 0)
        conceal_non_ref_pics(p_Dpb, 0);
#endif
    // mark all frames unused
    for (int i = 0; i < p_Dpb->used_size; i++) {
#if MVC_EXTENSION_ENABLE
        assert(p_Dpb->fs[i]->view_id == p_Dpb->layer_id);
#endif
        unmark_for_reference(p_Dpb->fs[i]);
    }

    while (remove_unused_frame_from_dpb(p_Dpb));

    // output frames in POC order
    while (p_Dpb->used_size && output_one_frame_from_dpb(p_Dpb));

    p_Dpb->last_output_poc = INT_MIN;
}

#if (MVC_EXTENSION_ENABLE)
void flush_dpbs(dpb_t **p_Dpb_layers, int nLayers)
{
    VideoParameters *p_Vid = p_Dpb_layers[0]->p_Vid;
    dpb_t *p_Dpb;

#if (DISABLE_ERC == 0)
    if (p_Vid->conceal_mode != 0)
        conceal_non_ref_pics(p_Dpb_layers[0], 0);
#endif
    // mark all frames unused
    for (int j = 0; j < nLayers; j++) {
        p_Dpb = p_Dpb_layers[j];
        if (p_Dpb->init_done) {
            for (int i = 0; i < (int)p_Dpb->used_size; i++)
                unmark_for_reference (p_Dpb->fs[i]);
            while (remove_unused_frame_from_dpb(p_Dpb));
        }
    }
    // output frames in POC order
    int used_size = p_Dpb_layers[0]->used_size;
    for (int j = 1; j < nLayers; j++) {
        p_Dpb = p_Dpb_layers[j];
        if (p_Dpb->init_done) {
            if (p_Dpb->used_size && p_Dpb->used_size != used_size)
                assert(!"DPB used_size is not equal!");
            if ((int)p_Dpb->used_size > used_size)
                used_size = (int)p_Dpb->used_size;
        }
    }
    while (used_size) {
        for (int j = 0; j < nLayers; j++) {
            p_Dpb = p_Dpb_layers[j];
            if (p_Dpb->used_size)
                output_one_frame_from_dpb(p_Dpb);
        }
        used_size--;
    }
    for (int j = 0; j < nLayers; j++) {
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
void dpb_split_field(VideoParameters *p_Vid, frame_store *fs)
{
  int i, j, ii, jj, jj4;
  int idiv,jdiv;
  int currentmb;
  int twosz16 = 2 * (fs->frame->size_x >> 4);
  storable_picture *fs_top = NULL, *fs_btm = NULL; 
  storable_picture *frame = fs->frame;

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
      pic_motion_params_old* frm_motion = &frame->motion;
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
void dpb_combine_field_yuv(VideoParameters *p_Vid, frame_store *fs)
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
  storable_picture *picture = NULL;
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
    picture->frame_poc  = currSlice->PicOrderCnt;
    picture->poc        = currSlice->PicOrderCnt;

    store_picture_in_dpb(currSlice->p_Dpb, picture);

    picture=NULL;
    p_Vid->pre_frame_num = UnusedShortTermFrameNum;
    UnusedShortTermFrameNum = (UnusedShortTermFrameNum + 1) % sps->MaxFrameNum;
  }
  currSlice->delta_pic_order_cnt[0] = tmp1;
  currSlice->delta_pic_order_cnt[1] = tmp2;
  currSlice->frame_num = CurrFrameNum;
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

