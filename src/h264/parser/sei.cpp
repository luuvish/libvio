#include <math.h>

#include "global.h"
#include "memalloc.h"
#include "sei.h"
#include "data_partition.h"
#include "slice.h"
#include "dpb.h"
#include "parset.h"


typedef enum {
    SEI_BUFFERING_PERIOD = 0,
    SEI_PIC_TIMING,
    SEI_PAN_SCAN_RECT,
    SEI_FILLER_PAYLOAD,
    SEI_USER_DATA_REGISTERED_ITU_T_T35,
    SEI_USER_DATA_UNREGISTERED,
    SEI_RECOVERY_POINT,
    SEI_DEC_REF_PIC_MARKING_REPETITION,
    SEI_SPARE_PIC,
    SEI_SCENE_INFO,
    SEI_SUB_SEQ_INFO,
    SEI_SUB_SEQ_LAYER_CHARACTERISTICS,
    SEI_SUB_SEQ_CHARACTERISTICS,
    SEI_FULL_FRAME_FREEZE,
    SEI_FULL_FRAME_FREEZE_RELEASE,
    SEI_FULL_FRAME_SNAPSHOT,
    SEI_PROGRESSIVE_REFINEMENT_SEGMENT_START,
    SEI_PROGRESSIVE_REFINEMENT_SEGMENT_END,
    SEI_MOTION_CONSTRAINED_SLICE_GROUP_SET,
    SEI_FILM_GRAIN_CHARACTERISTICS,
    SEI_DEBLOCKING_FILTER_DISPLAY_PREFERENCE,
    SEI_STEREO_VIDEO_INFO,
    SEI_POST_FILTER_HINTS,
    SEI_TONE_MAPPING,
    SEI_SCALABILITY_INFO,
    SEI_SUB_PIC_SCALABLE_LAYER,
    SEI_NON_REQUIRED_LAYER_REP,
    SEI_PRIORITY_LAYER_INFO,
    SEI_LAYERS_NOT_PRESENT,
    SEI_LAYER_DEPENDENCY_CHANGE,
    SEI_SCALABLE_NESTING,
    SEI_BASE_LAYER_TEMPORAL_HRD,
    SEI_QUALITY_LAYER_INTEGRITY_CHECK,
    SEI_REDUNDANT_PIC_PROPERTY,
    SEI_TL0_DEP_REP_INDEX,
    SEI_TL_SWITCHING_POINT,
    SEI_PARALLEL_DECODING_INFO,
    SEI_MVC_SCALABLE_NESTING,
    SEI_VIEW_SCALABILITY_INFO,
    SEI_MULTIVIEW_SCENE_INFO,
    SEI_MULTIVIEW_ACQUISITION_INFO,
    SEI_NON_REQUIRED_VIEW_COMPONENT,
    SEI_VIEW_DEPENDENCY_CHANGE,
    SEI_OPERATION_POINTS_NOT_PRESENT,
    SEI_BASE_VIEW_TEMPORAL_HRD,
    SEI_FRAME_PACKING_ARRANGEMENT,

    SEI_MAX_ELEMENTS  //!< number of maximum syntax elements
} SEI_type;

//! Frame packing arrangement Information
typedef struct {
    unsigned int  frame_packing_arrangement_id;
    Boolean       frame_packing_arrangement_cancel_flag;
    unsigned char frame_packing_arrangement_type;
    Boolean       quincunx_sampling_flag;
    unsigned char content_interpretation_type;
    Boolean       spatial_flipping_flag;
    Boolean       frame0_flipped_flag;
    Boolean       field_views_flag;
    Boolean       current_frame_is_frame0_flag;
    Boolean       frame0_self_contained_flag;
    Boolean       frame1_self_contained_flag;
    unsigned char frame0_grid_position_x;
    unsigned char frame0_grid_position_y;
    unsigned char frame1_grid_position_x;
    unsigned char frame1_grid_position_y;
    unsigned char frame_packing_arrangement_reserved_byte;
    unsigned int  frame_packing_arrangement_repetition_period;
    Boolean       frame_packing_arrangement_extension_flag;
} frame_packing_arrangement_information_struct;

/*!
************************************************************************
*  \brief
*     Interpret the spare picture SEI message
*  \param payload
*     a pointer that point to the sei payload
*  \param size
*     the size of the sei message
*  \param p_Vid
*     the image pointer
*
************************************************************************
*/
void interpret_spare_pic( byte* payload, int size, VideoParameters *p_Vid )
{
  int i,x,y;
  data_partition_t* buf;
  int bit0, bit1, bitc, no_bit0;
  int target_frame_num = 0;
  int num_spare_pics;
  int delta_spare_frame_num, CandidateSpareFrameNum, SpareFrameNum = 0;
  int ref_area_indicator;

  sps_t *sps = p_Vid->active_sps;

  int m, n, left, right, top, bottom,directx, directy;
  byte ***map;

  assert( payload!=NULL);
  assert( p_Vid!=NULL);

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  target_frame_num = buf->ue("SEI: target_frame_num");

  num_spare_pics = 1 + buf->ue("SEI: num_spare_pics_minus1");

  get_mem3D(&map, num_spare_pics, sps->FrameHeightInMbs, sps->PicWidthInMbs);

  for (i=0; i<num_spare_pics; i++)
  {
    if (i==0)
    {
      CandidateSpareFrameNum = target_frame_num - 1;
      if ( CandidateSpareFrameNum < 0 ) CandidateSpareFrameNum = MAX_FN - 1;
    }
    else
      CandidateSpareFrameNum = SpareFrameNum;

    delta_spare_frame_num = buf->ue("SEI: delta_spare_frame_num");

    SpareFrameNum = CandidateSpareFrameNum - delta_spare_frame_num;
    if( SpareFrameNum < 0 )
      SpareFrameNum = MAX_FN + SpareFrameNum;

    ref_area_indicator = buf->ue("SEI: ref_area_indicator");

    switch ( ref_area_indicator )
    {
    case 0:   // The whole frame can serve as spare picture
      for (y=0; y<sps->FrameHeightInMbs; y++)
        for (x=0; x<sps->PicWidthInMbs; x++)
          map[i][y][x] = 0;
      break;
    case 1:   // The map is not compressed
      for (y=0; y<sps->FrameHeightInMbs; y++)
        for (x=0; x<sps->PicWidthInMbs; x++)
        {
          map[i][y][x] = (byte) buf->u(1, "SEI: ref_mb_indicator");
        }
      break;
    case 2:   // The map is compressed
              //!KS: could not check this function, description is unclear (as stated in Ed. Note)
      bit0 = 0;
      bit1 = 1;
      bitc = bit0;
      no_bit0 = -1;

      x = ( sps->PicWidthInMbs - 1 ) / 2;
      y = ( sps->FrameHeightInMbs - 1 ) / 2;
      left = right = x;
      top = bottom = y;
      directx = 0;
      directy = 1;

      for (m=0; m<sps->FrameHeightInMbs; m++)
        for (n=0; n<sps->PicWidthInMbs; n++)
        {

          if (no_bit0<0)
          {
            no_bit0 = buf->ue("SEI: zero_run_length");
          }
          if (no_bit0>0) 
            map[i][y][x] = (byte) bit0;
          else 
            map[i][y][x] = (byte) bit1;
          no_bit0--;

          // go to the next mb:
          if ( directx == -1 && directy == 0 )
          {
            if (x > left) x--;
            else if (x == 0)
            {
              y = bottom + 1;
              bottom++;
              directx = 1;
              directy = 0;
            }
            else if (x == left)
            {
              x--;
              left--;
              directx = 0;
              directy = 1;
            }
          }
          else if ( directx == 1 && directy == 0 )
          {
            if (x < right) x++;
            else if (x == sps->PicWidthInMbs - 1)
            {
              y = top - 1;
              top--;
              directx = -1;
              directy = 0;
            }
            else if (x == right)
            {
              x++;
              right++;
              directx = 0;
              directy = -1;
            }
          }
          else if ( directx == 0 && directy == -1 )
          {
            if ( y > top) y--;
            else if (y == 0)
            {
              x = left - 1;
              left--;
              directx = 0;
              directy = 1;
            }
            else if (y == top)
            {
              y--;
              top--;
              directx = -1;
              directy = 0;
            }
          }
          else if ( directx == 0 && directy == 1 )
          {
            if (y < bottom) y++;
            else if (y == sps->FrameHeightInMbs - 1)
            {
              x = right+1;
              right++;
              directx = 0;
              directy = -1;
            }
            else if (y == bottom)
            {
              y++;
              bottom++;
              directx = 1;
              directy = 0;
            }
          }


        }
      break;
    default:
      printf( "Wrong ref_area_indicator %d!\n", ref_area_indicator );
      exit(0);
      break;
    }

  } // end of num_spare_pics

  free_mem3D( map );

  free(buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Sub-sequence information SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_subsequence_info( byte* payload, int size, VideoParameters *p_Vid )
{
  data_partition_t* buf;
  int sub_seq_layer_num, sub_seq_id, first_ref_pic_flag, leading_non_ref_pic_flag, last_pic_flag,
    sub_seq_frame_num_flag, sub_seq_frame_num;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  sub_seq_layer_num        = buf->ue("SEI: sub_seq_layer_num");
  sub_seq_id               = buf->ue("SEI: sub_seq_id");
  first_ref_pic_flag       = buf->u(1, "SEI: first_ref_pic_flag");
  leading_non_ref_pic_flag = buf->u(1, "SEI: leading_non_ref_pic_flag");
  last_pic_flag            = buf->u(1, "SEI: last_pic_flag");
  sub_seq_frame_num_flag   = buf->u(1, "SEI: sub_seq_frame_num_flag");
  if (sub_seq_frame_num_flag)
  {
    sub_seq_frame_num      = buf->ue("SEI: sub_seq_frame_num");
  }

  free(buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the Sub-sequence layer characteristics SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_subsequence_layer_characteristics_info( byte* payload, int size, VideoParameters *p_Vid )
{
  data_partition_t* buf;
  long num_sub_layers, accurate_statistics_flag, average_bit_rate, average_frame_rate;
  int i;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  num_sub_layers = 1 + buf->ue("SEI: num_sub_layers_minus1");

  for (i=0; i<num_sub_layers; i++)
  {
    accurate_statistics_flag = buf->u(1, "SEI: accurate_statistics_flag");
    average_bit_rate         = buf->u(16, "SEI: average_bit_rate");
    average_frame_rate       = buf->u(16, "SEI: average_frame_rate");
  }
  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Sub-sequence characteristics SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_subsequence_characteristics_info( byte* payload, int size, VideoParameters *p_Vid )
{
  data_partition_t* buf;
  int i;
  int sub_seq_layer_num, sub_seq_id, duration_flag, average_rate_flag, accurate_statistics_flag;
  unsigned long sub_seq_duration, average_bit_rate, average_frame_rate;
  int num_referenced_subseqs, ref_sub_seq_layer_num, ref_sub_seq_id, ref_sub_seq_direction;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  sub_seq_layer_num = buf->ue("SEI: sub_seq_layer_num");
  sub_seq_id        = buf->ue("SEI: sub_seq_id");
  duration_flag     = buf->u(1, "SEI: duration_flag");

  if ( duration_flag )
  {
    sub_seq_duration = buf->u(32, "SEI: duration_flag");
  }

  average_rate_flag = buf->u(1, "SEI: average_rate_flag");

  if ( average_rate_flag )
  {
    accurate_statistics_flag = buf->u(1, "SEI: accurate_statistics_flag");
    average_bit_rate         = buf->u(16, "SEI: average_bit_rate");
    average_frame_rate       = buf->u(16, "SEI: average_frame_rate");
  }

  num_referenced_subseqs  = buf->ue("SEI: num_referenced_subseqs");

  for (i=0; i<num_referenced_subseqs; i++)
  {
    ref_sub_seq_layer_num  = buf->ue("SEI: ref_sub_seq_layer_num");
    ref_sub_seq_id         = buf->ue("SEI: ref_sub_seq_id");
    ref_sub_seq_direction  = buf->u(1, "SEI: ref_sub_seq_direction");
  }

  free( buf );
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Scene information SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_scene_information( byte* payload, int size, VideoParameters *p_Vid )
{
  data_partition_t* buf;
  int scene_id, scene_transition_type, second_scene_id;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  scene_id              = buf->ue("SEI: scene_id");
  scene_transition_type = buf->ue("SEI: scene_transition_type");
  if ( scene_transition_type > 3 )
  {
    second_scene_id     = buf->ue("SEI: scene_transition_type");
  }

  free( buf );
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Filler payload SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_filler_payload_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int payload_cnt = 0;

  while (payload_cnt<size)
  {
    if (payload[payload_cnt] == 0xFF)
    {
       payload_cnt++;
    }
  }
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the User data unregistered SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_user_data_unregistered_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int offset = 0;
  byte payload_byte;

  assert (size>=16);

  for (offset = 0; offset < 16; offset++)
  {
  }

  while (offset < size)
  {
    payload_byte = payload[offset];
    offset ++;
  }
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the User data registered by ITU-T T.35 SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_user_data_registered_itu_t_t35_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int offset = 0;
  byte itu_t_t35_country_code, itu_t_t35_country_code_extension_byte, payload_byte;

  itu_t_t35_country_code = payload[offset];
  offset++;
  if(itu_t_t35_country_code == 0xFF)
  {
    itu_t_t35_country_code_extension_byte = payload[offset];
    offset++;
  }
  while (offset < size)
  {
    payload_byte = payload[offset];
    offset ++;
  }
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Pan scan rectangle SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_pan_scan_rect_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int pan_scan_rect_cancel_flag;
  int pan_scan_cnt_minus1, i;
  int pan_scan_rect_repetition_period;
  int pan_scan_rect_id, pan_scan_rect_left_offset, pan_scan_rect_right_offset;
  int pan_scan_rect_top_offset, pan_scan_rect_bottom_offset;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  pan_scan_rect_id = buf->ue("SEI: pan_scan_rect_id");

  pan_scan_rect_cancel_flag = buf->u(1, "SEI: pan_scan_rect_cancel_flag");
  if (!pan_scan_rect_cancel_flag) 
  {
    pan_scan_cnt_minus1 = buf->ue("SEI: pan_scan_cnt_minus1");
    for (i = 0; i <= pan_scan_cnt_minus1; i++) 
    {
      pan_scan_rect_left_offset   = buf->se("SEI: pan_scan_rect_left_offset");
      pan_scan_rect_right_offset  = buf->se("SEI: pan_scan_rect_right_offset");
      pan_scan_rect_top_offset    = buf->se("SEI: pan_scan_rect_top_offset");
      pan_scan_rect_bottom_offset = buf->se("SEI: pan_scan_rect_bottom_offset");
    }
    pan_scan_rect_repetition_period = buf->ue("SEI: pan_scan_rect_repetition_period");
  }

  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Random access point SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_recovery_point_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int recovery_frame_cnt, exact_match_flag, broken_link_flag, changing_slice_group_idc;


  data_partition_t* buf;


  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  recovery_frame_cnt       = buf->ue("SEI: recovery_frame_cnt");
  exact_match_flag         = buf->u(1, "SEI: exact_match_flag");
  broken_link_flag         = buf->u(1, "SEI: broken_link_flag");
  changing_slice_group_idc = buf->u(2, "SEI: changing_slice_group_idc");

  p_Vid->recovery_point = 1;
  p_Vid->recovery_frame_cnt = recovery_frame_cnt;

  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Decoded Picture Buffer Management Repetition SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_dec_ref_pic_marking_repetition_info( byte* payload, int size, VideoParameters *p_Vid, slice_t *pSlice )
{
  int original_idr_flag, original_frame_num;
  int original_field_pic_flag, original_bottom_field_flag;

  DecRefPicMarking_t *tmp_drpm;
  DecRefPicMarking_t *old_drpm;
  int old_idr_flag, old_no_output_of_prior_pics_flag, old_long_term_reference_flag , old_adaptive_ref_pic_buffering_flag;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  original_idr_flag     = buf->u(1, "SEI: original_idr_flag");
  original_frame_num    = buf->ue("SEI: original_frame_num");

  if ( !p_Vid->active_sps->frame_mbs_only_flag )
  {
    original_field_pic_flag = buf->u(1, "SEI: original_field_pic_flag");
    if ( original_field_pic_flag )
    {
      original_bottom_field_flag = buf->u(1, "SEI: original_bottom_field_flag");
    }
  }

  // we need to save everything that is probably overwritten in dec_ref_pic_marking()
  old_drpm = pSlice->dec_ref_pic_marking_buffer;
  old_idr_flag = pSlice->idr_flag; //p_Vid->idr_flag;

  old_no_output_of_prior_pics_flag = pSlice->no_output_of_prior_pics_flag; //p_Vid->no_output_of_prior_pics_flag;
  old_long_term_reference_flag = pSlice->long_term_reference_flag;
  old_adaptive_ref_pic_buffering_flag = pSlice->adaptive_ref_pic_marking_mode_flag;

  // set new initial values
  //p_Vid->idr_flag = original_idr_flag;
  pSlice->idr_flag = original_idr_flag;
  pSlice->dec_ref_pic_marking_buffer = NULL;

  dec_ref_pic_marking(p_Vid, buf, pSlice);

  while (pSlice->dec_ref_pic_marking_buffer)
  {
    tmp_drpm=pSlice->dec_ref_pic_marking_buffer;

    pSlice->dec_ref_pic_marking_buffer=tmp_drpm->Next;
    free (tmp_drpm);
  }

  // restore old values in p_Vid
  pSlice->dec_ref_pic_marking_buffer = old_drpm;
  pSlice->idr_flag = old_idr_flag;
  pSlice->no_output_of_prior_pics_flag = old_no_output_of_prior_pics_flag;
  p_Vid->no_output_of_prior_pics_flag = pSlice->no_output_of_prior_pics_flag;
  pSlice->long_term_reference_flag = old_long_term_reference_flag;
  pSlice->adaptive_ref_pic_marking_mode_flag = old_adaptive_ref_pic_buffering_flag;

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the Full-frame freeze SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_full_frame_freeze_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int full_frame_freeze_repetition_period;
  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  full_frame_freeze_repetition_period  = buf->ue("SEI: full_frame_freeze_repetition_period");

  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Full-frame freeze release SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_full_frame_freeze_release_info( byte* payload, int size, VideoParameters *p_Vid )
{
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the Full-frame snapshot SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_full_frame_snapshot_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int snapshot_id;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  snapshot_id = buf->ue("SEI: snapshot_id");

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the Progressive refinement segment start SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_progressive_refinement_start_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int progressive_refinement_id, num_refinement_steps_minus1;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  progressive_refinement_id   = buf->ue("SEI: progressive_refinement_id");
  num_refinement_steps_minus1 = buf->ue("SEI: num_refinement_steps_minus1");

  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Progressive refinement segment end SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_progressive_refinement_end_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int progressive_refinement_id;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  progressive_refinement_id   = buf->ue("SEI: progressive_refinement_id");

  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Motion-constrained slice group set SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_motion_constrained_slice_group_set_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int num_slice_groups_minus1, slice_group_id, exact_match_flag, pan_scan_rect_flag, pan_scan_rect_id;
  int i;
  int sliceGroupSize;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  num_slice_groups_minus1   = buf->ue("SEI: num_slice_groups_minus1");
  sliceGroupSize = ceil(log2(num_slice_groups_minus1 + 1));

  for (i=0; i<=num_slice_groups_minus1;i++)
  {

    slice_group_id   = buf->u(sliceGroupSize, "SEI: slice_group_id");
  }

  exact_match_flag   = buf->u(1, "SEI: exact_match_flag");
  pan_scan_rect_flag = buf->u(1, "SEI: pan_scan_rect_flag");

  if (pan_scan_rect_flag)
  {
    pan_scan_rect_id = buf->ue("SEI: pan_scan_rect_id");
  }

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the film grain characteristics SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_film_grain_characteristics_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int film_grain_characteristics_cancel_flag;
  int model_id, separate_colour_description_present_flag;
  int film_grain_bit_depth_luma_minus8, film_grain_bit_depth_chroma_minus8, film_grain_full_range_flag, film_grain_colour_primaries, film_grain_transfer_characteristics, film_grain_matrix_coefficients;
  int blending_mode_id, log2_scale_factor, comp_model_present_flag[3];
  int num_intensity_intervals_minus1, num_model_values_minus1;
  int intensity_interval_lower_bound, intensity_interval_upper_bound;
  int comp_model_value;
  int film_grain_characteristics_repetition_period;

  int c, i, j;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  film_grain_characteristics_cancel_flag = buf->u(1, "SEI: film_grain_characteristics_cancel_flag");
  if(!film_grain_characteristics_cancel_flag)
  {

    model_id                                    = buf->u(2, "SEI: model_id");
    separate_colour_description_present_flag    = buf->u(1, "SEI: separate_colour_description_present_flag");
    if (separate_colour_description_present_flag)
    {
      film_grain_bit_depth_luma_minus8          = buf->u(3, "SEI: film_grain_bit_depth_luma_minus8");
      film_grain_bit_depth_chroma_minus8        = buf->u(3, "SEI: film_grain_bit_depth_chroma_minus8");
      film_grain_full_range_flag                = buf->u(1, "SEI: film_grain_full_range_flag");
      film_grain_colour_primaries               = buf->u(8, "SEI: film_grain_colour_primaries");
      film_grain_transfer_characteristics       = buf->u(8, "SEI: film_grain_transfer_characteristics");
      film_grain_matrix_coefficients            = buf->u(8, "SEI: film_grain_matrix_coefficients");
    }
    blending_mode_id                            = buf->u(2, "SEI: blending_mode_id");
    log2_scale_factor                           = buf->u(4, "SEI: log2_scale_factor");
    for (c = 0; c < 3; c ++)
    {
      comp_model_present_flag[c]                = buf->u(1, "SEI: comp_model_present_flag");
    }
    for (c = 0; c < 3; c ++)
      if (comp_model_present_flag[c])
      {
        num_intensity_intervals_minus1          = buf->u(8, "SEI: num_intensity_intervals_minus1");
        num_model_values_minus1                 = buf->u(3, "SEI: num_model_values_minus1");
        for (i = 0; i <= num_intensity_intervals_minus1; i ++)
        {
          intensity_interval_lower_bound        = buf->u(8, "SEI: intensity_interval_lower_bound");
          intensity_interval_upper_bound        = buf->u(8, "SEI: intensity_interval_upper_bound");
          for (j = 0; j <= num_model_values_minus1; j++)
          {
            comp_model_value                    = buf->se("SEI: comp_model_value");
          }
        }
      }
    film_grain_characteristics_repetition_period = buf->ue("SEI: film_grain_characteristics_repetition_period");
  }

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the deblocking filter display preference SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_deblocking_filter_display_preference_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int deblocking_display_preference_cancel_flag;
  int display_prior_to_deblocking_preferred_flag, dec_frame_buffering_constraint_flag, deblocking_display_preference_repetition_period;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  deblocking_display_preference_cancel_flag             = buf->u(1, "SEI: deblocking_display_preference_cancel_flag");
  if(!deblocking_display_preference_cancel_flag)
  {
    display_prior_to_deblocking_preferred_flag            = buf->u(1, "SEI: display_prior_to_deblocking_preferred_flag");
    dec_frame_buffering_constraint_flag                   = buf->u(1, "SEI: dec_frame_buffering_constraint_flag");
    deblocking_display_preference_repetition_period       = buf->ue("SEI: deblocking_display_preference_repetition_period");
  }

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the stereo video info SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_stereo_video_info_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int field_views_flags;
  int top_field_is_left_view_flag, current_frame_is_left_view_flag, next_frame_is_second_view_flag;
  int left_view_self_contained_flag;
  int right_view_self_contained_flag;

  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  field_views_flags = buf->u(1, "SEI: field_views_flags");
  if (field_views_flags)
  {
    top_field_is_left_view_flag         = buf->u(1, "SEI: top_field_is_left_view_flag");
  }
  else
  {
    current_frame_is_left_view_flag     = buf->u(1, "SEI: current_frame_is_left_view_flag");
    next_frame_is_second_view_flag      = buf->u(1, "SEI: next_frame_is_second_view_flag");
  }

  left_view_self_contained_flag         = buf->u(1, "SEI: left_view_self_contained_flag");
  right_view_self_contained_flag        = buf->u(1, "SEI: right_view_self_contained_flag");

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the Reserved SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_reserved_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int offset = 0;
  byte payload_byte;

  while (offset < size)
  {
    payload_byte = payload[offset];
    offset ++;
  }
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Buffering period SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_buffering_period_info( byte* payload, int size, VideoParameters *p_Vid )
{
  int seq_parameter_set_id, initial_cpb_removal_delay, initial_cpb_removal_delay_offset;
  unsigned int k;

  data_partition_t* buf;
  sps_t *sps;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  seq_parameter_set_id   = buf->ue("SEI: seq_parameter_set_id");

  sps = &p_Vid->SeqParSet[seq_parameter_set_id];

  activate_sps(p_Vid, sps);

  // Note: NalHrdBpPresentFlag and CpbDpbDelaysPresentFlag can also be set "by some means not specified in this Recommendation | International Standard"
  if (sps->vui_parameters_present_flag)
  {

    if (sps->vui_parameters.nal_hrd_parameters_present_flag)
    {
      for (k=0; k<sps->vui_parameters.nal_hrd_parameters.cpb_cnt_minus1+1; k++)
      {
        initial_cpb_removal_delay        = buf->u(sps->vui_parameters.nal_hrd_parameters.initial_cpb_removal_delay_length_minus1+1, "SEI: initial_cpb_removal_delay");
        initial_cpb_removal_delay_offset = buf->u(sps->vui_parameters.nal_hrd_parameters.initial_cpb_removal_delay_length_minus1+1, "SEI: initial_cpb_removal_delay_offset");

      }
    }

    if (sps->vui_parameters.vcl_hrd_parameters_present_flag)
    {
      for (k=0; k<sps->vui_parameters.vcl_hrd_parameters.cpb_cnt_minus1+1; k++)
      {
        initial_cpb_removal_delay        = buf->u(sps->vui_parameters.vcl_hrd_parameters.initial_cpb_removal_delay_length_minus1+1, "SEI: initial_cpb_removal_delay");
        initial_cpb_removal_delay_offset = buf->u(sps->vui_parameters.vcl_hrd_parameters.initial_cpb_removal_delay_length_minus1+1, "SEI: initial_cpb_removal_delay_offset");
      }
    }
  }

  free (buf);
}


/*!
 ************************************************************************
 *  \brief
 *     Interpret the Picture timing SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void interpret_picture_timing_info( byte* payload, int size, VideoParameters *p_Vid )
{
  sps_t *active_sps = p_Vid->active_sps;

  int cpb_removal_delay, dpb_output_delay, pic_struct_present_flag, pic_struct;
  int clock_timestamp_flag;
  int ct_type, nuit_field_based_flag, counting_type, full_timestamp_flag, discontinuity_flag, cnt_dropped_flag, nframes;
  int seconds_value, minutes_value, hours_value, seconds_flag, minutes_flag, hours_flag, time_offset;
  int NumClockTs = 0;
  int i;

  int cpb_removal_len = 24;
  int dpb_output_len  = 24;

  Boolean CpbDpbDelaysPresentFlag;

  data_partition_t* buf;

  if (NULL==active_sps)
  {
    fprintf (stderr, "Warning: no active SPS, timing SEI cannot be parsed\n");
    return;
  }

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  // CpbDpbDelaysPresentFlag can also be set "by some means not specified in this Recommendation | International Standard"
  CpbDpbDelaysPresentFlag =  (Boolean) (active_sps->vui_parameters_present_flag
                              && (   (active_sps->vui_parameters.nal_hrd_parameters_present_flag != 0)
                                   ||(active_sps->vui_parameters.vcl_hrd_parameters_present_flag != 0)));

  if (CpbDpbDelaysPresentFlag )
  {
    if (active_sps->vui_parameters_present_flag)
    {
      if (active_sps->vui_parameters.nal_hrd_parameters_present_flag)
      {
        cpb_removal_len = active_sps->vui_parameters.nal_hrd_parameters.cpb_removal_delay_length_minus1 + 1;
        dpb_output_len  = active_sps->vui_parameters.nal_hrd_parameters.dpb_output_delay_length_minus1  + 1;
      }
      else if (active_sps->vui_parameters.vcl_hrd_parameters_present_flag)
      {
        cpb_removal_len = active_sps->vui_parameters.vcl_hrd_parameters.cpb_removal_delay_length_minus1 + 1;
        dpb_output_len  = active_sps->vui_parameters.vcl_hrd_parameters.dpb_output_delay_length_minus1  + 1;
      }
    }

    if ((active_sps->vui_parameters.nal_hrd_parameters_present_flag)||
      (active_sps->vui_parameters.vcl_hrd_parameters_present_flag))
    {
      cpb_removal_delay = buf->u(cpb_removal_len, "SEI: cpb_removal_delay");
      dpb_output_delay  = buf->u(dpb_output_len,  "SEI: dpb_output_delay");
    }
  }

  if (!active_sps->vui_parameters_present_flag)
  {
    pic_struct_present_flag = 0;
  }
  else
  {
    pic_struct_present_flag  =  active_sps->vui_parameters.pic_struct_present_flag;
  }

  if (pic_struct_present_flag)
  {
    pic_struct = buf->u(4, "SEI: pic_struct");
    switch (pic_struct)
    {
    case 0:
    case 1:
    case 2:
      NumClockTs = 1;
      break;
    case 3:
    case 4:
    case 7:
      NumClockTs = 2;
      break;
    case 5:
    case 6:
    case 8:
      NumClockTs = 3;
      break;
    default:
      error("reserved pic_struct used (can't determine NumClockTs)", 500);
    }
    for (i=0; i<NumClockTs; i++)
    {
      clock_timestamp_flag = buf->u(1, "SEI: clock_timestamp_flag");
      if (clock_timestamp_flag)
      {
        ct_type               = buf->u(2, "SEI: ct_type");
        nuit_field_based_flag = buf->u(1, "SEI: nuit_field_based_flag");
        counting_type         = buf->u(5, "SEI: counting_type");
        full_timestamp_flag   = buf->u(1, "SEI: full_timestamp_flag");
        discontinuity_flag    = buf->u(1, "SEI: discontinuity_flag");
        cnt_dropped_flag      = buf->u(1, "SEI: cnt_dropped_flag");
        nframes               = buf->u(8, "SEI: nframes");

        if (full_timestamp_flag)
        {
          seconds_value         = buf->u(6, "SEI: seconds_value");
          minutes_value         = buf->u(6, "SEI: minutes_value");
          hours_value           = buf->u(5, "SEI: hours_value");
        }
        else
        {
          seconds_flag          = buf->u(1, "SEI: seconds_flag");
          if (seconds_flag)
          {
            seconds_value         = buf->u(6, "SEI: seconds_value");
            minutes_flag          = buf->u(1, "SEI: minutes_flag");
            if(minutes_flag)
            {
              minutes_value         = buf->u(6, "SEI: minutes_value");
              hours_flag            = buf->u(1, "SEI: hours_flag");
              if(hours_flag)
              {
                hours_value           = buf->u(5, "SEI: hours_value");
              }
            }
          }
        }
        {
          int time_offset_length;
          if (active_sps->vui_parameters.vcl_hrd_parameters_present_flag)
            time_offset_length = active_sps->vui_parameters.vcl_hrd_parameters.time_offset_length;
          else if (active_sps->vui_parameters.nal_hrd_parameters_present_flag)
            time_offset_length = active_sps->vui_parameters.nal_hrd_parameters.time_offset_length;
          else
            time_offset_length = 24;
          if (time_offset_length)
            time_offset = buf->i(time_offset_length, "SEI: time_offset");
          else
            time_offset = 0;
        }
      }
    }
  }

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the Frame Packing Arrangement SEI message
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 ************************************************************************
 */
void interpret_frame_packing_arrangement_info( byte* payload, int size, VideoParameters *p_Vid )
{
  frame_packing_arrangement_information_struct seiFramePackingArrangement;
  data_partition_t* buf;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  seiFramePackingArrangement.frame_packing_arrangement_id = (unsigned int) buf->ue( "SEI: frame_packing_arrangement_id");
  seiFramePackingArrangement.frame_packing_arrangement_cancel_flag = buf->u(1, "SEI: frame_packing_arrangement_cancel_flag");
  if ( seiFramePackingArrangement.frame_packing_arrangement_cancel_flag == FALSE )
  {
    seiFramePackingArrangement.frame_packing_arrangement_type = (unsigned char) buf->u(7, "SEI: frame_packing_arrangement_type");
    seiFramePackingArrangement.quincunx_sampling_flag         = buf->u(1, "SEI: quincunx_sampling_flag");
    seiFramePackingArrangement.content_interpretation_type    = (unsigned char) buf->u(6, "SEI: content_interpretation_type");
    seiFramePackingArrangement.spatial_flipping_flag          = buf->u(1, "SEI: spatial_flipping_flag");
    seiFramePackingArrangement.frame0_flipped_flag            = buf->u(1, "SEI: frame0_flipped_flag");
    seiFramePackingArrangement.field_views_flag               = buf->u(1, "SEI: field_views_flag");
    seiFramePackingArrangement.current_frame_is_frame0_flag   = buf->u(1, "SEI: current_frame_is_frame0_flag");
    seiFramePackingArrangement.frame0_self_contained_flag     = buf->u(1, "SEI: frame0_self_contained_flag");
    seiFramePackingArrangement.frame1_self_contained_flag     = buf->u(1, "SEI: frame1_self_contained_flag");
    if ( seiFramePackingArrangement.quincunx_sampling_flag == FALSE && seiFramePackingArrangement.frame_packing_arrangement_type != 5 )
    {
      seiFramePackingArrangement.frame0_grid_position_x = (unsigned char)buf->u(4, "SEI: frame0_grid_position_x");
      seiFramePackingArrangement.frame0_grid_position_y = (unsigned char)buf->u(4, "SEI: frame0_grid_position_y");
      seiFramePackingArrangement.frame1_grid_position_x = (unsigned char)buf->u(4, "SEI: frame1_grid_position_x");
      seiFramePackingArrangement.frame1_grid_position_y = (unsigned char)buf->u(4, "SEI: frame1_grid_position_y");
    }
    seiFramePackingArrangement.frame_packing_arrangement_reserved_byte = (unsigned char) buf->u(8, "SEI: frame_packing_arrangement_reserved_byte");
    seiFramePackingArrangement.frame_packing_arrangement_repetition_period = (unsigned int) buf->ue("SEI: frame_packing_arrangement_repetition_period");
  }
  seiFramePackingArrangement.frame_packing_arrangement_extension_flag = buf->u(1, "SEI: frame_packing_arrangement_extension_flag");

  free (buf);
}

/*!
 ************************************************************************
 *  \brief
 *     Interpret the HDR tone-mapping SEI message (JVT-T060)
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
typedef struct
{
  unsigned int  tone_map_id;
  unsigned char tone_map_cancel_flag;
  unsigned int  tone_map_repetition_period;
  unsigned char coded_data_bit_depth;
  unsigned char sei_bit_depth;
  unsigned int  model_id;
  // variables for model 0
  int  min_value;
  int  max_value;
  // variables for model 1
  int  sigmoid_midpoint;
  int  sigmoid_width;
  // variables for model 2
  int start_of_coded_interval[1<<MAX_SEI_BIT_DEPTH];
  // variables for model 3
  int num_pivots;
  int coded_pivot_value[MAX_NUM_PIVOTS];
  int sei_pivot_value[MAX_NUM_PIVOTS];
} tone_mapping_struct_tmp;

void interpret_tone_mapping( byte* payload, int size, VideoParameters *p_Vid )
{
  tone_mapping_struct_tmp seiToneMappingTmp;
  data_partition_t* buf;
  int i = 0, max_coded_num, max_output_num;

  memset (&seiToneMappingTmp, 0, sizeof (tone_mapping_struct_tmp));

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  seiToneMappingTmp.tone_map_id = buf->ue("SEI: tone_map_id");
  seiToneMappingTmp.tone_map_cancel_flag = (unsigned char) buf->u(1, "SEI: tone_map_cancel_flag");

  if (!seiToneMappingTmp.tone_map_cancel_flag) 
  {
    seiToneMappingTmp.tone_map_repetition_period  = buf->ue("SEI: tone_map_repetition_period");
    seiToneMappingTmp.coded_data_bit_depth        = (unsigned char) buf->u(8, "SEI: coded_data_bit_depth");
    seiToneMappingTmp.sei_bit_depth               = (unsigned char) buf->u(8, "SEI: sei_bit_depth");
    seiToneMappingTmp.model_id                    = buf->ue("SEI: model_id");

    max_coded_num  = 1<<seiToneMappingTmp.coded_data_bit_depth;
    max_output_num = 1<<seiToneMappingTmp.sei_bit_depth;

    if (seiToneMappingTmp.model_id == 0) 
    { // linear mapping with clipping
      seiToneMappingTmp.min_value   = buf->u(32, "SEI: min_value");
      seiToneMappingTmp.max_value   = buf->u(32, "SEI: min_value");
    }
    else if (seiToneMappingTmp.model_id == 1) 
    { // sigmoidal mapping
      seiToneMappingTmp.sigmoid_midpoint = buf->u(32, "SEI: sigmoid_midpoint");
      seiToneMappingTmp.sigmoid_width    = buf->u(32, "SEI: sigmoid_width");
    }
    else if (seiToneMappingTmp.model_id == 2) 
    { // user defined table mapping
      for (i=0; i<max_output_num; i++) 
      {
        seiToneMappingTmp.start_of_coded_interval[i] = buf->u((((seiToneMappingTmp.coded_data_bit_depth+7)>>3)<<3), "SEI: start_of_coded_interval");
      }
    }
    else if (seiToneMappingTmp.model_id == 3) 
    {  // piece-wise linear mapping
      seiToneMappingTmp.num_pivots = buf->u(16, "SEI: num_pivots");
      seiToneMappingTmp.coded_pivot_value[0] = 0;
      seiToneMappingTmp.sei_pivot_value[0] = 0;
      seiToneMappingTmp.coded_pivot_value[seiToneMappingTmp.num_pivots+1] = max_coded_num-1;
      seiToneMappingTmp.sei_pivot_value[seiToneMappingTmp.num_pivots+1] = max_output_num-1;

      for (i=1; i < seiToneMappingTmp.num_pivots+1; i++) 
      {
        seiToneMappingTmp.coded_pivot_value[i] = buf->u( (((seiToneMappingTmp.coded_data_bit_depth+7)>>3)<<3), "SEI: coded_pivot_value");
        seiToneMappingTmp.sei_pivot_value[i] = buf->u( (((seiToneMappingTmp.sei_bit_depth+7)>>3)<<3), "SEI: sei_pivot_value");
      }
    }

#if (ENABLE_OUTPUT_TONEMAPPING)
    // Currently, only when the map_id == 0, the tone-mapping is actually applied.
    if (seiToneMappingTmp.tone_map_id== 0) 
    {
      int j;
      p_Vid->seiToneMapping->seiHasTone_mapping = TRUE;
      p_Vid->seiToneMapping->tone_map_repetition_period = seiToneMappingTmp.tone_map_repetition_period;
      p_Vid->seiToneMapping->coded_data_bit_depth = seiToneMappingTmp.coded_data_bit_depth;
      p_Vid->seiToneMapping->sei_bit_depth = seiToneMappingTmp.sei_bit_depth;
      p_Vid->seiToneMapping->model_id = seiToneMappingTmp.model_id;
      p_Vid->seiToneMapping->count = 0;

      // generate the look up table of tone mapping
      switch(seiToneMappingTmp.model_id)
      {
      case 0:            // linear mapping with clipping
        for (i=0; i<=seiToneMappingTmp.min_value; i++)
          p_Vid->seiToneMapping->lut[i] = 0;

        for (i=seiToneMappingTmp.min_value+1; i < seiToneMappingTmp.max_value; i++)
          p_Vid->seiToneMapping->lut[i] = (imgpel) ((i-seiToneMappingTmp.min_value) * (max_output_num-1)/(seiToneMappingTmp.max_value- seiToneMappingTmp.min_value));

        for (i=seiToneMappingTmp.max_value; i<max_coded_num; i++)
          p_Vid->seiToneMapping->lut[i] = (imgpel) (max_output_num - 1);
        break;
      case 1: // sigmoid mapping

        for (i=0; i < max_coded_num; i++) 
        {
          double tmp = 1.0 + exp( -6*(double)(i-seiToneMappingTmp.sigmoid_midpoint)/seiToneMappingTmp.sigmoid_width);
          p_Vid->seiToneMapping->lut[i] = (imgpel)( (double)(max_output_num-1)/ tmp + 0.5);
        }
        break;
      case 2: // user defined table
        if (0 < max_output_num-1)
        {
          for (j=0; j<max_output_num-1; j++) 
          {
            for (i=seiToneMappingTmp.start_of_coded_interval[j]; i<seiToneMappingTmp.start_of_coded_interval[j+1]; i++) 
            {
              p_Vid->seiToneMapping->lut[i] = (imgpel) j;
            }
          }
          p_Vid->seiToneMapping->lut[i] = (imgpel) (max_output_num - 1);
        }
        break;
      case 3: // piecewise linear mapping
        for (j=0; j<seiToneMappingTmp.num_pivots+1; j++) 
        {
          double slope = (double)(seiToneMappingTmp.sei_pivot_value[j+1] - seiToneMappingTmp.sei_pivot_value[j])/(seiToneMappingTmp.coded_pivot_value[j+1]-seiToneMappingTmp.coded_pivot_value[j]);
          for (i=seiToneMappingTmp.coded_pivot_value[j]; i <= seiToneMappingTmp.coded_pivot_value[j+1]; i++) 
          {
            p_Vid->seiToneMapping->lut[i] = (imgpel) (seiToneMappingTmp.sei_pivot_value[j] + (int)(( (i - seiToneMappingTmp.coded_pivot_value[j]) * slope)));
          }
        }
        break;

      default:
        break;
      } // end switch
    }
#endif
  } // end !tone_map_cancel_flag
  free (buf);
}

#if (ENABLE_OUTPUT_TONEMAPPING)
// tone map using the look-up-table generated according to SEI tone mapping message
void tone_map (imgpel** imgX, imgpel* lut, int size_x, int size_y)
{
  int i, j;

  for(i=0;i<size_y;i++)
  {
    for(j=0;j<size_x;j++)
    {
      imgX[i][j] = (imgpel)lut[imgX[i][j]];
    }
  }
}

void init_tone_mapping_sei(ToneMappingSEI *seiToneMapping) 
{
  seiToneMapping->seiHasTone_mapping = FALSE;
  seiToneMapping->count = 0;
}

void update_tone_mapping_sei(ToneMappingSEI *seiToneMapping) 
{

  if(seiToneMapping->tone_map_repetition_period == 0)
  {
    seiToneMapping->seiHasTone_mapping = FALSE;
    seiToneMapping->count = 0;
  }
  else if (seiToneMapping->tone_map_repetition_period>1)
  {
    seiToneMapping->count++;
    if (seiToneMapping->count>=seiToneMapping->tone_map_repetition_period) 
    {
      seiToneMapping->seiHasTone_mapping = FALSE;
      seiToneMapping->count = 0;
    }
  }
}
#endif

/*!
 ************************************************************************
 *  \brief
 *     Interpret the post filter hints SEI message (JVT-U035)
 *  \param payload
 *     a pointer that point to the sei payload
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *    
 ************************************************************************
 */
void interpret_post_filter_hints_info( byte* payload, int size, VideoParameters *p_Vid )
{
  data_partition_t* buf;
  unsigned int filter_hint_size_y, filter_hint_size_x, filter_hint_type, color_component, cx, cy, additional_extension_flag;
  int ***filter_hint;

  buf = new data_partition_t;
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  filter_hint_size_y = buf->ue("SEI: filter_hint_size_y"); // interpret post-filter hint SEI here
  filter_hint_size_x = buf->ue("SEI: filter_hint_size_x"); // interpret post-filter hint SEI here
  filter_hint_type   = buf->u(2, "SEI: filter_hint_type"); // interpret post-filter hint SEI here

  get_mem3Dint (&filter_hint, 3, filter_hint_size_y, filter_hint_size_x);

  for (color_component = 0; color_component < 3; color_component ++)
    for (cy = 0; cy < filter_hint_size_y; cy ++)
      for (cx = 0; cx < filter_hint_size_x; cx ++)
        filter_hint[color_component][cy][cx] = buf->se("SEI: filter_hint"); // interpret post-filter hint SEI here

  additional_extension_flag = buf->u(1, "SEI: additional_extension_flag"); // interpret post-filter hint SEI here

  free_mem3Dint (filter_hint);
  free( buf );
}



/*!
 ************************************************************************
 *  \brief
 *     Interpret the SEI rbsp
 *  \param msg
 *     a pointer that point to the sei message.
 *  \param size
 *     the size of the sei message
 *  \param p_Vid
 *     the image pointer
 *
 ************************************************************************
 */
void parse_sei(byte *msg, int size, VideoParameters *p_Vid, slice_t *pSlice)
{
    int payload_type = 0;
    int payload_size = 0;
    int offset = 1;
    byte tmp_byte;

    do {
        // sei_message();
        payload_type = 0;
        tmp_byte = msg[offset++];
        while (tmp_byte == 0xFF) {
            payload_type += 255;
            tmp_byte = msg[offset++];
        }
        payload_type += tmp_byte;   // this is the last byte

        payload_size = 0;
        tmp_byte = msg[offset++];
        while (tmp_byte == 0xFF) {
            payload_size += 255;
            tmp_byte = msg[offset++];
        }
        payload_size += tmp_byte;   // this is the last byte

        switch (payload_type) {    // sei_payload( type, size );

        case SEI_BUFFERING_PERIOD:
            interpret_buffering_period_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PIC_TIMING:
            interpret_picture_timing_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PAN_SCAN_RECT:
            interpret_pan_scan_rect_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FILLER_PAYLOAD:
            interpret_filler_payload_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_USER_DATA_REGISTERED_ITU_T_T35:
            interpret_user_data_registered_itu_t_t35_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_USER_DATA_UNREGISTERED:
            interpret_user_data_unregistered_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_RECOVERY_POINT:
            interpret_recovery_point_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_DEC_REF_PIC_MARKING_REPETITION:
            interpret_dec_ref_pic_marking_repetition_info( msg+offset, payload_size, p_Vid, pSlice );
            break;
        case SEI_SPARE_PIC:
            interpret_spare_pic( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SCENE_INFO:
            interpret_scene_information( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SUB_SEQ_INFO:
            interpret_subsequence_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SUB_SEQ_LAYER_CHARACTERISTICS:
            interpret_subsequence_layer_characteristics_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_SUB_SEQ_CHARACTERISTICS:
            interpret_subsequence_characteristics_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FULL_FRAME_FREEZE:
            interpret_full_frame_freeze_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FULL_FRAME_FREEZE_RELEASE:
            interpret_full_frame_freeze_release_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FULL_FRAME_SNAPSHOT:
            interpret_full_frame_snapshot_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PROGRESSIVE_REFINEMENT_SEGMENT_START:
            interpret_progressive_refinement_start_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_PROGRESSIVE_REFINEMENT_SEGMENT_END:
            interpret_progressive_refinement_end_info( msg+offset, payload_size, p_Vid );
            break;
        case SEI_MOTION_CONSTRAINED_SLICE_GROUP_SET:
            interpret_motion_constrained_slice_group_set_info( msg+offset, payload_size, p_Vid );
        case SEI_FILM_GRAIN_CHARACTERISTICS:
            interpret_film_grain_characteristics_info ( msg+offset, payload_size, p_Vid );
            break;
        case SEI_DEBLOCKING_FILTER_DISPLAY_PREFERENCE:
            interpret_deblocking_filter_display_preference_info ( msg+offset, payload_size, p_Vid );
            break;
        case SEI_STEREO_VIDEO_INFO:
            interpret_stereo_video_info_info ( msg+offset, payload_size, p_Vid );
            break;
        case SEI_TONE_MAPPING:
            interpret_tone_mapping( msg+offset, payload_size, p_Vid );
            break;
        case SEI_POST_FILTER_HINTS:
            interpret_post_filter_hints_info ( msg+offset, payload_size, p_Vid );
            break;
        case SEI_FRAME_PACKING_ARRANGEMENT:
            interpret_frame_packing_arrangement_info( msg+offset, payload_size, p_Vid );
            break;
        default:
            interpret_reserved_info( msg+offset, payload_size, p_Vid );
            break;    
        }
        offset += payload_size;

    } while (msg[offset] != 0x80);    // more_rbsp_data()  msg[offset] != 0x80

    // ignore the trailing bits rbsp_trailing_bits();
    assert(msg[offset] == 0x80);      // this is the trailing bits
    assert(offset + 1 == size);
}
