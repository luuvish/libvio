#ifndef _MBUFFERDEC_H_
#define _MBUFFERDEC_H_

#include "global.h"
#include "bitstream_cabac.h"

#define MAX_LIST_SIZE 33
//! definition of pic motion parameters
typedef struct pic_motion_params_old
{
  byte *      mb_field_decoding_flag;      //!< field macroblock indicator
} PicMotionParamsOld;

//! definition of pic motion parameters
typedef struct pic_motion_params
{
  struct storable_picture *ref_pic[2];  //!< referrence picture pointer
  MotionVector             mv[2];       //!< motion vector  
  char                     ref_idx[2];  //!< reference picture   [list][subblock_y][subblock_x]
  //byte                   mb_field_decoding_flag;    //!< field macroblock indicator
  byte                     slice_no;
} PicMotionParams;

struct DecRefPicMarking_s;

//! definition a picture (field or frame)
typedef struct storable_picture
{
  PictureStructure structure;

  int         poc;
  int         top_poc;
  int         bottom_poc;
  int         frame_poc;
  unsigned int  frame_num;
  unsigned int  recovery_frame;

  int         PicNum;
  int         LongTermPicNum;
  int         LongTermFrameIdx;

  byte        is_long_term;
  int         used_for_reference;
  int         is_output;
  int         non_existing;
  int         separate_colour_plane_flag;

  int         size_x, size_y, size_x_cr, size_y_cr;
  int         size_x_m1, size_y_m1, size_x_cr_m1, size_y_cr_m1;
  int         coded_frame;
  int         mb_aff_frame_flag;
  unsigned    PicWidthInMbs;
  unsigned    PicSizeInMbs;
  int         iChromaPadX;
  int         iChromaPadY;


  imgpel **     imgY;         //!< Y picture component
  imgpel ***    imgUV;        //!< U and V picture components

  struct pic_motion_params **mv_info;          //!< Motion info
  struct pic_motion_params **JVmv_info[MAX_PLANE];          //!< Motion info

  struct pic_motion_params_old  motion;              //!< Motion info  
  struct pic_motion_params_old  JVmotion[MAX_PLANE]; //!< Motion info for 4:4:4 independent mode decoding

  struct storable_picture *top_field;     // for mb aff, if frame for referencing the top field
  struct storable_picture *bottom_field;  // for mb aff, if frame for referencing the bottom field
  struct storable_picture *frame;         // for mb aff, if field for referencing the combined frame

  int         slice_type;
  int         idr_flag;
  int         no_output_of_prior_pics_flag;
  int         long_term_reference_flag;
  int         adaptive_ref_pic_buffering_flag;

  int         chroma_format_idc;
  int         frame_mbs_only_flag;
  int         frame_cropping_flag;
  int         frame_crop_left_offset;
  int         frame_crop_right_offset;
  int         frame_crop_top_offset;
  int         frame_crop_bottom_offset;
  struct DecRefPicMarking_s *dec_ref_pic_marking_buffer;                    //!< stores the memory management control operations

  // picture error concealment
  int         concealed_pic; //indicates if this is a concealed picture
  
  // variables for tone mapping
  int         seiHasTone_mapping;
  int         tone_mapping_model_id;
  int         tonemapped_bit_depth;  
  imgpel*     tone_mapping_lut;                //!< tone mapping look up table

  int         proc_flag;
#if (MVC_EXTENSION_ENABLE)
  int         view_id;
  int         inter_view_flag;
  int         anchor_pic_flag;
#endif
  int         iLumaStride;
  int         iChromaStride;
  int         iLumaExpandedHeight;
  int         iChromaExpandedHeight;
  imgpel **cur_imgY; // for more efficient get_block_luma
  int no_ref;
  int iCodingType;
  //
  char listXsize[MAX_NUM_SLICES][2];
  struct storable_picture **listX[MAX_NUM_SLICES][2];
  int         layer_id;
} StorablePicture;

typedef StorablePicture *StorablePicturePtr;

typedef struct frame_store {
    int       is_used;                //!< 0=empty; 1=top; 2=bottom; 3=both fields (or frame)
    int       is_reference;           //!< 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
    int       is_long_term;           //!< 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
    int       is_orig_reference;      //!< original marking by nal_ref_idc: 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used

    int       is_non_existent;

    uint32_t    FrameNum;
    int32_t     FrameNumWrap;
    unsigned  recovery_frame;

    int         LongTermFrameIdx;
    int       is_output;
    int       poc;

    // picture error concealment
    int concealment_reference;

    StorablePicture *frame;
    StorablePicture *top_field;
    StorablePicture *bottom_field;

#if (MVC_EXTENSION_ENABLE)
    int       view_id;
    int       inter_view_flag[2];
    int       anchor_pic_flag[2];
#endif
    int       layer_id;
} FrameStore;

typedef struct decoded_picture_buffer_t {
    VideoParameters *p_Vid;
    InputParameters *p_Inp;
    FrameStore  **fs;
    FrameStore  **fs_ref;
    FrameStore  **fs_ltref;
    FrameStore  **fs_ilref; // inter-layer reference (for multi-layered codecs)
    unsigned      size;
    unsigned      used_size;
    unsigned      ref_frames_in_buffer;
    unsigned      ltref_frames_in_buffer;
    int           last_output_poc;
#if (MVC_EXTENSION_ENABLE)
    int           last_output_view_id;
#endif
    int           max_long_term_pic_idx;  

    int           init_done;
    int           num_ref_frames;

    FrameStore   *last_picture;
    unsigned      used_size_il;
    int           layer_id;
} dpb_t;


extern void              init_dpb(VideoParameters *p_Vid, dpb_t *p_Dpb, int type);
extern void              re_init_dpb(VideoParameters *p_Vid, dpb_t *p_Dpb, int type);
extern void              free_dpb(dpb_t *p_Dpb);
extern FrameStore*       alloc_frame_store(void);
extern void              free_frame_store (FrameStore* f);
extern StorablePicture*  alloc_storable_picture(VideoParameters *p_Vid, PictureStructure type, int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output);
extern void              free_storable_picture (StorablePicture* p);

#if (MVC_EXTENSION_ENABLE)
extern void             idr_memory_management(dpb_t *p_Dpb, StorablePicture* p);
extern void             flush_dpbs(dpb_t **p_Dpb, int nLayers);
extern int              GetMaxDecFrameBuffering(VideoParameters *p_Vid);
extern void             append_interview_list(dpb_t *p_Dpb, 
                                              bool field_pic_flag, bool bottom_field_flag,
                                              int list_idx, FrameStore **list, int *listXsize, int currPOC, 
                                              int curr_view_id, int anchor_pic_flag);
#endif

struct slice_t;


extern void update_ref_list  (dpb_t *p_Dpb);
extern void update_ltref_list(dpb_t *p_Dpb);
extern void update_pic_num   (slice_t *currSlice);

extern void gen_pic_list_from_frame_list(bool bottom_field_flag, FrameStore **fs_list, int list_idx, StorablePicture **list, char *list_size, int long_term);
extern void init_lists_i_slice (slice_t *currSlice);
extern void init_lists_p_slice (slice_t *currSlice);
extern void init_lists_b_slice (slice_t *currSlice);
extern void init_lists         (slice_t *currSlice);

extern void reorder_ref_pic_list(struct slice_t *currSlice, int cur_list);
extern void reorder_lists   (slice_t *currSlice);
extern void init_mbaff_lists(VideoParameters *p_Vid, slice_t *currSlice);

extern void store_picture_in_dpb(dpb_t *p_Dpb, StorablePicture* p);




extern void unmark_for_reference(FrameStore* fs);
extern void unmark_for_long_term_reference(FrameStore* fs);
extern void insert_picture_in_dpb(VideoParameters *p_Vid, FrameStore* fs, StorablePicture* p);
extern void remove_frame_from_dpb(dpb_t *p_Dpb, int pos);
extern int  output_one_frame_from_dpb(dpb_t *p_Dpb);

extern void flush_dpb(dpb_t *p_Dpb);

extern void dpb_split_field      (VideoParameters *p_Vid, FrameStore *fs);
extern void dpb_combine_field    (VideoParameters *p_Vid, FrameStore *fs);
extern void dpb_combine_field_yuv(VideoParameters *p_Vid, FrameStore *fs);

extern void fill_frame_num_gap(VideoParameters *p_Vid, slice_t *pSlice);


extern int  init_img_data(VideoParameters *p_Vid, ImageData *p_ImgData, sps_t *sps);
extern void free_img_data(VideoParameters *p_Vid, ImageData *p_ImgData);
extern void pad_dec_picture(VideoParameters *p_Vid, StorablePicture *dec_picture);
extern void pad_buf(imgpel *pImgBuf, int iWidth, int iHeight, int iStride, int iPadX, int iPadY);
extern void process_picture_in_dpb_s(VideoParameters *p_Vid, StorablePicture *p_pic);
extern StorablePicture * clone_storable_picture( VideoParameters *p_Vid, StorablePicture *p_pic );
extern void store_proc_picture_in_dpb(dpb_t *p_Dpb, StorablePicture* p);


#if (MVC_EXTENSION_ENABLE)
extern void reorder_lists_mvc     (slice_t *currSlice, int currPOC);
extern void init_lists_p_slice_mvc(slice_t *currSlice);
extern void init_lists_b_slice_mvc(slice_t *currSlice);
extern void init_lists_i_slice_mvc(slice_t *currSlice);

extern void reorder_ref_pic_list_mvc(slice_t *currSlice, int cur_list, int **anchor_ref, int **non_anchor_ref,
                                                 int view_id, int anchor_pic_flag, int currPOC, int listidx);

extern void reorder_short_term(slice_t *currSlice, int cur_list, int num_ref_idx_lX_active_minus1, int picNumLX, int *refIdxLX, int currViewID);
extern void reorder_long_term(slice_t *currSlice, StorablePicture **RefPicListX, int num_ref_idx_lX_active_minus1, int LongTermPicNum, int *refIdxLX, int currViewID);
#endif



static inline int compare_pic_by_pic_num_desc( const void *arg1, const void *arg2 )
{
  int pic_num1 = (*(StorablePicture**)arg1)->PicNum;
  int pic_num2 = (*(StorablePicture**)arg2)->PicNum;

  if (pic_num1 < pic_num2)
    return 1;
  if (pic_num1 > pic_num2)
    return -1;
  else
    return 0;
}

static inline int compare_pic_by_lt_pic_num_asc( const void *arg1, const void *arg2 )
{
  int long_term_pic_num1 = (*(StorablePicture**)arg1)->LongTermPicNum;
  int long_term_pic_num2 = (*(StorablePicture**)arg2)->LongTermPicNum;

  if ( long_term_pic_num1 < long_term_pic_num2)
    return -1;
  if ( long_term_pic_num1 > long_term_pic_num2)
    return 1;
  else
    return 0;
}

static inline int compare_fs_by_frame_num_desc( const void *arg1, const void *arg2 )
{
  int frame_num_wrap1 = (*(FrameStore**)arg1)->FrameNumWrap;
  int frame_num_wrap2 = (*(FrameStore**)arg2)->FrameNumWrap;
  if ( frame_num_wrap1 < frame_num_wrap2)
    return 1;
  if ( frame_num_wrap1 > frame_num_wrap2)
    return -1;
  else
    return 0;
}

static inline int compare_fs_by_lt_pic_idx_asc( const void *arg1, const void *arg2 )
{
  int long_term_frame_idx1 = (*(FrameStore**)arg1)->LongTermFrameIdx;
  int long_term_frame_idx2 = (*(FrameStore**)arg2)->LongTermFrameIdx;

  if ( long_term_frame_idx1 < long_term_frame_idx2)
    return -1;
  else if ( long_term_frame_idx1 > long_term_frame_idx2)
    return 1;
  else
    return 0;
}

static inline int compare_pic_by_poc_asc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(StorablePicture**)arg1)->poc;
  int poc2 = (*(StorablePicture**)arg2)->poc;

  if ( poc1 < poc2)
    return -1;  
  else if ( poc1 > poc2)
    return 1;
  else
    return 0;
}

static inline int compare_pic_by_poc_desc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(StorablePicture**)arg1)->poc;
  int poc2 = (*(StorablePicture**)arg2)->poc;

  if (poc1 < poc2)
    return 1;
  else if (poc1 > poc2)
    return -1;
  else
    return 0;
}

static inline int compare_fs_by_poc_asc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(FrameStore**)arg1)->poc;
  int poc2 = (*(FrameStore**)arg2)->poc;

  if (poc1 < poc2)
    return -1;
  else if (poc1 > poc2)
    return 1;
  else
    return 0;
}

static inline int compare_fs_by_poc_desc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(FrameStore**)arg1)->poc;
  int poc2 = (*(FrameStore**)arg2)->poc;

  if (poc1 < poc2)
    return 1;
  else if (poc1 > poc2)
    return -1;
  else
    return 0;
}


extern void get_smallest_poc(dpb_t *p_Dpb, int *poc,int * pos);
extern int remove_unused_frame_from_dpb(dpb_t *p_Dpb);

#endif
