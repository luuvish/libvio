#ifndef _GLOBAL_H_
#define _GLOBAL_H_


#include <cstdint>

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

#include "defines.h"
#include "ifunctions.h"

#include "parset.h"
#include "bitstream.h"
#include "macroblock.h"

using vio::h264::bitstream_t;
using vio::h264::mb_t;


typedef enum {
    CF_UNKNOWN = -1,     //!< Unknown color format
    YUV400     =  0,     //!< Monochrome
    YUV420     =  1,     //!< 4:2:0
    YUV422     =  2,     //!< 4:2:2
    YUV444     =  3      //!< 4:4:4
} ColorFormat;


typedef struct image_data
{
  ColorFormat yuv_format;                    //!< YUV format (0=4:0:0, 1=4:2:0, 2=4:2:2, 3=4:4:4)
  // Standard data
  imgpel **frm_data[MAX_PLANE];     //!< Frame Data
  imgpel **top_data[MAX_PLANE];     //!< pointers to top field data
  imgpel **bot_data[MAX_PLANE];     //!< pointers to bottom field data

  imgpel **frm_data_buf[2][MAX_PLANE];     //!< Frame Data
  imgpel **top_data_buf[2][MAX_PLANE];     //!< pointers to top field data
  imgpel **bot_data_buf[2][MAX_PLANE];     //!< pointers to bottom field data
  
  //! Optional data (could also add uint8 data in case imgpel is of type uint16_t)
  //! These can be useful for enabling input/conversion of content of different types
  //! while keeping optimal processing size.
  uint16_t **frm_uint16[MAX_PLANE];   //!< optional frame Data for uint16_t
  uint16_t **top_uint16[MAX_PLANE];   //!< optional pointers to top field data
  uint16_t **bot_uint16[MAX_PLANE];   //!< optional pointers to bottom field data

  int frm_stride[MAX_PLANE];
  int top_stride[MAX_PLANE];
  int bot_stride[MAX_PLANE];
} ImageData;

/***********************************************************************
 * T y p e    d e f i n i t i o n s    f o r    T M L
 ***********************************************************************
 */

typedef enum {
    PLANE_Y = 0,  // PLANE_Y
    PLANE_U = 1,  // PLANE_Cb
    PLANE_V = 2,  // PLANE_Cr
} ColorPlane;

enum {
    LIST_0 = 0,
    LIST_1 = 1
};

//! Field Coding Types
typedef enum {
    FRAME_CODING         = 0,
    FIELD_CODING         = 1,
    ADAPTIVE_CODING      = 2,
    FRAME_MB_PAIR_CODING = 3
} CodingType;


typedef enum {
    FRAME,
    TOP_FIELD,
    BOTTOM_FIELD
} PictureStructure;           //!< New enum for field processing


#define ET_SIZE 300      //!< size of error text buffer
extern char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

struct pic_motion_params_old;
struct pic_motion_params;

struct slice_t;
struct macroblock_t;

struct motion_vector_t {
    int16_t mv_x;
    int16_t mv_y;
};

using mv_t = motion_vector_t;

inline bool operator == (const mv_t& l, const mv_t& r)
{
    return (l.mv_x == r.mv_x) && (l.mv_y == r.mv_y);
}

inline mv_t operator + (const mv_t& l, const mv_t& r)
{
    return {(int16_t)(l.mv_x + r.mv_x), (int16_t)(l.mv_y + r.mv_y)};
}

inline mv_t operator + (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x + r), (int16_t)(l.mv_y + r)};
}

inline mv_t operator - (const mv_t& l, const mv_t& r)
{
    return {(int16_t)(l.mv_x - r.mv_x), (int16_t)(l.mv_y - r.mv_y)};
}

inline mv_t operator - (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x - r), (int16_t)(l.mv_y - r)};
}

inline mv_t operator * (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x * r), (int16_t)(l.mv_y * r)};
}

inline mv_t operator * (int l, const mv_t& r)
{
    return {(int16_t)(l * r.mv_x), (int16_t)(l * r.mv_y)};
}

inline mv_t operator >> (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x >> r), (int16_t)(l.mv_y >> r)};
}



struct CodingParameters {
    int layer_id;

    //padding info;
    void (*img2buf)(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);

    mb_t* mb_data;
    mb_t* mb_data_JV[MAX_PLANE];
};

struct decoded_picture_buffer_t;

struct InputParameters;
struct VideoParameters;
struct DecodedPicList;

struct LayerParameters {
    int               layer_id;
    VideoParameters*  p_Vid;
    CodingParameters* p_Cps;
    sps_t*            p_SPS;
    decoded_picture_buffer_t* p_Dpb;
};

struct OldSliceParams {
    unsigned    field_pic_flag;   
    unsigned    frame_num;
    int         nal_ref_idc;
    unsigned    pic_oder_cnt_lsb;
    int         delta_pic_oder_cnt_bottom;
    int         delta_pic_order_cnt[2];
    bool        bottom_field_flag;
    bool        idr_flag;
    int         idr_pic_id;
    int         pps_id;
#if (MVC_EXTENSION_ENABLE)
    int         view_id;
    int         inter_view_flag;
    int         anchor_pic_flag;
#endif
    int         layer_id;
};

struct SNRParameters;

struct tone_mapping_struct_s;
struct frame_store;
struct sei_params;
struct concealment_node;
struct object_buffer;
struct ercVariables_s;
struct storable_picture;

// video parameters
struct VideoParameters {
    InputParameters* p_Inp;

    int         p_out; // FILE handle for YUV output file
#if (MVC_EXTENSION_ENABLE)
    int         p_out_mvc[MAX_VIEW_NUM];
#endif
    int         p_ref;

    bitstream_t bitstream;



    OldSliceParams* old_slice;
    SNRParameters*  snr;

    decoded_picture_buffer_t* p_Dpb_layer[MAX_NUM_DPB_LAYERS];
    CodingParameters* p_EncodePar[MAX_NUM_DPB_LAYERS];
    LayerParameters*  p_LayerPar [MAX_NUM_DPB_LAYERS];

    bool        global_init_done[2];

#if (ENABLE_OUTPUT_TONEMAPPING)
    tone_mapping_struct_s *seiToneMapping;
#endif

    int         iNumOfSlicesAllocated;
    int         iSliceNumOfCurrPic;
    slice_t**   ppSliceList;
    slice_t*    pNextSlice;
    int         newframe;

    nalu_t*     nalu;

    DecodedPicList* pDecOuputPic;
    pps_t*      pNextPPS;
    bool        first_sps; // use only for print first sps


    frame_store* out_buffer;


    int32_t     prevPicOrderCntMsb;
    uint32_t    prevPicOrderCntLsb;
    uint32_t    prevFrameNum;
    uint32_t    prevFrameNumOffset;
    int         last_has_mmco_5;
    int         last_pic_bottom_field;

    // FMO
    int*        MbToSliceGroupMap;
    int*        MapUnitToSliceGroupMap;

    pps_t*      active_pps;
    sps_t*      active_sps;
    sps_t       SeqParSet[MAXSPS];
    pps_t       PicParSet[MAXPPS];

#if (MVC_EXTENSION_ENABLE)
    sub_sps_t*  active_subset_sps;
    sub_sps_t   SubsetSeqParSet[MAXSPS];
#endif

    sei_params* p_SEI;

    int         number;                                 //!< frame number

    //current picture property;
    unsigned    num_dec_mb;

    int         type;                                   //!< image type INTER/INTRA
    int         structure;                     //!< Identify picture structure type

    // global picture format dependent buffers, memory allocation in decod.c
    mb_t*       mb_data;               //!< array containing all MBs of a whole frame
    mb_t*       mb_data_JV[MAX_PLANE]; //!< mb_data to be used for 4:4:4 independent mode


    // picture error concealment
    // concealment_head points to first node in list, concealment_end points to
    // last node in list. Initialize both to NULL, meaning no nodes in list yet
    concealment_node* concealment_head;
    concealment_node* concealment_end;

    unsigned    pre_frame_num;           //!< store the frame_num in the last decoded slice. For detecting gap in frame_num.
    int         non_conforming_stream;

    // for POC mode 1:
    int         ThisPOC;

    int         no_output_of_prior_pics_flag;

    // picture error concealment
    int         last_ref_pic_poc;
    int         conceal_mode;
    int         earlier_missing_poc;
    unsigned    frame_to_conceal;
    int         IDR_concealment_flag;
    int         conceal_slice_type;

    // random access point decoding
    int         recovery_point;
    int         recovery_point_found;
    int         recovery_frame_cnt;
    int         recovery_frame_num;
    int         recovery_poc;

    // Redundant slices. Should be moved to another structure and allocated only if extended profile
    unsigned    previous_frame_num; //!< frame number of previous slice
    //!< non-zero: i-th previous frame is correct
    int         Is_primary_correct;          //!< if primary frame is correct, 0: incorrect
    int         Is_redundant_correct;        //!< if redundant frame is correct, 0:incorrect

    frame_store* last_out_fs;
    int         pocs_in_dpb[100];

    storable_picture* dec_picture;
    storable_picture* dec_picture_JV[MAX_PLANE];  //!< dec_picture to be used during 4:4:4 independent mode decoding
    storable_picture* no_reference_picture; //!< dummy storable picture for recovery point

    // Error parameters
    object_buffer*  erc_object_list;
    ercVariables_s* erc_errorVar;

    int         erc_mvperMB;
    VideoParameters* erc_img;

    int         recovery_flag;

    void (*img2buf)(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);

    ImageData   tempData3;
    //control;
    int         last_dec_layer_id;
    int         dpb_layer_id;

    int         profile_idc;

    VideoParameters();
    ~VideoParameters();

    void calculate_frame_no(storable_picture *p);
    void status(storable_picture** dec_picture);
    void report();
};


// prototypes
extern void error(const char *text, int code);

// dynamic mem allocation
extern void free_layer_buffers( VideoParameters *p_Vid, int layer_id );
extern void OpenOutputFiles(VideoParameters *p_Vid, int view0_id, int view1_id);

static inline int is_FREXT_profile(unsigned int profile_idc) 
{
    return ( profile_idc >= FREXT_HP || profile_idc == FREXT_CAVLC444 );
}


#endif
