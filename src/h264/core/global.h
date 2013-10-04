#ifndef _GLOBAL_H_
#define _GLOBAL_H_


#include <chrono>

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
  CM_UNKNOWN = -1,
  CM_YUV     =  0,
  CM_RGB     =  1,
  CM_XYZ     =  2
} ColorModel;

typedef enum {
  CF_UNKNOWN = -1,     //!< Unknown color format
  YUV400     =  0,     //!< Monochrome
  YUV420     =  1,     //!< 4:2:0
  YUV422     =  2,     //!< 4:2:2
  YUV444     =  3      //!< 4:4:4
} ColorFormat;

typedef enum {
  PF_UNKNOWN = -1,     //!< Unknown color ordering
  UYVY       =  0,     //!< UYVY
  YUY2       =  1,     //!< YUY2
  YUYV       =  1,     //!< YUYV
  YVYU       =  2,     //!< YVYU
  BGR        =  3,     //!< BGR
  V210       =  4      //!< Video Clarity 422 format (10 bits)
} PixelFormat;


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

typedef enum
{
  // YUV
  PLANE_Y = 0,  // PLANE_Y
  PLANE_U = 1,  // PLANE_Cb
  PLANE_V = 2,  // PLANE_Cr
  // RGB
  PLANE_G = 0,
  PLANE_B = 1,
  PLANE_R = 2
} ColorPlane;

enum {
  LIST_0 = 0,
  LIST_1 = 1,
  BI_PRED = 2,
  BI_PRED_L0 = 3,
  BI_PRED_L1 = 4
};

//! Field Coding Types
typedef enum
{
  FRAME_CODING         = 0,
  FIELD_CODING         = 1,
  ADAPTIVE_CODING      = 2,
  FRAME_MB_PAIR_CODING = 3
} CodingType;


typedef enum
{
  FRAME,
  TOP_FIELD,
  BOTTOM_FIELD
} PictureStructure;           //!< New enum for field processing

typedef enum
{
  P_SLICE = 0,
  B_SLICE = 1,
  I_SLICE = 2,
  SP_SLICE = 3,
  SI_SLICE = 4,
  NUM_SLICE_TYPES = 5
} SliceType;


#define ET_SIZE 300      //!< size of error text buffer
extern char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

struct pic_motion_params_old;
struct pic_motion_params;

struct slice_t;
struct macroblock_t;

/***********************************************************************
 * T y p e    d e f i n i t i o n s    f o r    J M
 ***********************************************************************
 */
typedef enum {
    LumaComp = 0,
    CrComp = 1,
    CbComp = 2
} Color_Component;

// Motion Vector structure
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



/***********************************************************************
 * N e w   D a t a    t y p e s   f o r    T M L
 ***********************************************************************
 */

typedef struct decodedpic_t {
    int                  bValid;                 //0: invalid, 1: valid, 3: valid for 3D output;
    int                  iViewId;                //-1: single view, >=0 multiview[VIEW1|VIEW0];
    int                  iPOC;
    int                  iYUVFormat;             //0: 4:0:0, 1: 4:2:0, 2: 4:2:2, 3: 4:4:4
    int                  iYUVStorageFormat;      //0: YUV seperate; 1: YUV interleaved; 2: 3D output;
    int                  iBitDepth;
    byte                *pY;                     //if iPictureFormat is 1, [0]: top; [1] bottom;
    byte                *pU;
    byte                *pV;
    int                  iWidth;                 //frame width;              
    int                  iHeight;                //frame height;
    int                  iYBufStride;            //stride of pY[0/1] buffer in bytes;
    int                  iUVBufStride;           //stride of pU[0/1] and pV[0/1] buffer in bytes;
    int                  iSkipPicNum;
    int                  iBufSize;
    struct decodedpic_t *pNext;
} DecodedPicList;

//****************************** ~DM ***********************************
typedef struct coding_par {
    int layer_id;

    //padding info;
    void (*img2buf)(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);

    mb_t *mb_data;               //!< array containing all MBs of a whole frame
    mb_t *mb_data_JV[MAX_PLANE]; //!< mb_data to be used for 4:4:4 independent mode
} CodingParameters;

struct decoded_picture_buffer_t;

typedef struct layer_par {
    int               layer_id;
    struct video_par* p_Vid;
    CodingParameters* p_Cps;
    sps_t*            p_SPS;
    decoded_picture_buffer_t* p_Dpb;
} LayerParameters;



// input parameters from configuration file
typedef struct inp_par
{
    char infile[FILE_NAME_SIZE];                       //!< H.264 inputfile
    char outfile[FILE_NAME_SIZE];                      //!< Decoded YUV 4:2:0 output
    char reffile[FILE_NAME_SIZE];                      //!< Optional YUV 4:2:0 reference file for SNR measurement

    int FileFormat;                         //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
    int ref_offset;
    int poc_scale;
    int write_uv;
    int silent;

    // Input/output sequence format related variables

#if (MVC_EXTENSION_ENABLE)
    int  DecodeAllLayers;
#endif

    // picture error concealment
    int conceal_mode;
    int ref_poc_gap;
    int poc_gap;
  
    int iDecFrmNum;

    int bDisplayDecParams;
    int dpb_plus[2];
} InputParameters;

typedef struct old_slice_par
{
    unsigned field_pic_flag;   
    unsigned frame_num;
    int      nal_ref_idc;
    unsigned pic_oder_cnt_lsb;
    int      delta_pic_oder_cnt_bottom;
    int      delta_pic_order_cnt[2];
    byte     bottom_field_flag;
    byte     idr_flag;
    int      idr_pic_id;
    int      pps_id;
#if (MVC_EXTENSION_ENABLE)
    int      view_id;
    int      inter_view_flag;
    int      anchor_pic_flag;
#endif
    int      layer_id;
} OldSliceParams;

// signal to noise ratio parameters
typedef struct snr_par
{
    int   frame_ctr;
    float snr[3];                                //!< current SNR (component)
    float snr1[3];                               //!< SNR (dB) first frame (component)
    float snra[3];                               //!< Average component SNR (dB) remaining frames
    float sse[3];                                //!< component SSE 
    float msse[3];                                //!< Average component SSE 
} SNRParameters;

struct tone_mapping_struct_s;
struct frame_store;
struct sei_params;
struct concealment_node;
struct object_buffer;
struct ercVariables_s;
struct storable_picture;

// video parameters
typedef struct video_par {
    InputParameters      *p_Inp;

    int p_out; // FILE handle for YUV output file
#if (MVC_EXTENSION_ENABLE)
    int p_out_mvc[MAX_VIEW_NUM];
#endif
    int p_ref;

    bitstream_t bitstream;



    OldSliceParams               *old_slice;
    SNRParameters                *snr;

    decoded_picture_buffer_t* p_Dpb_layer[MAX_NUM_DPB_LAYERS];
    CodingParameters* p_EncodePar[MAX_NUM_DPB_LAYERS];
    LayerParameters*  p_LayerPar [MAX_NUM_DPB_LAYERS];

    bool global_init_done[2];

#if (ENABLE_OUTPUT_TONEMAPPING)
    tone_mapping_struct_s *seiToneMapping;
#endif

    int              iNumOfSlicesAllocated;
    int              iSliceNumOfCurrPic;
    slice_t **ppSliceList;
    slice_t  *pNextSlice;
    int       newframe;

    nalu_t *nalu;

    DecodedPicList *pDecOuputPic;
    pps_t *pNextPPS;
    bool first_sps; // use only for print first sps


    frame_store *out_buffer;

    // Timing related variables
    std::chrono::system_clock::time_point start_time;
    std::chrono::system_clock::time_point end_time;
    int64_t                               tot_time;

    // B pictures
    int  Bframe_ctr;
    int  frame_no;
    int  g_nFrame;


    int32_t     prevPicOrderCntMsb;
    uint32_t    prevPicOrderCntLsb;
    uint32_t    prevFrameNum;
    uint32_t    prevFrameNumOffset;
    int         last_has_mmco_5;
    int         last_pic_bottom_field;

    // FMO
    int        *MbToSliceGroupMap;
    int        *MapUnitToSliceGroupMap;

    pps_t      *active_pps;
    sps_t      *active_sps;
    sps_t       SeqParSet[MAXSPS];
    pps_t       PicParSet[MAXPPS];

#if (MVC_EXTENSION_ENABLE)
    sub_sps_t  *active_subset_sps;
    sub_sps_t   SubsetSeqParSet[MAXSPS];
#endif

  sei_params        *p_SEI;

  int number;                                 //!< frame number

  //current picture property;
  unsigned int num_dec_mb;

  int type;                                   //!< image type INTER/INTRA
  int structure;                     //!< Identify picture structure type

    // global picture format dependent buffers, memory allocation in decod.c
    mb_t *mb_data;               //!< array containing all MBs of a whole frame
    mb_t *mb_data_JV[MAX_PLANE]; //!< mb_data to be used for 4:4:4 independent mode


  // picture error concealment
  // concealment_head points to first node in list, concealment_end points to
  // last node in list. Initialize both to NULL, meaning no nodes in list yet
  concealment_node *concealment_head;
  concealment_node *concealment_end;

  unsigned int pre_frame_num;           //!< store the frame_num in the last decoded slice. For detecting gap in frame_num.
  int non_conforming_stream;

  // for POC mode 1:
  int ThisPOC;
  // /////////////////////////

  int no_output_of_prior_pics_flag;

  int idr_psnr_number;
  int psnr_number;


  // picture error concealment
  int last_ref_pic_poc;
  int ref_poc_gap;
  int poc_gap;
  int conceal_mode;
  int earlier_missing_poc;
  unsigned int frame_to_conceal;
  int IDR_concealment_flag;
  int conceal_slice_type;

  // random access point decoding
  int recovery_point;
  int recovery_point_found;
  int recovery_frame_cnt;
  int recovery_frame_num;
  int recovery_poc;

  // Redundant slices. Should be moved to another structure and allocated only if extended profile
  unsigned int previous_frame_num; //!< frame number of previous slice
  //!< non-zero: i-th previous frame is correct
  int Is_primary_correct;          //!< if primary frame is correct, 0: incorrect
  int Is_redundant_correct;        //!< if redundant frame is correct, 0:incorrect


  frame_store *last_out_fs;
  int pocs_in_dpb[100];

  storable_picture *dec_picture;
  storable_picture *dec_picture_JV[MAX_PLANE];  //!< dec_picture to be used during 4:4:4 independent mode decoding
  storable_picture *no_reference_picture; //!< dummy storable picture for recovery point

  // Error parameters
  object_buffer  *erc_object_list;
  ercVariables_s *erc_errorVar;

  int erc_mvperMB;
  video_par *erc_img;

  storable_picture *pending_output;
  int    recovery_flag;

  void (*buf2img)(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int o_size_x, int o_size_y, int symbol_size_in_bytes, int bitshift);
  void (*img2buf)(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);

  ImageData tempData3;
  //control;
  int last_dec_layer_id;
  int dpb_layer_id;

/******************* deprecative variables; ***************************************/
  // Fidelity Range Extensions Stuff

  int profile_idc;
} VideoParameters;


typedef struct decoder_params {
    InputParameters *p_Inp;
    VideoParameters *p_Vid;
} DecoderParams;

extern DecoderParams  *p_Dec;

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
