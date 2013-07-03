
/*!
 ************************************************************************
 *  \file
 *     global.h
 *  \brief
 *     global definitions for H.264 decoder.
 *  \author
 *     Copyright (C) 1999  Telenor Satellite Services,Norway
 *                         Ericsson Radio Systems, Sweden
 *
 *     Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *
 *     Telenor Satellite Services
 *     Keysers gt.13                       tel.:   +47 23 13 86 98
 *     N-0130 Oslo,Norway                  fax.:   +47 22 77 79 80
 *
 *     Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *
 *     Ericsson Radio Systems
 *     KI/ERA/T/VV
 *     164 80 Stockholm, Sweden
 *
 ************************************************************************
 */
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/timeb.h>

#include "win32.h"
#include "defines.h"
#include "ifunctions.h"

#include "parset.h"

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


typedef struct frame_format
{  
  ColorFormat yuv_format;                    //!< YUV format (0=4:0:0, 1=4:2:0, 2=4:2:2, 3=4:4:4)
  ColorModel  color_model;                   //!< 4:4:4 format (0: YUV, 1: RGB, 2: XYZ)
  PixelFormat pixel_format;                  //!< pixel format support for certain interleaved yuv sources
  double      frame_rate;                    //!< frame rate
  int         width[3];                      //!< component frame width
  int         height[3];                     //!< component frame height    
  int         auto_crop_right;               //!< luma component auto crop right
  int         auto_crop_bottom;              //!< luma component auto crop bottom
  int         auto_crop_right_cr;            //!< chroma component auto crop right
  int         auto_crop_bottom_cr;           //!< chroma component auto crop bottom
  int         width_crop;                    //!< width after cropping consideration
  int         height_crop;                   //!< height after cropping consideration
  int         mb_width;                      //!< luma component frame width
  int         mb_height;                     //!< luma component frame height    
  int         size_cmp[3];                   //!< component sizes (width * height)
  int         size;                          //!< total image size (sum of size_cmp)
  int         bit_depth[3];                  //!< component bit depth  
  int         max_value[3];                  //!< component max value
  int         max_value_sq[3];               //!< component max value squared
  int         pic_unit_size_on_disk;         //!< picture sample unit size on storage medium
  int         pic_unit_size_shift3;          //!< pic_unit_size_on_disk >> 3
} FrameFormat;

typedef struct image_data
{
  FrameFormat format;               //!< image format
  // Standard data
  imgpel **frm_data[MAX_PLANE];     //!< Frame Data
  imgpel **top_data[MAX_PLANE];     //!< pointers to top field data
  imgpel **bot_data[MAX_PLANE];     //!< pointers to bottom field data

  imgpel **frm_data_buf[2][MAX_PLANE];     //!< Frame Data
  imgpel **top_data_buf[2][MAX_PLANE];     //!< pointers to top field data
  imgpel **bot_data_buf[2][MAX_PLANE];     //!< pointers to bottom field data
  
  //! Optional data (could also add uint8 data in case imgpel is of type uint16)
  //! These can be useful for enabling input/conversion of content of different types
  //! while keeping optimal processing size.
  uint16 **frm_uint16[MAX_PLANE];   //!< optional frame Data for uint16
  uint16 **top_uint16[MAX_PLANE];   //!< optional pointers to top field data
  uint16 **bot_uint16[MAX_PLANE];   //!< optional pointers to bottom field data

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
  CAVLC,
  CABAC
} SymbolMode;


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

typedef enum
{
  IS_LUMA = 0,
  IS_CHROMA = 1
} Component_Type;


#define ET_SIZE 300      //!< size of error text buffer
extern char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

struct pic_motion_params_old;
struct pic_motion_params;

struct slice_t;
struct macroblock_dec;

/***********************************************************************
 * T y p e    d e f i n i t i o n s    f o r    J M
 ***********************************************************************
 */
typedef enum {
    LumaComp = 0,
    CrComp = 1,
    CbComp = 2
} Color_Component;

typedef struct pix_pos {
    int   available;
    int   mb_addr;
    short x;
    short y;
    short pos_x;
    short pos_y;
} PixelPos;

// Motion Vector structure
typedef struct {
    short mv_x;
    short mv_y;
} MotionVector;

static const MotionVector zero_mv = {0, 0};

typedef struct {
    short x;
    short y;
} BlockPos;

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
    int profile_idc;
    int width;
    int height;
    int width_cr;                               //!< width chroma  
    int height_cr;                              //!< height chroma

    int pic_unit_bitsize_on_disk;
    short bitdepth_luma;
    short bitdepth_chroma;
    int bitdepth_scale[2];
    int bitdepth_luma_qp_scale;
    int bitdepth_chroma_qp_scale;
    unsigned int dc_pred_value_comp[MAX_PLANE]; //!< component value for DC prediction (depends on component pel bit depth)
    int max_pel_value_comp[MAX_PLANE];       //!< max value that one picture element (pixel) can take (depends on pic_unit_bitdepth)

    int yuv_format;
    int lossless_qpprime_flag;
    int num_blk8x8_uv;
    int num_uv_blocks;
    int num_cdc_coeff;
    int mb_cr_size_x;
    int mb_cr_size_y;
    int mb_cr_size_x_blk;
    int mb_cr_size_y_blk;
    int mb_cr_size;
    int mb_size[3][2];                         //!< component macroblock dimensions
    int mb_size_blk[3][2];                     //!< component macroblock dimensions 
    int mb_size_shift[3][2];

    int max_vmv_r;                             //!< maximum vertical motion vector range in luma quarter frame pixel units for the current level_idc
    int separate_colour_plane_flag;
    int ChromaArrayType;
    int max_frame_num;
    unsigned int PicWidthInMbs;
    unsigned int PicHeightInMapUnits;
    unsigned int FrameHeightInMbs;
    unsigned int FrameSizeInMbs;
    int iLumaPadX;
    int iLumaPadY;
    int iChromaPadX;
    int iChromaPadY;

    int subpel_x;
    int subpel_y;
    int shiftpel_x;
    int shiftpel_y;
    int total_scale;
    unsigned int oldFrameSizeInMbs;

    //padding info;
    void (*img2buf)(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);
    int rgb_output;

    imgpel **imgY_ref;                              //!< reference frame find snr
    imgpel ***imgUV_ref;
    struct macroblock_dec *mb_data;               //!< array containing all MBs of a whole frame
    struct macroblock_dec *mb_data_JV[MAX_PLANE]; //!< mb_data to be used for 4:4:4 independent mode
    char  *intra_block;
    char  *intra_block_JV[MAX_PLANE];
    BlockPos *PicPos;  
    byte **ipredmode;                  //!< prediction type [90][74]
    byte **ipredmode_JV[MAX_PLANE];
    byte ****nz_coeff;
    int **siblock;
    int **siblock_JV[MAX_PLANE];
    int *qp_per_matrix;
    int *qp_rem_matrix;
} CodingParameters;

typedef struct layer_par {
  int                            layer_id;
  struct video_par              *p_Vid;
  CodingParameters              *p_Cps;
  sps_t      *p_SPS;
  struct decoded_picture_buffer *p_Dpb;
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
    int intra_profile_deblocking;               //!< Loop filter usage determined by flags and parameters in bitstream 

    // Input/output sequence format related variables
    FrameFormat source;                   //!< source related information
    FrameFormat output;                   //!< output related information

#if (MVC_EXTENSION_ENABLE)
    int  DecodeAllLayers;
#endif

    // picture error concealment
    int conceal_mode;
    int ref_poc_gap;
    int poc_gap;

    // Needed to allow compilation for decoder. May be used later for distortion computation operations
    int stdRange;                         //!< 1 - standard range, 0 - full range
    int videoCode;                        //!< 1 - 709, 3 - 601:  See VideoCode in io_tiff.
    int export_views;
  
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

// video parameters
typedef struct video_par {
    InputParameters      *p_Inp;

    int p_out; // FILE handle for YUV output file
#if (MVC_EXTENSION_ENABLE)
    int p_out_mvc[MAX_VIEW_NUM];
#endif
    int p_ref;

    struct bitstream_t *bitstream;



    OldSliceParams               *old_slice;
    SNRParameters                *snr;

    struct decoded_picture_buffer *p_Dpb_layer[MAX_NUM_DPB_LAYERS];
    CodingParameters              *p_EncodePar[MAX_NUM_DPB_LAYERS];
    LayerParameters               *p_LayerPar [MAX_NUM_DPB_LAYERS];

    bool global_init_done[2];

#if (ENABLE_OUTPUT_TONEMAPPING)
    struct tone_mapping_struct_s *seiToneMapping;
#endif

    int              iNumOfSlicesAllocated;
    int              iSliceNumOfCurrPic;
    struct slice_t **ppSliceList;
    struct slice_t  *pNextSlice;
    int              newframe;

    struct nalu_t *nalu;

    DecodedPicList *pDecOuputPic;
    pps_t *pNextPPS;
    bool first_sps; // use only for print first sps


    struct frame_store *out_buffer;




  pps_t *active_pps;
  sps_t *active_sps;
  sps_t SeqParSet[MAXSPS];
  pps_t PicParSet[MAXPPS];


#if (MVC_EXTENSION_ENABLE)
    sub_sps_t *active_subset_sps;
    sub_sps_t SubsetSeqParSet[MAXSPS];
    int last_pic_width_in_mbs_minus1;
    int last_pic_height_in_map_units_minus1;
    int last_max_dec_frame_buffering;
    int last_profile_idc;
#endif

  struct sei_params        *p_SEI;

  int number;                                 //!< frame number

  //current picture property;
  unsigned int num_dec_mb;
    char  *intra_block;
    char  *intra_block_JV[MAX_PLANE];

  int type;                                   //!< image type INTER/INTRA

  byte **ipredmode;                  //!< prediction type [90][74]
  byte **ipredmode_JV[MAX_PLANE];
  byte ****nz_coeff;
  int **siblock;
  int **siblock_JV[MAX_PLANE];
  BlockPos *PicPos;

  int structure;                     //!< Identify picture structure type

  struct macroblock_dec *mb_data;               //!< array containing all MBs of a whole frame
  struct macroblock_dec *mb_data_JV[MAX_PLANE]; //!< mb_data to be used for 4:4:4 independent mode
  int ChromaArrayType;

  // picture error concealment
  // concealment_head points to first node in list, concealment_end points to
  // last node in list. Initialize both to NULL, meaning no nodes in list yet
  struct concealment_node *concealment_head;
  struct concealment_node *concealment_end;

  unsigned int pre_frame_num;           //!< store the frame_num in the last decoded slice. For detecting gap in frame_num.
  int non_conforming_stream;

  // ////////////////////////
  // for POC mode 0:
  signed   int PrevPicOrderCntMsb;
  unsigned int PrevPicOrderCntLsb;

  // for POC mode 1:
  signed int ExpectedPicOrderCnt, PicOrderCntCycleCnt, FrameNumInPicOrderCntCycle;
  unsigned int PreviousFrameNum, FrameNumOffset;
  int ExpectedDeltaPerPicOrderCntCycle;
  int ThisPOC;
  int PreviousFrameNumOffset;
  // /////////////////////////

  unsigned int PicHeightInMbs;
  unsigned int PicSizeInMbs;

  int no_output_of_prior_pics_flag;

  int last_has_mmco_5;
  int last_pic_bottom_field;

  int idr_psnr_number;
  int psnr_number;

  // Timing related variables
  TIME_T start_time;
  TIME_T end_time;

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

  byte *buf;
  byte *ibuf;

  // Redundant slices. Should be moved to another structure and allocated only if extended profile
  unsigned int previous_frame_num; //!< frame number of previous slice
  //!< non-zero: i-th previous frame is correct
  int Is_primary_correct;          //!< if primary frame is correct, 0: incorrect
  int Is_redundant_correct;        //!< if redundant frame is correct, 0:incorrect

  // Time 
  int64 tot_time;

  // B pictures
  int  Bframe_ctr;
  int  frame_no;

  int  g_nFrame;

  // global picture format dependent buffers, memory allocation in decod.c
  imgpel **imgY_ref;                              //!< reference frame find snr
  imgpel ***imgUV_ref;

  int *qp_per_matrix;
  int *qp_rem_matrix;

  struct frame_store *last_out_fs;
  int pocs_in_dpb[100];

  struct storable_picture *dec_picture;
  struct storable_picture *dec_picture_JV[MAX_PLANE];  //!< dec_picture to be used during 4:4:4 independent mode decoding
  struct storable_picture *no_reference_picture; //!< dummy storable picture for recovery point

  // Error parameters
  struct object_buffer  *erc_object_list;
  struct ercVariables_s *erc_errorVar;

  int erc_mvperMB;
  struct video_par *erc_img;

  struct storable_picture *pending_output;
  int    pending_output_state;
  int    recovery_flag;



  // report
  char cslice_type[9];  
  // FMO
  int *MbToSliceGroupMap;
  int *MapUnitToSliceGroupMap;
  int  NumberOfSliceGroups;    // the number of slice groups -1 (0 == scan order, 7 == maximum)

  void (*buf2img)          (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int o_size_x, int o_size_y, int symbol_size_in_bytes, int bitshift);
  void (*getNeighbour)     (struct macroblock_dec *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
  void (*get_mb_block_pos) (BlockPos *PicPos, int mb_addr, short *x, short *y);
  void (*GetStrengthVer)   (struct macroblock_dec *MbQ, int edge, int mvlimit, struct storable_picture *p);
  void (*GetStrengthHor)   (struct macroblock_dec *MbQ, int edge, int mvlimit, struct storable_picture *p);
  void (*EdgeLoopLumaVer)  (ColorPlane pl, imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, struct storable_picture *p);
  void (*EdgeLoopLumaHor)  (ColorPlane pl, imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, struct storable_picture *p);
  void (*EdgeLoopChromaVer)(imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, int uv, struct storable_picture *p);
  void (*EdgeLoopChromaHor)(imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, int uv, struct storable_picture *p);
  void (*img2buf)          (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);

  ImageData tempData3;
  int iDeblockMode;  //0: deblock in picture, 1: deblock in slice;
  int iLumaPadX;
  int iLumaPadY;
  int iChromaPadX;
  int iChromaPadY;
  //control;
  int bDeblockEnable;
  int iPostProcess;
  int bFrameInit;
  int last_dec_poc;
  int last_dec_view_id;
  int last_dec_layer_id;
  int dpb_layer_id;

/******************* deprecative variables; ***************************************/
  int width;
  int height;
  int width_cr;                               //!< width chroma  
  int height_cr;                              //!< height chroma
  // Fidelity Range Extensions Stuff
  int pic_unit_bitsize_on_disk;
  short bitdepth_luma;
  short bitdepth_chroma;
  int bitdepth_scale[2];
  int bitdepth_luma_qp_scale;
  int bitdepth_chroma_qp_scale;
  unsigned int dc_pred_value_comp[MAX_PLANE]; //!< component value for DC prediction (depends on component pel bit depth)
  int max_pel_value_comp[MAX_PLANE];       //!< max value that one picture element (pixel) can take (depends on pic_unit_bitdepth)

  int separate_colour_plane_flag;
  int pic_unit_size_on_disk;

  int profile_idc;
  int yuv_format;
  int lossless_qpprime_flag;
  int num_blk8x8_uv;
  int num_uv_blocks;
  int num_cdc_coeff;
  int mb_cr_size_x;
  int mb_cr_size_y;
  int mb_cr_size_x_blk;
  int mb_cr_size_y_blk;
  int mb_cr_size;
  int mb_size[3][2];                         //!< component macroblock dimensions
  int mb_size_blk[3][2];                     //!< component macroblock dimensions 
  int mb_size_shift[3][2];
  int subpel_x;
  int subpel_y;
  int shiftpel_x;
  int shiftpel_y;
  int total_scale;
  int max_frame_num;

  unsigned int PicWidthInMbs;
  unsigned int PicHeightInMapUnits;
  unsigned int FrameHeightInMbs;
  unsigned int FrameSizeInMbs;
  unsigned int oldFrameSizeInMbs;
  int max_vmv_r;                             //!< maximum vertical motion vector range in luma quarter frame pixel units for the current level_idc
  //int max_mb_vmv_r;                        //!< maximum vertical motion vector range in luma quarter pixel units for the current level_idc
/******************* end deprecative variables; ***************************************/

  struct dec_stat_parameters *dec_stats;
} VideoParameters;


typedef struct decoder_params
{
  InputParameters   *p_Inp;          //!< Input Parameters
  VideoParameters   *p_Vid;          //!< Image Parameters
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


#ifdef __cplusplus
}
#endif

#endif

