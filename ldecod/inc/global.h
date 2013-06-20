
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
#include "types.h"
#include "io_image.h"
#include "frame.h"

#include "parser.h"
#include "parset.h"


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
    DEC_OPENED = 0,
    DEC_STOPPED
} DecoderStatus_e;

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

#if (MVC_EXTENSION_ENABLE)
typedef struct nalunitheadermvcext_tag {
    unsigned int non_idr_flag;
    unsigned int priority_id;
    unsigned int view_id;
    unsigned int temporal_id;
    unsigned int anchor_pic_flag;
    unsigned int inter_view_flag;
    unsigned int reserved_one_bit;
    unsigned int iPrefixNALU;
} NALUnitHeaderMVCExt_t;
#endif

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
  seq_parameter_set_rbsp_t      *p_SPS;
  struct decoded_picture_buffer *p_Dpb;
} LayerParameters;

// video parameters
typedef struct video_par {
  struct inp_par      *p_Inp;
  pic_parameter_set_rbsp_t *active_pps;
  seq_parameter_set_rbsp_t *active_sps;
  seq_parameter_set_rbsp_t SeqParSet[MAXSPS];
  pic_parameter_set_rbsp_t PicParSet[MAXPPS];
  struct decoded_picture_buffer *p_Dpb_layer[MAX_NUM_DPB_LAYERS];
  CodingParameters *p_EncodePar[MAX_NUM_DPB_LAYERS];
  LayerParameters *p_LayerPar[MAX_NUM_DPB_LAYERS];

#if (MVC_EXTENSION_ENABLE)
  subset_seq_parameter_set_rbsp_t *active_subset_sps;
  //int svc_extension_flag;
  subset_seq_parameter_set_rbsp_t SubsetSeqParSet[MAXSPS];
  int last_pic_width_in_mbs_minus1;
  int last_pic_height_in_map_units_minus1;
  int last_max_dec_frame_buffering;
  int last_profile_idc;
#endif

  struct sei_params        *p_SEI;

  struct old_slice_par *old_slice;
  struct snr_par       *snr;
  int number;                                 //!< frame number
  
  //current picture property;
  unsigned int num_dec_mb;
  int iSliceNumOfCurrPic;
  int iNumOfSlicesAllocated;
  int iNumOfSlicesDecoded;
  struct slice_t **ppSliceList;
  char  *intra_block;
  char  *intra_block_JV[MAX_PLANE];
  //int qp;                                     //!< quant for the current frame

  //int sp_switch;                              //!< 1 for switching sp, 0 for normal sp  
  int type;                                   //!< image type INTER/INTRA

  byte **ipredmode;                  //!< prediction type [90][74]
  byte **ipredmode_JV[MAX_PLANE];
  byte ****nz_coeff;
  int **siblock;
  int **siblock_JV[MAX_PLANE];
  BlockPos *PicPos;

  int newframe;
  int structure;                     //!< Identify picture structure type

  //Slice      *currentSlice;          //!< pointer to current Slice data struct
  struct slice_t      *pNextSlice;             //!< pointer to first Slice of next picture;
  struct macroblock_dec *mb_data;               //!< array containing all MBs of a whole frame
  struct macroblock_dec *mb_data_JV[MAX_PLANE]; //!< mb_data to be used for 4:4:4 independent mode
  //int colour_plane_id;               //!< colour_plane_id of the current coded slice
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

  Boolean first_sps;
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

  // files
  int p_out;                       //!< file descriptor to output YUV file
#if (MVC_EXTENSION_ENABLE)
  int p_out_mvc[MAX_VIEW_NUM];     //!< file descriptor to output YUV file for MVC
#endif
  int p_ref;                       //!< pointer to input original reference YUV file file

  // B pictures
  int  Bframe_ctr;
  int  frame_no;

  int  g_nFrame;
  Boolean global_init_done[2];

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

  struct frame_store *out_buffer;

  struct storable_picture *pending_output;
  int    pending_output_state;
  int    recovery_flag;


    struct bitstream_t *bitstream;


  // report
  char cslice_type[9];  
  // FMO
  int *MbToSliceGroupMap;
  int *MapUnitToSliceGroupMap;
  int  NumberOfSliceGroups;    // the number of slice groups -1 (0 == scan order, 7 == maximum)

#if (ENABLE_OUTPUT_TONEMAPPING)
  struct tone_mapping_struct_s *seiToneMapping;
#endif

  void (*buf2img)          (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int o_size_x, int o_size_y, int symbol_size_in_bytes, int bitshift);
  void (*getNeighbour)     (struct macroblock_dec *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
  void (*get_mb_block_pos) (BlockPos *PicPos, int mb_addr, short *x, short *y);
  void (*GetStrengthVer)   (struct macroblock_dec *MbQ, int edge, int mvlimit, struct storable_picture *p);
  void (*GetStrengthHor)   (struct macroblock_dec *MbQ, int edge, int mvlimit, struct storable_picture *p);
  void (*EdgeLoopLumaVer)  (ColorPlane pl, imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge);
  void (*EdgeLoopLumaHor)  (ColorPlane pl, imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, struct storable_picture *p);
  void (*EdgeLoopChromaVer)(imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, int uv, struct storable_picture *p);
  void (*EdgeLoopChromaHor)(imgpel** Img, byte *Strength, struct macroblock_dec *MbQ, int edge, int uv, struct storable_picture *p);
  void (*img2buf)          (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride);

  ImageData tempData3;
  DecodedPicList *pDecOuputPic;
  int iDeblockMode;  //0: deblock in picture, 1: deblock in slice;
  struct nalu_t *nalu;
  int iLumaPadX;
  int iLumaPadY;
  int iChromaPadX;
  int iChromaPadY;
  //control;
  int bDeblockEnable;
  int iPostProcess;
  int bFrameInit;
  pic_parameter_set_rbsp_t *pNextPPS;
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

  int  ProcessInput;
  int  enable_32_pulldown;
#if (MVC_EXTENSION_ENABLE)
  int  DecodeAllLayers;
#endif

  // picture error concealment
  int conceal_mode;
  int ref_poc_gap;
  int poc_gap;


  // dummy for encoder
  int start_frame;

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

typedef struct decoder_params
{
  InputParameters   *p_Inp;          //!< Input Parameters
  VideoParameters   *p_Vid;          //!< Image Parameters
  int64              bufferSize;     //!< buffersize for tiff reads (not currently supported)
  int                UsedBits;      // for internal statistics, is adjusted by read_se_v, read_ue_v, read_u_1
  FILE              *p_trace;        //!< Trace file
  int                bitcounter;
} DecoderParams;

extern DecoderParams  *p_Dec;

// prototypes
extern void error(const char *text, int code);

// dynamic mem allocation
extern int  init_global_buffers( VideoParameters *p_Vid, int layer_id );
extern void free_global_buffers( VideoParameters *p_Vid);
extern void free_layer_buffers( VideoParameters *p_Vid, int layer_id );

extern void tracebits (const char *trace_str, int len, int info, int value1);
extern void tracebits2(const char *trace_str, int len, int info);

extern unsigned CeilLog2   ( unsigned uiVal);
extern unsigned CeilLog2_sf( unsigned uiVal);


extern void FreeDecPicList ( DecodedPicList *pDecPicList );
extern void ClearDecPicList( VideoParameters *p_Vid );
extern DecodedPicList *get_one_avail_dec_pic_from_list(DecodedPicList *pDecPicList, int b3D, int view_id);
extern struct slice_t *malloc_slice( InputParameters *p_Inp, VideoParameters *p_Vid );
extern void copy_slice_info (struct slice_t *currSlice, OldSliceParams *p_old_slice );
extern void OpenOutputFiles(VideoParameters *p_Vid, int view0_id, int view1_id);
extern void set_global_coding_par(VideoParameters *p_Vid, CodingParameters *cps);

static inline int is_FREXT_profile(unsigned int profile_idc) 
{
  return ( profile_idc >= FREXT_HP || profile_idc == FREXT_CAVLC444 );
}
static inline int HI_intra_only_profile(unsigned int profile_idc, Boolean constraint_set3_flag)
{
  return ( ((profile_idc >= FREXT_Hi10P) && constraint_set3_flag) || (profile_idc == FREXT_CAVLC444) );
}
static inline int is_BL_profile(unsigned int profile_idc) 
{
  return ( profile_idc == FREXT_CAVLC444 || profile_idc == BASELINE || profile_idc == MAIN || profile_idc == EXTENDED ||
           profile_idc == FREXT_HP || profile_idc == FREXT_Hi10P || profile_idc == FREXT_Hi422 || profile_idc == FREXT_Hi444);
}
static inline int is_EL_profile(unsigned int profile_idc) 
{
  return ( (profile_idc == MVC_HIGH) || (profile_idc == STEREO_HIGH)
           );
}

static inline int is_MVC_profile(unsigned int profile_idc)
{
  return ( (0)
#if (MVC_EXTENSION_ENABLE)
  || (profile_idc == MVC_HIGH) || (profile_idc == STEREO_HIGH)
#endif
  );
}

#ifdef __cplusplus
}
#endif

#endif

