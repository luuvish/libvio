#ifndef _GLOBAL_H_
#define _GLOBAL_H_


#include <cstdint>
#include <vector>

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#include "defines.h"

#include "parset.h"
#include "bitstream.h"
#include "macroblock.h"

#include "image_data.h"


using vio::h264::bitstream_t;
using vio::h264::mb_t;


enum {
    EOS = 1, //!< End Of Sequence
    SOP = 2, //!< Start Of Picture
    SOS = 3, //!< Start Of slice_t
};

enum ColorPlane {
    PLANE_Y = 0,
    PLANE_U = 1,
    PLANE_V = 2
};

enum ColorFormat {
    YUV400 = 0,
    YUV420 = 1,
    YUV422 = 2,
    YUV444 = 3
};

enum PictureStructure {
    FRAME,
    TOP_FIELD,
    BOTTOM_FIELD
};


struct pic_motion_params_old;
struct pic_motion_params;

struct slice_t;
struct macroblock_t;

struct decoded_picture_buffer_t;

struct InputParameters;
struct VideoParameters;
struct DecodedPicList;

struct SNRParameters;

struct tone_mapping_struct_s;
struct picture_t;
using pic_t = picture_t;
struct sei_params;
struct concealment_node;
struct object_buffer;
struct ercVariables_s;
struct storable_picture;

struct CodingParameters {
    int layer_id;

    mb_t* mb_data;
    mb_t* mb_data_JV[3];
};

struct LayerParameters {
    int               layer_id;
    VideoParameters*  p_Vid;
    CodingParameters* p_Cps;
    sps_t*            p_SPS;
    decoded_picture_buffer_t* p_Dpb;
};

#define MAX_NUM_DPB_LAYERS 2

// video parameters
struct VideoParameters {
    InputParameters* p_Inp;

    int         p_out; // FILE handle for YUV output file
#if (MVC_EXTENSION_ENABLE)
    int         p_out_mvc[MAX_VIEW_NUM];
#endif
    int         p_ref;

    bitstream_t bitstream;



    SNRParameters*  snr;

    decoded_picture_buffer_t* p_Dpb_layer[MAX_NUM_DPB_LAYERS];
    CodingParameters* p_EncodePar[MAX_NUM_DPB_LAYERS];
    LayerParameters*  p_LayerPar [MAX_NUM_DPB_LAYERS];

    bool        global_init_done[2];

    tone_mapping_struct_s *seiToneMapping;

    int         iSliceNumOfCurrPic;
    std::vector<slice_t*> ppSliceList;
    slice_t*    pNextSlice;
    int         newframe;

    nal_unit_t*     nalu;

    struct {
        uint8_t* pY;
        uint8_t* pU;
        uint8_t* pV;
    } pDecOuputPic;
    pps_t*      pNextPPS;


    pic_t*      out_buffer;


    int32_t     prevPicOrderCntMsb;
    uint32_t    prevPicOrderCntLsb;
    uint32_t    prevFrameNum;
    uint32_t    prevFrameNumOffset;
    int         last_has_mmco_5;
    int         last_pic_bottom_field;

    pps_t*      active_pps;
    sps_t*      active_sps;
    sps_t       SeqParSet[MAXSPS];
    pps_t       PicParSet[MAXPPS];

#if (MVC_EXTENSION_ENABLE)
    sub_sps_t*  active_subset_sps;
    sub_sps_t   SubsetSeqParSet[MAXSPS];
#endif

    int         number;                                 //!< frame number

    //current picture property;
    unsigned    num_dec_mb;

    int         profile_idc;
    int         type;                                   //!< image type INTER/INTRA
    int         structure;                     //!< Identify picture structure type

    // global picture format dependent buffers, memory allocation in decod.c
    mb_t*       mb_data;               //!< array containing all MBs of a whole frame
    mb_t*       mb_data_JV[3]; //!< mb_data to be used for 4:4:4 independent mode


    // picture error concealment
    // concealment_head points to first node in list, concealment_end points to
    // last node in list. Initialize both to NULL, meaning no nodes in list yet
    concealment_node* concealment_head;
    concealment_node* concealment_end;

    // for POC mode 1:
    int         PicOrderCnt;

    int         no_output_of_prior_pics_flag;

    // picture error concealment
    int         last_ref_pic_poc;
    int         conceal_mode;
    int         earlier_missing_poc;
    unsigned    frame_to_conceal;
    int         IDR_concealment_flag;
    int         conceal_slice_type;

    // SEI
    bool        recovery_point;
    uint32_t    recovery_frame_cnt;

    // random access point decoding
    bool        recovery_flag;
    bool        recovery_point_found;
    int         recovery_frame_num;
    int         recovery_poc;
    bool        non_conforming_stream;

    unsigned    pre_frame_num;           //!< store the frame_num in the last decoded slice. For detecting gap in frame_num.

    // Redundant slices. Should be moved to another structure and allocated only if extended profile
    unsigned    previous_frame_num; //!< frame number of previous slice
    //!< non-zero: i-th previous frame is correct
    int         Is_primary_correct;          //!< if primary frame is correct, 0: incorrect
    int         Is_redundant_correct;        //!< if redundant frame is correct, 0:incorrect

    pic_t*      last_out_fs;
    int         pocs_in_dpb[100];

    storable_picture* dec_picture;
    storable_picture* dec_picture_JV[3];  //!< dec_picture to be used during 4:4:4 independent mode decoding
    storable_picture* no_reference_picture; //!< dummy storable picture for recovery point

    // Error parameters
    object_buffer*  erc_object_list;
    ercVariables_s* erc_errorVar;

    int         erc_mvperMB;
    VideoParameters* erc_img;

    ImageData   tempData3;
    //control;
    int         last_dec_layer_id;
    int         dpb_layer_id;

    VideoParameters();
    ~VideoParameters();

    void OpenOutputFiles(int view0_id, int view1_id);

    void calculate_frame_no(storable_picture *p);
    void status(storable_picture** dec_picture);
    void report();
};

extern void error(int code, const char* format, ...);

// dynamic mem allocation
extern void free_layer_buffers(VideoParameters *p_Vid, int layer_id);


void init_picture(slice_t* currSlice);
void init_picture_decoding(VideoParameters *p_Vid);
void exit_picture(VideoParameters *p_Vid);

#if (MVC_EXTENSION_ENABLE)
extern int GetVOIdx(VideoParameters *p_Vid, int iViewId);
#endif


#endif
