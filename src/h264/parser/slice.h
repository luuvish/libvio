#ifndef _SLICE_H_
#define _SLICE_H_


#include "global.h"
#include "sets.h"
#include "dpb.h"
#include "macroblock.h"

#include "neighbour.h"
#include "interpret.h"
#include "decoder.h"

#define MAX_NUM_DECSLICES      16

#define MAX_REFERENCE_PICTURES 32               //!< H.264 allows 32 fields

#define MAX_NUM_REF_IDX 32

namespace vio { namespace h264 {
struct cabac_contexts_t;
}}

using vio::h264::cabac_contexts_t;
using vio::h264::mb_t;

using vio::h264::Neighbour;
using vio::h264::Parser;
using vio::h264::Decoder;


enum {
    P_slice  = 0, // or 5
    B_slice  = 1, // or 6
    I_slice  = 2, // or 7
    SP_slice = 3, // or 8
    SI_slice = 4  // or 9
};


struct decoded_reference_picture_marking_t;
using drpm_t = decoded_reference_picture_marking_t;

struct decoded_reference_picture_marking_t {
    uint32_t    memory_management_control_operation;                  // ue(v)
    uint32_t    difference_of_pic_nums_minus1;                        // ue(v)
    uint32_t    long_term_pic_num;                                    // ue(v)
    uint32_t    long_term_frame_idx;                                  // ue(v)
    uint32_t    max_long_term_frame_idx_plus1;                        // ue(v)
    drpm_t*     Next;
};

struct VideoParameters;

struct slice_header_t {
    uint32_t    first_mb_in_slice;                                    // ue(v)
    uint8_t     slice_type;                                           // ue(v)
    uint8_t     pic_parameter_set_id;                                 // ue(v)
    uint8_t     colour_plane_id;                                      // u(2)
    uint32_t    frame_num;                                            // u(v)
    bool        field_pic_flag;                                       // u(1)
    bool        bottom_field_flag;                                    // u(1)
    uint16_t    idr_pic_id;                                           // ue(v)
    uint32_t    pic_order_cnt_lsb;                                    // u(v)
    int32_t     delta_pic_order_cnt_bottom;                           // se(v)
    int32_t     delta_pic_order_cnt[2];                               // se(v)
    uint8_t     redundant_pic_cnt;                                    // ue(v)
    bool        direct_spatial_mv_pred_flag;                          // u(1)
    bool        num_ref_idx_active_override_flag;                     // u(1)
    uint8_t     num_ref_idx_l0_active_minus1;                         // ue(v)
    uint8_t     num_ref_idx_l1_active_minus1;                         // ue(v)

    bool        ref_pic_list_modification_flag_l0;                    // u(1)
    bool        ref_pic_list_modification_flag_l1;                    // u(1)
    uint8_t     modification_of_pic_nums_idc[2][32+1];                // ue(v)
    uint32_t    abs_diff_pic_num_minus1     [2][32+1];                // ue(v)
    uint32_t    long_term_pic_num           [2][32+1];                // ue(v)
#if (MVC_EXTENSION_ENABLE)
    uint32_t    abs_diff_view_idx_minus1    [2][32+1];                // ue(v)
#endif

    uint8_t     luma_log2_weight_denom;                               // ue(v)
    uint8_t     chroma_log2_weight_denom;                             // ue(v)
    bool        luma_weight_l0_flag  [MAX_NUM_REF_IDX];               // u(1)
    int8_t      luma_weight_l0       [MAX_NUM_REF_IDX];               // se(v)
    int8_t      luma_offset_l0       [MAX_NUM_REF_IDX];               // se(v)
    bool        chroma_weight_l0_flag[MAX_NUM_REF_IDX];               // u(1)
    int8_t      chroma_weight_l0     [MAX_NUM_REF_IDX][2];            // se(v)
    int8_t      chroma_offset_l0     [MAX_NUM_REF_IDX][2];            // se(v)
    bool        luma_weight_l1_flag  [MAX_NUM_REF_IDX];               // u(1)
    int8_t      luma_weight_l1       [MAX_NUM_REF_IDX];               // se(v)
    int8_t      luma_offset_l1       [MAX_NUM_REF_IDX];               // se(v)
    bool        chroma_weight_l1_flag[MAX_NUM_REF_IDX];               // u(1)
    int8_t      chroma_weight_l1     [MAX_NUM_REF_IDX][2];            // se(v)
    int8_t      chroma_offset_l1     [MAX_NUM_REF_IDX][2];            // se(v)

    bool        no_output_of_prior_pics_flag;                         // u(1)
    bool        long_term_reference_flag;                             // u(1)
    bool        adaptive_ref_pic_marking_mode_flag;                   // u(1)
    //struct decoded_reference_picture_marking_t;
    //using drpm_t = decoded_reference_picture_marking_t;
    //struct decoded_reference_picture_marking_t {
    //    uint32_t    memory_management_control_operation;              // ue(v)
    //    uint32_t    difference_of_pic_nums_minus1;                    // ue(v)
    //    uint32_t    long_term_pic_num;                                // ue(v)
    //    uint32_t    long_term_frame_idx;                              // ue(v)
    //    uint32_t    max_long_term_frame_idx_plus1;                    // ue(v)
    //    drpm_t*     Next;
    //};
    drpm_t*     dec_ref_pic_marking_buffer;

    uint8_t     cabac_init_idc;                                       // ue(v)
    int8_t      slice_qp_delta;                                       // se(v)
    bool        sp_for_switch_flag;                                   // u(1)
    int8_t      slice_qs_delta;                                       // se(v)
    uint8_t     disable_deblocking_filter_idc;                        // ue(v)
    int8_t      slice_alpha_c0_offset_div2;                           // se(v)
    int8_t      slice_beta_offset_div2;                               // se(v)
    uint32_t    slice_group_change_cycle;                             // u(v)

    PictureStructure    structure;
    bool        MbaffFrameFlag;
    uint32_t    PicHeightInMbs;
    uint32_t    PicHeightInSampleL;
    uint32_t    PicHeightInSampleC;
    uint32_t    PicSizeInMbs;
    uint32_t    MaxPicNum;
    uint32_t    CurrPicNum;
    int8_t      SliceQpY;
    int8_t      QsY;
    int8_t      FilterOffsetA;
    int8_t      FilterOffsetB;
    uint32_t    MapUnitsInSliceGroup0;
    uint8_t*    MbToSliceGroupMap;
    uint8_t*    MapUnitToSliceGroupMap;

    int32_t     PicOrderCntMsb;
    int32_t     FrameNumOffset;
    int32_t     TopFieldOrderCnt;
    int32_t     BottomFieldOrderCnt;
    int32_t     PicOrderCnt;
};

using shr_t = slice_header_t;


#if (MVC_EXTENSION_ENABLE)
struct NALUnitHeaderMVCExt_t {
    unsigned int non_idr_flag;
    unsigned int priority_id;
    unsigned int view_id;
    unsigned int temporal_id;
    unsigned int anchor_pic_flag;
    unsigned int inter_view_flag;
    unsigned int reserved_one_bit;
    unsigned int iPrefixNALU;
};
#endif


struct slice_t {
    VideoParameters* p_Vid;
    pps_t*      active_pps;
    sps_t*      active_sps;
    int         svc_extension_flag;

    // dpb pointer
    dpb_t*      p_Dpb;

    //slice property;
    bool        idr_flag;
    int         nal_ref_idc;
    uint8_t     nal_unit_type;

    shr_t       header;

    unsigned    num_dec_mb;
    short       current_slice_nr;

    int         layer_id;
#if (MVC_EXTENSION_ENABLE)
    int         view_id;
    int         inter_view_flag;
    int         anchor_pic_flag;

    NALUnitHeaderMVCExt_t NaluHeaderMVCExt;
#endif

    //slice header information;
    int               ref_flag[17]; //!< 0: i-th previous frame is incorrect
    char              RefPicSize[2];
    storable_picture* RefPicList[2][33];

#if (MVC_EXTENSION_ENABLE)
    int         listinterviewidx0;
    int         listinterviewidx1;
    pic_t**     fs_listinterview0;
    pic_t**     fs_listinterview1;
#endif


    int         dpB_NotPresent;    //!< non-zero, if data partition B is lost
    int         dpC_NotPresent;    //!< non-zero, if data partition C is lost

    px_t***     mb_pred; // IntraPrediction()

    Neighbour   neighbour;
    Parser      parser;
    Decoder     decoder;

    // for signalling to the neighbour logic that this is a deblocker call

    int           erc_mvperMB;
    storable_picture* dec_picture;

    slice_t();
    ~slice_t();

    void        update_pic_num();
    void        init_lists    ();
    void        init_ref_lists();

    void        fmo_init();
    void        fmo_close();
    int         FmoGetNextMBNr(int CurrentMbNr);

    void        FmoGenerateType0MapUnitMap();
    void        FmoGenerateType1MapUnitMap();
    void        FmoGenerateType2MapUnitMap();
    void        FmoGenerateType3MapUnitMap();
    void        FmoGenerateType4MapUnitMap();
    void        FmoGenerateType5MapUnitMap();
    void        FmoGenerateType6MapUnitMap();
    void        FmoGenerateMapUnitToSliceGroupMap();
    void        FmoGenerateMbToSliceGroupMap();

    void        init();
    void        decode();

    bool        operator!=(const slice_t& slice);
};


void decode_poc(VideoParameters *p_Vid, slice_t *pSlice);

#endif /* _SLICE_H_ */
