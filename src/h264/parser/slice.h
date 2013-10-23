#ifndef _SLICE_H_
#define _SLICE_H_


#include "global.h"
#include "sets.h"
#include "dpb.h"
#include "macroblock.h"

#include "neighbour.h"
#include "interpret.h"
#include "decoder.h"


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
    struct ref_pic_list_modification_t {
        uint8_t     modification_of_pic_nums_idc;                     // ue(v)
        uint32_t    abs_diff_pic_num_minus1;                          // ue(v)
        uint32_t    long_term_pic_num;                                // ue(v)
        uint32_t    abs_diff_view_idx_minus1;                         // ue(v)
    };
    using ref_pic_list_modification_v = std::vector<ref_pic_list_modification_t>;
    ref_pic_list_modification_v ref_pic_list_modifications[2];

    uint8_t     luma_log2_weight_denom;                               // ue(v)
    uint8_t     chroma_log2_weight_denom;                             // ue(v)
    struct pred_weight_t {
        bool        weight_flag;                                      // u(1)
        int8_t      weight;                                           // se(v)
        int8_t      offset;                                           // se(v)
    };
    using pred_weight_v = std::vector<pred_weight_t>;
    pred_weight_v pred_weight_l[2][3];

    bool        no_output_of_prior_pics_flag;                         // u(1)
    bool        long_term_reference_flag;                             // u(1)
    bool        adaptive_ref_pic_marking_mode_flag;                   // u(1)
    struct adaptive_ref_pic_marking_t {
        uint32_t    memory_management_control_operation;              // ue(v)
        uint32_t    difference_of_pic_nums_minus1;                    // ue(v)
        uint32_t    long_term_pic_num;                                // ue(v)
        uint32_t    long_term_frame_idx;                              // ue(v)
        uint32_t    max_long_term_frame_idx_plus1;                    // ue(v)
    };
    using adaptive_ref_pic_marking_v = std::vector<adaptive_ref_pic_marking_t>;
    adaptive_ref_pic_marking_v adaptive_ref_pic_markings;

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

    int32_t     PicOrderCntMsb;
    int32_t     FrameNumOffset;
    int32_t     TopFieldOrderCnt;
    int32_t     BottomFieldOrderCnt;
    int32_t     PicOrderCnt;
};

using shr_t = slice_header_t;


struct slice_t {
    VideoParameters* p_Vid;
    dpb_t*      p_Dpb;
    pps_t*      active_pps;
    sps_t*      active_sps;

    bool        forbidden_zero_bit;                                   // f(1)
    uint8_t     nal_ref_idc;                                          // u(2)
    uint8_t     nal_unit_type;                                        // u(5)

    bool        mvc_extension_flag;
    bool        svc_extension_flag;                                   // u(1)
    bool        non_idr_flag;                                         // u(1)
    uint8_t     priority_id;                                          // u(6)
//    uint16_t    view_id;                                              // u(10)
    uint8_t     temporal_id;                                          // u(3)
    bool        anchor_pic_flag;                                      // u(1)
    bool        inter_view_flag;                                      // u(1)
    bool        reserved_one_bit;                                     // u(1)

    bool        IdrPicFlag;
    shr_t       header;

    using uint8_v = std::vector<uint8_t>;
    uint8_v     MbToSliceGroupMap;
    //slice header information;
    int               ref_flag[17]; //!< 0: i-th previous frame is incorrect
    char              RefPicSize[2];
    storable_picture* RefPicList[2][33];

    unsigned    num_dec_mb;
    short       current_slice_nr;


    int         layer_id;
    int         view_id;


    int         listinterviewidx0;
    int         listinterviewidx1;
    pic_t**     fs_listinterview0;
    pic_t**     fs_listinterview1;

    int         dpB_NotPresent;    //!< non-zero, if data partition B is lost
    int         dpC_NotPresent;    //!< non-zero, if data partition C is lost

    px_t***     mb_pred; // IntraPrediction()

    Neighbour   neighbour;
    Parser      parser;
    Decoder     decoder;

    // for signalling to the neighbour logic that this is a deblocker call

    storable_picture* dec_picture;

    slice_t();
    ~slice_t();

    void        init_lists    ();
    void        init_ref_lists();

    void        decode_poc();
    void        init_slice_group_map();

    int         NextMbAddress(int n);
    void        init();
    void        decode();

    bool        operator!=(const slice_t& slice);

protected:
    void        FmoGenerateType0MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
    void        FmoGenerateType1MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
    void        FmoGenerateType2MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
    void        FmoGenerateType3MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
    void        FmoGenerateType4MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
    void        FmoGenerateType5MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
    void        FmoGenerateType6MapUnitMap(uint8_v& mapUnitToSliceGroupMap);
};


#endif /* _SLICE_H_ */
