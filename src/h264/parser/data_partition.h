/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : data_partition.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _DATA_PARTITION_H_
#define _DATA_PARTITION_H_

#include <cstdint>


namespace vio  {
namespace h264 {


struct nal_unit_t {
    static const uint32_t MAX_NAL_UNIT_SIZE = 8000000;

    enum {
        NALU_TYPE_SLICE    =  1,
        NALU_TYPE_DPA      =  2,
        NALU_TYPE_DPB      =  3,
        NALU_TYPE_DPC      =  4,
        NALU_TYPE_IDR      =  5,
        NALU_TYPE_SEI      =  6,
        NALU_TYPE_SPS      =  7,
        NALU_TYPE_PPS      =  8,
        NALU_TYPE_AUD      =  9,
        NALU_TYPE_EOSEQ    = 10,
        NALU_TYPE_EOSTREAM = 11,
        NALU_TYPE_FILL     = 12,
        NALU_TYPE_SPS_EXT  = 13,
#if (MVC_EXTENSION_ENABLE)
        NALU_TYPE_PREFIX   = 14,
        NALU_TYPE_SUB_SPS  = 15,
        NALU_TYPE_SLC_EXT  = 20,
        NALU_TYPE_VDRD     = 24
#endif
    };

    uint16_t    lost_packets;
    uint32_t    max_size;

    uint32_t    num_bytes_in_nal_unit;
    uint32_t    num_bytes_in_rbsp;
    uint8_t*    rbsp_byte;

    bool        forbidden_zero_bit;                                   // f(1)
    uint8_t     nal_ref_idc;                                          // u(2)
    uint8_t     nal_unit_type;                                        // u(5)

#if (MVC_EXTENSION_ENABLE)
    bool        svc_extension_flag;                                   // u(1)
    bool        non_idr_flag;                                         // u(1)
    uint8_t     priority_id;                                          // u(6)
    uint16_t    view_id;                                              // u(10)
    uint8_t     temporal_id;                                          // u(3)
    bool        anchor_pic_flag;                                      // u(1)
    bool        inter_view_flag;                                      // u(1)
    bool        reserved_one_bit;                                     // u(1)
#endif

    nal_unit_t(uint32_t size=MAX_NAL_UNIT_SIZE) :
        max_size { size }, rbsp_byte { new uint8_t[size] } {}

    ~nal_unit_t() {
        if (this->rbsp_byte)
            delete []this->rbsp_byte;
    }
};



struct data_partition_t : public nal_unit_t {
    int         frame_bitoffset;

                data_partition_t(uint32_t size=MAX_NAL_UNIT_SIZE);
                data_partition_t(const nal_unit_t& nal);

    data_partition_t& operator=(const nal_unit_t& nal);

    bool        byte_aligned            (void);
    bool        more_data_in_byte_stream(void);
    bool        more_rbsp_data          (void);
    bool        more_rbsp_trailing_data (void);

    uint32_t    next_bits(uint8_t n);
    uint32_t    read_bits(uint8_t n);

    uint32_t    u (uint8_t n,   const char* name="");
    int32_t     i (uint8_t n,   const char* name="");
    uint32_t    f (uint8_t n,   const char* name="");
    uint32_t    b (uint8_t n=8, const char* name="");

    uint32_t    ue(const char* name="");
    int32_t     se(const char* name="");
    uint32_t    ae(const char* name="");
    uint32_t    ce(const char* name="");
    uint32_t    me(const char* name="");
    uint32_t    te(const char* name="");


    void        vui_parameters(vui_t& vui);
    void        hrd_parameters(hrd_t& hrd);

    void        seq_parameter_set_rbsp(sps_t& sps);
    void        seq_parameter_set_data(sps_t& sps);
    void        scaling_list(int* scalingList, int sizeOfScalingList, bool* useDefaultScalingMatrixFlag);
    void        seq_parameter_set_extension_rbsp(sps_ext_t& sps_ext);
    void        seq_parameter_set_svc_extension (sps_svc_t& sps_svc);
    void        svc_vui_parameters_extension    (svc_vui_t& svc_vui);
    void        seq_parameter_set_mvc_extension (sps_mvc_t& sps_mvc);
    void        mvc_vui_parameters_extension    (mvc_vui_t& mvc_vui);
    void        seq_parameter_set_mvcd_extension(sps_mvcd_t& sps_mvcd);
    void        subset_seq_parameter_set_rbsp   (sub_sps_t& sub_sps);

    void        rbsp_trailing_bits(void);
};


struct cabac_context_t;

struct cabac_engine_t {
    data_partition_t* dp;
    uint16_t    codIRange;
    uint16_t    codIOffset;

    void        init(data_partition_t* dp);

    bool        decode_decision (cabac_context_t* ctx);
    bool        decode_bypass   ();
    bool        decode_terminate();
    void        renormD();

    uint32_t    u  (cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx);
    uint32_t    tu (cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx, uint32_t cMax);
    int32_t     ueg(cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx, uint32_t cMax, uint8_t k);
    uint32_t    fl (cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx, uint32_t cMax);
};


}
}


#endif /* _DATA_PARTITION_H_ */
