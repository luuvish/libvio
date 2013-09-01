#ifndef _NALU_H_
#define _NALU_H_


#include "memalloc.h"

//! Data Partitioning Modes
typedef enum {
    PAR_DP_1,   //!< no data partitioning is supported
    PAR_DP_3    //!< data partitioning with 3 partitions
} PAR_DP_TYPE;

//! values for nal_unit_type
typedef enum {
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
#if (MVC_EXTENSION_ENABLE)
    NALU_TYPE_PREFIX   = 14,
    NALU_TYPE_SUB_SPS  = 15,
    NALU_TYPE_SLC_EXT  = 20,
    NALU_TYPE_VDRD     = 24  // View and Dependency Representation Delimiter NAL Unit
#endif
} NaluType;

//! NAL unit structure
struct nalu_t {
    int32_t     startcodeprefix_len;   //!< 4 for parameter sets and first slice in picture, 3 for everything else (suggested)
    uint32_t    len;                   //!< Length of the NAL unit (Excluding the start code, which does not belong to the NALU)
    uint32_t    max_size;              //!< NAL Unit Buffer size
    int32_t     forbidden_bit;         //!< should be always FALSE
    NaluType    nal_unit_type;         //!< NALU_TYPE_xxxx
    uint8_t     nal_ref_idc;           //!< NALU_PRIORITY_xxxx  
    uint8_t*    buf;                   //!< contains the first byte followed by the EBSP
    uint16_t    lost_packets;          //!< true, if packet loss is detected
#if (MVC_EXTENSION_ENABLE)
    int32_t     svc_extension_flag;    //!< should be always 0, for MVC
    int32_t     non_idr_flag;          //!< 0 = current is IDR
    int32_t     priority_id;           //!< a lower value of priority_id specifies a higher priority
    int32_t     view_id;               //!< view identifier for the NAL unit
    int32_t     temporal_id;           //!< temporal identifier for the NAL unit
    int32_t     anchor_pic_flag;       //!< anchor access unit
    int32_t     inter_view_flag;       //!< inter-view prediction enable
    int32_t     reserved_one_bit;      //!< shall be equal to 1
#endif

    nalu_t(int buffersize) {
        this->max_size = buffersize;
        this->buf = new uint8_t[buffersize]; 
        if (!this->buf)
            no_mem_exit("AllocNALU: n->buf");
    }

    ~nalu_t() {
        if (this->buf)
            delete this->buf;
    }
};


struct cabac_context_t;
struct data_partition_t;

struct cabac_engine_t {
    data_partition_t* dp;
    uint16_t    codIRange;
    uint16_t    codIOffset;

    void        init(data_partition_t* dp);

    bool        decode_decision (cabac_context_t* ctx);
    bool        decode_bypass   ();
    bool        decode_terminate();
    void        renormD();

    uint32_t    u (cabac_context_t* ctx, int ctx_offset);
    uint32_t    tu(cabac_context_t* ctx, int ctx_offset, unsigned int max_symbol);
};

struct data_partition_t {
    int         frame_bitoffset;    //!< actual position in the codebuffer, bit-oriented, CAVLC only
    int         bitstream_length;   //!< over codebuffer lnegth, byte oriented, CAVLC only
    uint8_t*    streamBuffer;      //!< actual codebuffer for read bytes

    cabac_engine_t de_cabac;

    uint32_t    nal_unit_header_bytes;
    uint8_t     nal_ref_idc;
    uint8_t     nal_unit_type;
    bool        svc_extension_flag;
    uint32_t    num_bytes_in_rbsp;
    uint8_t*    rbsp_byte;

                data_partition_t();
                ~data_partition_t();

    void        init(nalu_t* nalu);

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
};


#endif /* _NALU_H_ */
