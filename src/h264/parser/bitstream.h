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
 *  File      : bitstream.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _BITSTREAM_H_
#define _BITSTREAM_H_

#ifdef __cplusplus
extern "C" {
#endif

struct annex_b_t;
struct nalu_t;
struct macroblock_t;
struct datapartition_dec;
struct cabac_context_t;

//! Data Partitioning Modes
typedef enum {
    PAR_DP_1,   //!< no data partitioning is supported
    PAR_DP_3    //!< data partitioning with 3 partitions
} PAR_DP_TYPE;

//! Output File Types
typedef enum {
    PAR_OF_ANNEXB,    //!< Annex B byte stream format
    PAR_OF_RTP       //!< RTP packets in outfile
} PAR_OF_TYPE;

struct bitstream_t {
    int        FileFormat; //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
    int        BitStreamFile;
    annex_b_t* annex_b;
    int        LastAccessUnitExists;
    int        NALUCount;
};


//! struct to characterize the state of the arithmetic coding engine
struct cabac_engine_t {
    uint16_t codIRange;
    uint32_t Dvalue;
    int      DbitsLeft;
    byte    *Dcodestrm;
    int     *Dcodestrm_len;

    void     init(unsigned char *code_buffer, int firstbyte, int *code_len);

    bool     decode_decision(cabac_context_t *ctx);
    bool     decode_bypass();
    bool     decode_terminate();
    void     renormD();

    uint32_t u (cabac_context_t* ctx, int ctx_offset);
    uint32_t tu(cabac_context_t* ctx, int ctx_offset, unsigned int max_symbol);
};

//! Bitstream
typedef struct bit_stream_dec {
    // CABAC Decoding
    int           read_len;           //!< actual position in the codebuffer, CABAC only
    int           code_len;           //!< overall codebuffer length, CABAC only
    // CAVLC Decoding
    int           frame_bitoffset;    //!< actual position in the codebuffer, bit-oriented, CAVLC only
    int           bitstream_length;   //!< over codebuffer lnegth, byte oriented, CAVLC only
    // ErrorConcealment
    byte          *streamBuffer;      //!< actual codebuffer for read bytes

    cabac_engine_t de_cabac;

    uint32_t nal_unit_header_bytes;
    uint8_t  nal_ref_idc;
    uint8_t  nal_unit_type;
    bool     svc_extension_flag;
    uint32_t num_bytes_in_rbsp;
    uint8_t *rbsp_byte;

    bool     byte_aligned            (void);
    bool     more_data_in_byte_stream(void);
    bool     more_rbsp_data          (void);
    bool     more_rbsp_trailing_data (void);

    uint32_t next_bits(uint8_t n);
    uint32_t read_bits(uint8_t n);

    uint32_t u (uint8_t n,   const char *name="");
    int32_t  i (uint8_t n,   const char *name="");
    uint32_t f (uint8_t n,   const char *name="");
    uint32_t b (uint8_t n=8, const char *name="");

    uint32_t ue(const char *name="");
    int32_t  se(const char *name="");
    uint32_t ae(const char *name="");
    uint32_t ce(const char *name="");
    uint32_t me(const char *name="");
    uint32_t te(const char *name="");
} Bitstream;

typedef struct datapartition_dec {
    Bitstream *bitstream;
} DataPartition;




void open_bitstream(struct bitstream_t **bitstream,
                    char *name, int format, unsigned max_size);
void close_bitstream(struct bitstream_t *bitstream);
void reset_bitstream(struct bitstream_t *bitstream);

DataPartition *AllocPartition(int n);
void FreePartition(DataPartition *dp, int n);
Bitstream *InitPartition(DataPartition *dp, struct nalu_t *nalu);


struct syntax_element_t {
    int type;    //!< type of syntax element for data part.
    int value1;  //!< numerical value of syntax element
    int value2;  //!< for blocked symbols, e.g. run/level
    int len;     //!< length of code
    int inf;     //!< info part of CAVLC code
    int context; //!< CABAC context
};


// CAVLC mapping
int GetBits(byte buffer[],int totbitoffset,int *info, int bitcount, int numbits);



#ifdef __cplusplus
}
#endif

#endif /* _BITSTREAM_H_ */
