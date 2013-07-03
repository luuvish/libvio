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

// CAVLC block types
typedef enum {
    LUMA              =  0,
    LUMA_INTRA16x16DC =  1,
    LUMA_INTRA16x16AC =  2,
    CB                =  3,
    CB_INTRA16x16DC   =  4,
    CB_INTRA16x16AC   =  5,
    CR                =  8,
    CR_INTRA16x16DC   =  9,
    CR_INTRA16x16AC   = 10
} CAVLCBlockTypes;

// CABAC block types
typedef enum {
    LUMA_16DC     =   0,
    LUMA_16AC     =   1,
    LUMA_8x8      =   2,
    LUMA_8x4      =   3,
    LUMA_4x8      =   4,
    LUMA_4x4      =   5,
    CHROMA_DC     =   6,
    CHROMA_AC     =   7,
    CHROMA_DC_2x4 =   8,
    CHROMA_DC_4x4 =   9,
    CB_16DC       =  10,
    CB_16AC       =  11,
    CB_8x8        =  12,
    CB_8x4        =  13,
    CB_4x8        =  14,
    CB_4x4        =  15,
    CR_16DC       =  16,
    CR_16AC       =  17,
    CR_8x8        =  18,
    CR_8x4        =  19,
    CR_4x8        =  20,
    CR_4x4        =  21
} CABACBlockTypes;

struct annex_b_struct;
struct nalu_t;
struct macroblock_dec;
struct syntaxelement_dec;
struct datapartition_dec;
struct syntaxelement_dec;

//! Data Partitioning Modes
typedef enum
{
  PAR_DP_1,   //!< no data partitioning is supported
  PAR_DP_3    //!< data partitioning with 3 partitions
} PAR_DP_TYPE;

//! Output File Types
typedef enum
{
  PAR_OF_ANNEXB,    //!< Annex B byte stream format
  PAR_OF_RTP       //!< RTP packets in outfile
} PAR_OF_TYPE;

typedef struct bitstream_t {
    int                    FileFormat; //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
    int                    BitStreamFile;
    struct annex_b_struct *annex_b;
    int                    LastAccessUnitExists;
    int                    NALUCount;
} bitstream_t;


//! struct to characterize the state of the arithmetic coding engine
typedef struct {
    unsigned int  Drange;
    unsigned int  Dvalue;
    int           DbitsLeft;
    byte         *Dcodestrm;
    int          *Dcodestrm_len;
} DecodingEnvironment;

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
    int           ei_flag;            //!< error indication, 0: no error, else unspecified error

    DecodingEnvironment de_cabac;

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

//! Syntaxelement
typedef struct syntaxelement_dec {
    int           type;                  //!< type of syntax element for data part.
    int           value1;                //!< numerical value of syntax element
    int           value2;                //!< for blocked symbols, e.g. run/level
    int           len;                   //!< length of code
    int           inf;                   //!< info part of CAVLC code
    unsigned int  bitpattern;            //!< CAVLC bitpattern
    int           context;               //!< CABAC context
    int           k;                     //!< CABAC context for coeff_count,uv

    //! for mapping of CAVLC to syntaxElement
    void (*mapping)(int len, int info, int *value1, int *value2);
    //! used for CABAC: refers to actual coding method of each individual syntax element type
    void (*reading)(struct macroblock_dec *currMB, struct syntaxelement_dec *, DecodingEnvironment *);
} SyntaxElement;

//! DataPartition
typedef struct datapartition_dec {
    Bitstream *bitstream;

    int (*readSyntaxElement)(struct macroblock_dec *currMB, struct syntaxelement_dec *, struct datapartition_dec *);
          /*!< virtual function;
               actual method depends on chosen data partition and
               entropy coding method  */
} DataPartition;




void open_bitstream(struct bitstream_t **bitstream,
                    char *name, int format, unsigned max_size);
void close_bitstream(struct bitstream_t *bitstream);
void reset_bitstream(struct bitstream_t *bitstream);

DataPartition *AllocPartition(int n);
void FreePartition(DataPartition *dp, int n);
Bitstream *InitPartition(DataPartition *dp, struct nalu_t *nalu);



// CAVLC mapping
void linfo_ue(int len, int info, int *value1, int *dummy);
void linfo_se(int len, int info, int *value1, int *dummy);

int  readSyntaxElement_VLC (struct syntaxelement_dec *sym, Bitstream *currStream);
int  readSyntaxElement_UVLC(struct macroblock_dec *currMB, struct syntaxelement_dec *sym, struct datapartition_dec *dp);
int  GetBits  (byte buffer[],int totbitoffset,int *info, int bitcount, int numbits);
int  ShowBits (byte buffer[],int totbitoffset,int bitcount, int numbits);



#ifdef __cplusplus
}
#endif

#endif /* _BITSTREAM_H_ */
