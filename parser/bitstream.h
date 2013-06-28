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
} Bitstream;

typedef struct bitstream_t {
    int                    FileFormat; //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
    int                    BitStreamFile;
    struct annex_b_struct *annex_b;
    int                    LastAccessUnitExists;
    int                    NALUCount;
} bitstream_t;

void open_bitstream(struct bitstream_t **bitstream,
                    char *name, int format, unsigned max_size);
void close_bitstream(struct bitstream_t *bitstream);
void reset_bitstream(struct bitstream_t *bitstream);

struct datapartition_dec;
struct syntaxelement_dec;
struct macroblock_dec;

int read_se_v (const char *tracestring, Bitstream *bitstream, int *used_bits);
int read_ue_v (const char *tracestring, Bitstream *bitstream, int *used_bits);
Boolean read_u_1 (const char *tracestring, Bitstream *bitstream, int *used_bits);
int read_u_v (int LenInBits, const char *tracestring, Bitstream *bitstream, int *used_bits);
int read_i_v (int LenInBits, const char *tracestring, Bitstream *bitstream, int *used_bits);

// CAVLC mapping
void linfo_ue(int len, int info, int *value1, int *dummy);
void linfo_se(int len, int info, int *value1, int *dummy);

int  readSyntaxElement_VLC (struct syntaxelement_dec *sym, Bitstream *currStream);
int  readSyntaxElement_UVLC(struct macroblock_dec *currMB, struct syntaxelement_dec *sym, struct datapartition_dec *dp);
int  readSyntaxElement_FLC(struct syntaxelement_dec *sym, Bitstream *currStream);
int  GetBits  (byte buffer[],int totbitoffset,int *info, int bitcount, int numbits);
int  ShowBits (byte buffer[],int totbitoffset,int bitcount, int numbits);

int  more_rbsp_data (byte buffer[],int totbitoffset,int bytecount);

#ifdef __cplusplus
}
#endif

#endif /* _BITSTREAM_H_ */
