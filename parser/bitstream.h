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
