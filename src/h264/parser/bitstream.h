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


struct annex_b_t;
struct nalu_t;

struct bitstream_t {
    enum class type { ANNEX_B, RTP };

    type        FileFormat; //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
    int         BitStreamFile;
    annex_b_t*  annex_b;
    int         LastAccessUnitExists;
    int         NALUCount;

    void        open (const char* name, type format, unsigned max_size);
    void        close();
    void        reset();

    int         read_next_nalu     (nalu_t* nalu);
    void        CheckZeroByteNonVCL(nalu_t* nalu);
    void        CheckZeroByteVCL   (nalu_t* nalu);
};


void open_rtp         (const char* fn, int* p_BitStreamFile);
void close_rtp        (int* p_BitStreamFile);
int  get_nalu_from_rtp(nalu_t* nalu, int BitStreamFile);


#endif /* _BITSTREAM_H_ */
