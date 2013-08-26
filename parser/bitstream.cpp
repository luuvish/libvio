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
 *  File      : bitstream.cc
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "bitstream.h"
#include "bitstream_annex_b.h"
#include "bitstream_rtp.h"
#include "bitstream_nal.h"


void open_bitstream(struct bitstream_t **bitstream, char *name, int format, unsigned max_size)
{
    *bitstream = (bitstream_t *)calloc(1, sizeof(bitstream_t));

    (*bitstream)->FileFormat = format;

    (*bitstream)->LastAccessUnitExists  = 0;
    (*bitstream)->NALUCount = 0;

    switch (format) {
    case PAR_OF_RTP:
        open_rtp(name, &(*bitstream)->BitStreamFile);
        break;
    case PAR_OF_ANNEXB:
    default:
        malloc_annex_b(max_size, &(*bitstream)->annex_b);
        open_annex_b(name, (*bitstream)->annex_b);
        break;
    }
}

void close_bitstream(struct bitstream_t *bitstream)
{
    switch (bitstream->FileFormat) {
    case PAR_OF_RTP:
        close_rtp(&bitstream->BitStreamFile);
        break;   
    case PAR_OF_ANNEXB:
    default:
        close_annex_b(bitstream->annex_b);
        free_annex_b(&bitstream->annex_b);
        break;
    }

    free(bitstream);
}

void reset_bitstream(struct bitstream_t *bitstream)
{
    if (bitstream->FileFormat == PAR_OF_ANNEXB)
        reset_annex_b(bitstream->annex_b); 
}



DataPartition *AllocPartition(int n)
{
    DataPartition *partArr, *dataPart;
    int i;

    partArr = (DataPartition *)calloc(n, sizeof(DataPartition));
    if (partArr == NULL) {
        snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Data Partition failed");
        error(errortext, 100);
    }

    for (i = 0; i < n; ++i) { // loop over all data partitions
        dataPart = &(partArr[i]);
        dataPart->bitstream = (Bitstream *)calloc(1, sizeof(Bitstream));
        if (dataPart->bitstream == NULL) {
            snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Bitstream failed");
            error(errortext, 100);
        }
        dataPart->bitstream->streamBuffer = (byte *)calloc(MAX_CODED_FRAME_SIZE, sizeof(byte));
        if (dataPart->bitstream->streamBuffer == NULL) {
            snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for streamBuffer failed");
            error(errortext, 100);
        }
    }
    return partArr;
}

void FreePartition(DataPartition *dp, int n)
{
    int i;

    assert(dp != NULL);
    assert(dp->bitstream != NULL);
    assert(dp->bitstream->streamBuffer != NULL);
    for (i = 0; i < n; ++i) {
        free(dp[i].bitstream->streamBuffer);
        free(dp[i].bitstream);
    }
    free(dp);
}

Bitstream *InitPartition(DataPartition *dp, struct nalu_t *nalu)
{
    Bitstream *currStream = dp->bitstream;
    currStream->ei_flag    = 0;
    currStream->frame_bitoffset = currStream->read_len = 0;

    memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
    currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

    return currStream;
}
