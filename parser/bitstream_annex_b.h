
/*!
 *************************************************************************************
 * \file annexb.h
 *
 * \brief
 *    Annex B byte stream buffer handling.
 *
 *************************************************************************************
 */

#ifndef _ANNEXB_H_
#define _ANNEXB_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bitstream_nal.h"

typedef struct annex_b_struct {
  	int   BitStreamFile;                //!< the bit stream file
  	byte *iobuffer;
  	byte *iobufferread;
  	int   bytesinbuffer;
  	int   is_eof;
  	int   iIOBufferSize;

  	int   IsFirstByteStreamNALU;
  	int   nextstartcodebytes;
  	byte *Buf;  
} ANNEXB_t;

int  get_nalu_from_annex_b(NALU_t *nalu, ANNEXB_t *annex_b);

void open_annex_b  (char *fn, ANNEXB_t *annex_b);
void close_annex_b (ANNEXB_t *annex_b);
void malloc_annex_b(unsigned max_size, ANNEXB_t **p_annex_b);
void free_annex_b  (ANNEXB_t **p_annex_b);
void init_annex_b  (ANNEXB_t *annex_b);
void reset_annex_b (ANNEXB_t *annex_b);

#ifdef __cplusplus
}
#endif

#endif /* _ANNEXB_H_ */
