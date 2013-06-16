
/*!
 *************************************************************************************
 * \file rtp.h
 *
 * \brief
 *    Prototypes for rtp.c
 *************************************************************************************
 */

#ifndef _RTP_H_
#define _RTP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bitstream_nal.h"

void open_rtp         (char *fn, int *p_BitStreamFile);
void close_rtp        (int *p_BitStreamFile);
int  get_nalu_from_rtp(NALU_t *nalu, int BitStreamFile);

#ifdef __cplusplus
}
#endif

#endif /* _RTP_H_ */
