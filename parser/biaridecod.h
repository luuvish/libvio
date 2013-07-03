
/*!
 ***************************************************************************
 * \file
 *    biaridecod.h
 *
 * \brief
 *    Headerfile for binary arithmetic decoder routines
 *
 * \author
 *    Detlev Marpe,
 *    Gabi Blättermann
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000
 **************************************************************************
 */

#ifndef _BIARIDECOD_H_
#define _BIARIDECOD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bitstream.h"

/************************************************************************
 * D e f i n i t i o n s
 ***********************************************************************
 */

void arideco_start_decoding(DecodingEnvironment *eep, unsigned char *code_buffer, int firstbyte, int *code_len);
int  arideco_bits_read(DecodingEnvironment *dep);
void arideco_done_decoding(DecodingEnvironment *dep);
void biari_init_context (int qp, BiContextTypePtr ctx, const char* ini);
unsigned int biari_decode_symbol(DecodingEnvironment *dep, BiContextType *bi_ct );
unsigned int biari_decode_symbol_eq_prob(DecodingEnvironment *dep);
unsigned int biari_decode_final(DecodingEnvironment *dep);

#ifdef __cplusplus
}
#endif

#endif /* _BIARIDECOD_H_ */
