/*!
 ************************************************************************
 * \file macroblock.h
 *
 * \brief
 *    Arrays for macroblock encoding
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    Copyright (C) 1999 Telenor Satellite Services, Norway
 ************************************************************************
 */

#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "mbuffer.h"
#include "block.h"

struct slice_t;

extern void setup_slice_methods(struct slice_t *currSlice);

extern void start_macroblock     (struct slice_t *currSlice, Macroblock **currMB);
extern int  decode_one_macroblock(Macroblock *currMB, StorablePicture *dec_picture);
extern Boolean  exit_macroblock  (struct slice_t *currSlice, int eos_bit);
extern void update_qp            (Macroblock *currMB, int qp);


#ifdef __cplusplus
}
#endif

#endif
