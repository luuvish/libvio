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

#ifndef _MACROB_H_
#define _MACROB_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "dpb.h"

void decode_one_macroblock(Macroblock *currMB, StorablePicture *dec_picture);

#ifdef __cplusplus
}
#endif

#endif
