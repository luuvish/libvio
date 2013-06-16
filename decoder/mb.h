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
#include "mbuffer.h"
#include "block.h"

int decode_one_component_i_slice (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int decode_one_component_p_slice (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int decode_one_component_b_slice (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int decode_one_component_sp_slice(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);

#ifdef __cplusplus
}
#endif

#endif
