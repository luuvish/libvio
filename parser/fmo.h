
/*!
 ***************************************************************************
 *
 * \file fmo.h
 *
 * \brief
 *    Support for Flexilble Macroblock Ordering (FMO)
 *
 * \date
 *    19 June, 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifndef _FMO_H_
#define _FMO_H_

#ifdef __cplusplus
extern "C" {
#endif

struct slice_t;

int fmo_init(VideoParameters *p_Vid, struct slice_t *pSlice);
int FmoFinit(VideoParameters *p_Vid);

int FmoGetNextMBNr(VideoParameters *p_Vid, int CurrentMbNr);

#ifdef __cplusplus
}
#endif

#endif
