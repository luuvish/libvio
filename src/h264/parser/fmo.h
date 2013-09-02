#ifndef _FMO_H_
#define _FMO_H_

struct slice_t;

int fmo_init(VideoParameters *p_Vid, slice_t *pSlice);
int FmoFinit(VideoParameters *p_Vid);

int FmoGetNextMBNr(VideoParameters *p_Vid, int CurrentMbNr);

#endif
