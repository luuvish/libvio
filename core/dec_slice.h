

#ifndef _DEC_SLICE_H_
#define _DEC_SLICE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "slice.h"

bool init_slice(Slice *currSlice);
void decode_one_slice(Slice *currSlice);

#ifdef __cplusplus
}
#endif

#endif /* _DEC_SLICE_H_ */
