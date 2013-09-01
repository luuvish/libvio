#ifndef _H264DECODER_H_
#define _H264DECODER_H_


#include "global.h"

typedef enum {
  	DEC_GEN_NOERR = 0,
  	DEC_OPEN_NOERR = 0,
  	DEC_CLOSE_NOERR = 0,  
  	DEC_SUCCEED = 0,
  	DEC_EOS =1,
  	DEC_NEED_DATA = 2,
  	DEC_INVALID_PARAM = 3,
  	DEC_ERRMASK = 0x8000
} DecErrCode;

int OpenDecoder(InputParameters *p_Inp);
int DecodeOneFrame(DecodedPicList **ppDecPic);
int FinitDecoder(DecodedPicList **ppDecPicList);
int CloseDecoder();


#endif
