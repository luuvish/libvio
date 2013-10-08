#ifndef _H264DECODER_H_
#define _H264DECODER_H_


#include <cstdint>


enum {
    DEC_SUCCEED = 0,
    DEC_EOS     = 1,
    DEC_ERRMASK = 0x8000
};

struct InputParameters;
struct VideoParameters;

struct DecoderParams {
    InputParameters* p_Inp;
    VideoParameters* p_Vid;

    void OpenDecoder(InputParameters* p_Inp);
    int  DecodeOneFrame();
    void FinitDecoder();
    void CloseDecoder();

    int  decode_one_frame();
};


#endif
