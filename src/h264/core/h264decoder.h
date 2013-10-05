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

struct DecodedPicList {
    int         bValid;                 //0: invalid, 1: valid, 3: valid for 3D output;
    int         iViewId;                //-1: single view, >=0 multiview[VIEW1|VIEW0];
    int         iPOC;
    int         iYUVFormat;             //0: 4:0:0, 1: 4:2:0, 2: 4:2:2, 3: 4:4:4
    int         iYUVStorageFormat;      //0: YUV seperate; 1: YUV interleaved; 2: 3D output;
    int         iBitDepth;
    uint8_t*    pY;                     //if iPictureFormat is 1, [0]: top; [1] bottom;
    uint8_t*    pU;
    uint8_t*    pV;
    int         iWidth;                 //frame width;              
    int         iHeight;                //frame height;
    int         iYBufStride;            //stride of pY[0/1] buffer in bytes;
    int         iUVBufStride;           //stride of pU[0/1] and pV[0/1] buffer in bytes;
    int         iSkipPicNum;
    int         iBufSize;
    DecodedPicList* pNext;
};

struct DecoderParams {
    InputParameters *p_Inp;
    VideoParameters *p_Vid;

    void OpenDecoder(InputParameters* p_Inp);
    int  DecodeOneFrame(DecodedPicList** ppDecPic);
    void FinitDecoder(DecodedPicList** ppDecPicList);
    void CloseDecoder();
};


#endif
