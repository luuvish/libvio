#ifndef _IMAGE_DATA_H_
#define _IMAGE_DATA_H_

struct ImageData {
    int yuv_format;

    imgpel** frm_data[3];
    imgpel** top_data[3];
    imgpel** bot_data[3];

    int frm_stride[3];
    int top_stride[3];
    int bot_stride[3];
};

struct VideoParameters;
struct storable_picture;
struct slice_t;

extern void free_img_data(VideoParameters *p_Vid, ImageData *p_ImgData);
extern void picture_in_dpb(slice_t* currSlice, VideoParameters *p_Vid, storable_picture *p_pic);

#endif // _IMAGE_DATA_H_
