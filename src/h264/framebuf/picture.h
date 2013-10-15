#ifndef _FRAME_BUFFER_H_
#define _FRAME_BUFFER_H_


#include <cstdint>
#include <vector>
#include "global.h"


struct motion_vector_t {
    int16_t mv_x;
    int16_t mv_y;
};

using mv_t = motion_vector_t;

inline bool operator == (const mv_t& l, const mv_t& r)
{
    return (l.mv_x == r.mv_x) && (l.mv_y == r.mv_y);
}

inline mv_t operator + (const mv_t& l, const mv_t& r)
{
    return {(int16_t)(l.mv_x + r.mv_x), (int16_t)(l.mv_y + r.mv_y)};
}

inline mv_t operator + (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x + r), (int16_t)(l.mv_y + r)};
}

inline mv_t operator - (const mv_t& l, const mv_t& r)
{
    return {(int16_t)(l.mv_x - r.mv_x), (int16_t)(l.mv_y - r.mv_y)};
}

inline mv_t operator - (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x - r), (int16_t)(l.mv_y - r)};
}

inline mv_t operator * (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x * r), (int16_t)(l.mv_y * r)};
}

inline mv_t operator * (int l, const mv_t& r)
{
    return {(int16_t)(l * r.mv_x), (int16_t)(l * r.mv_y)};
}

inline mv_t operator >> (const mv_t& l, int r)
{
    return {(int16_t)(l.mv_x >> r), (int16_t)(l.mv_y >> r)};
}

#include "slice.h"


struct storable_picture;

struct pic_motion_params_old {
    bool*       mb_field_decoding_flag;
};

struct pic_motion_params {
    storable_picture* ref_pic[2];
    mv_t        mv[2];
    char        ref_idx[2];
    uint8_t     slice_no;
};

struct decoded_reference_picture_marking_t;
using drpm_t = decoded_reference_picture_marking_t;

struct slice_header_t;
struct macroblock_t;


struct storable_picture {
    sps_t*                       sps;
    pps_t*                       pps;
    std::vector<slice_t*>        slice_headers;
    std::vector<macroblock_t*>   mbs;
    pic_motion_params**          mv_infos;
    imgpel**                     pixels[3];

    int         poc;
    int         top_poc;
    int         bottom_poc;
    int         frame_poc;
    unsigned    frame_num;
    unsigned    recovery_frame;

    int         PicNum;
    int         LongTermPicNum;
    int         LongTermFrameIdx;

    int         is_long_term;
    int         used_for_reference;
    int         is_output;
    int         no_ref;
    int         non_existing;

    struct {
        int         layer_id;
#if (MVC_EXTENSION_ENABLE)
        int         view_id;
        int         inter_view_flag;
        int         anchor_pic_flag;
#endif

        PictureStructure structure;
        int         iCodingType;
        int         idr_flag;
        int         slice_type;
    } slice;

    int         size_x, size_y, size_x_cr, size_y_cr;

    int         iLumaStride;
    int         iChromaStride;
    int         iChromaPadX;
    int         iChromaPadY;
    imgpel**    imgY;
    imgpel**    imgUV[2];

    pic_motion_params** mv_info;
    pic_motion_params** JVmv_info[3];

    pic_motion_params_old motion;
    pic_motion_params_old JVmotion[3];

    storable_picture* top_field;
    storable_picture* bottom_field;
    storable_picture* frame;

    // picture error concealment
    int         concealed_pic; //indicates if this is a concealed picture

    // variables for tone mapping
    int         seiHasTone_mapping;
    int         tone_mapping_model_id;
    int         tonemapped_bit_depth;  
    imgpel*     tone_mapping_lut;                //!< tone mapping look up table

    bool        is_short_ref();
    bool        is_long_ref();

    void        clear();

    storable_picture(VideoParameters *p_Vid, PictureStructure type, int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output);
    ~storable_picture();
};

struct picture_t {
    //sps_t*                       sps;
    //pps_t*                       pps;
    //std::vector<slice_header_t*> slice_headers;
    //std::vector<macroblock_t*>   mbs;
    //pic_motion_params**          mv_infos;
    //imgpel**                     pixels[3];

    //picture_t*  pic[2][3] = {{nullptr}};

    int         is_used;                //!< 0=empty; 1=top; 2=bottom; 3=both fields (or frame)
    int         is_reference;           //!< 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
    int         is_long_term;           //!< 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
    int         is_orig_reference;      //!< original marking by nal_ref_idc: 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used

    int         is_non_existent;

    uint32_t    FrameNum;
    int32_t     FrameNumWrap;
    unsigned    recovery_frame;

    int         LongTermFrameIdx;
    int         is_output;
    int         poc;

    // picture error concealment
    int         concealment_reference;

    storable_picture* frame;
    storable_picture* top_field;
    storable_picture* bottom_field;

#if (MVC_EXTENSION_ENABLE)
    int         view_id;
    int         inter_view_flag[2];
    int         anchor_pic_flag[2];
#endif
    int         layer_id;

    picture_t() =default;
    ~picture_t();

    void        unmark_for_reference();
    void        unmark_for_long_term_reference();

    bool        is_short_term_reference();
    bool        is_long_term_reference();
    bool        is_used_for_reference();

    void        insert_picture(VideoParameters* p_Vid, storable_picture* p);

    void        dpb_combine_field_yuv(VideoParameters* p_Vid);

protected:
    void        dpb_split_field  (VideoParameters* p_Vid);
    void        dpb_combine_field(VideoParameters* p_Vid);
};

using pic_t = picture_t;


#endif // _FRAME_BUFFER_H_
