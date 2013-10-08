#ifndef _MBUFFERDEC_H_
#define _MBUFFERDEC_H_

#include "global.h"
#include "bitstream_cabac.h"

#define MAX_NUM_SLICES 50

#define MAX_LIST_SIZE 33

//! Field Coding Types
enum {
    FRAME_CODING         = 0,
    FIELD_CODING         = 1,
    ADAPTIVE_CODING      = 2,
    FRAME_MB_PAIR_CODING = 3
};

enum {
    LIST_0 = 0,
    LIST_1 = 1
};


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

//! definition a picture (field or frame)
struct storable_picture {
    PictureStructure structure;

    int         poc;
    int         top_poc;
    int         bottom_poc;
    int         frame_poc;
    unsigned    frame_num;
    unsigned    recovery_frame;

    int         PicNum;
    int         LongTermPicNum;
    int         LongTermFrameIdx;

    byte        is_long_term;
    int         used_for_reference;
    int         is_output;
    int         non_existing;
    int         separate_colour_plane_flag;

    int         size_x, size_y, size_x_cr, size_y_cr;
    int         coded_frame;
    int         mb_aff_frame_flag;
    unsigned    PicWidthInMbs;
    unsigned    PicSizeInMbs;
    int         iChromaPadX;
    int         iChromaPadY;


    imgpel**    imgY;
    imgpel***   imgUV;

    pic_motion_params** mv_info;
    pic_motion_params** JVmv_info[3];

    pic_motion_params_old motion;
    pic_motion_params_old JVmotion[3];

    storable_picture* top_field;
    storable_picture* bottom_field;
    storable_picture* frame;

    int         slice_type;
    int         idr_flag;
    int         no_output_of_prior_pics_flag;
    int         long_term_reference_flag;
    int         adaptive_ref_pic_buffering_flag;

    drpm_t*     dec_ref_pic_marking_buffer;

    // picture error concealment
    int         concealed_pic; //indicates if this is a concealed picture

    // variables for tone mapping
    int         seiHasTone_mapping;
    int         tone_mapping_model_id;
    int         tonemapped_bit_depth;  
    imgpel*     tone_mapping_lut;                //!< tone mapping look up table

#if (MVC_EXTENSION_ENABLE)
    int         view_id;
    int         inter_view_flag;
    int         anchor_pic_flag;
#endif
    int         iLumaStride;
    int         iChromaStride;
    imgpel**    cur_imgY; // for more efficient get_block_luma
    int         no_ref;
    int         iCodingType;
    //
    char              listXsize[MAX_NUM_SLICES][2];
    storable_picture* listX[MAX_NUM_SLICES][2][33];
    int         layer_id;

    int         get_pic_num_x(int difference_of_pic_nums_minus1);

    bool        is_short_ref();
    bool        is_long_ref();
};


extern storable_picture* alloc_storable_picture(VideoParameters *p_Vid, PictureStructure type, int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output);
extern void              free_storable_picture (storable_picture* p);


struct frame_store {
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

    frame_store() =default;
    ~frame_store();

    void        unmark_for_reference();
    void        unmark_for_long_term_reference();

    bool        is_short_term_reference();
    bool        is_long_term_reference();
    bool        is_used_for_reference();
};

struct decoded_picture_buffer_t {
    VideoParameters* p_Vid;
    InputParameters* p_Inp;
    frame_store** fs;
    frame_store** fs_ref;
    frame_store** fs_ltref;
    frame_store** fs_ilref; // inter-layer reference (for multi-layered codecs)
    unsigned    size;
    unsigned    used_size;
    unsigned    ref_frames_in_buffer;
    unsigned    ltref_frames_in_buffer;
    int         last_output_poc;
#if (MVC_EXTENSION_ENABLE)
    int         last_output_view_id;
#endif
    int         max_long_term_pic_idx;  

    int         init_done;
    int         num_ref_frames;

    frame_store* last_picture;
    unsigned    used_size_il;
    int         layer_id;

    void        init(VideoParameters* p_Vid, int type);
    void        free();

    void        store_picture(storable_picture* p);

#if (MVC_EXTENSION_ENABLE)
    void        idr_memory_management(storable_picture* p);
#endif
    void        flush();

    void        update_ref_list();
    void        update_ltref_list();

protected:
    void        get_smallest_poc(int* poc, int* pos);

    bool        remove_unused_frame();
    void        remove_frame(int pos);
    bool        output_one_frame();

    void        check_num_ref();
    void        mm_unmark_short_term_for_reference(storable_picture* p, int difference_of_pic_nums_minus1);
    void        mm_unmark_long_term_for_reference(storable_picture* p, int long_term_pic_num);
    void        unmark_long_term_frame_for_reference_by_frame_idx(int long_term_frame_idx);
    void        unmark_long_term_field_for_reference_by_frame_idx(PictureStructure structure, int long_term_frame_idx, int mark_current, unsigned curr_frame_num, int curr_pic_num);
    void        mark_pic_long_term(storable_picture* p, int long_term_frame_idx, int picNumX);
    void        mm_assign_long_term_frame_idx(storable_picture* p, int difference_of_pic_nums_minus1, int long_term_frame_idx);
    void        mm_update_max_long_term_frame_idx(int max_long_term_frame_idx_plus1);
    void        mm_unmark_all_short_term_for_reference();
    void        mm_unmark_all_long_term_for_reference();
    void        mm_mark_current_picture_long_term(storable_picture* p, int long_term_frame_idx);
    void        adaptive_memory_management(storable_picture* p);
    void        sliding_window_memory_management(storable_picture* p);

    void        conceal_non_ref_pics(int diff);
    void        sliding_window_poc_management(storable_picture* p);
    void        write_lost_non_ref_pic(int poc, int p_out);
    void        write_lost_ref_after_idr(int pos);
};

using dpb_t = decoded_picture_buffer_t;


#endif
