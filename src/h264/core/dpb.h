#ifndef _DPB_H_
#define _DPB_H_

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


struct frame_store;

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

public:
    void        init(VideoParameters* p_Vid, int type);
    void        free();

    void        store_picture(storable_picture* p);

#if (MVC_EXTENSION_ENABLE)
    void        idr_memory_management(storable_picture* p);
    void        store_proc_picture(storable_picture* p);
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


extern void fill_frame_num_gap(VideoParameters *p_Vid, slice_t *pSlice);
extern void pad_dec_picture(VideoParameters *p_Vid, storable_picture *dec_picture);
extern void pad_buf(imgpel *pImgBuf, int iWidth, int iHeight, int iStride, int iPadX, int iPadY);


storable_picture* get_ref_pic(mb_t& mb, storable_picture** RefPicListX, int ref_idx);


#endif // _DPB_H_
