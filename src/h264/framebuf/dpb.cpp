#include <limits.h>

#include "global.h"
#include "input_parameters.h"

#include "slice.h"
#include "dpb.h"
#include "memalloc.h"
#include "output.h"

using vio::h264::mb_t;

#include "erc_api.h"

static inline int RSD(int x)
{
    return ((x&2)?(x|1):(x&(~1)));
}


static inline int is_FREXT_profile(unsigned int profile_idc) 
{
    return ( profile_idc >= FREXT_HP || profile_idc == FREXT_CAVLC444 );
}



void decoded_picture_buffer_t::get_smallest_poc(int* poc, int* pos)
{
    if (this->used_size < 1)
        error("Cannot determine smallest POC, DPB empty.", 150);

    *pos = -1;
    *poc = INT_MAX;
    for (int i = 0; i < this->used_size; i++) {
        if (*poc > this->fs[i]->poc && !this->fs[i]->is_output) {
            *poc = this->fs[i]->poc;
            *pos = i;
        }
    }
}

bool decoded_picture_buffer_t::remove_unused_frame()
{
    for (int i = 0; i < this->used_size; i++) {
        if (this->fs[i]->is_output && !this->fs[i]->is_used_for_reference()) {
            this->remove_frame(i);
            return true;
        }
    }
    return false;
}


int getDpbSize(VideoParameters* p_Vid, sps_t* active_sps)
{
    int pic_size = (active_sps->pic_width_in_mbs_minus1 + 1) * (active_sps->pic_height_in_map_units_minus1 + 1) * (active_sps->frame_mbs_only_flag ? 1 : 2) * 384;

    int size = 0;

    switch (active_sps->level_idc) {
    case 9:
    case 10:
        size = 152064;
        break;
    case 11:
        if (!is_FREXT_profile(active_sps->profile_idc) && (active_sps->constraint_set3_flag == 1))
            size = 152064;
        else
            size = 345600;
        break;
    case 12:
    case 13:
    case 20:
        size = 912384;
        break;
    case 21:
        size = 1824768;
        break;
    case 22:
    case 30:
        size = 3110400;
        break;
    case 31:
        size = 6912000;
        break;
    case 32:
        size = 7864320;
        break;
    case 40:
    case 41:
        size = 12582912;
        break;
    case 42:
        size = 13369344;
        break;
    case 50:
        size = 42393600;
        break;
    case 51:
    case 52:
        size = 70778880;
        break;
    default:
        error("undefined level", 500);
        break;
    }

    size /= pic_size;
#if MVC_EXTENSION_ENABLE
    if (p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH) {
        int num_views = p_Vid->active_subset_sps->num_views_minus1+1;
        size = min(2 * size, max<int>(1, round(log2(num_views))) * 16) / num_views;
    } else
#endif
        size = min(size, 16);

    if (active_sps->vui_parameters_present_flag && active_sps->vui_parameters.bitstream_restriction_flag) {
        int size_vui;
        if ((int)active_sps->vui_parameters.max_dec_frame_buffering > size)
            error("max_dec_frame_buffering larger than MaxDpbSize", 500);
        size_vui = max<int>(1, active_sps->vui_parameters.max_dec_frame_buffering);
#ifdef _DEBUG
        if (size_vui < size)
            printf("Warning: max_dec_frame_buffering(%d) is less than DPB size(%d) calculated from Profile/Level.\n", size_vui, size);
#endif
        size = size_vui;    
    }

    return size;
}


void decoded_picture_buffer_t::init(VideoParameters* p_Vid, int type)
{
    sps_t& sps = *p_Vid->active_sps;

    this->p_Vid = p_Vid;
    if (this->init_done)
        this->free();

    this->size = getDpbSize(p_Vid, &sps) + p_Vid->p_Inp->dpb_plus[type == 2 ? 1 : 0];
    this->num_ref_frames = sps.max_num_ref_frames; 

#if (MVC_EXTENSION_ENABLE)
    if (sps.max_dec_frame_buffering < sps.max_num_ref_frames)
#else
    if (this->size < sps.max_num_ref_frames)
#endif
        error ("DPB size at specified level is smaller than the specified number of reference frames. This is not allowed.\n", 1000);

    this->used_size    = 0;
    this->last_picture = nullptr;

    this->ref_frames_in_buffer   = 0;
    this->ltref_frames_in_buffer = 0;

    this->fs       = new frame_store*[this->size];
    this->fs_ref   = new frame_store*[this->size];
    this->fs_ltref = new frame_store*[this->size];
#if (MVC_EXTENSION_ENABLE)
    this->fs_ilref = new frame_store*[1];
#endif

    for (int i = 0; i < this->size; i++) {
        this->fs[i]       = new frame_store {};
        this->fs_ref[i]   = nullptr;
        this->fs_ltref[i] = nullptr;
        this->fs[i]->layer_id = -1;
#if (MVC_EXTENSION_ENABLE)
        this->fs[i]->view_id = -1;
        this->fs[i]->inter_view_flag[0] = this->fs[i]->inter_view_flag[1] = 0;
        this->fs[i]->anchor_pic_flag[0] = this->fs[i]->anchor_pic_flag[1] = 0;
#endif
    }

#if (MVC_EXTENSION_ENABLE)
    if (type == 2) {
        this->fs_ilref[0] = new frame_store {};
        // These may need some cleanups
        this->fs_ilref[0]->view_id = -1;
        this->fs_ilref[0]->inter_view_flag[0] = this->fs_ilref[0]->inter_view_flag[1] = 0;
        this->fs_ilref[0]->anchor_pic_flag[0] = this->fs_ilref[0]->anchor_pic_flag[1] = 0;
    } else
        this->fs_ilref[0] = nullptr;
#endif

    this->init_done = 1;
    this->last_output_poc = INT_MIN;
#if (MVC_EXTENSION_ENABLE)
    this->last_output_view_id = -1;
#endif

    /* allocate a dummy storable picture */
    if (!p_Vid->no_reference_picture) {
        p_Vid->no_reference_picture = new storable_picture(p_Vid, FRAME,
            sps.PicWidthInMbs * 16, sps.FrameHeightInMbs * 16,
            sps.PicWidthInMbs * sps.MbWidthC, sps.FrameHeightInMbs * sps.MbHeightC, 1);
        p_Vid->no_reference_picture->top_field    = p_Vid->no_reference_picture;
        p_Vid->no_reference_picture->bottom_field = p_Vid->no_reference_picture;
        p_Vid->no_reference_picture->frame        = p_Vid->no_reference_picture;
    }

    p_Vid->last_has_mmco_5 = 0;
    // picture error concealment
    if (p_Vid->conceal_mode != 0 && !p_Vid->last_out_fs)
        p_Vid->last_out_fs = new frame_store {};
}

void decoded_picture_buffer_t::free()
{
    VideoParameters* p_Vid = this->p_Vid;

    if (this->fs) {
        for (int i = 0; i < this->size; i++)
            delete this->fs[i];
        delete []this->fs;
        this->fs = nullptr;
    }

    if (this->fs_ref)
        delete []this->fs_ref;
    if (this->fs_ltref)
        delete []this->fs_ltref;

#if (MVC_EXTENSION_ENABLE)
    if (this->fs_ilref) {
        for (int i = 0; i < 1; i++)
            delete this->fs_ilref[i];
        delete []this->fs_ilref;
        this->fs_ilref = nullptr;
    }

    this->last_output_view_id = -1;
#endif
    this->last_output_poc = INT_MIN;
    this->init_done = 0;

    // picture error concealment
    if (p_Vid->conceal_mode != 0 || p_Vid->last_out_fs)
        delete p_Vid->last_out_fs;

    if (p_Vid->no_reference_picture) {
        delete p_Vid->no_reference_picture;
        p_Vid->no_reference_picture = nullptr;
    }
}

void decoded_picture_buffer_t::idr_memory_management(storable_picture* p)
{
    if (p->slice.no_output_of_prior_pics_flag) {
        // free all stored pictures
        for (int i = 0; i < this->used_size; i++) {
            // reset all reference settings
            delete this->fs[i];
            this->fs[i] = new frame_store {};
        }
        for (int i = 0; i < this->ref_frames_in_buffer; i++)
            this->fs_ref[i] = nullptr;
        for (int i = 0; i < this->ltref_frames_in_buffer; i++)
            this->fs_ltref[i] = nullptr;
        this->used_size = 0;
    } else
        this->flush();
    this->last_picture = nullptr;

    this->update_ref_list();
    this->update_ltref_list();
    this->last_output_poc = INT_MIN;

    if (p->slice.long_term_reference_flag) {
        this->max_long_term_pic_idx = 0;
        p->is_long_term             = 1;
        p->LongTermFrameIdx         = 0;
    } else {
        this->max_long_term_pic_idx = -1;
        p->is_long_term             = 0;
    }

#if (MVC_EXTENSION_ENABLE)
    this->last_output_view_id = -1;
#endif
}

void decoded_picture_buffer_t::store_proc_picture(storable_picture* p)
{
    VideoParameters* p_Vid = this->p_Vid;
    frame_store* fs = this->fs_ilref[0];
    if (this->used_size_il > 0 && fs->is_used == 3) {
        if (fs->frame) {
            delete fs->frame;
            fs->frame = nullptr;
        }
        if (fs->top_field) {
            delete fs->top_field;
            fs->top_field = nullptr;
        }
        if (fs->bottom_field) {
            delete fs->bottom_field;
            fs->bottom_field = nullptr;
        }
        fs->is_used = 0;
        fs->is_reference = 0;
        this->used_size_il--;   
    }

    fs->insert_picture(p_Vid, p);
    if ((p->slice.structure == FRAME && fs->is_used == 3) ||
        (p->slice.structure != FRAME && fs->is_used && fs->is_used < 3))
        this->used_size_il++;  
}

void decoded_picture_buffer_t::remove_frame(int pos)
{
    frame_store* fs = this->fs[pos];

    switch (fs->is_used) {
    case 3:
        delete fs->frame;
        delete fs->top_field;
        delete fs->bottom_field;
        fs->frame        = nullptr;
        fs->top_field    = nullptr;
        fs->bottom_field = nullptr;
        break;
    case 2:
        delete fs->bottom_field;
        fs->bottom_field = nullptr;
        break;
    case 1:
        delete fs->top_field;
        fs->top_field = nullptr;
        break;
    case 0:
        break;
    default:
        error("invalid frame store type", 500);
    }
    fs->is_used           = 0;
    fs->is_long_term      = 0;
    fs->is_reference      = 0;
    fs->is_orig_reference = 0;

    // move empty framestore to end of buffer
    frame_store* tmp = this->fs[pos];

    for (int i = pos; i < this->used_size - 1; i++)
        this->fs[i] = this->fs[i + 1];
    this->fs[this->used_size - 1] = tmp;
    this->used_size--;
}


bool decoded_picture_buffer_t::output_one_frame()
{
    VideoParameters* p_Vid = this->p_Vid;
    //diagnostics
    if (this->used_size < 1)
        error("Cannot output frame, DPB empty.",150);

    // find smallest POC
    int poc, pos;
    this->get_smallest_poc(&poc, &pos);
    if (pos == -1)
        return false;

#if (DISABLE_ERC == 0)
    // picture error concealment
    if (p_Vid->conceal_mode != 0) {
        if (this->last_output_poc == 0)
            this->write_lost_ref_after_idr(pos);
#if (MVC_EXTENSION_ENABLE)
        this->write_lost_non_ref_pic(poc, p_Vid->p_out_mvc[this->layer_id]);
#else
        this->write_lost_non_ref_pic(poc, p_Vid->p_out);
#endif
    }
#endif
// JVT-P072 ends

#if (MVC_EXTENSION_ENABLE)
    write_stored_frame(p_Vid, this->fs[pos], p_Vid->p_out_mvc[this->layer_id]);
#else
    write_stored_frame(p_Vid, this->fs[pos], p_Vid->p_out);
#endif

    // picture error concealment
    if (p_Vid->conceal_mode == 0) {
        if (this->last_output_poc >= poc)
            error("output POC must be in ascending order", 150);
    }

    this->last_output_poc = poc;

    // free frame store and move empty store to end of buffer
    if (!this->fs[pos]->is_used_for_reference())
        this->remove_frame(pos);
    return true;
}


void decoded_picture_buffer_t::flush()
{
    VideoParameters *p_Vid = this->p_Vid;

    if (!this->init_done)
        return;
#if (DISABLE_ERC == 0)
    if (p_Vid->conceal_mode != 0)
        this->conceal_non_ref_pics(0);
#endif
    // mark all frames unused
    for (int i = 0; i < this->used_size; i++) {
#if MVC_EXTENSION_ENABLE
        assert(this->fs[i]->view_id == this->layer_id);
#endif
        this->fs[i]->unmark_for_reference();
    }

    while (this->remove_unused_frame());

    // output frames in POC order
    while (this->used_size && this->output_one_frame());

    this->last_output_poc = INT_MIN;
}

void decoded_picture_buffer_t::update_ref_list()
{
    int i, j;
    for (i = 0, j = 0; i < this->used_size; i++) {
        if (this->fs[i]->is_short_term_reference())
            this->fs_ref[j++] = this->fs[i];
    }

    this->ref_frames_in_buffer = j;

    while (j < this->size)
        this->fs_ref[j++] = nullptr;
}

void decoded_picture_buffer_t::update_ltref_list()
{
    int i, j;
    for (i = 0, j = 0; i < this->used_size; i++) {
        if (this->fs[i]->is_long_term_reference())
            this->fs_ltref[j++] = this->fs[i];
    }

    this->ltref_frames_in_buffer = j;

    while (j < this->size)
        this->fs_ltref[j++] = nullptr;
}

void decoded_picture_buffer_t::check_num_ref()
{
    if (this->ltref_frames_in_buffer + this->ref_frames_in_buffer > max(1, this->num_ref_frames))
        error("Max. number of reference frames exceeded. Invalid stream.", 500);
}

void decoded_picture_buffer_t::mm_unmark_short_term_for_reference(storable_picture* p, int difference_of_pic_nums_minus1)
{
    int currPicNum = (p->slice.structure == FRAME) ? p->frame_num : 2 * p->frame_num + 1;
    int picNumX = currPicNum - (difference_of_pic_nums_minus1 + 1);

    for (int i = 0; i < this->ref_frames_in_buffer; i++) {
        if (p->slice.structure == FRAME) {
            if (this->fs_ref[i]->is_reference == 3 && this->fs_ref[i]->is_long_term == 0) {
                if (this->fs_ref[i]->frame->PicNum == picNumX)
                    return this->fs_ref[i]->unmark_for_reference();
            }
        } else {
            if ((this->fs_ref[i]->is_reference & 1) && !(this->fs_ref[i]->is_long_term & 1)) {
                if (this->fs_ref[i]->top_field->PicNum == picNumX) {
                    this->fs_ref[i]->top_field->used_for_reference = 0;
                    this->fs_ref[i]->is_reference &= 2;
                    if (this->fs_ref[i]->is_used == 3)
                        this->fs_ref[i]->frame->used_for_reference = 0;
                    return;
                }
            }
            if ((this->fs_ref[i]->is_reference & 2) && !(this->fs_ref[i]->is_long_term & 2)) {
                if (this->fs_ref[i]->bottom_field->PicNum == picNumX) {
                    this->fs_ref[i]->bottom_field->used_for_reference = 0;
                    this->fs_ref[i]->is_reference &= 1;
                    if (this->fs_ref[i]->is_used == 3)
                        this->fs_ref[i]->frame->used_for_reference = 0;
                    return;
                }
            }
        }
    }
}

void decoded_picture_buffer_t::mm_unmark_long_term_for_reference(storable_picture* p, int long_term_pic_num)
{
    for (int i = 0; i < this->ltref_frames_in_buffer; i++) {
        if (p->slice.structure == FRAME) {
            if (this->fs_ltref[i]->is_reference == 3 && this->fs_ltref[i]->is_long_term == 3) {
                if (this->fs_ltref[i]->frame->LongTermPicNum == long_term_pic_num)
                    this->fs_ltref[i]->unmark_for_long_term_reference();
            }
        } else {
            if ((this->fs_ltref[i]->is_reference & 1) && (this->fs_ltref[i]->is_long_term & 1)) {
                if (this->fs_ltref[i]->top_field->LongTermPicNum == long_term_pic_num) {
                    this->fs_ltref[i]->top_field->used_for_reference = 0;
                    this->fs_ltref[i]->top_field->is_long_term = 0;
                    this->fs_ltref[i]->is_reference &= 2;
                    this->fs_ltref[i]->is_long_term &= 2;
                    if (this->fs_ltref[i]->is_used == 3) {
                        this->fs_ltref[i]->frame->used_for_reference = 0;
                        this->fs_ltref[i]->frame->is_long_term = 0;
                    }
                    return;
                }
            }
            if ((this->fs_ltref[i]->is_reference & 2) && (this->fs_ltref[i]->is_long_term & 2)) {
                if (this->fs_ltref[i]->bottom_field->LongTermPicNum == long_term_pic_num) {
                    this->fs_ltref[i]->bottom_field->used_for_reference = 0;
                    this->fs_ltref[i]->bottom_field->is_long_term = 0;
                    this->fs_ltref[i]->is_reference &= 1;
                    this->fs_ltref[i]->is_long_term &= 1;
                    if (this->fs_ltref[i]->is_used == 3) {
                        this->fs_ltref[i]->frame->used_for_reference = 0;
                        this->fs_ltref[i]->frame->is_long_term = 0;
                    }
                    return;
                }
            }
        }
    }
}

void decoded_picture_buffer_t::unmark_long_term_frame_for_reference_by_frame_idx(int long_term_frame_idx)
{
    for (int i = 0; i < this->ltref_frames_in_buffer; i++) {
        if (this->fs_ltref[i]->LongTermFrameIdx == long_term_frame_idx)
            this->fs_ltref[i]->unmark_for_long_term_reference();
    }
}

void decoded_picture_buffer_t::unmark_long_term_field_for_reference_by_frame_idx(PictureStructure structure, int long_term_frame_idx, int mark_current, unsigned curr_frame_num, int curr_pic_num)
{
    VideoParameters *p_Vid = this->p_Vid;

    assert(structure != FRAME);
    if (curr_pic_num < 0)
        curr_pic_num += (2 * p_Vid->active_sps->MaxFrameNum);

    for (int i = 0; i < this->ltref_frames_in_buffer; i++) {
        if (this->fs_ltref[i]->LongTermFrameIdx == long_term_frame_idx) {
            if (structure == TOP_FIELD) {
                if (this->fs_ltref[i]->is_long_term == 3)
                    this->fs_ltref[i]->unmark_for_long_term_reference();
                else if (this->fs_ltref[i]->is_long_term == 1)
                    this->fs_ltref[i]->unmark_for_long_term_reference();
                else if (mark_current) {
                    if (this->last_picture) {
                        if (this->last_picture != this->fs_ltref[i] || this->last_picture->FrameNum != curr_frame_num)
                            this->fs_ltref[i]->unmark_for_long_term_reference();
                    } else
                        this->fs_ltref[i]->unmark_for_long_term_reference();
                } else if (this->fs_ltref[i]->FrameNum != (unsigned)(curr_pic_num >> 1))
                    this->fs_ltref[i]->unmark_for_long_term_reference();
            }
            if (structure == BOTTOM_FIELD) {
                if (this->fs_ltref[i]->is_long_term == 3)
                    this->fs_ltref[i]->unmark_for_long_term_reference();
                else if (this->fs_ltref[i]->is_long_term == 2)
                    this->fs_ltref[i]->unmark_for_long_term_reference();
                else if (mark_current) {
                    if (this->last_picture) {
                        if (this->last_picture != this->fs_ltref[i] || this->last_picture->FrameNum != curr_frame_num)
                            this->fs_ltref[i]->unmark_for_long_term_reference();
                    } else
                        this->fs_ltref[i]->unmark_for_long_term_reference();
                } else if (this->fs_ltref[i]->FrameNum != (unsigned)(curr_pic_num >> 1))
                    this->fs_ltref[i]->unmark_for_long_term_reference();
            }
        }
    }
}

void decoded_picture_buffer_t::mark_pic_long_term(storable_picture* p, int long_term_frame_idx, int picNumX)
{
    int add_top, add_bottom;

    if (p->slice.structure == FRAME) {
        for (int i = 0; i < this->ref_frames_in_buffer; i++) {
            if (this->fs_ref[i]->is_reference == 3) {
                if (!this->fs_ref[i]->frame->is_long_term &&
                     this->fs_ref[i]->frame->PicNum == picNumX) {
                    this->fs_ref[i]->LongTermFrameIdx =
                    this->fs_ref[i]->frame->LongTermFrameIdx = long_term_frame_idx;
                    this->fs_ref[i]->frame->LongTermPicNum = long_term_frame_idx;
                    this->fs_ref[i]->frame->is_long_term = 1;

                    if (this->fs_ref[i]->top_field && this->fs_ref[i]->bottom_field) {
                        this->fs_ref[i]->top_field->LongTermFrameIdx =
                        this->fs_ref[i]->bottom_field->LongTermFrameIdx = long_term_frame_idx;
                        this->fs_ref[i]->top_field->LongTermPicNum = long_term_frame_idx;
                        this->fs_ref[i]->bottom_field->LongTermPicNum = long_term_frame_idx;
                        this->fs_ref[i]->top_field->is_long_term =
                        this->fs_ref[i]->bottom_field->is_long_term = 1;
                    }
                    this->fs_ref[i]->is_long_term = 3;
                    return;
                }
            }
        }
        printf("Warning: reference frame for long term marking not found\n");
    } else {
        if (p->slice.structure == TOP_FIELD) {
            add_top    = 1;
            add_bottom = 0;
        } else {
            add_top    = 0;
            add_bottom = 1;
        }
        for (int i = 0; i < this->ref_frames_in_buffer; i++) {
            if (this->fs_ref[i]->is_reference & 1) {
                if (!this->fs_ref[i]->top_field->is_long_term &&
                     this->fs_ref[i]->top_field->PicNum == picNumX) {
                    if (this->fs_ref[i]->is_long_term &&
                        this->fs_ref[i]->LongTermFrameIdx != long_term_frame_idx)
                        printf("Warning: assigning long_term_frame_idx different from other field\n");

                    this->fs_ref[i]->LongTermFrameIdx =
                    this->fs_ref[i]->top_field->LongTermFrameIdx = long_term_frame_idx;
                    this->fs_ref[i]->top_field->LongTermPicNum = 2 * long_term_frame_idx + add_top;
                    this->fs_ref[i]->top_field->is_long_term = 1;
                    this->fs_ref[i]->is_long_term |= 1;
                    if (this->fs_ref[i]->is_long_term == 3) {
                        this->fs_ref[i]->frame->is_long_term = 1;
                        this->fs_ref[i]->frame->LongTermFrameIdx =
                        this->fs_ref[i]->frame->LongTermPicNum = long_term_frame_idx;
                    }
                    return;
                }
            }
            if (this->fs_ref[i]->is_reference & 2) {
                if (!this->fs_ref[i]->bottom_field->is_long_term &&
                     this->fs_ref[i]->bottom_field->PicNum == picNumX) {
                    if (this->fs_ref[i]->is_long_term &&
                        this->fs_ref[i]->LongTermFrameIdx != long_term_frame_idx)
                        printf("Warning: assigning long_term_frame_idx different from other field\n");

                    this->fs_ref[i]->LongTermFrameIdx =
                    this->fs_ref[i]->bottom_field->LongTermFrameIdx = long_term_frame_idx;
                    this->fs_ref[i]->bottom_field->LongTermPicNum = 2 * long_term_frame_idx + add_bottom;
                    this->fs_ref[i]->bottom_field->is_long_term = 1;
                    this->fs_ref[i]->is_long_term |= 2;
                    if (this->fs_ref[i]->is_long_term == 3) {
                        this->fs_ref[i]->frame->is_long_term = 1;
                        this->fs_ref[i]->frame->LongTermFrameIdx =
                        this->fs_ref[i]->frame->LongTermPicNum = long_term_frame_idx;
                    }
                    return;
                }
            }
        }
        printf("Warning: reference field for long term marking not found\n");
    }
}

void decoded_picture_buffer_t::mm_assign_long_term_frame_idx(storable_picture* p, int difference_of_pic_nums_minus1, int long_term_frame_idx)
{
    int currPicNum = (p->slice.structure == FRAME) ? p->frame_num : 2 * p->frame_num + 1;
    int picNumX = currPicNum - (difference_of_pic_nums_minus1 + 1);

    // remove frames/fields with same long_term_frame_idx
    if (p->slice.structure == FRAME)
        this->unmark_long_term_frame_for_reference_by_frame_idx(long_term_frame_idx);
    else {
        PictureStructure structure = FRAME;

        for (int i = 0; i < this->ref_frames_in_buffer; i++) {
            if (this->fs_ref[i]->is_reference & 1) {
                if (this->fs_ref[i]->top_field->PicNum == picNumX) {
                    structure = TOP_FIELD;
                    break;
                }
            }
            if (this->fs_ref[i]->is_reference & 2) {
                if (this->fs_ref[i]->bottom_field->PicNum == picNumX) {
                    structure = BOTTOM_FIELD;
                    break;
                }
            }
        }
        if (structure == FRAME)
            error("field for long term marking not found", 200);

        this->unmark_long_term_field_for_reference_by_frame_idx(structure, long_term_frame_idx, 0, 0, picNumX);
    }

    this->mark_pic_long_term(p, long_term_frame_idx, picNumX);
}

void decoded_picture_buffer_t::mm_update_max_long_term_frame_idx(int max_long_term_frame_idx_plus1)
{
    this->max_long_term_pic_idx = max_long_term_frame_idx_plus1 - 1;

    // check for invalid frames
    for (int i = 0; i < this->ltref_frames_in_buffer; i++) {
        if (this->fs_ltref[i]->LongTermFrameIdx > this->max_long_term_pic_idx)
            this->fs_ltref[i]->unmark_for_long_term_reference();
    }
}

void decoded_picture_buffer_t::mm_unmark_all_short_term_for_reference()
{
    for (int i = 0; i < this->ref_frames_in_buffer; i++)
        this->fs_ref[i]->unmark_for_reference();
    this->update_ref_list();
}

void decoded_picture_buffer_t::mm_unmark_all_long_term_for_reference()
{
    this->mm_update_max_long_term_frame_idx(0);
}

void decoded_picture_buffer_t::mm_mark_current_picture_long_term(storable_picture* p, int long_term_frame_idx)
{
    // remove long term pictures with same long_term_frame_idx
    if (p->slice.structure == FRAME)
        this->unmark_long_term_frame_for_reference_by_frame_idx(long_term_frame_idx);
    else
        this->unmark_long_term_field_for_reference_by_frame_idx(p->slice.structure, long_term_frame_idx, 1, p->PicNum, 0);

    p->is_long_term = 1;
    p->LongTermFrameIdx = long_term_frame_idx;
}

void decoded_picture_buffer_t::adaptive_memory_management(storable_picture* p)
{
    drpm_t* tmp_drpm;
    VideoParameters *p_Vid = this->p_Vid;

    p_Vid->last_has_mmco_5 = 0;

    assert(!p->slice.idr_flag);
    assert(p->slice.adaptive_ref_pic_buffering_flag);

    while (p->slice.dec_ref_pic_marking_buffer) {
        tmp_drpm = p->slice.dec_ref_pic_marking_buffer;
        switch (tmp_drpm->memory_management_control_operation) {
        case 0:
            if (tmp_drpm->Next != NULL)
                error("memory_management_control_operation = 0 not last operation in buffer", 500);
            break;
        case 1:
            this->mm_unmark_short_term_for_reference(p, tmp_drpm->difference_of_pic_nums_minus1);
            this->update_ref_list();
            break;
        case 2:
            this->mm_unmark_long_term_for_reference(p, tmp_drpm->long_term_pic_num);
            this->update_ltref_list();
            break;
        case 3:
            this->mm_assign_long_term_frame_idx(p, tmp_drpm->difference_of_pic_nums_minus1, tmp_drpm->long_term_frame_idx);
            this->update_ref_list();
            this->update_ltref_list();
            break;
        case 4:
            this->mm_update_max_long_term_frame_idx(tmp_drpm->max_long_term_frame_idx_plus1);
            this->update_ltref_list();
            break;
        case 5:
            this->mm_unmark_all_short_term_for_reference();
            this->mm_unmark_all_long_term_for_reference();
            p_Vid->last_has_mmco_5 = 1;
            break;
        case 6:
            this->mm_mark_current_picture_long_term(p, tmp_drpm->long_term_frame_idx);
            this->check_num_ref();
            break;
        default:
            error("invalid memory_management_control_operation in buffer", 500);
        }
        p->slice.dec_ref_pic_marking_buffer = tmp_drpm->Next;
        delete tmp_drpm;
    }
    if (p_Vid->last_has_mmco_5) {
        p->PicNum = p->frame_num = 0;

        switch (p->slice.structure) {
        case TOP_FIELD:
            p->poc = p->top_poc = 0;
            break;
        case BOTTOM_FIELD:
            p->poc = p->bottom_poc = 0;
            break;
        case FRAME:
            p->top_poc    -= p->poc;
            p->bottom_poc -= p->poc;
            p->poc = min(p->top_poc, p->bottom_poc);
            break;
        }
#if (MVC_EXTENSION_ENABLE)
        if (p->slice.view_id == 0) {
            p_Vid->p_Dpb_layer[0]->flush();
            p_Vid->p_Dpb_layer[1]->flush();
        } else
            this->flush();
#endif
    }
}

void decoded_picture_buffer_t::sliding_window_memory_management(storable_picture* p)
{
    assert(!p->slice.idr_flag);

    // if this is a reference pic with sliding window, unmark first ref frame
    if (this->ref_frames_in_buffer == max(1, this->num_ref_frames) - this->ltref_frames_in_buffer) {
        for (int i = 0; i < this->used_size; i++) {
            if (this->fs[i]->is_reference && !this->fs[i]->is_long_term) {
                this->fs[i]->unmark_for_reference();
                this->update_ref_list();
                break;
            }
        }
    }

    p->is_long_term = 0;
}

void decoded_picture_buffer_t::store_picture(storable_picture* p)
{
    VideoParameters *p_Vid = this->p_Vid;
    unsigned i;
    int poc, pos;
    // picture error concealment
  
    // diagnostics
    // if frame, check for new store,
    assert(p != NULL);

    p_Vid->last_has_mmco_5       = 0;
    p_Vid->last_pic_bottom_field = p->slice.structure == BOTTOM_FIELD;

    if (p->slice.idr_flag) {
        this->idr_memory_management(p);
        // picture error concealment
        memset(p_Vid->pocs_in_dpb, 0, sizeof(int)*100);
    } else {
        // adaptive memory management
        if (p->used_for_reference && p->slice.adaptive_ref_pic_buffering_flag)
            this->adaptive_memory_management(p);
    }

    if (p->slice.structure == TOP_FIELD || p->slice.structure == BOTTOM_FIELD) {
        // check for frame store with same pic_number
        if (this->last_picture) {
            if ((int)this->last_picture->FrameNum == p->PicNum) {
                if ((p->slice.structure == TOP_FIELD && this->last_picture->is_used == 2) ||
                    (p->slice.structure == BOTTOM_FIELD && this->last_picture->is_used == 1)) {
                    if (( p->used_for_reference && this->last_picture->is_orig_reference != 0) ||
                        (!p->used_for_reference && this->last_picture->is_orig_reference == 0)) {
                        this->last_picture->insert_picture(p_Vid, p);
                        this->update_ref_list();
                        this->update_ltref_list();
                        this->last_picture = NULL;
                        return;
                    }
                }
            }
        }
    }

    // this is a frame or a field which has no stored complementary field

    // sliding window, if necessary
    if (!p->slice.idr_flag && p->used_for_reference && !p->slice.adaptive_ref_pic_buffering_flag)
        this->sliding_window_memory_management(p);

    // picture error concealment
    if (p_Vid->conceal_mode != 0) {
        for (i = 0; i < this->size; i++)
            if (this->fs[i]->is_reference)
                this->fs[i]->concealment_reference = 1;
    }

    // first try to remove unused frames
    if (this->used_size == this->size) {
#if (DISABLE_ERC == 0)
        // picture error concealment
        if (p_Vid->conceal_mode != 0)
            this->conceal_non_ref_pics(2);
#endif
        this->remove_unused_frame();

#if (DISABLE_ERC == 0)
        if (p_Vid->conceal_mode != 0)
            this->sliding_window_poc_management(p);
#endif
    }

    // then output frames until one can be removed
    while (this->used_size == this->size) {
        // non-reference frames may be output directly
        if (!p->used_for_reference) {
            this->get_smallest_poc(&poc, &pos);
            if (-1 == pos || p->poc < poc) {
#if (MVC_EXTENSION_ENABLE)
                direct_output(p_Vid, p, p_Vid->p_out_mvc[this->layer_id]);
#endif
                return;
            }
        }
        // flush a frame
        this->output_one_frame();
    }

    // check for duplicate frame number in short term reference buffer
    if (p->used_for_reference && !p->is_long_term) {
        for (i = 0; i < this->ref_frames_in_buffer; i++) {
            if (this->fs_ref[i]->FrameNum == p->frame_num)
                error("duplicate frame_num in short-term reference picture buffer", 500);
        }
    }
    // store at end of buffer
    this->fs[this->used_size]->insert_picture(p_Vid, p);

    // picture error concealment
    if (p->slice.idr_flag)
        p_Vid->earlier_missing_poc = 0;

    if (p->slice.structure != FRAME)
        this->last_picture = this->fs[this->used_size];
    else
        this->last_picture = NULL;

    this->used_size++;

    if (p_Vid->conceal_mode != 0)
        p_Vid->pocs_in_dpb[this->used_size-1] = p->poc;

    this->update_ref_list();
    this->update_ltref_list();

    this->check_num_ref();
}


void fill_frame_num_gap(VideoParameters *p_Vid, slice_t *currSlice)
{
    sps_t& sps = *p_Vid->active_sps;
    slice_t& slice = *currSlice;
    shr_t& shr = slice.header;
  
    int tmp1 = shr.delta_pic_order_cnt[0];
    int tmp2 = shr.delta_pic_order_cnt[1];
    shr.delta_pic_order_cnt[0] =
    shr.delta_pic_order_cnt[1] = 0;

    printf("A gap in frame number is found, try to fill it.\n");

    int CurrFrameNum            = shr.frame_num;
    int UnusedShortTermFrameNum = (p_Vid->pre_frame_num + 1) % sps.MaxFrameNum;

    while (CurrFrameNum != UnusedShortTermFrameNum) {
        storable_picture* picture = new storable_picture(p_Vid, FRAME,
            sps.PicWidthInMbs * 16, sps.FrameHeightInMbs * 16,
            sps.PicWidthInMbs * sps.MbWidthC, sps.FrameHeightInMbs * sps.MbHeightC, 1);
        picture->slice.coded_frame                     = 1;
        picture->slice.adaptive_ref_pic_buffering_flag = 0;
        picture->PicNum                                = UnusedShortTermFrameNum;
        picture->frame_num                             = UnusedShortTermFrameNum;
        picture->non_existing                          = 1;
        picture->is_output                             = 1;
        picture->used_for_reference                    = 1;
#if (MVC_EXTENSION_ENABLE)
        picture->slice.view_id                         = slice.view_id;
#endif

        shr.frame_num = UnusedShortTermFrameNum;
        if (sps.pic_order_cnt_type != 0)
            decode_poc(p_Vid, p_Vid->ppSliceList[0]);
        picture->top_poc    = shr.TopFieldOrderCnt;
        picture->bottom_poc = shr.BottomFieldOrderCnt;
        picture->frame_poc  = shr.PicOrderCnt;
        picture->poc        = shr.PicOrderCnt;

        slice.p_Dpb->store_picture(picture);

        p_Vid->pre_frame_num = UnusedShortTermFrameNum;
        UnusedShortTermFrameNum = (UnusedShortTermFrameNum + 1) % sps.MaxFrameNum;
    }

    shr.delta_pic_order_cnt[0] = tmp1;
    shr.delta_pic_order_cnt[1] = tmp2;
    shr.frame_num = CurrFrameNum;
}


storable_picture* get_ref_pic(mb_t& mb, storable_picture** RefPicListX, int ref_idx)
{
    slice_t& slice = *mb.p_Slice;
    shr_t& shr = slice.header;

    return (ref_idx < 0) ? nullptr :
            !shr.MbaffFrameFlag || !mb.mb_field_decoding_flag ? RefPicListX[ref_idx] :
            mb.mbAddrX % 2 == ref_idx % 2 ? 
            RefPicListX[ref_idx / 2]->top_field : RefPicListX[ref_idx / 2]->bottom_field;
}
