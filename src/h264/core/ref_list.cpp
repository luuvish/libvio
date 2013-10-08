#include <functional>

#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "dpb.h"
#include "ref_list.h"
#include "memalloc.h"
#include "output.h"

using vio::h264::mb_t;

#include "erc_api.h"



void update_pic_num(slice_t *currSlice)
{
    dpb_t *dpb = currSlice->p_Dpb;
    sps_t *sps = currSlice->active_sps;
    int i;

    if (!currSlice->field_pic_flag) {
        for (i = 0; i < dpb->ref_frames_in_buffer; i++) {
            frame_store *fs = dpb->fs_ref[i];
            if (fs->is_used == 3) {
                if (fs->frame->used_for_reference && !fs->frame->is_long_term) {
                    if (fs->FrameNum > currSlice->frame_num)
                        fs->FrameNumWrap = fs->FrameNum - sps->MaxFrameNum;
                    else
                        fs->FrameNumWrap = fs->FrameNum;
                    fs->frame->PicNum = fs->FrameNumWrap;
                }
            }
        }
        for (i = 0; i < dpb->ltref_frames_in_buffer; i++) {
            frame_store *fs = dpb->fs_ltref[i];
            if (fs->is_used == 3) {
                if (fs->frame->is_long_term)
                    fs->frame->LongTermPicNum = fs->frame->LongTermFrameIdx;
            }
        }
    } else {
        int add_top    = !currSlice->bottom_field_flag ? 1 : 0;
        int add_bottom = !currSlice->bottom_field_flag ? 0 : 1;

        for (i = 0; i < dpb->ref_frames_in_buffer; i++) {
            frame_store *fs = dpb->fs_ref[i];
            if (fs->is_reference) {
                if (fs->FrameNum > currSlice->frame_num)
                    fs->FrameNumWrap = fs->FrameNum - sps->MaxFrameNum;
                else
                    fs->FrameNumWrap = fs->FrameNum;
                if (fs->is_reference & 1)
                    fs->top_field->PicNum = 2 * fs->FrameNumWrap + add_top;
                if (fs->is_reference & 2)
                    fs->bottom_field->PicNum = 2 * fs->FrameNumWrap + add_bottom;
            }
        }
        for (i = 0; i < dpb->ltref_frames_in_buffer; i++) {
            frame_store *fs = dpb->fs_ltref[i];
            if (fs->is_long_term & 1)
                fs->top_field->LongTermPicNum = 2 * fs->top_field->LongTermFrameIdx + add_top;
            if (fs->is_long_term & 2)
                fs->bottom_field->LongTermPicNum = 2 * fs->bottom_field->LongTermFrameIdx + add_bottom;
        }
    }
}



void gen_pic_list_from_frame_list(bool bottom_field_flag, frame_store **fs_list, int list_idx, storable_picture **list, char *list_size, int long_term)
{
    int top_idx = 0;
    int bot_idx = 0;

    auto is_ref = long_term ? std::mem_fn(&storable_picture::is_long_ref) : std::mem_fn(&storable_picture::is_short_ref);

    if (!bottom_field_flag) {
        while (top_idx < list_idx || bot_idx<list_idx) {
            for (; top_idx < list_idx; top_idx++) {
                if (fs_list[top_idx]->is_used & 1) {
                    if (is_ref(fs_list[top_idx]->top_field)) {
                        // short term ref pic
                        list[(short) *list_size] = fs_list[top_idx]->top_field;
                        (*list_size)++;
                        top_idx++;
                        break;
                    }
                }
            }
            for (; bot_idx < list_idx; bot_idx++) {
                if (fs_list[bot_idx]->is_used & 2) {
                    if (is_ref(fs_list[bot_idx]->bottom_field)) {
                        // short term ref pic
                        list[(short) *list_size] = fs_list[bot_idx]->bottom_field;
                        (*list_size)++;
                        bot_idx++;
                        break;
                    }
                }
            }
        }
    } else {
        while (top_idx < list_idx || bot_idx<list_idx) {
            for (; bot_idx < list_idx; bot_idx++) {
                if (fs_list[bot_idx]->is_used & 2) {
                    if (is_ref(fs_list[bot_idx]->bottom_field)) {
                        // short term ref pic
                        list[(short) *list_size] = fs_list[bot_idx]->bottom_field;
                        (*list_size)++;
                        bot_idx++;
                        break;
                    }
                }
            }
            for (; top_idx < list_idx; top_idx++) {
                if (fs_list[top_idx]->is_used & 1) {
                    if (is_ref(fs_list[top_idx]->top_field)) {
                        // short term ref pic
                        list[(short) *list_size] = fs_list[top_idx]->top_field;
                        (*list_size)++;
                        top_idx++;
                        break;
                    }
                }
            }
        }
    }
}

static void init_lists_i_slice(slice_t *currSlice)
{
#if (MVC_EXTENSION_ENABLE)
    currSlice->listinterviewidx0 = 0;
    currSlice->listinterviewidx1 = 0;
#endif
    currSlice->listXsize[0] = 0;
    currSlice->listXsize[1] = 0;
}

static void init_lists_p_slice(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    dpb_t *p_Dpb = currSlice->p_Dpb;

#if (MVC_EXTENSION_ENABLE)
    currSlice->listinterviewidx0 = 0;
    currSlice->listinterviewidx1 = 0;
#endif

    if (!currSlice->field_pic_flag) {
        int list0idx = 0;
        for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used == 3) {
                if (p_Dpb->fs_ref[i]->frame->used_for_reference && !p_Dpb->fs_ref[i]->frame->is_long_term)
                    currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
        }

        // order list 0 by PicNum
        std::qsort(currSlice->listX[0], list0idx, sizeof(storable_picture*), [](const void* a, const void* b) {
            int pic_num1 = (*(storable_picture**)a)->PicNum;
            int pic_num2 = (*(storable_picture**)b)->PicNum;
            //int pic_num1 = (*(reinterpret_cast<const storable_picture**>(a)))->PicNum;
            //int pic_num2 = (*(reinterpret_cast<const storable_picture**>(b)))->PicNum;
            return (pic_num1 < pic_num2) ? 1 : (pic_num1 > pic_num2) ? -1 : 0;
        });
        currSlice->listXsize[0] = (char) list0idx;

        // long term handling
        for (int i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ltref[i]->is_used == 3) {
                if (p_Dpb->fs_ltref[i]->frame->is_long_term)
                    currSlice->listX[0][list0idx++] = p_Dpb->fs_ltref[i]->frame;
            }
        }

        std::qsort(&currSlice->listX[0][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0],
                   sizeof(storable_picture*), [](const void* a, const void* b) {
            int pic_num1 = (*(storable_picture**)a)->LongTermPicNum;
            int pic_num2 = (*(storable_picture**)b)->LongTermPicNum;
            return (pic_num1 < pic_num2) ? -1 : (pic_num1 > pic_num2) ? 1 : 0;
        });
        currSlice->listXsize[0] = (char) list0idx;
    } else {
        frame_store** fs_list0  = new frame_store*[p_Dpb->size];
        frame_store** fs_listlt = new frame_store*[p_Dpb->size];

        int list0idx = 0;
        for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_reference)
                fs_list0[list0idx++] = p_Dpb->fs_ref[i];
        }

        std::qsort(fs_list0, list0idx, sizeof(frame_store*), [](const void* a, const void* b) {
            int pic_num1 = (*(frame_store**)a)->FrameNumWrap;
            int pic_num2 = (*(frame_store**)b)->FrameNumWrap;
            return (pic_num1 < pic_num2) ? 1 : (pic_num1 > pic_num2) ? -1 : 0;
        });

        currSlice->listXsize[0] = 0;
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);

        // long term handling
        int listltidx = 0;
        for (int i = 0; i < p_Dpb->ltref_frames_in_buffer; i++)
            fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];

        std::qsort(fs_listlt, listltidx, sizeof(frame_store*), [](const void* a, const void* b) {
            int idx1 = (*(frame_store**)a)->LongTermFrameIdx;
            int idx2 = (*(frame_store**)b)->LongTermFrameIdx;
            return (idx1 < idx2) ? -1 : (idx1 > idx2) ? 1 : 0;
        });

        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);

        delete []fs_list0;
        delete []fs_listlt;
    }

    // set max size
    currSlice->listXsize[0] = (char) min<int>(currSlice->listXsize[0], currSlice->num_ref_idx_l0_active_minus1 + 1);
    currSlice->listXsize[1] = 0;

    // set the unused list entries to NULL
    for (int i = currSlice->listXsize[0]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[0][i] = p_Vid->no_reference_picture;
    for (int i = currSlice->listXsize[1]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[1][i] = p_Vid->no_reference_picture;
}

static void init_lists_b_slice(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    dpb_t *p_Dpb = currSlice->p_Dpb;

#if (MVC_EXTENSION_ENABLE)
    currSlice->listinterviewidx0 = 0;
    currSlice->listinterviewidx1 = 0;
#endif

    // B-slice_t
    if (!currSlice->field_pic_flag) {
        int list0idx = 0;
        for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used == 3) {
                if (p_Dpb->fs_ref[i]->frame->used_for_reference && !p_Dpb->fs_ref[i]->frame->is_long_term) {
                    if (currSlice->PicOrderCnt >= p_Dpb->fs_ref[i]->frame->poc) //!KS use >= for error concealment
                        currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
                }
            }
        }
        std::qsort(currSlice->listX[0], list0idx, sizeof(storable_picture*), [](const void* a, const void* b) {
            int poc1 = (*(storable_picture**)a)->poc;
            int poc2 = (*(storable_picture**)b)->poc;
            return (poc1 < poc2) ? 1 : (poc1 > poc2) ? -1 : 0;
        });

        //get the backward reference picture (POC>current POC) in list0;
        int list0idx_1 = list0idx;
        for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used == 3) {
                if (p_Dpb->fs_ref[i]->frame->used_for_reference && !p_Dpb->fs_ref[i]->frame->is_long_term) {
                    if (currSlice->PicOrderCnt < p_Dpb->fs_ref[i]->frame->poc)
                        currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
                }
            }
        }
        std::qsort(&currSlice->listX[0][list0idx_1], list0idx-list0idx_1, sizeof(storable_picture*), [](const void* a, const void* b) {
            int poc1 = (*(storable_picture**)a)->poc;
            int poc2 = (*(storable_picture**)b)->poc;
            return (poc1 < poc2) ? -1 : (poc1 > poc2) ? 1 : 0;
        });

        for (int j = 0; j < list0idx_1; j++)
            currSlice->listX[1][list0idx - list0idx_1 + j] = currSlice->listX[0][j];
        for (int j = list0idx_1; j < list0idx; j++)
            currSlice->listX[1][j - list0idx_1] = currSlice->listX[0][j];

        currSlice->listXsize[0] = currSlice->listXsize[1] = list0idx;

        // long term handling
        for (int i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ltref[i]->is_used == 3) {
                if (p_Dpb->fs_ltref[i]->frame->is_long_term) {
                    currSlice->listX[0][list0idx]   = p_Dpb->fs_ltref[i]->frame;
                    currSlice->listX[1][list0idx++] = p_Dpb->fs_ltref[i]->frame;
                }
            }
        }
        std::qsort(&currSlice->listX[0][(int)currSlice->listXsize[0]], list0idx - currSlice->listXsize[0],
                   sizeof(storable_picture*), [](const void* a, const void* b) {
            int pic_num1 = (*(storable_picture**)a)->LongTermPicNum;
            int pic_num2 = (*(storable_picture**)b)->LongTermPicNum;
            return (pic_num1 < pic_num2) ? -1 : (pic_num1 > pic_num2) ? 1 : 0;
        });
        std::qsort(&currSlice->listX[1][(int)currSlice->listXsize[0]], list0idx - currSlice->listXsize[0],
                   sizeof(storable_picture*), [](const void* a, const void* b) {
            int pic_num1 = (*(storable_picture**)a)->LongTermPicNum;
            int pic_num2 = (*(storable_picture**)b)->LongTermPicNum;
            return (pic_num1 < pic_num2) ? -1 : (pic_num1 > pic_num2) ? 1 : 0;
        });
        currSlice->listXsize[0] = currSlice->listXsize[1] = (char) list0idx;
    } else {
        frame_store** fs_list0  = new frame_store*[p_Dpb->size];
        frame_store** fs_list1  = new frame_store*[p_Dpb->size];
        frame_store** fs_listlt = new frame_store*[p_Dpb->size];

        currSlice->listXsize[0] = 0;
        currSlice->listXsize[1] = 1;

        int list0idx = 0;
        for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used) {
                if (currSlice->PicOrderCnt >= p_Dpb->fs_ref[i]->poc)
                    fs_list0[list0idx++] = p_Dpb->fs_ref[i];
            }
        }
        std::qsort(fs_list0, list0idx, sizeof(frame_store*), [](const void* a, const void* b) {
            int poc1 = (*(frame_store**)a)->poc;
            int poc2 = (*(frame_store**)b)->poc;
            return (poc1 < poc2) ? 1 : (poc1 > poc2) ? -1 : 0;
        });

        int list0idx_1 = list0idx;
        for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used) {
                if (currSlice->PicOrderCnt < p_Dpb->fs_ref[i]->poc)
                    fs_list0[list0idx++] = p_Dpb->fs_ref[i];
            }
        }
        std::qsort(&fs_list0[list0idx_1], list0idx-list0idx_1, sizeof(frame_store*), [](const void* a, const void* b) {
            int poc1 = (*(frame_store**)a)->poc;
            int poc2 = (*(frame_store**)b)->poc;
            return (poc1 < poc2) ? -1 : (poc1 > poc2) ? 1 : 0;
        });

        for (int j = 0; j < list0idx_1; j++)
            fs_list1[list0idx - list0idx_1 + j] = fs_list0[j];
        for (int j = list0idx_1; j < list0idx; j++)
            fs_list1[j - list0idx_1] = fs_list0[j];

        currSlice->listXsize[0] = 0;
        currSlice->listXsize[1] = 0;
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list1, list0idx, currSlice->listX[1], &currSlice->listXsize[1], 0);

        // long term handling
        int listltidx = 0;
        for (int i = 0; i < p_Dpb->ltref_frames_in_buffer; i++)
            fs_listlt[listltidx++] = p_Dpb->fs_ltref[i];

        std::qsort(fs_listlt, listltidx, sizeof(frame_store*), [](const void* a, const void* b) {
            int idx1 = (*(frame_store**)a)->LongTermFrameIdx;
            int idx2 = (*(frame_store**)b)->LongTermFrameIdx;
            return (idx1 < idx2) ? -1 : (idx1 > idx2) ? 1 : 0;
        });

        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[1], &currSlice->listXsize[1], 1);

        delete []fs_list0;
        delete []fs_list1;
        delete []fs_listlt;
    }

    if ((currSlice->listXsize[0] == currSlice->listXsize[1]) && (currSlice->listXsize[0] > 1)) {
        // check if lists are identical, if yes swap first two elements of currSlice->listX[1]
        int diff = 0;
        for (int j = 0; j < currSlice->listXsize[0]; j++) {
            if (currSlice->listX[0][j] != currSlice->listX[1][j]) {
                diff = 1;
                break;
            }
        }
        if (!diff) {
            storable_picture *tmp_s = currSlice->listX[1][0];
            currSlice->listX[1][0] = currSlice->listX[1][1];
            currSlice->listX[1][1] = tmp_s;
        }
    }

    // set max size
    currSlice->listXsize[0] = min<int>(currSlice->listXsize[0], currSlice->num_ref_idx_l0_active_minus1 + 1);
    currSlice->listXsize[1] = min<int>(currSlice->listXsize[1], currSlice->num_ref_idx_l1_active_minus1 + 1);

    // set the unused list entries to NULL
    for (int i = currSlice->listXsize[0]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[0][i] = p_Vid->no_reference_picture;
    for (int i = currSlice->listXsize[1]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[1][i] = p_Vid->no_reference_picture;
}

void init_lists(slice_t *currSlice)
{
#if (MVC_EXTENSION_ENABLE)
    if (currSlice->view_id) {
        switch (currSlice->slice_type) {
        case P_slice: 
        case SP_slice:
            init_lists_p_slice_mvc(currSlice);
            return;
        case B_slice:
            init_lists_b_slice_mvc(currSlice);
            return;
        case I_slice: 
        case SI_slice: 
            init_lists_i_slice_mvc(currSlice);
            return;
        default:
            printf("Unsupported slice type\n");
            break;
        }
    } else
#endif
    {
        switch (currSlice->slice_type) {
        case P_slice:
        case SP_slice:
            init_lists_p_slice(currSlice);
            return;
        case B_slice:
            init_lists_b_slice(currSlice);
            return;
        case I_slice:
        case SI_slice:
            init_lists_i_slice(currSlice);
            return;
        default:
            printf("Unsupported slice type\n");
            break;
        }
    }
}


static storable_picture *get_short_term_pic(slice_t *currSlice, dpb_t *p_Dpb, int picNum)
{
    unsigned i;

    for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
        if (!currSlice->field_pic_flag) {
            if (p_Dpb->fs_ref[i]->is_reference == 3)
                if (!p_Dpb->fs_ref[i]->frame->is_long_term && p_Dpb->fs_ref[i]->frame->PicNum == picNum)
                    return p_Dpb->fs_ref[i]->frame;
        } else {
            if (p_Dpb->fs_ref[i]->is_reference & 1)
                if (!p_Dpb->fs_ref[i]->top_field->is_long_term && p_Dpb->fs_ref[i]->top_field->PicNum == picNum)
                    return p_Dpb->fs_ref[i]->top_field;
            if (p_Dpb->fs_ref[i]->is_reference & 2)
                if (!p_Dpb->fs_ref[i]->bottom_field->is_long_term && p_Dpb->fs_ref[i]->bottom_field->PicNum == picNum)
                    return p_Dpb->fs_ref[i]->bottom_field;
        }
    }

    return currSlice->p_Vid->no_reference_picture;
}

static storable_picture *get_long_term_pic(slice_t *currSlice, dpb_t *p_Dpb, int LongtermPicNum)
{
    uint32_t i;

    for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
        if (!currSlice->field_pic_flag) {
            if (p_Dpb->fs_ltref[i]->is_reference == 3)
                if (p_Dpb->fs_ltref[i]->frame->is_long_term && p_Dpb->fs_ltref[i]->frame->LongTermPicNum == LongtermPicNum)
                    return p_Dpb->fs_ltref[i]->frame;
        } else {
            if (p_Dpb->fs_ltref[i]->is_reference & 1)
                if (p_Dpb->fs_ltref[i]->top_field->is_long_term && p_Dpb->fs_ltref[i]->top_field->LongTermPicNum == LongtermPicNum)
                    return p_Dpb->fs_ltref[i]->top_field;
            if (p_Dpb->fs_ltref[i]->is_reference & 2)
                if (p_Dpb->fs_ltref[i]->bottom_field->is_long_term && p_Dpb->fs_ltref[i]->bottom_field->LongTermPicNum == LongtermPicNum)
                    return p_Dpb->fs_ltref[i]->bottom_field;
        }
    }

    return NULL;
}

#if (MVC_EXTENSION_ENABLE)
void reorder_short_term(slice_t *currSlice, int cur_list, int num_ref_idx_lX_active_minus1, int picNumLX, int *refIdxLX, int currViewID)
{
    storable_picture **RefPicListX = currSlice->listX[cur_list]; 
    int cIdx, nIdx;

    for (cIdx = num_ref_idx_lX_active_minus1 + 1; cIdx > *refIdxLX; cIdx--)
        RefPicListX[cIdx] = RefPicListX[cIdx - 1];
    RefPicListX[(*refIdxLX)++] = get_short_term_pic(currSlice, currSlice->p_Dpb, picNumLX);
    nIdx = *refIdxLX;
    for (cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1 + 1; cIdx++) {
        if (RefPicListX[cIdx]) {
            if (RefPicListX[cIdx]->is_long_term ||
                RefPicListX[cIdx]->PicNum != picNumLX ||
                (currViewID != -1 && RefPicListX[cIdx]->layer_id != currViewID))
                RefPicListX[nIdx++] = RefPicListX[cIdx];
        }
    }
}

void reorder_long_term(slice_t *currSlice, storable_picture **RefPicListX, int num_ref_idx_lX_active_minus1, int LongTermPicNum, int *refIdxLX, int currViewID)
{
    int cIdx, nIdx;

    for (cIdx = num_ref_idx_lX_active_minus1 + 1; cIdx > *refIdxLX; cIdx--)
        RefPicListX[cIdx] = RefPicListX[cIdx - 1];
    RefPicListX[(*refIdxLX)++] = get_long_term_pic(currSlice, currSlice->p_Dpb, LongTermPicNum);
    nIdx = *refIdxLX;
    for (cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1 + 1; cIdx++) {
        if (RefPicListX[cIdx]) {
            if (!RefPicListX[cIdx]->is_long_term ||
                 RefPicListX[cIdx]->LongTermPicNum != LongTermPicNum ||
                (currViewID != -1 && RefPicListX[cIdx]->layer_id != currViewID))
                RefPicListX[nIdx++] = RefPicListX[cIdx];
        }
    }
}
#endif

static void reorder_ref_pic_list(slice_t *currSlice, int cur_list)
{
    uint8_t  *modification_of_pic_nums_idc = currSlice->modification_of_pic_nums_idc[cur_list];
    uint32_t *abs_diff_pic_num_minus1      = currSlice->abs_diff_pic_num_minus1     [cur_list];
    uint32_t *long_term_pic_num            = currSlice->long_term_pic_num           [cur_list];

    int num_ref_idx_lX_active_minus1 = 
        cur_list == 0 ? currSlice->num_ref_idx_l0_active_minus1
                      : currSlice->num_ref_idx_l1_active_minus1;

    int picNumLXNoWrap, picNumLXPred, picNumLX;
    int refIdxLX = 0;
    int i;

    picNumLXPred = currSlice->CurrPicNum;

    for (i = 0; modification_of_pic_nums_idc[i] != 3; i++) {
        if (modification_of_pic_nums_idc[i] > 3)
            error ("Invalid modification_of_pic_nums_idc command", 500);

        if (modification_of_pic_nums_idc[i] < 2) {
            if (modification_of_pic_nums_idc[i] == 0) {
                if (picNumLXPred < abs_diff_pic_num_minus1[i] + 1)
                    picNumLXNoWrap = picNumLXPred - (abs_diff_pic_num_minus1[i] + 1) + currSlice->MaxPicNum;
                else
                    picNumLXNoWrap = picNumLXPred - (abs_diff_pic_num_minus1[i] + 1);
            } else {
                if (picNumLXPred + (abs_diff_pic_num_minus1[i] + 1) >= currSlice->MaxPicNum)
                    picNumLXNoWrap = picNumLXPred + (abs_diff_pic_num_minus1[i] + 1) - currSlice->MaxPicNum;
                else
                    picNumLXNoWrap = picNumLXPred + (abs_diff_pic_num_minus1[i] + 1);
            }
            picNumLXPred = picNumLXNoWrap;

            if (picNumLXNoWrap > currSlice->CurrPicNum)
                picNumLX = picNumLXNoWrap - currSlice->MaxPicNum;
            else
                picNumLX = picNumLXNoWrap;

#if (MVC_EXTENSION_ENABLE)
            reorder_short_term(currSlice, cur_list, num_ref_idx_lX_active_minus1, picNumLX, &refIdxLX, -1);
#endif
        } else {
#if (MVC_EXTENSION_ENABLE)
            reorder_long_term(currSlice, currSlice->listX[cur_list], num_ref_idx_lX_active_minus1, long_term_pic_num[i], &refIdxLX, -1);
#endif
        }
    }
    // that's a definition
    currSlice->listXsize[cur_list] = num_ref_idx_lX_active_minus1 + 1;
}

void reorder_lists(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;

    if (currSlice->slice_type != I_slice && currSlice->slice_type != SI_slice) {
        if (currSlice->ref_pic_list_modification_flag_l0)
            reorder_ref_pic_list(currSlice, LIST_0);
        if (p_Vid->no_reference_picture == currSlice->listX[0][currSlice->num_ref_idx_l0_active_minus1]) {
            if (p_Vid->non_conforming_stream)
                printf("RefPicList0[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_l0_active_minus1);
            else
                error("RefPicList0[ num_ref_idx_l0_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
        }
        // that's a definition
        currSlice->listXsize[0] = (char) currSlice->num_ref_idx_l0_active_minus1 + 1;
    }

    if (currSlice->slice_type == B_slice) {
        if (currSlice->ref_pic_list_modification_flag_l1)
            reorder_ref_pic_list(currSlice, LIST_1);
        if (p_Vid->no_reference_picture == currSlice->listX[1][currSlice->num_ref_idx_l1_active_minus1]) {
            if (p_Vid->non_conforming_stream)
                printf("RefPicList1[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_l1_active_minus1);
            else
                error("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
        }
        // that's a definition
        currSlice->listXsize[1] = (char) currSlice->num_ref_idx_l1_active_minus1 + 1;
    }
}

/*!
 ************************************************************************
 * \brief
 *    Initialize listX[2..5] from lists 0 and 1
 *    listX[2]: list0 for current_field==top
 *    listX[3]: list1 for current_field==top
 *    listX[4]: list0 for current_field==bottom
 *    listX[5]: list1 for current_field==bottom
 *
 ************************************************************************
 */
void init_mbaff_lists(VideoParameters *p_Vid, slice_t *currSlice)
{
    unsigned j;
    int i;

    for (i = 2; i < 6; i++) {
        for (j = 0; j < MAX_LIST_SIZE; j++)
            currSlice->listX[i][j] = p_Vid->no_reference_picture;
        currSlice->listXsize[i]=0;
    }

    for (i = 0; i < currSlice->listXsize[0]; i++) {
        currSlice->listX[2][2*i  ] = currSlice->listX[0][i]->top_field;
        currSlice->listX[2][2*i+1] = currSlice->listX[0][i]->bottom_field;
        currSlice->listX[4][2*i  ] = currSlice->listX[0][i]->bottom_field;
        currSlice->listX[4][2*i+1] = currSlice->listX[0][i]->top_field;
    }
    currSlice->listXsize[2] = currSlice->listXsize[4] = currSlice->listXsize[0] * 2;

    for (i = 0; i < currSlice->listXsize[1]; i++) {
        currSlice->listX[3][2*i  ] = currSlice->listX[1][i]->top_field;
        currSlice->listX[3][2*i+1] = currSlice->listX[1][i]->bottom_field;
        currSlice->listX[5][2*i  ] = currSlice->listX[1][i]->bottom_field;
        currSlice->listX[5][2*i+1] = currSlice->listX[1][i]->top_field;
    }
    currSlice->listXsize[3] = currSlice->listXsize[5] = currSlice->listXsize[1] * 2;
}




void decoded_picture_buffer_t::check_num_ref()
{
    if (this->ltref_frames_in_buffer + this->ref_frames_in_buffer > max(1, this->num_ref_frames))
        error("Max. number of reference frames exceeded. Invalid stream.", 500);
}

void decoded_picture_buffer_t::mm_unmark_short_term_for_reference(storable_picture* p, int difference_of_pic_nums_minus1)
{
    int picNumX = p->get_pic_num_x(difference_of_pic_nums_minus1);

    for (int i = 0; i < this->ref_frames_in_buffer; i++) {
        if (p->structure == FRAME) {
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
        if (p->structure == FRAME) {
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

    if (p->structure == FRAME) {
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
        if (p->structure == TOP_FIELD) {
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
    int picNumX = p->get_pic_num_x(difference_of_pic_nums_minus1);

    // remove frames/fields with same long_term_frame_idx
    if (p->structure == FRAME)
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
    if (p->structure == FRAME)
        this->unmark_long_term_frame_for_reference_by_frame_idx(long_term_frame_idx);
    else
        this->unmark_long_term_field_for_reference_by_frame_idx(p->structure, long_term_frame_idx, 1, p->PicNum, 0);

    p->is_long_term = 1;
    p->LongTermFrameIdx = long_term_frame_idx;
}

void decoded_picture_buffer_t::adaptive_memory_management(storable_picture* p)
{
    drpm_t* tmp_drpm;
    VideoParameters *p_Vid = this->p_Vid;

    p_Vid->last_has_mmco_5 = 0;

    assert(!p->idr_flag);
    assert(p->adaptive_ref_pic_buffering_flag);

    while (p->dec_ref_pic_marking_buffer) {
        tmp_drpm = p->dec_ref_pic_marking_buffer;
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
        p->dec_ref_pic_marking_buffer = tmp_drpm->Next;
        ::free(tmp_drpm);
    }
    if (p_Vid->last_has_mmco_5) {
        p->PicNum = p->frame_num = 0;

        switch (p->structure) {
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
        if (p->view_id == 0) {
            p_Vid->p_Dpb_layer[0]->flush();
            p_Vid->p_Dpb_layer[1]->flush();
        } else
            this->flush();
#endif
    }
}

void decoded_picture_buffer_t::sliding_window_memory_management(storable_picture* p)
{
    assert(!p->idr_flag);

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
    p_Vid->last_pic_bottom_field = p->structure == BOTTOM_FIELD;

    if (p->idr_flag) {
        this->idr_memory_management(p);
        // picture error concealment
        memset(p_Vid->pocs_in_dpb, 0, sizeof(int)*100);
    } else {
        // adaptive memory management
        if (p->used_for_reference && p->adaptive_ref_pic_buffering_flag)
            this->adaptive_memory_management(p);
    }

    if (p->structure == TOP_FIELD || p->structure == BOTTOM_FIELD) {
        // check for frame store with same pic_number
        if (this->last_picture) {
            if ((int)this->last_picture->FrameNum == p->PicNum) {
                if ((p->structure == TOP_FIELD && this->last_picture->is_used == 2) ||
                    (p->structure == BOTTOM_FIELD && this->last_picture->is_used == 1)) {
                    if (( p->used_for_reference && this->last_picture->is_orig_reference != 0) ||
                        (!p->used_for_reference && this->last_picture->is_orig_reference == 0)) {
                        insert_picture_in_dpb(p_Vid, this->last_picture, p);
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
    if (!p->idr_flag && p->used_for_reference && !p->adaptive_ref_pic_buffering_flag)
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
    insert_picture_in_dpb(p_Vid, this->fs[this->used_size],p);

    // picture error concealment
    if (p->idr_flag)
        p_Vid->earlier_missing_poc = 0;

    if (p->structure != FRAME)
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
