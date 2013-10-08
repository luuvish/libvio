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


