#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "dpb.h"
#include "memalloc.h"
#include "output.h"

using vio::h264::mb_t;

#include "erc_api.h"

static inline int is_short_ref(storable_picture *s)
{
  	return ((s->used_for_reference) && (!(s->is_long_term)));
}

static inline int is_long_ref(storable_picture *s)
{
  	return ((s->used_for_reference) && (s->is_long_term));
}



void idr_memory_management(dpb_t *p_Dpb, storable_picture* p)
{
  	if (p->no_output_of_prior_pics_flag) {
    	// free all stored pictures
    	for (int i = 0; i < p_Dpb->used_size; i++) {
      		// reset all reference settings
      		delete p_Dpb->fs[i];
      		p_Dpb->fs[i] = new frame_store {};
    	}
    	for (int i = 0; i < p_Dpb->ref_frames_in_buffer; i++)
      		p_Dpb->fs_ref[i] = NULL;
    	for (int i = 0; i < p_Dpb->ltref_frames_in_buffer; i++)
      		p_Dpb->fs_ltref[i] = NULL;
    	p_Dpb->used_size = 0;
  	} else
    	flush_dpb(p_Dpb);
  	p_Dpb->last_picture = NULL;

  	update_ref_list(p_Dpb);
  	update_ltref_list(p_Dpb);
  	p_Dpb->last_output_poc = INT_MIN;

  	if (p->long_term_reference_flag) {
    	p_Dpb->max_long_term_pic_idx = 0;
    	p->is_long_term              = 1;
    	p->LongTermFrameIdx          = 0;
  	} else {
    	p_Dpb->max_long_term_pic_idx = -1;
    	p->is_long_term              = 0;
  	}

#if (MVC_EXTENSION_ENABLE)
  	p_Dpb->last_output_view_id = -1;
#endif
}


static int is_short_term_reference(frame_store* fs)
{
    if (fs->is_used == 3) { // frame
        if (fs->frame->used_for_reference && !fs->frame->is_long_term)
            return 1;
    }

    if (fs->is_used & 1) { // top field
        if (fs->top_field) {
            if (fs->top_field->used_for_reference && !fs->top_field->is_long_term)
                return 1;
        }
    }

    if (fs->is_used & 2) { // bottom field
        if (fs->bottom_field) {
            if (fs->bottom_field->used_for_reference && !fs->bottom_field->is_long_term)
                return 1;
        }
    }
    return 0;
}

static int is_long_term_reference(frame_store* fs)
{
    if (fs->is_used == 3) { // frame
        if (fs->frame->used_for_reference && fs->frame->is_long_term)
            return 1;
    }

    if (fs->is_used & 1) { // top field
        if (fs->top_field) {
            if (fs->top_field->used_for_reference && fs->top_field->is_long_term)
                return 1;
        }
    }

    if (fs->is_used & 2) { // bottom field
        if (fs->bottom_field) {
            if (fs->bottom_field->used_for_reference && fs->bottom_field->is_long_term)
                return 1;
        }
    }
    return 0;
}

void update_ref_list(dpb_t *p_Dpb)
{
    int i, j;
    for (i = 0, j = 0; i < p_Dpb->used_size; i++) {
        if (is_short_term_reference(p_Dpb->fs[i]))
            p_Dpb->fs_ref[j++] = p_Dpb->fs[i];
    }

    p_Dpb->ref_frames_in_buffer = j;

    while (j < p_Dpb->size)
        p_Dpb->fs_ref[j++] = NULL;
}

void update_ltref_list(dpb_t *p_Dpb)
{
    int i, j;
    for (i = 0, j = 0; i < p_Dpb->used_size; i++) {
        if (is_long_term_reference(p_Dpb->fs[i]))
            p_Dpb->fs_ltref[j++] = p_Dpb->fs[i];
    }

    p_Dpb->ltref_frames_in_buffer = j;

    while (j < p_Dpb->size)
        p_Dpb->fs_ltref[j++] = NULL;
}


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

  	int (*is_ref)(storable_picture *s) = (long_term) ? is_long_ref : is_short_ref;

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

    unsigned int i;

    int list0idx = 0;
    int listltidx = 0;

    frame_store **fs_list0;
    frame_store **fs_listlt;

#if (MVC_EXTENSION_ENABLE)
    currSlice->listinterviewidx0 = 0;
    currSlice->listinterviewidx1 = 0;
#endif

    if (!currSlice->field_pic_flag) {
        for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used == 3) {
                if (p_Dpb->fs_ref[i]->frame->used_for_reference && !p_Dpb->fs_ref[i]->frame->is_long_term)
                    currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
        }
        // order list 0 by PicNum
        qsort((void *)currSlice->listX[0], list0idx, sizeof(storable_picture*), compare_pic_by_pic_num_desc);
        currSlice->listXsize[0] = (char) list0idx;

        // long term handling
        for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ltref[i]->is_used == 3) {
                if (p_Dpb->fs_ltref[i]->frame->is_long_term)
                    currSlice->listX[0][list0idx++] = p_Dpb->fs_ltref[i]->frame;
            }
        }
        qsort((void *)&currSlice->listX[0][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(storable_picture*), compare_pic_by_lt_pic_num_asc);
        currSlice->listXsize[0] = (char) list0idx;
    } else {
        fs_list0 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
        if (NULL == fs_list0)
            no_mem_exit("init_lists: fs_list0");
        fs_listlt = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
        if (NULL == fs_listlt)
            no_mem_exit("init_lists: fs_listlt");

        for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_reference)
                fs_list0[list0idx++] = p_Dpb->fs_ref[i];
        }

        qsort((void *)fs_list0, list0idx, sizeof(frame_store*), compare_fs_by_frame_num_desc);

        currSlice->listXsize[0] = 0;
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);

        // long term handling
        for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++)
            fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];

        qsort((void *)fs_listlt, listltidx, sizeof(frame_store*), compare_fs_by_lt_pic_idx_asc);

        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);

        free(fs_list0);
        free(fs_listlt);
    }
    currSlice->listXsize[1] = 0;

    // set max size
    currSlice->listXsize[0] = (char) min<int>(currSlice->listXsize[0], currSlice->num_ref_idx_l0_active_minus1 + 1);
    currSlice->listXsize[1] = (char) min<int>(currSlice->listXsize[1], currSlice->num_ref_idx_l1_active_minus1 + 1);

    // set the unused list entries to NULL
    for (i = currSlice->listXsize[0]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[0][i] = p_Vid->no_reference_picture;
    for (i = currSlice->listXsize[1]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[1][i] = p_Vid->no_reference_picture;
}

static void init_lists_b_slice(slice_t *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    dpb_t *p_Dpb = currSlice->p_Dpb;

    int list0idx = 0;
    int list0idx_1 = 0;
    int listltidx = 0;
    int i, j;

    frame_store **fs_list0;
    frame_store **fs_list1;
    frame_store **fs_listlt;

#if (MVC_EXTENSION_ENABLE)
    currSlice->listinterviewidx0 = 0;
    currSlice->listinterviewidx1 = 0;
#endif

    // B-slice_t
    if (!currSlice->field_pic_flag) {
        for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used == 3) {
                if (p_Dpb->fs_ref[i]->frame->used_for_reference && !p_Dpb->fs_ref[i]->frame->is_long_term) {
                    if (currSlice->framepoc >= p_Dpb->fs_ref[i]->frame->poc) //!KS use >= for error concealment
                        currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
                }
            }
        }
        qsort((void *)currSlice->listX[0], list0idx, sizeof(storable_picture*), compare_pic_by_poc_desc);

        //get the backward reference picture (POC>current POC) in list0;
        list0idx_1 = list0idx;
        for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used == 3) {
                if (p_Dpb->fs_ref[i]->frame->used_for_reference && !p_Dpb->fs_ref[i]->frame->is_long_term) {
                    if (currSlice->framepoc < p_Dpb->fs_ref[i]->frame->poc)
                        currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
                }
            }
        }
        qsort((void *)&currSlice->listX[0][list0idx_1], list0idx-list0idx_1, sizeof(storable_picture*), compare_pic_by_poc_asc);

        for (j = 0; j < list0idx_1; j++)
            currSlice->listX[1][list0idx - list0idx_1 + j] = currSlice->listX[0][j];
        for (j = list0idx_1; j < list0idx; j++)
            currSlice->listX[1][j - list0idx_1] = currSlice->listX[0][j];

        currSlice->listXsize[0] = currSlice->listXsize[1] = list0idx;

        // long term handling
        for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ltref[i]->is_used == 3) {
                if (p_Dpb->fs_ltref[i]->frame->is_long_term) {
                    currSlice->listX[0][list0idx]   = p_Dpb->fs_ltref[i]->frame;
                    currSlice->listX[1][list0idx++] = p_Dpb->fs_ltref[i]->frame;
                }
            }
        }
        qsort((void *)&currSlice->listX[0][(int)currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(storable_picture*), compare_pic_by_lt_pic_num_asc);
        qsort((void *)&currSlice->listX[1][(int)currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(storable_picture*), compare_pic_by_lt_pic_num_asc);
        currSlice->listXsize[0] = currSlice->listXsize[1] = (char) list0idx;
    } else {
        fs_list0 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
        if (NULL == fs_list0)
            no_mem_exit("init_lists: fs_list0");
        fs_list1 = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
        if (NULL == fs_list1)
            no_mem_exit("init_lists: fs_list1");
        fs_listlt = (frame_store **)calloc(p_Dpb->size, sizeof (frame_store*));
        if (NULL == fs_listlt)
            no_mem_exit("init_lists: fs_listlt");

        currSlice->listXsize[0] = 0;
        currSlice->listXsize[1] = 1;

        for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used) {
                if (currSlice->ThisPOC >= p_Dpb->fs_ref[i]->poc)
                    fs_list0[list0idx++] = p_Dpb->fs_ref[i];
            }
        }
        qsort((void *)fs_list0, list0idx, sizeof(frame_store*), compare_fs_by_poc_desc);
        list0idx_1 = list0idx;
        for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
            if (p_Dpb->fs_ref[i]->is_used) {
                if (currSlice->ThisPOC < p_Dpb->fs_ref[i]->poc)
                    fs_list0[list0idx++] = p_Dpb->fs_ref[i];
            }
        }
        qsort((void *)&fs_list0[list0idx_1], list0idx-list0idx_1, sizeof(frame_store*), compare_fs_by_poc_asc);

        for (j = 0; j < list0idx_1; j++)
            fs_list1[list0idx - list0idx_1 + j] = fs_list0[j];
        for (j = list0idx_1; j < list0idx; j++)
            fs_list1[j - list0idx_1] = fs_list0[j];

        currSlice->listXsize[0] = 0;
        currSlice->listXsize[1] = 0;
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_list1, list0idx, currSlice->listX[1], &currSlice->listXsize[1], 0);

        // long term handling
        for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++)
            fs_listlt[listltidx++] = p_Dpb->fs_ltref[i];

        qsort((void *)fs_listlt, listltidx, sizeof(frame_store*), compare_fs_by_lt_pic_idx_asc);

        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);
        gen_pic_list_from_frame_list(currSlice->bottom_field_flag, fs_listlt, listltidx, currSlice->listX[1], &currSlice->listXsize[1], 1);

        free(fs_list0);
        free(fs_list1);
        free(fs_listlt);
    }

    if ((currSlice->listXsize[0] == currSlice->listXsize[1]) && (currSlice->listXsize[0] > 1)) {
        // check if lists are identical, if yes swap first two elements of currSlice->listX[1]
        int diff = 0;
        for (j = 0; j < currSlice->listXsize[0]; j++) {
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
    for (i = currSlice->listXsize[0]; i < MAX_LIST_SIZE; i++)
        currSlice->listX[0][i] = p_Vid->no_reference_picture;
    for (i = currSlice->listXsize[1]; i < MAX_LIST_SIZE; i++)
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




static void check_num_ref(dpb_t *p_Dpb)
{
  	if ((int)(p_Dpb->ltref_frames_in_buffer +  p_Dpb->ref_frames_in_buffer ) > max(1, p_Dpb->num_ref_frames))
    	error ("Max. number of reference frames exceeded. Invalid stream.", 500);
}

static int get_pic_num_x(storable_picture *p, int difference_of_pic_nums_minus1)
{
  	int currPicNum;

  	if (p->structure == FRAME)
    	currPicNum = p->frame_num;
  	else
    	currPicNum = 2 * p->frame_num + 1;

  	return currPicNum - (difference_of_pic_nums_minus1 + 1);
}

static void mm_unmark_short_term_for_reference(dpb_t *p_Dpb, storable_picture *p, int difference_of_pic_nums_minus1)
{
  	int picNumX = get_pic_num_x(p, difference_of_pic_nums_minus1);
  	uint32_t i;

  	for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
    	if (p->structure == FRAME) {
      		if (p_Dpb->fs_ref[i]->is_reference == 3 && p_Dpb->fs_ref[i]->is_long_term == 0) {
        		if (p_Dpb->fs_ref[i]->frame->PicNum == picNumX) {
          			unmark_for_reference(p_Dpb->fs_ref[i]);
          			return;
        		}
      		}
    	} else {
      		if ((p_Dpb->fs_ref[i]->is_reference & 1) && !(p_Dpb->fs_ref[i]->is_long_term & 1)) {
        		if (p_Dpb->fs_ref[i]->top_field->PicNum == picNumX) {
          			p_Dpb->fs_ref[i]->top_field->used_for_reference = 0;
          			p_Dpb->fs_ref[i]->is_reference &= 2;
          			if (p_Dpb->fs_ref[i]->is_used == 3)
            			p_Dpb->fs_ref[i]->frame->used_for_reference = 0;
          			return;
        		}
      		}
      		if ((p_Dpb->fs_ref[i]->is_reference & 2) && !(p_Dpb->fs_ref[i]->is_long_term & 2)) {
        		if (p_Dpb->fs_ref[i]->bottom_field->PicNum == picNumX) {
          			p_Dpb->fs_ref[i]->bottom_field->used_for_reference = 0;
          			p_Dpb->fs_ref[i]->is_reference &= 1;
          			if (p_Dpb->fs_ref[i]->is_used == 3)
            			p_Dpb->fs_ref[i]->frame->used_for_reference = 0;
          			return;
        		}
      		}
    	}
  	}
}

static void mm_unmark_long_term_for_reference(dpb_t *p_Dpb, storable_picture *p, int long_term_pic_num)
{
  	uint32_t i;
  	for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
    	if (p->structure == FRAME) {
      		if (p_Dpb->fs_ltref[i]->is_reference == 3 && p_Dpb->fs_ltref[i]->is_long_term == 3) {
        		if (p_Dpb->fs_ltref[i]->frame->LongTermPicNum == long_term_pic_num)
          			unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      		}
    	} else {
      		if ((p_Dpb->fs_ltref[i]->is_reference & 1) && (p_Dpb->fs_ltref[i]->is_long_term & 1)) {
        		if (p_Dpb->fs_ltref[i]->top_field->LongTermPicNum == long_term_pic_num) {
          			p_Dpb->fs_ltref[i]->top_field->used_for_reference = 0;
          			p_Dpb->fs_ltref[i]->top_field->is_long_term = 0;
          			p_Dpb->fs_ltref[i]->is_reference &= 2;
          			p_Dpb->fs_ltref[i]->is_long_term &= 2;
          			if (p_Dpb->fs_ltref[i]->is_used == 3) {
          			  	p_Dpb->fs_ltref[i]->frame->used_for_reference = 0;
          			  	p_Dpb->fs_ltref[i]->frame->is_long_term = 0;
          			}
          			return;
        		}
      		}
      		if ((p_Dpb->fs_ltref[i]->is_reference & 2) && (p_Dpb->fs_ltref[i]->is_long_term & 2)) {
        		if (p_Dpb->fs_ltref[i]->bottom_field->LongTermPicNum == long_term_pic_num) {
          			p_Dpb->fs_ltref[i]->bottom_field->used_for_reference = 0;
          			p_Dpb->fs_ltref[i]->bottom_field->is_long_term = 0;
          			p_Dpb->fs_ltref[i]->is_reference &= 1;
          			p_Dpb->fs_ltref[i]->is_long_term &= 1;
          			if (p_Dpb->fs_ltref[i]->is_used == 3) {
          			  	p_Dpb->fs_ltref[i]->frame->used_for_reference = 0;
          			  	p_Dpb->fs_ltref[i]->frame->is_long_term = 0;
          			}
          			return;
        		}
      		}
    	}
  	}
}

static void unmark_long_term_frame_for_reference_by_frame_idx(dpb_t *p_Dpb, int long_term_frame_idx)
{
  	uint32_t i;
  	for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
    	if (p_Dpb->fs_ltref[i]->LongTermFrameIdx == long_term_frame_idx)
      		unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
  	}
}

static void unmark_long_term_field_for_reference_by_frame_idx(dpb_t *p_Dpb, PictureStructure structure, int long_term_frame_idx, int mark_current, unsigned curr_frame_num, int curr_pic_num)
{
  	VideoParameters *p_Vid = p_Dpb->p_Vid;
  	unsigned i;

  	assert(structure != FRAME);
  	if (curr_pic_num < 0)
    	curr_pic_num += (2 * p_Vid->active_sps->MaxFrameNum);

  	for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
    	if (p_Dpb->fs_ltref[i]->LongTermFrameIdx == long_term_frame_idx) {
      		if (structure == TOP_FIELD) {
        		if (p_Dpb->fs_ltref[i]->is_long_term == 3)
          			unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
        		else if (p_Dpb->fs_ltref[i]->is_long_term == 1)
        			unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      			else if (mark_current) {
      				if (p_Dpb->last_picture) {
        				if (p_Dpb->last_picture != p_Dpb->fs_ltref[i] ||
        				 	p_Dpb->last_picture->FrameNum != curr_frame_num)
          					unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      				} else
        				unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
    			} else if (p_Dpb->fs_ltref[i]->FrameNum != (unsigned)(curr_pic_num >> 1))
        			unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      		}
      		if (structure == BOTTOM_FIELD) {
        		if (p_Dpb->fs_ltref[i]->is_long_term == 3)
          			unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
        		else if (p_Dpb->fs_ltref[i]->is_long_term == 2)
            		unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      			else if (mark_current) {
      				if (p_Dpb->last_picture) {
        				if (p_Dpb->last_picture != p_Dpb->fs_ltref[i] ||
        				 	p_Dpb->last_picture->FrameNum != curr_frame_num)
          					unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      				} else
        				unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
    			} else if (p_Dpb->fs_ltref[i]->FrameNum != (unsigned)(curr_pic_num >> 1))
        			unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
      		}
    	}
  	}
}

static void mark_pic_long_term(dpb_t *p_Dpb, storable_picture* p, int long_term_frame_idx, int picNumX)
{
  	uint32_t i;
  	int add_top, add_bottom;

  	if (p->structure == FRAME) {
    	for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
      		if (p_Dpb->fs_ref[i]->is_reference == 3) {
        		if (!p_Dpb->fs_ref[i]->frame->is_long_term &&
        			 p_Dpb->fs_ref[i]->frame->PicNum == picNumX) {
          			p_Dpb->fs_ref[i]->LongTermFrameIdx =
          			p_Dpb->fs_ref[i]->frame->LongTermFrameIdx = long_term_frame_idx;
          			p_Dpb->fs_ref[i]->frame->LongTermPicNum = long_term_frame_idx;
          			p_Dpb->fs_ref[i]->frame->is_long_term = 1;

          			if (p_Dpb->fs_ref[i]->top_field && p_Dpb->fs_ref[i]->bottom_field) {
            			p_Dpb->fs_ref[i]->top_field->LongTermFrameIdx =
            			p_Dpb->fs_ref[i]->bottom_field->LongTermFrameIdx = long_term_frame_idx;
            			p_Dpb->fs_ref[i]->top_field->LongTermPicNum = long_term_frame_idx;
            			p_Dpb->fs_ref[i]->bottom_field->LongTermPicNum = long_term_frame_idx;
            			p_Dpb->fs_ref[i]->top_field->is_long_term =
            			p_Dpb->fs_ref[i]->bottom_field->is_long_term = 1;
          			}
          			p_Dpb->fs_ref[i]->is_long_term = 3;
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
    	for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
      		if (p_Dpb->fs_ref[i]->is_reference & 1) {
        		if (!p_Dpb->fs_ref[i]->top_field->is_long_term &&
        			 p_Dpb->fs_ref[i]->top_field->PicNum == picNumX) {
          			if (p_Dpb->fs_ref[i]->is_long_term &&
          			 	p_Dpb->fs_ref[i]->LongTermFrameIdx != long_term_frame_idx)
              			printf("Warning: assigning long_term_frame_idx different from other field\n");

          			p_Dpb->fs_ref[i]->LongTermFrameIdx =
          			p_Dpb->fs_ref[i]->top_field->LongTermFrameIdx = long_term_frame_idx;
          			p_Dpb->fs_ref[i]->top_field->LongTermPicNum = 2 * long_term_frame_idx + add_top;
          			p_Dpb->fs_ref[i]->top_field->is_long_term = 1;
          			p_Dpb->fs_ref[i]->is_long_term |= 1;
          			if (p_Dpb->fs_ref[i]->is_long_term == 3) {
            			p_Dpb->fs_ref[i]->frame->is_long_term = 1;
            			p_Dpb->fs_ref[i]->frame->LongTermFrameIdx =
            			p_Dpb->fs_ref[i]->frame->LongTermPicNum = long_term_frame_idx;
          			}
          			return;
        		}
      		}
      		if (p_Dpb->fs_ref[i]->is_reference & 2) {
        		if (!p_Dpb->fs_ref[i]->bottom_field->is_long_term &&
        			 p_Dpb->fs_ref[i]->bottom_field->PicNum == picNumX) {
          			if (p_Dpb->fs_ref[i]->is_long_term &&
          			 	p_Dpb->fs_ref[i]->LongTermFrameIdx != long_term_frame_idx)
              			printf("Warning: assigning long_term_frame_idx different from other field\n");

          			p_Dpb->fs_ref[i]->LongTermFrameIdx =
          			p_Dpb->fs_ref[i]->bottom_field->LongTermFrameIdx = long_term_frame_idx;
          			p_Dpb->fs_ref[i]->bottom_field->LongTermPicNum = 2 * long_term_frame_idx + add_bottom;
          			p_Dpb->fs_ref[i]->bottom_field->is_long_term = 1;
          			p_Dpb->fs_ref[i]->is_long_term |= 2;
          			if (p_Dpb->fs_ref[i]->is_long_term == 3) {
          			  	p_Dpb->fs_ref[i]->frame->is_long_term = 1;
          			  	p_Dpb->fs_ref[i]->frame->LongTermFrameIdx =
          			  	p_Dpb->fs_ref[i]->frame->LongTermPicNum = long_term_frame_idx;
          			}
          			return;
        		}
      		}
    	}
    	printf("Warning: reference field for long term marking not found\n");
  	}
}

static void mm_assign_long_term_frame_idx(dpb_t *p_Dpb, storable_picture* p, int difference_of_pic_nums_minus1, int long_term_frame_idx)
{
  	int picNumX = get_pic_num_x(p, difference_of_pic_nums_minus1);

  	// remove frames/fields with same long_term_frame_idx
  	if (p->structure == FRAME)
    	unmark_long_term_frame_for_reference_by_frame_idx(p_Dpb, long_term_frame_idx);
  	else {
    	unsigned i;
    	PictureStructure structure = FRAME;

    	for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
      		if (p_Dpb->fs_ref[i]->is_reference & 1) {
        		if (p_Dpb->fs_ref[i]->top_field->PicNum == picNumX) {
          			structure = TOP_FIELD;
          			break;
        		}
      		}
      		if (p_Dpb->fs_ref[i]->is_reference & 2) {
        		if (p_Dpb->fs_ref[i]->bottom_field->PicNum == picNumX) {
        		  	structure = BOTTOM_FIELD;
        		  	break;
        		}
      		}
    	}
    	if (structure == FRAME)
      		error("field for long term marking not found",200);

    	unmark_long_term_field_for_reference_by_frame_idx(p_Dpb, structure, long_term_frame_idx, 0, 0, picNumX);
  	}

  	mark_pic_long_term(p_Dpb, p, long_term_frame_idx, picNumX);
}

static void mm_update_max_long_term_frame_idx(dpb_t *p_Dpb, int max_long_term_frame_idx_plus1)
{
  	uint32_t i;

  	p_Dpb->max_long_term_pic_idx = max_long_term_frame_idx_plus1 - 1;

  	// check for invalid frames
  	for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++) {
    	if (p_Dpb->fs_ltref[i]->LongTermFrameIdx > p_Dpb->max_long_term_pic_idx)
      	unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
  	}
}

static void mm_unmark_all_short_term_for_reference (dpb_t *p_Dpb)
{
  	int i;
  	for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++)
    	unmark_for_reference(p_Dpb->fs_ref[i]);
  	update_ref_list(p_Dpb);
}

static void mm_unmark_all_long_term_for_reference (dpb_t *p_Dpb)
{
  	mm_update_max_long_term_frame_idx(p_Dpb, 0);
}

static void mm_mark_current_picture_long_term(dpb_t *p_Dpb, storable_picture *p, int long_term_frame_idx)
{
  	// remove long term pictures with same long_term_frame_idx
  	if (p->structure == FRAME)
    	unmark_long_term_frame_for_reference_by_frame_idx(p_Dpb, long_term_frame_idx);
  	else
    	unmark_long_term_field_for_reference_by_frame_idx(p_Dpb, p->structure, long_term_frame_idx, 1, p->PicNum, 0);

  	p->is_long_term = 1;
  	p->LongTermFrameIdx = long_term_frame_idx;
}

static void adaptive_memory_management(dpb_t *p_Dpb, storable_picture* p)
{
  	DecRefPicMarking_t *tmp_drpm;
  	VideoParameters *p_Vid = p_Dpb->p_Vid;

  	p_Vid->last_has_mmco_5 = 0;

  	assert(!p->idr_flag);
  	assert(p->adaptive_ref_pic_buffering_flag);

  	while (p->dec_ref_pic_marking_buffer) {
    	tmp_drpm = p->dec_ref_pic_marking_buffer;
    	switch (tmp_drpm->memory_management_control_operation) {
      	case 0:
        	if (tmp_drpm->Next != NULL)
          		error ("memory_management_control_operation = 0 not last operation in buffer", 500);
        	break;
      	case 1:
        	mm_unmark_short_term_for_reference(p_Dpb, p, tmp_drpm->difference_of_pic_nums_minus1);
        	update_ref_list(p_Dpb);
        	break;
      	case 2:
        	mm_unmark_long_term_for_reference(p_Dpb, p, tmp_drpm->long_term_pic_num);
        	update_ltref_list(p_Dpb);
        	break;
      	case 3:
        	mm_assign_long_term_frame_idx(p_Dpb, p, tmp_drpm->difference_of_pic_nums_minus1, tmp_drpm->long_term_frame_idx);
        	update_ref_list(p_Dpb);
        	update_ltref_list(p_Dpb);
        	break;
      	case 4:
	        mm_update_max_long_term_frame_idx (p_Dpb, tmp_drpm->max_long_term_frame_idx_plus1);
    	    update_ltref_list(p_Dpb);
        	break;
      	case 5:
        	mm_unmark_all_short_term_for_reference(p_Dpb);
        	mm_unmark_all_long_term_for_reference(p_Dpb);
        	p_Vid->last_has_mmco_5 = 1;
        	break;
      	case 6:
        	mm_mark_current_picture_long_term(p_Dpb, p, tmp_drpm->long_term_frame_idx);
        	check_num_ref(p_Dpb);
        	break;
      	default:
        	error("invalid memory_management_control_operation in buffer", 500);
    	}
    	p->dec_ref_pic_marking_buffer = tmp_drpm->Next;
    	free(tmp_drpm);
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
	        p->poc = min (p->top_poc, p->bottom_poc);
	        break;
    	}
#if (MVC_EXTENSION_ENABLE)
    	if (p->view_id == 0) {
      		flush_dpb(p_Vid->p_Dpb_layer[0]);
      		flush_dpb(p_Vid->p_Dpb_layer[1]);
    	} else
      		flush_dpb(p_Dpb);
#endif
  	}
}

static void sliding_window_memory_management(dpb_t *p_Dpb, storable_picture* p)
{
  	uint32_t i;

  	assert(!p->idr_flag);

  	// if this is a reference pic with sliding window, unmark first ref frame
  	if (p_Dpb->ref_frames_in_buffer == max(1, p_Dpb->num_ref_frames) - p_Dpb->ltref_frames_in_buffer) {
    	for (i = 0; i < p_Dpb->used_size; i++) {
      		if (p_Dpb->fs[i]->is_reference && !p_Dpb->fs[i]->is_long_term) {
        		unmark_for_reference(p_Dpb->fs[i]);
        		update_ref_list(p_Dpb);
        		break;
      		}
    	}
  	}

  	p->is_long_term = 0;
}

void store_picture_in_dpb(dpb_t *p_Dpb, storable_picture* p)
{
  	VideoParameters *p_Vid = p_Dpb->p_Vid;
  	unsigned i;
  	int poc, pos;
  	// picture error concealment
  
  	// diagnostics
  	// if frame, check for new store,
  	assert(p != NULL);

  	p_Vid->last_has_mmco_5       = 0;
  	p_Vid->last_pic_bottom_field = p->structure == BOTTOM_FIELD;

  	if (p->idr_flag) {
    	idr_memory_management(p_Dpb, p);
  		// picture error concealment
    	memset(p_Vid->pocs_in_dpb, 0, sizeof(int)*100);
  	} else {
    	// adaptive memory management
    	if (p->used_for_reference && p->adaptive_ref_pic_buffering_flag)
      		adaptive_memory_management(p_Dpb, p);
  	}

  	if (p->structure == TOP_FIELD || p->structure == BOTTOM_FIELD) {
    	// check for frame store with same pic_number
    	if (p_Dpb->last_picture) {
      		if ((int)p_Dpb->last_picture->FrameNum == p->PicNum) {
        		if ((p->structure == TOP_FIELD && p_Dpb->last_picture->is_used == 2) ||
        			(p->structure == BOTTOM_FIELD && p_Dpb->last_picture->is_used == 1)) {
          			if (( p->used_for_reference && p_Dpb->last_picture->is_orig_reference != 0) ||
              			(!p->used_for_reference && p_Dpb->last_picture->is_orig_reference == 0)) {
            			insert_picture_in_dpb(p_Vid, p_Dpb->last_picture, p);
            			update_ref_list(p_Dpb);
            			update_ltref_list(p_Dpb);
            			p_Dpb->last_picture = NULL;
            			return;
          			}
        		}
      		}
    	}
  	}

  	// this is a frame or a field which has no stored complementary field

  	// sliding window, if necessary
  	if (!p->idr_flag && p->used_for_reference && !p->adaptive_ref_pic_buffering_flag)
    	sliding_window_memory_management(p_Dpb, p);

  	// picture error concealment
  	if (p_Vid->conceal_mode != 0) {
    	for (i = 0; i < p_Dpb->size; i++)
      		if (p_Dpb->fs[i]->is_reference)
        		p_Dpb->fs[i]->concealment_reference = 1;
  	}

  	// first try to remove unused frames
  	if (p_Dpb->used_size == p_Dpb->size) {
#if (DISABLE_ERC == 0)
    	// picture error concealment
    	if (p_Vid->conceal_mode != 0)
      		conceal_non_ref_pics(p_Dpb, 2);
#endif
    	remove_unused_frame_from_dpb(p_Dpb);

#if (DISABLE_ERC == 0)
    	if(p_Vid->conceal_mode != 0)
      		sliding_window_poc_management(p_Dpb, p);
#endif
  	}

  	// then output frames until one can be removed
  	while (p_Dpb->used_size == p_Dpb->size) {
    	// non-reference frames may be output directly
    	if (!p->used_for_reference) {
      		get_smallest_poc(p_Dpb, &poc, &pos);
      		if (-1 == pos || p->poc < poc) {
#if (MVC_EXTENSION_ENABLE)
        		direct_output(p_Vid, p, p_Vid->p_out_mvc[p_Dpb->layer_id]);
#endif
        		return;
      		}
    	}
    	// flush a frame
    	output_one_frame_from_dpb(p_Dpb);
  	}

  	// check for duplicate frame number in short term reference buffer
  	if (p->used_for_reference && !p->is_long_term) {
    	for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++) {
      		if (p_Dpb->fs_ref[i]->FrameNum == p->frame_num)
        		error("duplicate frame_num in short-term reference picture buffer", 500);
    	}
  	}
  	// store at end of buffer
  	insert_picture_in_dpb(p_Vid, p_Dpb->fs[p_Dpb->used_size],p);

  	// picture error concealment
  	if (p->idr_flag)
    	p_Vid->earlier_missing_poc = 0;

  	if (p->structure != FRAME)
    	p_Dpb->last_picture = p_Dpb->fs[p_Dpb->used_size];
  	else
    	p_Dpb->last_picture = NULL;

  	p_Dpb->used_size++;

  	if (p_Vid->conceal_mode != 0)
    	p_Vid->pocs_in_dpb[p_Dpb->used_size-1] = p->poc;

  	update_ref_list(p_Dpb);
  	update_ltref_list(p_Dpb);

  	check_num_ref(p_Dpb);
}
