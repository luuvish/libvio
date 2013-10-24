#include "global.h"
#include "input_parameters.h"
#include "slice.h"
#include "dpb.h"
#include "memalloc.h"
#include "macroblock.h"
#include "output.h"


/* Thomson APIs for concealing entire frame loss */

struct concealment_node {
    storable_picture* picture;
    int missingpocs;
    concealment_node* next;
};


//! look up tables for FRExt_chroma support
static const uint8_t subblk_offset_x[3][8][4] = {
  { {  0,  4,  0,  4 },
    {  0,  4,  0,  4 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 } },
  { {  0,  4,  0,  4 },
    {  0,  4,  0,  4 },
    {  0,  4,  0,  4 },
    {  0,  4,  0,  4 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 } },
  { {  0,  4,  0,  4 },
    {  8, 12,  8, 12 },
    {  0,  4,  0,  4 },
    {  8, 12,  8, 12 },
    {  0,  4,  0,  4 },
    {  8, 12,  8, 12 },
    {  0,  4,  0,  4 },
    {  8, 12,  8, 12 } }
};

static const uint8_t subblk_offset_y[3][8][4] = {
  { {  0,  0,  4,  4 },
    {  0,  0,  4,  4 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 } },
  { {  0,  0,  4,  4 },
    {  8,  8, 12, 12 },
    {  0,  0,  4,  4 },
    {  8,  8, 12, 12 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 },
    {  0,  0,  0,  0 } },
  { {  0,  0,  4,  4 },
    {  0,  0,  4,  4 },
    {  8,  8, 12, 12 },
    {  8,  8, 12, 12 },
    {  0,  0,  4,  4 },
    {  0,  0,  4,  4 },
    {  8,  8, 12, 12 },
    {  8,  8, 12, 12 } }
};


static void buildPredblockRegionYUV(VideoParameters* p_Vid, int* mv,
                                    int x, int y, px_t* predMB, int list, int current_mb_nr)
{
    storable_picture* dec_picture = p_Vid->dec_picture;

    int ref_frame = mv[2];
    int mb_nr = current_mb_nr;
  
    mb_t& mb = p_Vid->mb_data[mb_nr];   // intialization code deleted, see below, StW  
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;
    int yuv = sps.chroma_format_idc - 1;

    auto ref_pic = slice.RefPicList[list][ref_frame];

    /* Update coordinates of the current concealed macroblock */
    mb.mb.x = (short)(x / 4);
    mb.mb.y = (short)(y / 4);

    int mv_mul = 4;

    // luma *******************************************************
    px_t tmp_block[16][16];
    int vec1_x = x * mv_mul + mv[0];
    int vec1_y = y * mv_mul + mv[1];
    slice.decoder.get_block_luma(ref_pic,
        vec1_x, vec1_y, 4, 4, tmp_block,
        dec_picture->iLumaStride, dec_picture->size_x - 1,
        mb.mb_field_decoding_flag ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y - 1, PLANE_Y, &mb);

    for (int jj = 0; jj < 16/4; ++jj) {
        for (int ii = 0; ii < 4; ++ii)
            slice.mb_pred[PLANE_Y][jj][ii] = tmp_block[jj][ii];
    }

    px_t* pMB = predMB;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i)
            pMB[j * 4 + i] = slice.mb_pred[PLANE_Y][j][i];
    }
    pMB += 16;

    if (sps.chroma_format_idc != YUV400) {
        // chroma *******************************************************
        int f1_x = 64 / sps.MbWidthC;
        int f2_x = f1_x - 1;

        int f1_y = 64 / sps.MbHeightC;
        int f2_y = f1_y - 1;

        int f3 = f1_x * f1_y;
        int f4 = f3 >> 1;

        for (int uv = 0; uv < 2; ++uv) {
            int joff = subblk_offset_y[yuv][0][0];
            int j4   = mb.mb.y * sps.MbHeightC / 4 + joff;
            int ioff = subblk_offset_x[yuv][0][0];
            int i4   = mb.mb.x * sps.MbWidthC / 4 + ioff;

            for (int jj = 0; jj < 2; ++jj) {
                for (int ii = 0; ii < 2; ++ii) {
                    int i1 = (i4 + ii) * f1_x + mv[0];
                    int j1 = (j4 + jj) * f1_y + mv[1];

                    int ii0 = clip3(0, dec_picture->size_x_cr - 1, i1 / f1_x);
                    int jj0 = clip3(0, dec_picture->size_y_cr - 1, j1 / f1_y);
                    int ii1 = clip3(0, dec_picture->size_x_cr - 1, (i1 + f2_x) / f1_x);
                    int jj1 = clip3(0, dec_picture->size_y_cr - 1, (j1 + f2_y) / f1_y);

                    int if1 = (i1 & f2_x);
                    int jf1 = (j1 & f2_y);
                    int if0 = (f1_x - if1);
                    int jf0 = (f1_y - jf1);

                    slice.mb_pred[uv + 1][jj][ii] = (px_t)
                        ((if0 * jf0 * ref_pic->imgUV[uv][jj0][ii0] +
                          if1 * jf0 * ref_pic->imgUV[uv][jj0][ii1] +
                          if0 * jf1 * ref_pic->imgUV[uv][jj1][ii0] +
                          if1 * jf1 * ref_pic->imgUV[uv][jj1][ii1] + f4) / f3);
                }
            }

            for (int j = 0; j < 2; ++j) {
                for (int i = 0; i < 2; ++i)
                    pMB[j * 2 + i] = slice.mb_pred[uv + 1][j][i];
            }
            pMB += 4;
        }
    }
}


static void CopyImgData(px_t** inputY, px_t*** inputUV, px_t** outputY, px_t*** outputUV, 
                        int img_width, int img_height, int img_width_cr, int img_height_cr)
{
    for (int y = 0; y < img_height; ++y) {
        for (int x = 0; x < img_width; ++x)
            outputY[y][x] = inputY[y][x];
    }

    for (int y = 0; y < img_height_cr; ++y) {
        for (int x = 0; x < img_width_cr; ++x) {
            outputUV[0][y][x] = inputUV[0][y][x];
            outputUV[1][y][x] = inputUV[1][y][x];
        }
    }
}

static void copy_to_conceal(storable_picture *src, storable_picture *dst, VideoParameters *p_Vid)
{
    sps_t* sps = p_Vid->active_sps;

    int current_mb_nr = 0;

    dst->sps = src->sps;
    dst->pps = src->pps;
    dst->slice_headers = src->slice_headers;

    dst->slice.slice_type = src->slice.slice_type = p_Vid->conceal_slice_type;
    dst->slice.idr_flag = 0; //since we do not want to clears the ref list

    // Conceals the missing frame by frame copy concealment
    if (p_Vid->conceal_mode == 1) {
        // We need these initializations for using deblocking filter for frame copy
        // concealment as well.

        CopyImgData(src->imgY, src->imgUV, dst->imgY, dst->imgUV,
                    sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
                    sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC);
    }

    // Conceals the missing frame by motion vector copy concealment
    if (p_Vid->conceal_mode == 2) {
        px_t* storeYUV;
        if (sps->chroma_format_idc != YUV400)
            storeYUV = new px_t[16 + sps->MbWidthC * sps->MbHeightC * 2 / 16];
        else
            storeYUV = new px_t[16];

        int mb_width  = sps->PicWidthInMbs;
        int mb_height = sps->FrameHeightInMbs;
        int scale = (p_Vid->conceal_slice_type == B_slice) ? 2 : 1;

        if (p_Vid->conceal_slice_type == B_slice) {
            p_Vid->p_Dpb_layer[0]->init_lists_for_non_reference_loss(
                dst->slice.slice_type, p_Vid->ppSliceList[0]->header.field_pic_flag);
        } else
            p_Vid->ppSliceList[0]->init_lists();

        int multiplier = 4;

        for (int i = 0; i < mb_height * 4; ++i) {
            int mm = i * 4;
            for (int j = 0; j < mb_width * 4; ++j) {
                int nn = j * 4;

                int mv[3];
                mv[0] = src->mv_info[i][j].mv[LIST_0].mv_x / scale;
                mv[1] = src->mv_info[i][j].mv[LIST_0].mv_y / scale;
                mv[2] = src->mv_info[i][j].ref_idx[LIST_0];

                if (mv[2] < 0)
                    mv[2] = 0;

                dst->mv_info[i][j].mv[LIST_0].mv_x = (short) mv[0];
                dst->mv_info[i][j].mv[LIST_0].mv_y = (short) mv[1];
                dst->mv_info[i][j].ref_idx[LIST_0] = (char)  mv[2];

                int x = j * multiplier;
                int y = i * multiplier;

                if (mm % 16 == 0 && nn % 16 == 0)
                    current_mb_nr++;

                buildPredblockRegionYUV(p_Vid, mv, x, y, storeYUV, LIST_0, current_mb_nr);

                px_t* predMB = storeYUV;

                for (int ii = 0; ii < multiplier; ++ii) {
                    for (int jj = 0; jj < multiplier; ++jj)
                        dst->imgY[i * multiplier + ii][j * multiplier + jj] = predMB[ii * multiplier + jj];
                }

                predMB += (multiplier * multiplier);

                if (sps->chroma_format_idc != YUV400) {
                    for (int uv = 0; uv < 2; ++uv) {
                        for (int ii = 0; ii < multiplier / 2; ++ii) {
                            for (int jj = 0; jj < multiplier / 2; ++jj)
                                dst->imgUV[uv][i * multiplier / 2 + ii][j * multiplier / 2 + jj] =
                                    predMB[ii * (multiplier / 2) + jj];
                        }
                        predMB += (multiplier * multiplier / 4);
                    }
                }
            }
        }
        delete []storeYUV;
    }
}


storable_picture* decoded_picture_buffer_t::get_last_ref_pic_from_dpb()
{
    int used_size = this->used_size - 1;

    for (int i = used_size; i >= 0; --i) {
        if (this->fs[i]->is_used == 3) {
            if (this->fs[i]->frame->used_for_reference && !this->fs[i]->frame->is_long_term)
                return this->fs[i]->frame;
        }
    }

    return nullptr;
}

void decoded_picture_buffer_t::copy_prev_pic_to_concealed_pic(storable_picture* picture)
{
    VideoParameters* p_Vid = this->p_Vid;
    /* get the last ref pic in dpb */
    storable_picture* ref_pic = this->get_last_ref_pic_from_dpb();

    assert(ref_pic);

    /* copy all the struc from this to current concealment pic */
    p_Vid->conceal_slice_type = P_slice;
    copy_to_conceal(ref_pic, picture, p_Vid);
}


void decoded_picture_buffer_t::update_ref_list_for_concealment()
{
    VideoParameters* p_Vid = this->p_Vid;
    int j = 0;

    for (int i = 0; i < this->used_size; ++i) {
        if (this->fs[i]->concealment_reference)
            this->fs_ref[j++] = this->fs[i];
    }

    this->ref_frames_in_buffer = p_Vid->active_pps->num_ref_idx_l0_default_active_minus1;
}

void decoded_picture_buffer_t::init_lists_for_non_reference_loss(int currSliceType, bool field_pic_flag)
{
    VideoParameters *p_Vid = this->p_Vid;
    sps_t *active_sps = p_Vid->active_sps;

    unsigned i;
    int j;
    int diff;

    int list0idx = 0;
    int list0idx_1 = 0;

    storable_picture *tmp_s;

    if (!field_pic_flag) {
        for (i = 0; i < this->ref_frames_in_buffer; i++) {
            if (this->fs[i]->concealment_reference == 1) {
                if (this->fs[i]->FrameNum > p_Vid->frame_to_conceal)
                    this->fs_ref[i]->FrameNumWrap = this->fs[i]->FrameNum - active_sps->MaxFrameNum;
                else
                    this->fs_ref[i]->FrameNumWrap = this->fs[i]->FrameNum;
                this->fs_ref[i]->frame->PicNum = this->fs_ref[i]->FrameNumWrap;
            }
        }
    }

    if (currSliceType == P_slice) {
    // Calculate FrameNumWrap and PicNum
        if (!field_pic_flag) {
            for (i = 0; i < this->used_size; i++) {
                if (this->fs[i]->concealment_reference == 1)
                    p_Vid->ppSliceList[0]->RefPicList[0][list0idx++] = this->fs[i]->frame;
            }
            // order list 0 by PicNum
            std::qsort(p_Vid->ppSliceList[0]->RefPicList[0], list0idx, sizeof(storable_picture*), [](const void *arg1, const void *arg2) {
                int pic_num1 = (*(storable_picture**)arg1)->PicNum;
                int pic_num2 = (*(storable_picture**)arg2)->PicNum;
                return (pic_num1 < pic_num2) ? 1 : (pic_num1 > pic_num2) ? -1 : 0;
            });

            p_Vid->ppSliceList[0]->RefPicSize[0] = (char) list0idx;
        }
    }

    if (currSliceType == B_slice) {
        if (!field_pic_flag) {
            for (i = 0; i < this->used_size; i++) {
                if (this->fs[i]->concealment_reference == 1) {
                    if (p_Vid->earlier_missing_poc > this->fs[i]->frame->poc)
                        p_Vid->ppSliceList[0]->RefPicList[0][list0idx++] = this->fs[i]->frame;
                }
            }

            std::qsort(p_Vid->ppSliceList[0]->RefPicList[0], list0idx, sizeof(storable_picture*), [](const void *arg1, const void *arg2) {
                int poc1 = (*(storable_picture**)arg1)->poc;
                int poc2 = (*(storable_picture**)arg2)->poc;
                return (poc1 < poc2) ? 1 : (poc1 > poc2) ? -1 : 0;
            });

            list0idx_1 = list0idx;

            for (i = 0; i < this->used_size; i++) {
                if (this->fs[i]->concealment_reference == 1) {
                    if (p_Vid->earlier_missing_poc < this->fs[i]->frame->poc)
                        p_Vid->ppSliceList[0]->RefPicList[0][list0idx++] = this->fs[i]->frame;
                }
            }

            std::qsort(&p_Vid->ppSliceList[0]->RefPicList[0][list0idx_1], list0idx-list0idx_1,
                       sizeof(storable_picture*), [](const void *arg1, const void *arg2) {
                int poc1 = (*(storable_picture**)arg1)->poc;
                int poc2 = (*(storable_picture**)arg2)->poc;
                return (poc1 < poc2) ? -1 : (poc1 > poc2) ? 1 : 0;
            });

            for (j = 0; j < list0idx_1; j++)
                p_Vid->ppSliceList[0]->RefPicList[1][list0idx-list0idx_1+j]=p_Vid->ppSliceList[0]->RefPicList[0][j];
            for (j = list0idx_1; j < list0idx; j++)
                p_Vid->ppSliceList[0]->RefPicList[1][j-list0idx_1]=p_Vid->ppSliceList[0]->RefPicList[0][j];

            p_Vid->ppSliceList[0]->RefPicSize[0] = p_Vid->ppSliceList[0]->RefPicSize[1] = (char) list0idx;

            std::qsort(&p_Vid->ppSliceList[0]->RefPicList[0][(short) p_Vid->ppSliceList[0]->RefPicSize[0]], list0idx-p_Vid->ppSliceList[0]->RefPicSize[0],
                       sizeof(storable_picture*), [](const void *arg1, const void *arg2) {
                int long_term_pic_num1 = (*(storable_picture**)arg1)->LongTermPicNum;
                int long_term_pic_num2 = (*(storable_picture**)arg2)->LongTermPicNum;
                return (long_term_pic_num1 < long_term_pic_num2) ? -1 : (long_term_pic_num1 > long_term_pic_num2) ? 1 : 0;
            });
            std::qsort(&p_Vid->ppSliceList[0]->RefPicList[1][(short) p_Vid->ppSliceList[0]->RefPicSize[0]], list0idx-p_Vid->ppSliceList[0]->RefPicSize[0],
                       sizeof(storable_picture*), [](const void *arg1, const void *arg2) {
                int long_term_pic_num1 = (*(storable_picture**)arg1)->LongTermPicNum;
                int long_term_pic_num2 = (*(storable_picture**)arg2)->LongTermPicNum;
                return (long_term_pic_num1 < long_term_pic_num2) ? -1 : (long_term_pic_num1 > long_term_pic_num2) ? 1 : 0;
            });
            p_Vid->ppSliceList[0]->RefPicSize[0] = p_Vid->ppSliceList[0]->RefPicSize[1] = (char) list0idx;
        }
    }

    if ((p_Vid->ppSliceList[0]->RefPicSize[0] == p_Vid->ppSliceList[0]->RefPicSize[1]) &&
        (p_Vid->ppSliceList[0]->RefPicSize[0] > 1)) {
        // check if lists are identical, if yes swap first two elements of listX[1]
        diff = 0;
        for (j = 0; j < p_Vid->ppSliceList[0]->RefPicSize[0]; j++) {
            if (p_Vid->ppSliceList[0]->RefPicList[0][j]!=p_Vid->ppSliceList[0]->RefPicList[1][j])
                diff = 1;
        }
        if (!diff) {
            tmp_s = p_Vid->ppSliceList[0]->RefPicList[1][0];
            p_Vid->ppSliceList[0]->RefPicList[1][0]=p_Vid->ppSliceList[0]->RefPicList[1][1];
            p_Vid->ppSliceList[0]->RefPicList[1][1]=tmp_s;
        }
    }

    // set max size
    p_Vid->ppSliceList[0]->RefPicSize[0] = (char) min<int>(p_Vid->ppSliceList[0]->RefPicSize[0], (int)active_sps->max_num_ref_frames);
    p_Vid->ppSliceList[0]->RefPicSize[1] = 0;

    // set the unused list entries to NULL
    for (i = p_Vid->ppSliceList[0]->RefPicSize[0]; i < MAX_LIST_SIZE; i++)
        p_Vid->ppSliceList[0]->RefPicList[0][i] = NULL;
    for (i = p_Vid->ppSliceList[0]->RefPicSize[1]; i < MAX_LIST_SIZE; i++)
        p_Vid->ppSliceList[0]->RefPicList[1][i] = NULL;
}


storable_picture* decoded_picture_buffer_t::get_pic_from_dpb(int missingpoc, unsigned int* pos)
{
    VideoParameters* p_Vid = this->p_Vid;
    int used_size = this->used_size - 1;

    int concealfrom = 0;
    if (p_Vid->conceal_mode == 1)
        concealfrom = missingpoc - p_Vid->p_Inp->poc_gap;
    else if (p_Vid->conceal_mode == 2)
        concealfrom = missingpoc + p_Vid->p_Inp->poc_gap;

    for (int i = used_size; i >= 0; --i) {
        if (this->fs[i]->poc == concealfrom) {
            *pos = i;
            return this->fs[i]->frame;
        }
    }

    return nullptr;
}

static concealment_node* init_node(storable_picture* picture, int missingpoc)
{
    concealment_node* ptr = new concealment_node {};

    if (!ptr)
        return nullptr;
    else {
        ptr->picture = picture;
        ptr->missingpocs = missingpoc;
        ptr->next = nullptr;
        return ptr;
    }
}

static void add_node(VideoParameters* p_Vid, concealment_node* concealment_new)
{
    if (!p_Vid->concealment_head) {
        p_Vid->concealment_end = p_Vid->concealment_head = concealment_new;
        return;
    }
    p_Vid->concealment_end->next = concealment_new;
    p_Vid->concealment_end = concealment_new;
}

static void delete_node(VideoParameters* p_Vid, concealment_node* ptr)
{
    // We only need to delete the first node in the linked list
    if (ptr == p_Vid->concealment_head) {
        p_Vid->concealment_head = p_Vid->concealment_head->next;
        if (p_Vid->concealment_end == ptr)
            p_Vid->concealment_end = p_Vid->concealment_end->next;
        delete ptr;
    }
}



void decoded_picture_buffer_t::conceal_lost_frames(slice_t *pSlice)
{
    VideoParameters *p_Vid = this->p_Vid;
    sps_t *sps = p_Vid->active_sps;
    int tmp1 = pSlice->header.delta_pic_order_cnt[0];
    int tmp2 = pSlice->header.delta_pic_order_cnt[1];

    pSlice->header.delta_pic_order_cnt[0] = pSlice->header.delta_pic_order_cnt[1] = 0;

    // printf("A gap in frame number is found, try to fill it.\n");

    int UnusedShortTermFrameNum;
    if (p_Vid->IDR_concealment_flag == 1) {
        // Conceals an IDR frame loss. Uses the reference frame in the previous
        // GOP for concealment.
        UnusedShortTermFrameNum = 0;
        p_Vid->last_ref_pic_poc = -p_Vid->p_Inp->poc_gap;
        p_Vid->earlier_missing_poc = 0;
    } else
        UnusedShortTermFrameNum = (p_Vid->PrevRefFrameNum + 1) % p_Vid->active_sps->MaxFrameNum;

    int CurrFrameNum = pSlice->header.frame_num;

    while (CurrFrameNum != UnusedShortTermFrameNum) {
        storable_picture* picture = new storable_picture(p_Vid, FRAME,
            sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
            sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);

        picture->PicNum             = UnusedShortTermFrameNum;
        picture->frame_num          = UnusedShortTermFrameNum;
        picture->non_existing       = 0;
        picture->is_output          = 0;
        picture->used_for_reference = 1;
        picture->concealed_pic      = 1;

        pSlice->header.frame_num = UnusedShortTermFrameNum;

        picture->top_poc        = p_Vid->last_ref_pic_poc + p_Vid->p_Inp->ref_poc_gap;
        picture->bottom_poc     = picture->top_poc;
        picture->frame_poc      = picture->top_poc;
        picture->poc            = picture->top_poc;
        p_Vid->last_ref_pic_poc = picture->poc;

        this->copy_prev_pic_to_concealed_pic(picture);

        if (p_Vid->IDR_concealment_flag == 1) {
            picture->slice.slice_type = I_slice;
            picture->slice.idr_flag = 1;
            this->flush();
            picture->top_poc        = 0;
            picture->bottom_poc     = picture->top_poc;
            picture->frame_poc      = picture->top_poc;
            picture->poc            = picture->top_poc;
            p_Vid->last_ref_pic_poc = picture->poc;
        }

        p_Vid->p_Dpb_layer[0]->store_picture(picture);

        p_Vid->PrevRefFrameNum = UnusedShortTermFrameNum;
        UnusedShortTermFrameNum = (UnusedShortTermFrameNum + 1) % p_Vid->active_sps->MaxFrameNum;

        // update reference flags and set current flag.
        for (int i = 16; i > 0; --i)
            pSlice->ref_flag[i] = pSlice->ref_flag[i - 1];
        pSlice->ref_flag[0] = 0;
    }
    pSlice->header.delta_pic_order_cnt[0] = tmp1;
    pSlice->header.delta_pic_order_cnt[1] = tmp2;
    pSlice->header.frame_num = CurrFrameNum;
}

void decoded_picture_buffer_t::conceal_non_ref_pics(int diff)
{
    VideoParameters* p_Vid = this->p_Vid;
    sps_t* sps = p_Vid->active_sps;
    int temp_used_size = this->used_size;

    if (this->used_size == 0)
        return;

    std::qsort(p_Vid->pocs_in_dpb, this->size, sizeof(int), [](const void* i, const void* j) {
        return *(int*)i - *(int*)j;
    });

    for (int i = 0; i < this->size - diff; ++i) {
        this->used_size = this->size;

        if (p_Vid->pocs_in_dpb[i + 1] - p_Vid->pocs_in_dpb[i] > p_Vid->p_Inp->poc_gap) {
            storable_picture* conceal_to_picture = new storable_picture(p_Vid, FRAME,
                sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
                sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);

            int missingpoc = p_Vid->pocs_in_dpb[i] + p_Vid->p_Inp->poc_gap;

            if (missingpoc > p_Vid->earlier_missing_poc) {
                p_Vid->earlier_missing_poc     = missingpoc;
                conceal_to_picture->top_poc    = missingpoc;
                conceal_to_picture->bottom_poc = missingpoc;
                conceal_to_picture->frame_poc  = missingpoc;
                conceal_to_picture->poc        = missingpoc;

                unsigned pos = 0;
                storable_picture* conceal_from_picture = this->get_pic_from_dpb(missingpoc, &pos);

                assert(conceal_from_picture);

                this->used_size = pos + 1;

                p_Vid->frame_to_conceal = conceal_from_picture->frame_num + 1;

                this->update_ref_list_for_concealment();
                p_Vid->conceal_slice_type = B_slice;
                copy_to_conceal(conceal_from_picture, conceal_to_picture, p_Vid);
                concealment_node* concealment_ptr = init_node(conceal_to_picture, missingpoc);
                add_node(p_Vid, concealment_ptr);
            }
        }
    }

    this->used_size = temp_used_size;
}


void decoded_picture_buffer_t::sliding_window_poc_management(storable_picture* p)
{
    if (this->used_size == this->size) {
        VideoParameters* p_Vid = this->p_Vid;
        for (int i = 0; i < this->size - 1; ++i)
            p_Vid->pocs_in_dpb[i] = p_Vid->pocs_in_dpb[i + 1];
    }
}

void decoded_picture_buffer_t::write_lost_non_ref_pic(int poc, int p_out)
{
    VideoParameters* p_Vid = this->p_Vid;
    pic_t concealment_fs;
    if (poc > 0) {
        if (poc - this->last_output_poc > p_Vid->p_Inp->poc_gap) {
            concealment_fs.frame = p_Vid->concealment_head->picture;
            concealment_fs.is_output = 0;
            concealment_fs.is_reference = 0;
            concealment_fs.is_used = 3;

            write_stored_frame(p_Vid, &concealment_fs, p_out);
            delete_node(p_Vid, p_Vid->concealment_head);
        }
    }
}

void decoded_picture_buffer_t::write_lost_ref_after_idr(int pos)
{
    VideoParameters *p_Vid = this->p_Vid;
    sps_t *sps = p_Vid->active_sps;

    int temp = 1;

    if (!p_Vid->last_out_fs->frame) {
        p_Vid->last_out_fs->frame = new storable_picture(p_Vid, FRAME,
            sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
            sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 1);
        p_Vid->last_out_fs->is_used = 3;
    }

    if (p_Vid->conceal_mode == 2) {
        temp = 2;
        p_Vid->conceal_mode = 1;
    }
    copy_to_conceal(this->fs[pos]->frame, p_Vid->last_out_fs->frame, p_Vid);

    p_Vid->conceal_mode = temp;
}
