
#include "picture.h"
#include "memalloc.h"


static inline int RSD(int x)
{
    return ((x&2)?(x|1):(x&(~1)));
}


storable_picture::storable_picture(VideoParameters *p_Vid, PictureStructure structure,
                                   int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output)
{
    sps_t* sps = p_Vid->active_sps;  

    if (structure != FRAME) {
        size_y    /= 2;
        size_y_cr /= 2;
    }

    this->iChromaPadX = sps->chroma_format_idc == YUV444 ? MCBUF_LUMA_PAD_X : MCBUF_CHROMA_PAD_X;
    this->iChromaPadY = sps->chroma_format_idc == YUV444 ? MCBUF_LUMA_PAD_Y :
                        sps->chroma_format_idc == YUV422 ? MCBUF_CHROMA_PAD_Y * 2 : MCBUF_CHROMA_PAD_Y;
    this->iLumaStride   = size_x    + 2 * MCBUF_LUMA_PAD_X;
    this->iChromaStride = size_x_cr + 2 * this->iChromaPadX;

    get_mem2Dpel_pad(&this->imgY, size_y, size_x, MCBUF_LUMA_PAD_Y, MCBUF_LUMA_PAD_X);
    if (sps->chroma_format_idc != YUV400) {
        get_mem2Dpel_pad(&this->imgUV[0], size_y_cr, size_x_cr, this->iChromaPadY, this->iChromaPadX);
        get_mem2Dpel_pad(&this->imgUV[1], size_y_cr, size_x_cr, this->iChromaPadY, this->iChromaPadX);
    }

    get_mem2Dmp(&this->mv_info, (size_y / 4), (size_x / 4));
    this->motion.mb_field_decoding_flag = new bool[(size_x / 4) * (size_y / 4)];

    if (sps->separate_colour_plane_flag) {
        for (int nplane = 0; nplane < 3; ++nplane) {
            get_mem2Dmp(&this->JVmv_info[nplane], (size_y / 4), (size_x / 4));
            this->JVmotion[nplane].mb_field_decoding_flag = new bool[(size_x / 4) * (size_y / 4)];
        }
    }

    this->sps = sps;

    this->top_poc = this->bottom_poc = this->poc = 0;
    this->PicNum             = 0;
    this->frame_num          = 0;
    this->LongTermFrameIdx   = 0;
    this->LongTermPicNum     = 0;
    this->used_for_reference = 0;
    this->is_long_term       = 0;
    this->non_existing       = 0;
    this->is_output          = 0;
    this->no_ref             = 0;
    this->recovery_frame     = 0;
    this->slice.layer_id        = 0;
#if (MVC_EXTENSION_ENABLE)
    this->slice.view_id         = -1;
    this->slice.inter_view_flag = 0;
    this->slice.anchor_pic_flag = 0;
#endif

    this->size_x       = size_x;
    this->size_y       = size_y;
    this->size_x_cr    = size_x_cr;
    this->size_y_cr    = size_y_cr;

    this->top_field    = p_Vid->no_reference_picture;
    this->bottom_field = p_Vid->no_reference_picture;
    this->frame        = p_Vid->no_reference_picture;

    this->slice.structure   = structure;
    this->slice.iCodingType = 0;
    this->slice.idr_flag    = 0;

    this->seiHasTone_mapping = 0;
}

storable_picture::~storable_picture()
{
    if (this->mv_info) {
        free_mem2Dmp(this->mv_info);
        this->mv_info = nullptr;
    }
    delete []this->motion.mb_field_decoding_flag;

    if (this->sps->separate_colour_plane_flag) {
        for (int nplane = 0; nplane < 3; nplane++) {
            if (this->JVmv_info[nplane]) {
                free_mem2Dmp(this->JVmv_info[nplane]);
                this->JVmv_info[nplane] = nullptr;
            }
            delete []this->JVmotion[nplane].mb_field_decoding_flag;
        }
    }

    if (this->imgY) {
        free_mem2Dpel_pad(this->imgY, MCBUF_LUMA_PAD_Y, MCBUF_LUMA_PAD_X);
        this->imgY = nullptr;
    }

    if (this->imgUV[0]) {
        free_mem2Dpel_pad(this->imgUV[0], this->iChromaPadY, this->iChromaPadX);
        free_mem2Dpel_pad(this->imgUV[1], this->iChromaPadY, this->iChromaPadX);
        this->imgUV[0] = nullptr;
        this->imgUV[1] = nullptr;
    }

    if (this->seiHasTone_mapping)
        delete this->tone_mapping_lut;
}

bool storable_picture::is_short_ref()
{
    return this->used_for_reference && !this->is_long_term;
}

bool storable_picture::is_long_ref()
{
    return this->used_for_reference && this->is_long_term;
}

void storable_picture::clear()
{
    for (int i = 0; i < this->size_y; i++) {
        for (int j = 0; j < this->size_x; j++)
            this->imgY[i][j] = (imgpel) (1 << (this->sps->BitDepthY - 1));
    }
    for (int i = 0; i < this->size_y_cr; i++) {
        for (int j = 0; j < this->size_x_cr; j++)
            this->imgUV[0][i][j] = (imgpel) (1 << (this->sps->BitDepthC - 1));
    }
    for (int i = 0; i < this->size_y_cr; i++) {
        for (int j = 0; j < this->size_x_cr; j++)
            this->imgUV[1][i][j] = (imgpel) (1 << (this->sps->BitDepthC - 1));
    }
}


picture_t::~picture_t()
{
    if (this->frame)
        delete this->frame;
    if (this->top_field)
        delete this->top_field;
    if (this->bottom_field)
        delete this->bottom_field;
}

void picture_t::unmark_for_reference()
{
    if ((this->is_used & 1) && this->top_field)
        this->top_field->used_for_reference = 0;
    if ((this->is_used & 2) && this->bottom_field)
        this->bottom_field->used_for_reference = 0;
    if (this->is_used == 3) {
        if (this->top_field && this->bottom_field) {
            this->top_field->used_for_reference = 0;
            this->bottom_field->used_for_reference = 0;
        }
        this->frame->used_for_reference = 0;
    }
    this->is_reference = 0;

    if (this->frame) {
        delete []this->frame->motion.mb_field_decoding_flag;
        this->frame->motion.mb_field_decoding_flag = nullptr;
    }

    if (this->top_field) {
        delete []this->top_field->motion.mb_field_decoding_flag;
        this->top_field->motion.mb_field_decoding_flag = nullptr;
    }

    if (this->bottom_field) {
        delete []this->bottom_field->motion.mb_field_decoding_flag;
        this->bottom_field->motion.mb_field_decoding_flag = nullptr;
    }
}

void picture_t::unmark_for_long_term_reference()
{
    if ((this->is_used & 1) && this->top_field) {
        this->top_field->used_for_reference = 0;
        this->top_field->is_long_term = 0;
    }
    if ((this->is_used & 2) && this->bottom_field) {
        this->bottom_field->used_for_reference = 0;
        this->bottom_field->is_long_term = 0;
    }
    if (this->is_used == 3) {
        if (this->top_field && this->bottom_field) {
            this->top_field->used_for_reference = 0;
            this->top_field->is_long_term = 0;
            this->bottom_field->used_for_reference = 0;
            this->bottom_field->is_long_term = 0;
        }
        this->frame->used_for_reference = 0;
        this->frame->is_long_term = 0;
    }

    this->is_reference = 0;
    this->is_long_term = 0;
}

bool picture_t::is_short_term_reference()
{
    if (this->is_used == 3) {
        storable_picture* pic = this->frame;
        if (pic->used_for_reference && !pic->is_long_term)
            return true;
    }

    if ((this->is_used & 1) && this->top_field) {
        storable_picture* pic = this->top_field;
        if (pic->used_for_reference && !pic->is_long_term)
            return true;
    }

    if ((this->is_used & 2) && this->bottom_field) {
        storable_picture* pic = this->bottom_field;
        if (pic->used_for_reference && !pic->is_long_term)
            return true;
    }

    return false;
}

bool picture_t::is_long_term_reference()
{
    if (this->is_used == 3) {
        storable_picture* pic = this->frame;
        if (pic->used_for_reference && pic->is_long_term)
            return true;
    }

    if ((this->is_used & 1) && this->top_field) {
        storable_picture* pic = this->top_field;
        if (pic->used_for_reference && pic->is_long_term)
            return true;
    }

    if ((this->is_used & 2) && this->bottom_field) {
        storable_picture* pic = this->bottom_field;
        if (pic->used_for_reference && pic->is_long_term)
            return true;
    }

    return false;
}

bool picture_t::is_used_for_reference()
{
    if (this->is_reference)
        return true;

    if (this->is_used == 3) {
        if (this->frame->used_for_reference)
            return true;
    }

    if ((this->is_used & 1) && this->top_field) {
        if (this->top_field->used_for_reference)
            return true;
    }

    if ((this->is_used & 2) && this->bottom_field) {
        if (this->bottom_field->used_for_reference)
            return true;
    }

    return false;
}

void picture_t::dpb_split_field(VideoParameters* p_Vid)
{
    int currentmb;
    int twosz16 = 2 * (this->frame->size_x >> 4);
    storable_picture* fs_top = NULL, *fs_btm = NULL; 
    storable_picture* frame = this->frame;

    this->poc = frame->poc;

    if (!p_Vid->active_sps->frame_mbs_only_flag) {
        fs_top = this->top_field    = new storable_picture(p_Vid, TOP_FIELD,
            frame->size_x, frame->size_y, frame->size_x_cr, frame->size_y_cr, 1);
        fs_btm = this->bottom_field = new storable_picture(p_Vid, BOTTOM_FIELD,
            frame->size_x, frame->size_y, frame->size_x_cr, frame->size_y_cr, 1);

        for (int i = 0; i < (frame->size_y >> 1); i++)
            memcpy(fs_top->imgY[i], frame->imgY[i * 2], frame->size_x * sizeof(imgpel));

        for (int i = 0; i < (frame->size_y_cr >> 1); i++) {
            memcpy(fs_top->imgUV[0][i], frame->imgUV[0][i * 2], frame->size_x_cr * sizeof(imgpel));
            memcpy(fs_top->imgUV[1][i], frame->imgUV[1][i * 2], frame->size_x_cr * sizeof(imgpel));
        }

        for (int i = 0; i < (frame->size_y >> 1); i++)
            memcpy(fs_btm->imgY[i], frame->imgY[i * 2 + 1], frame->size_x * sizeof(imgpel));

        for (int i = 0; i < (frame->size_y_cr >> 1); i++) {
            memcpy(fs_btm->imgUV[0][i], frame->imgUV[0][i * 2 + 1], frame->size_x_cr * sizeof(imgpel));
            memcpy(fs_btm->imgUV[1][i], frame->imgUV[1][i * 2 + 1], frame->size_x_cr * sizeof(imgpel));
        }

        fs_top->sps = frame->sps;
        fs_top->pps = frame->pps;
        fs_top->slice_headers = frame->slice_headers;
        fs_btm->sps = frame->sps;
        fs_btm->pps = frame->pps;
        fs_btm->slice_headers = frame->slice_headers;

        fs_top->poc = frame->top_poc;
        fs_btm->poc = frame->bottom_poc;

#if (MVC_EXTENSION_ENABLE)
        fs_top->slice.view_id = frame->slice.view_id;
        fs_btm->slice.view_id = frame->slice.view_id;
#endif

        fs_top->frame_poc  = fs_btm->frame_poc  = frame->frame_poc;
        fs_top->top_poc    = fs_btm->top_poc    = frame->top_poc;
        fs_top->bottom_poc = fs_btm->bottom_poc = frame->bottom_poc;

        fs_top->used_for_reference = fs_btm->used_for_reference = frame->used_for_reference;
        fs_top->is_long_term = fs_btm->is_long_term = frame->is_long_term;
        this->LongTermFrameIdx = fs_top->LongTermFrameIdx = fs_btm->LongTermFrameIdx = frame->LongTermFrameIdx;

        fs_top->slice.iCodingType = fs_btm->slice.iCodingType = frame->slice.iCodingType;

        frame->top_field     = fs_top;
        frame->bottom_field  = fs_btm;
        frame->frame         = frame;
        fs_top->bottom_field = fs_btm;
        fs_top->frame        = frame;
        fs_top->top_field    = fs_top;
        fs_btm->top_field    = fs_top;
        fs_btm->frame        = frame;
        fs_btm->bottom_field = fs_btm;

#if (MVC_EXTENSION_ENABLE)
        fs_top->slice.view_id = fs_btm->slice.view_id = this->view_id;
        fs_top->slice.inter_view_flag = this->inter_view_flag[0];
        fs_btm->slice.inter_view_flag = this->inter_view_flag[1];
#endif

        if (frame->used_for_reference) {
            pad_dec_picture(p_Vid, fs_top);
            pad_dec_picture(p_Vid, fs_btm);
        }
    } else {
        this->top_field     = NULL;
        this->bottom_field  = NULL;
        frame->top_field    = NULL;
        frame->bottom_field = NULL;
        frame->frame        = frame;
    }

    if (!p_Vid->active_sps->frame_mbs_only_flag) {
        if (frame->slice_headers[0]->header.MbaffFrameFlag) {
            pic_motion_params_old* frm_motion = &frame->motion;
            for (int j = 0 ; j < (frame->size_y >> 3); j++) {
                int jj = (j >> 2)*8 + (j & 0x03);
                int jj4 = jj + 4;
                int jdiv = (j >> 1);
                for (int i = 0; i < (frame->size_x >> 2); i++) {
                    int idiv = (i >> 2);

                    currentmb = twosz16 * (jdiv >> 1) + idiv * 2 + (jdiv & 0x01);
                    // Assign field mvs attached to MB-Frame buffer to the proper buffer
                    if (frm_motion->mb_field_decoding_flag[currentmb]) {
                        auto RefPicList0 = p_Vid->ppSliceList[frame->mv_info[jj4][i].slice_no]->RefPicList[LIST_0];
                        auto RefPicList1 = p_Vid->ppSliceList[frame->mv_info[jj4][i].slice_no]->RefPicList[LIST_1];

                        fs_btm->mv_info[j][i].mv[LIST_0] = frame->mv_info[jj4][i].mv[LIST_0];
                        fs_btm->mv_info[j][i].mv[LIST_1] = frame->mv_info[jj4][i].mv[LIST_1];
                        fs_btm->mv_info[j][i].ref_idx[LIST_0] = frame->mv_info[jj4][i].ref_idx[LIST_0];
                        fs_btm->mv_info[j][i].ref_pic[LIST_0] = NULL;
                        if (fs_btm->mv_info[j][i].ref_idx[LIST_0] >= 0) {
                            int ref_idx = fs_btm->mv_info[j][i].ref_idx[LIST_0];
                            fs_btm->mv_info[j][i].ref_pic[LIST_0] = ref_idx % 2 ? RefPicList0[ref_idx / 2]->top_field : RefPicList0[ref_idx / 2]->bottom_field;
                        }
                        fs_btm->mv_info[j][i].ref_idx[LIST_1] = frame->mv_info[jj4][i].ref_idx[LIST_1];
                        fs_btm->mv_info[j][i].ref_pic[LIST_1] = NULL;
                        if (fs_btm->mv_info[j][i].ref_idx[LIST_1] >= 0) {
                            int ref_idx = fs_btm->mv_info[j][i].ref_idx[LIST_1];
                            fs_btm->mv_info[j][i].ref_pic[LIST_1] = ref_idx % 2 ? RefPicList1[ref_idx / 2]->top_field : RefPicList1[ref_idx / 2]->bottom_field;
                        }
          
                        fs_top->mv_info[j][i].mv[LIST_0] = frame->mv_info[jj][i].mv[LIST_0];
                        fs_top->mv_info[j][i].mv[LIST_1] = frame->mv_info[jj][i].mv[LIST_1];
                        fs_top->mv_info[j][i].ref_idx[LIST_0] = frame->mv_info[jj][i].ref_idx[LIST_0];
                        fs_top->mv_info[j][i].ref_pic[LIST_0] = NULL;
                        if (fs_top->mv_info[j][i].ref_idx[LIST_0] >= 0) {
                            int ref_idx = fs_top->mv_info[j][i].ref_idx[LIST_0];
                            fs_top->mv_info[j][i].ref_pic[LIST_0] = ref_idx % 2 ? RefPicList0[ref_idx / 2]->bottom_field : RefPicList0[ref_idx / 2]->top_field;
                        }
                        fs_top->mv_info[j][i].ref_idx[LIST_1] = frame->mv_info[jj][i].ref_idx[LIST_1];
                        fs_top->mv_info[j][i].ref_pic[LIST_1] = NULL;
                        if (fs_top->mv_info[j][i].ref_idx[LIST_1] >= 0) {
                            int ref_idx = fs_top->mv_info[j][i].ref_idx[LIST_1];
                            fs_top->mv_info[j][i].ref_pic[LIST_1] = ref_idx % 2 ? RefPicList1[ref_idx / 2]->bottom_field : RefPicList1[ref_idx / 2]->top_field;
                        }
                    }
                }
            }
        }

        //! Generate field MVs from Frame MVs
        for (int j = 0; j < (frame->size_y >> 3); j++) {
            int jj = 2* RSD(j);
            int jdiv = (j >> 1);
            for (int i = 0; i < (frame->size_x >> 2); i++) {
                int ii = RSD(i);
                int idiv = (i >> 2);

                currentmb = twosz16 * (jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);

                if (!frame->slice_headers[0]->header.MbaffFrameFlag || !frame->motion.mb_field_decoding_flag[currentmb]) {
                    fs_top->mv_info[j][i].mv[LIST_0] = fs_btm->mv_info[j][i].mv[LIST_0] = frame->mv_info[jj][ii].mv[LIST_0];
                    fs_top->mv_info[j][i].mv[LIST_1] = fs_btm->mv_info[j][i].mv[LIST_1] = frame->mv_info[jj][ii].mv[LIST_1];

                    // Scaling of references is done here since it will not affect spatial direct (2*0 =0)
                    if (frame->mv_info[jj][ii].ref_idx[LIST_0] == -1) {
                        fs_top->mv_info[j][i].ref_idx[LIST_0] = fs_btm->mv_info[j][i].ref_idx[LIST_0] = - 1;
                        fs_top->mv_info[j][i].ref_pic[LIST_0] = fs_btm->mv_info[j][i].ref_pic[LIST_0] = NULL;
                    } else {
                        fs_top->mv_info[j][i].ref_idx[LIST_0] = fs_btm->mv_info[j][i].ref_idx[LIST_0] = frame->mv_info[jj][ii].ref_idx[LIST_0];
                        fs_top->mv_info[j][i].ref_pic[LIST_0] = fs_btm->mv_info[j][i].ref_pic[LIST_0] = p_Vid->ppSliceList[frame->mv_info[jj][ii].slice_no]->RefPicList[LIST_0][(short) frame->mv_info[jj][ii].ref_idx[LIST_0]];
                    }

                    if (frame->mv_info[jj][ii].ref_idx[LIST_1] == -1) {
                        fs_top->mv_info[j][i].ref_idx[LIST_1] = fs_btm->mv_info[j][i].ref_idx[LIST_1] = - 1;
                        fs_top->mv_info[j][i].ref_pic[LIST_1] = fs_btm->mv_info[j][i].ref_pic[LIST_1] = NULL;
                    } else {
                        fs_top->mv_info[j][i].ref_idx[LIST_1] = fs_btm->mv_info[j][i].ref_idx[LIST_1] = frame->mv_info[jj][ii].ref_idx[LIST_1];
                        fs_top->mv_info[j][i].ref_pic[LIST_1] = fs_btm->mv_info[j][i].ref_pic[LIST_1] = p_Vid->ppSliceList[frame->mv_info[jj][ii].slice_no]->RefPicList[LIST_1][(short) frame->mv_info[jj][ii].ref_idx[LIST_1]];
                    }
                }
            }
        }
    }
}

void picture_t::dpb_combine_field_yuv(VideoParameters* p_Vid)
{
    if (!this->frame)
        this->frame = new storable_picture(p_Vid, FRAME,
            this->top_field->size_x, this->top_field->size_y * 2,
            this->top_field->size_x_cr, this->top_field->size_y_cr * 2, 1);

    for (int i = 0; i < this->top_field->size_y; i++) {
        memcpy(this->frame->imgY[i*2],     this->top_field->imgY[i]   , this->top_field->size_x * sizeof(imgpel));     // top field
        memcpy(this->frame->imgY[i*2 + 1], this->bottom_field->imgY[i], this->bottom_field->size_x * sizeof(imgpel)); // bottom field
    }

    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < this->top_field->size_y_cr; i++) {
            memcpy(this->frame->imgUV[j][i*2],     this->top_field->imgUV[j][i],    this->top_field->size_x_cr*sizeof(imgpel));
            memcpy(this->frame->imgUV[j][i*2 + 1], this->bottom_field->imgUV[j][i], this->bottom_field->size_x_cr*sizeof(imgpel));
        }
    }

    this->poc = this->frame->poc = this->frame->frame_poc = min(this->top_field->poc, this->bottom_field->poc);

    this->bottom_field->frame_poc = this->top_field->frame_poc = this->frame->poc;

    this->bottom_field->top_poc = this->frame->top_poc = this->top_field->poc;
    this->top_field->bottom_poc = this->frame->bottom_poc = this->bottom_field->poc;

    this->frame->used_for_reference = (this->top_field->used_for_reference && this->bottom_field->used_for_reference );
    this->frame->is_long_term = (this->top_field->is_long_term && this->bottom_field->is_long_term );

    if (this->frame->is_long_term)
        this->frame->LongTermFrameIdx = this->LongTermFrameIdx;

    this->frame->top_field    = this->top_field;
    this->frame->bottom_field = this->bottom_field;
    this->frame->frame        = this->frame;

    this->top_field->frame           = this->bottom_field->frame = this->frame;
    this->top_field->top_field       = this->top_field;
    this->top_field->bottom_field    = this->bottom_field;
    this->bottom_field->top_field    = this->top_field;
    this->bottom_field->bottom_field = this->bottom_field;
    if (this->top_field->used_for_reference || this->bottom_field->used_for_reference)
        pad_dec_picture(p_Vid, this->frame);
}

void picture_t::dpb_combine_field(VideoParameters* p_Vid)
{
    this->dpb_combine_field_yuv(p_Vid);

#if (MVC_EXTENSION_ENABLE)
    this->frame->slice.view_id = this->view_id;
#endif
    this->frame->slice.iCodingType = this->top_field->slice.iCodingType; //FIELD_CODING;
    //! Use inference flag to remap mvs/references

    //! Generate Frame parameters from field information.
    for (int j = 0; j < (this->top_field->size_y >> 2); j++) {
        int jj = (j << 1);
        int jj4 = jj + 1;
        for (int i = 0; i < (this->top_field->size_x >> 2); i++) {
            auto& mv_info_t = this->top_field->mv_info[j][i];
            this->frame->mv_info[jj][i].mv     [LIST_0] = mv_info_t.mv     [LIST_0];
            this->frame->mv_info[jj][i].mv     [LIST_1] = mv_info_t.mv     [LIST_1];
            this->frame->mv_info[jj][i].ref_idx[LIST_0] = mv_info_t.ref_idx[LIST_0];
            this->frame->mv_info[jj][i].ref_idx[LIST_1] = mv_info_t.ref_idx[LIST_1];

            /* bug: top field list doesnot exist.*/
            this->frame->mv_info[jj][i].ref_pic[LIST_0] = mv_info_t.ref_idx[LIST_0] >= 0 ?
                p_Vid->ppSliceList[mv_info_t.slice_no]->RefPicList[LIST_0][(int)mv_info_t.ref_idx[LIST_0]] : NULL;
            this->frame->mv_info[jj][i].ref_pic[LIST_1] = mv_info_t.ref_idx[LIST_1] >= 0 ?
                p_Vid->ppSliceList[mv_info_t.slice_no]->RefPicList[LIST_1][(int)mv_info_t.ref_idx[LIST_1]] : NULL;

            //! association with id already known for fields.
            auto& mv_info_b = this->bottom_field->mv_info[j][i];
            this->frame->mv_info[jj4][i].mv     [LIST_0] = mv_info_b.mv     [LIST_0];
            this->frame->mv_info[jj4][i].mv     [LIST_1] = mv_info_b.mv     [LIST_1];
            this->frame->mv_info[jj4][i].ref_idx[LIST_0] = mv_info_b.ref_idx[LIST_0];
            this->frame->mv_info[jj4][i].ref_idx[LIST_1] = mv_info_b.ref_idx[LIST_1];

            this->frame->mv_info[jj4][i].ref_pic[LIST_0] = mv_info_b.ref_idx[LIST_0] >= 0 ?
                p_Vid->ppSliceList[mv_info_b.slice_no]->RefPicList[LIST_0][(int)mv_info_b.ref_idx[LIST_0]] : NULL;
            this->frame->mv_info[jj4][i].ref_pic[LIST_1] = mv_info_b.ref_idx[LIST_1] >= 0 ?
                p_Vid->ppSliceList[mv_info_b.slice_no]->RefPicList[LIST_1][(int)mv_info_b.ref_idx[LIST_1]] : NULL;
        }
    }
}

void picture_t::insert_picture(VideoParameters* p_Vid, storable_picture* p)
{
    assert(p);

    switch (p->slice.structure) {
    case FRAME:
        this->frame = p;
        this->is_used = 3;
        if (p->used_for_reference) {
            this->is_reference = 3;
            this->is_orig_reference = 3;
            if (p->is_long_term) {
                this->is_long_term = 3;
                this->LongTermFrameIdx = p->LongTermFrameIdx;
            }
        }
        this->layer_id = p->slice.layer_id;
#if (MVC_EXTENSION_ENABLE)
        this->view_id = p->slice.view_id;
        this->inter_view_flag[0] = this->inter_view_flag[1] = p->slice.inter_view_flag;
        this->anchor_pic_flag[0] = this->anchor_pic_flag[1] = p->slice.anchor_pic_flag;
#endif
        // generate field views
        this->dpb_split_field(p_Vid);
        break;
    case TOP_FIELD:
        this->top_field = p;
        this->is_used |= 1;
        this->layer_id = p->slice.layer_id;
#if (MVC_EXTENSION_ENABLE)
        this->view_id = p->slice.view_id;
        this->inter_view_flag[0] = p->slice.inter_view_flag;
        this->anchor_pic_flag[0] = p->slice.anchor_pic_flag;
#endif
        if (p->used_for_reference) {
            this->is_reference |= 1;
            this->is_orig_reference |= 1;
            if (p->is_long_term) {
                this->is_long_term |= 1;
                this->LongTermFrameIdx = p->LongTermFrameIdx;
            }
        }
        if (this->is_used == 3) {
            // generate frame view
            this->dpb_combine_field(p_Vid);
        } else
            this->poc = p->poc;
        break;
    case BOTTOM_FIELD:
        this->bottom_field = p;
        this->is_used |= 2;
        this->layer_id = p->slice.layer_id;
#if (MVC_EXTENSION_ENABLE)
        this->view_id = p->slice.view_id;
        this->inter_view_flag[1] = p->slice.inter_view_flag;
        this->anchor_pic_flag[1] = p->slice.anchor_pic_flag;
#endif
        if (p->used_for_reference) {
            this->is_reference |= 2;
            this->is_orig_reference |= 2;
            if (p->is_long_term) {
                this->is_long_term |= 2;
                this->LongTermFrameIdx = p->LongTermFrameIdx;
            }
        }
        if (this->is_used == 3) {
            // generate frame view
            this->dpb_combine_field(p_Vid);
        } else
            this->poc = p->poc;
        break;
    }

    this->FrameNum = p->PicNum;
    this->recovery_frame = p->recovery_frame;
    this->is_output = p->is_output;

    if (this->is_used == 3)
        p_Vid->calculate_frame_no(p);
}
