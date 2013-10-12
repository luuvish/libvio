#include "global.h"
#include "memalloc.h"
#include "slice.h"


static int init_top_bot_planes(imgpel **imgFrame, int dim0, imgpel ***imgTopField, imgpel ***imgBotField)
{
    *imgTopField = new imgpel*[dim0 >> 1];
    *imgBotField = new imgpel*[dim0 >> 1];

    for (int i = 0; i < (dim0 >> 1); i++) {
        (*imgTopField)[i] = imgFrame[2 * i    ];
        (*imgBotField)[i] = imgFrame[2 * i + 1];
    }

    return dim0 * sizeof(imgpel*);
}

static void free_top_bot_planes(imgpel **imgTopField, imgpel **imgBotField)
{
    delete []imgTopField;
    delete []imgBotField;
}


static void init_img_data(VideoParameters *p_Vid, ImageData *p_ImgData, sps_t *sps)
{
    // allocate memory for reference frame buffers: p_ImgData->frm_data
    p_ImgData->yuv_format    = sps->chroma_format_idc;
    p_ImgData->frm_stride[0] = sps->PicWidthInMbs * 16;
    p_ImgData->frm_stride[1] =
    p_ImgData->frm_stride[2] = sps->PicWidthInMbs * sps->MbWidthC;
    p_ImgData->top_stride[0] =
    p_ImgData->bot_stride[0] = p_ImgData->frm_stride[0] << 1;
    p_ImgData->top_stride[1] =
    p_ImgData->top_stride[2] =
    p_ImgData->bot_stride[1] =
    p_ImgData->bot_stride[2] = p_ImgData->frm_stride[1] << 1;

    if (sps->separate_colour_plane_flag) {
        for (int nplane = 0; nplane < 3; nplane++)
            get_mem2Dpel(&(p_ImgData->frm_data[nplane]), sps->FrameHeightInMbs * 16, sps->PicWidthInMbs * 16);
    } else {
        get_mem2Dpel(&(p_ImgData->frm_data[0]), sps->FrameHeightInMbs * 16, sps->PicWidthInMbs * 16);

        if (sps->chroma_format_idc != YUV400) {
            get_mem2Dpel(&(p_ImgData->frm_data[1]), sps->FrameHeightInMbs * sps->MbHeightC, sps->PicWidthInMbs * sps->MbWidthC);
            get_mem2Dpel(&(p_ImgData->frm_data[2]), sps->FrameHeightInMbs * sps->MbHeightC, sps->PicWidthInMbs * sps->MbWidthC);

            if (sizeof(imgpel) == sizeof(unsigned char)) {
                for (int k = 1; k < 3; k++)
                    memset(p_ImgData->frm_data[k][0], 128, sps->FrameHeightInMbs * sps->MbHeightC * sps->PicWidthInMbs * sps->MbWidthC * sizeof(imgpel));
            } else {
                imgpel mean_val = (imgpel)sps->BitDepthC;
                for (int k = 1; k < 3; k++) {
                    for (int j = 0; j < sps->FrameHeightInMbs * sps->MbHeightC; j++) {
                        for (int i = 0; i < sps->PicWidthInMbs * sps->MbWidthC; i++)
                            p_ImgData->frm_data[k][j][i] = mean_val;
                    }
                }
            }
        }
    }

    if (!sps->frame_mbs_only_flag) {
        // allocate memory for field reference frame buffers
        init_top_bot_planes(p_ImgData->frm_data[0], sps->FrameHeightInMbs * 16, &(p_ImgData->top_data[0]), &(p_ImgData->bot_data[0]));

        if (sps->chroma_format_idc != YUV400) {
            init_top_bot_planes(p_ImgData->frm_data[1], sps->FrameHeightInMbs * sps->MbHeightC, &(p_ImgData->top_data[1]), &(p_ImgData->bot_data[1]));
            init_top_bot_planes(p_ImgData->frm_data[2], sps->FrameHeightInMbs * sps->MbHeightC, &(p_ImgData->top_data[2]), &(p_ImgData->bot_data[2]));
        }
    }
}

void free_img_data(VideoParameters *p_Vid, ImageData *p_ImgData)
{
    if (p_Vid->active_sps->separate_colour_plane_flag) {
        for (int nplane = 0; nplane < 3; nplane++) {
            if (p_ImgData->frm_data[nplane]) {
                free_mem2Dpel(p_ImgData->frm_data[nplane]);      // free ref frame buffers
                p_ImgData->frm_data[nplane] = NULL;
            }
        }
    } else {
        if (p_ImgData->frm_data[0]) {
            free_mem2Dpel(p_ImgData->frm_data[0]);      // free ref frame buffers
            p_ImgData->frm_data[0] = NULL;
        }
    
        if (p_ImgData->yuv_format != YUV400) {
            if (p_ImgData->frm_data[1]) {
                free_mem2Dpel(p_ImgData->frm_data[1]);
                p_ImgData->frm_data[1] = NULL;
            }
            if (p_ImgData->frm_data[2]) {
                free_mem2Dpel(p_ImgData->frm_data[2]);
                p_ImgData->frm_data[2] = NULL;
            }
        }
    }
  
    if (!p_Vid->active_sps->frame_mbs_only_flag) {
        free_top_bot_planes(p_ImgData->top_data[0], p_ImgData->bot_data[0]);

        if (p_ImgData->yuv_format != YUV400) {
            free_top_bot_planes(p_ImgData->top_data[1], p_ImgData->bot_data[1]);
            free_top_bot_planes(p_ImgData->top_data[2], p_ImgData->bot_data[2]);
        }
    }
}

static void copy_img_data(imgpel *out_img, imgpel *in_img, int ostride, int istride, unsigned int size_y, unsigned int size_x)
{
    for (int i = 0; i < size_y; i++) {
        memcpy(out_img, in_img, size_x);
        out_img += ostride;
        in_img += istride;
    }
}

static void process_picture_in_dpb_s(VideoParameters* p_Vid, storable_picture* p_pic)
{
    ImageData* p_img_out = &p_Vid->tempData3;

    if (!p_Vid->tempData3.frm_data[0])
        init_img_data(p_Vid, &p_Vid->tempData3, p_Vid->active_sps);

    imgpel*** d_img;
    if (p_pic->slice.structure == FRAME)
        d_img = p_img_out->frm_data;
    else { //If reference picture is a field, then frm_data will actually contain field data and therefore top/bottom stride is set accordingly.
        if (p_pic->slice.structure == TOP_FIELD)
            d_img = p_img_out->top_data;
        else
            d_img = p_img_out->bot_data;
    }

    for (int i = 0; i < p_pic->size_y; i++)
        memcpy(d_img[0][i], p_pic->imgY[i], p_pic->size_x*sizeof(imgpel));
    if (p_Vid->active_sps->chroma_format_idc != YUV400) {
        for (int i = 0; i < p_pic->size_y_cr; i++)
            memcpy(d_img[1][i], p_pic->imgUV[0][i], p_pic->size_x_cr * sizeof(imgpel));
        for (int i = 0; i < p_pic->size_y_cr; i++)
            memcpy(d_img[2][i], p_pic->imgUV[1][i], p_pic->size_x_cr * sizeof(imgpel));
    }
}

#if (MVC_EXTENSION_ENABLE)
static storable_picture* clone_storable_picture(VideoParameters* p_Vid, storable_picture* p_pic)
{
    int i, j;
    int nplane;
    int *istride = NULL;
    int ostride[2];
    imgpel ***img_in = NULL;

    sps_t *sps = p_Vid->active_sps;
    storable_picture *p_stored_pic = new storable_picture(p_Vid, (PictureStructure)p_Vid->structure,
        sps->PicWidthInMbs * 16, sps->FrameHeightInMbs * 16,
        sps->PicWidthInMbs * sps->MbWidthC, sps->FrameHeightInMbs * sps->MbHeightC, 0);

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (sps->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (sps->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

    p_stored_pic->PicNum                          = p_pic->PicNum;
    p_stored_pic->frame_num                       = p_pic->frame_num;
    p_stored_pic->LongTermFrameIdx                = p_pic->LongTermFrameIdx;
    p_stored_pic->LongTermPicNum                  = p_pic->LongTermPicNum;
    p_stored_pic->is_long_term                    = 0;
    p_stored_pic->non_existing                    = p_pic->non_existing;
    p_stored_pic->size_x                          = p_pic->size_x;
    p_stored_pic->size_y                          = p_pic->size_y;
    p_stored_pic->size_x_cr                       = p_pic->size_x_cr;
    p_stored_pic->size_y_cr                       = p_pic->size_y_cr;
  
    p_stored_pic->seiHasTone_mapping              = p_pic->seiHasTone_mapping;
    p_stored_pic->poc                             = p_pic->poc;
    p_stored_pic->top_poc                         = p_pic->top_poc;
    p_stored_pic->bottom_poc                      = p_pic->bottom_poc;
    p_stored_pic->frame_poc                       = p_pic->frame_poc;
    p_stored_pic->PicNum                          = p_pic->PicNum;
    p_stored_pic->frame_num                       = p_pic->frame_num;

    p_stored_pic->slice.structure                       = p_pic->slice.structure;
    p_stored_pic->slice.coded_frame                     = 1;
    p_stored_pic->slice.slice_type                      = p_pic->slice.slice_type;
    p_stored_pic->slice.idr_flag                        = p_pic->slice.idr_flag;
    p_stored_pic->slice.mb_aff_frame_flag               = p_pic->slice.mb_aff_frame_flag;
    p_stored_pic->slice.no_output_of_prior_pics_flag    = p_pic->slice.no_output_of_prior_pics_flag;
    p_stored_pic->slice.long_term_reference_flag        = 0;
    p_stored_pic->slice.adaptive_ref_pic_buffering_flag = 0;
    p_stored_pic->slice.dec_ref_pic_marking_buffer      = NULL;
    p_stored_pic->PicWidthInMbs                   = p_pic->PicWidthInMbs;
    p_stored_pic->recovery_frame                  = p_pic->recovery_frame;

    // store BL reconstruction

    ostride[0] = p_stored_pic->iLumaStride;
    ostride[1] = p_stored_pic->iChromaStride;
    if (p_stored_pic->slice.structure == FRAME) {
        istride = p_Vid->tempData3.frm_stride;
        img_in  = p_Vid->tempData3.frm_data;
    } else if (p_stored_pic->slice.structure == TOP_FIELD) {
        istride = p_Vid->tempData3.top_stride;
        img_in  = p_Vid->tempData3.top_data;
    } else {
        istride = p_Vid->tempData3.bot_stride;
        img_in  = p_Vid->tempData3.bot_data;
    }

    copy_img_data(&p_stored_pic->imgY[0][0], &img_in[0][0][0], ostride[0], istride[0], p_pic->size_y, p_pic->size_x * sizeof(imgpel)); 

    pad_buf(*p_stored_pic->imgY, p_stored_pic->size_x, p_stored_pic->size_y, p_stored_pic->iLumaStride, MCBUF_LUMA_PAD_X, MCBUF_LUMA_PAD_Y);

    if (p_Vid->active_sps->chroma_format_idc != YUV400) {    
        copy_img_data(&p_stored_pic->imgUV[0][0][0], &img_in[1][0][0], ostride[1], istride[1], p_pic->size_y_cr, p_pic->size_x_cr*sizeof(imgpel));
        pad_buf(*p_stored_pic->imgUV[0], p_stored_pic->size_x_cr, p_stored_pic->size_y_cr, p_stored_pic->iChromaStride, iChromaPadX, iChromaPadY);
        copy_img_data(&p_stored_pic->imgUV[1][0][0], &img_in[2][0][0], ostride[1], istride[2], p_pic->size_y_cr, p_pic->size_x_cr*sizeof(imgpel));
        pad_buf(*p_stored_pic->imgUV[1], p_stored_pic->size_x_cr, p_stored_pic->size_y_cr, p_stored_pic->iChromaStride, iChromaPadX, iChromaPadY);
    }

    for (j = 0; j < (p_pic->size_y / 4); j++) {
        char *ref_idx = p_stored_pic->mv_info[j][0].ref_idx;
        for (i = 0; i < (p_pic->size_x / 4); i++) {          
            *((short *) ref_idx) = -1;
            ref_idx += sizeof(pic_motion_params);
        }
    }

    if (p_Vid->active_sps->separate_colour_plane_flag != 0) {
        for (nplane = 0; nplane < 3; nplane++) {
            for (j = 0; j < (p_pic->size_y / 4); j++) {
                for (i = 0; i < (p_pic->size_x / 4); i++) {
                    p_stored_pic->JVmv_info[nplane][j][i].ref_idx[LIST_0] = -1;
                    p_stored_pic->JVmv_info[nplane][j][i].ref_idx[LIST_1] = -1;
                }
            }
        }
    }

    // MVC-related parameters
    p_stored_pic->slice.inter_view_flag = p_pic->slice.inter_view_flag;
    p_stored_pic->slice.anchor_pic_flag = 0;
    p_stored_pic->slice.view_id = 0;
    p_stored_pic->is_output = 1;
    p_stored_pic->used_for_reference = 1;

    return p_stored_pic;
}
#endif

void picture_in_dpb(slice_t* currSlice, VideoParameters *p_Vid, storable_picture *p_pic)
{
    process_picture_in_dpb_s(p_Vid, p_pic);
    currSlice->p_Dpb->store_proc_picture(clone_storable_picture(p_Vid, p_pic));
}
