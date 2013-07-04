
#include <math.h>
#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "fmo.h"
#include "bitstream_nal.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"

#include "sei.h"
#include "output.h"
#include "neighbour.h"
#include "memalloc.h"
#include "macroblock.h"
#include "mb.h"
#include "mb_read.h"

#include "intra_prediction.h"
#include "deblock.h"

#include "biaridecod.h"

#include "erc_api.h"
#include "dpb_common.h"
#include "dpb_mvc.h"


#include "dec_slice.h"


static void init_picture_decoding(VideoParameters *p_Vid)
{
    Slice *pSlice = p_Vid->ppSliceList[0];
    int j, i, iDeblockMode = 1;

    if (p_Vid->iSliceNumOfCurrPic >= MAX_NUM_SLICES)
        error ("Maximum number of supported slices exceeded. \nPlease recompile with increased value for MAX_NUM_SLICES", 200);

    if (p_Vid->pNextPPS->Valid && p_Vid->pNextPPS->pic_parameter_set_id == pSlice->pic_parameter_set_id) {
        pps_t tmpPPS;
        memcpy(&tmpPPS, &(p_Vid->PicParSet[pSlice->pic_parameter_set_id]), sizeof (pps_t));
        (p_Vid->PicParSet[pSlice->pic_parameter_set_id]).slice_group_id = NULL;
        MakePPSavailable(p_Vid, p_Vid->pNextPPS->pic_parameter_set_id, p_Vid->pNextPPS);
        memcpy(p_Vid->pNextPPS, &tmpPPS, sizeof (pps_t));
        tmpPPS.slice_group_id = NULL;
    }

    UseParameterSet(pSlice);
    if (pSlice->idr_flag)
        p_Vid->number = 0;

    p_Vid->PicHeightInMbs = p_Vid->FrameHeightInMbs / ( 1 + pSlice->field_pic_flag );
    p_Vid->PicSizeInMbs   = p_Vid->PicWidthInMbs * p_Vid->PicHeightInMbs;
    p_Vid->FrameSizeInMbs = p_Vid->PicWidthInMbs * p_Vid->FrameHeightInMbs;
    p_Vid->structure = pSlice->structure;

    fmo_init(p_Vid, pSlice);

#if (MVC_EXTENSION_ENABLE)
    if (pSlice->layer_id > 0 && pSlice->svc_extension_flag == 0 && pSlice->NaluHeaderMVCExt.non_idr_flag == 0)
        idr_memory_management(p_Vid->p_Dpb_layer[pSlice->layer_id], p_Vid->dec_picture);
    update_ref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_ltref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_pic_num(pSlice);
    i = pSlice->view_id;
#endif
    init_Deblock(p_Vid, pSlice->mb_aff_frame_flag);
    //init mb_data;
    for (j = 0; j < p_Vid->iSliceNumOfCurrPic; j++) {
        if (p_Vid->ppSliceList[j]->disable_deblocking_filter_idc != 1)
            iDeblockMode = 0;
#if (MVC_EXTENSION_ENABLE)
        assert(p_Vid->ppSliceList[j]->view_id == i);
#endif
    }
    p_Vid->iDeblockMode = iDeblockMode;
}


void decode_picture(VideoParameters *p_Vid)
{
    Slice **ppSliceList = p_Vid->ppSliceList;
    int iSliceNo;

    p_Vid->num_dec_mb = 0;

    init_picture_decoding(p_Vid);

    for (iSliceNo = 0; iSliceNo < p_Vid->iSliceNumOfCurrPic; iSliceNo++) {
        Slice *currSlice = ppSliceList[iSliceNo];

        assert(currSlice->current_header != EOS);
        assert(currSlice->current_slice_nr == iSliceNo);

        if (!init_slice(currSlice))
            continue;
        decode_one_slice(currSlice);

        p_Vid->num_dec_mb  += currSlice->num_dec_mb;
        p_Vid->erc_mvperMB += currSlice->erc_mvperMB;
    }

#if MVC_EXTENSION_ENABLE
    p_Vid->last_dec_view_id = p_Vid->dec_picture->view_id;
#endif
    if (p_Vid->dec_picture->structure == FRAME)
        p_Vid->last_dec_poc = p_Vid->dec_picture->frame_poc;
    else if (p_Vid->dec_picture->structure == TOP_FIELD)
        p_Vid->last_dec_poc = p_Vid->dec_picture->top_poc;
    else if (p_Vid->dec_picture->structure == BOTTOM_FIELD)
        p_Vid->last_dec_poc = p_Vid->dec_picture->bottom_poc;

    exit_picture(p_Vid, &p_Vid->dec_picture);

    p_Vid->previous_frame_num = ppSliceList[0]->frame_num;
}


static void update_mbaff_macroblock_data(imgpel **cur_img, imgpel (*temp)[16], int x0, int width, int height)
{
    imgpel (*temp_evn)[16] = temp;
    imgpel (*temp_odd)[16] = temp + height; 
    imgpel **temp_img = cur_img;
    int y;

    for (y = 0; y < 2 * height; ++y)
        memcpy(*temp++, (*temp_img++ + x0), width * sizeof(imgpel));

    for (y = 0; y < height; ++y) {
        memcpy((*cur_img++ + x0), *temp_evn++, width * sizeof(imgpel));
        memcpy((*cur_img++ + x0), *temp_odd++, width * sizeof(imgpel));
    }
}

static void MbAffPostProc(VideoParameters *p_Vid)
{
    imgpel temp_buffer[32][16];

    StorablePicture *dec_picture = p_Vid->dec_picture;
    imgpel ** imgY  = dec_picture->imgY;
    imgpel ***imgUV = dec_picture->imgUV;

    short i, x0, y0;

    for (i = 0; i < (int)dec_picture->PicSizeInMbs; i += 2) {
        if (dec_picture->motion.mb_field_decoding_flag[i]) {
            get_mb_pos(p_Vid, i, p_Vid->mb_size[IS_LUMA], &x0, &y0);
            update_mbaff_macroblock_data(imgY + y0, temp_buffer, x0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

            if (dec_picture->chroma_format_idc != YUV400) {
                x0 = (short) ((x0 * p_Vid->mb_cr_size_x) >> 4);
                y0 = (short) ((y0 * p_Vid->mb_cr_size_y) >> 4);

                update_mbaff_macroblock_data(imgUV[0] + y0, temp_buffer, x0, p_Vid->mb_cr_size_x, p_Vid->mb_cr_size_y);
                update_mbaff_macroblock_data(imgUV[1] + y0, temp_buffer, x0, p_Vid->mb_cr_size_x, p_Vid->mb_cr_size_y);
            }
        }
    }
}

void pad_buf(imgpel *pImgBuf, int iWidth, int iHeight, int iStride, int iPadX, int iPadY)
{
    int j;
    imgpel *pLine0 = pImgBuf - iPadX, *pLine;
    int i;
    for (i = -iPadX; i < 0; i++)
        pImgBuf[i] = *pImgBuf;
    for (i = 0; i < iPadX; i++)
        pImgBuf[i+iWidth] = *(pImgBuf+iWidth-1);

    for (j = -iPadY; j < 0; j++)
        memcpy(pLine0+j*iStride, pLine0, iStride*sizeof(imgpel));
    for (j = 1; j < iHeight; j++) {
        pLine = pLine0 + j*iStride;
        for (i = 0; i < iPadX; i++)
            pLine[i] = pLine[iPadX];
        pLine += iPadX+iWidth-1;
        for (i = 1; i < iPadX + 1; i++)
            pLine[i] = *pLine;
    }
    pLine = pLine0 + (iHeight - 1) * iStride;
    for (j = iHeight; j < iHeight + iPadY; j++)
        memcpy(pLine0+j*iStride,  pLine, iStride*sizeof(imgpel));
}

void pad_dec_picture(VideoParameters *p_Vid, StorablePicture *dec_picture)
{
    int iPadX = p_Vid->iLumaPadX;
    int iPadY = p_Vid->iLumaPadY;
    int iWidth = dec_picture->size_x;
    int iHeight = dec_picture->size_y;
    int iStride = dec_picture->iLumaStride;

    pad_buf(*dec_picture->imgY, iWidth, iHeight, iStride, iPadX, iPadY);

    if (dec_picture->chroma_format_idc != YUV400) {
        iPadX = p_Vid->iChromaPadX;
        iPadY = p_Vid->iChromaPadY;
        iWidth = dec_picture->size_x_cr;
        iHeight = dec_picture->size_y_cr;
        iStride = dec_picture->iChromaStride;
        pad_buf(*dec_picture->imgUV[0], iWidth, iHeight, iStride, iPadX, iPadY);
        pad_buf(*dec_picture->imgUV[1], iWidth, iHeight, iStride, iPadX, iPadY);
    }
}

static void status_picture(VideoParameters *p_Vid, StorablePicture **dec_picture)
{
    InputParameters *p_Inp = p_Vid->p_Inp;
    SNRParameters   *snr   = p_Vid->snr;

    char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};
    char yuvFormat[10];
    int64 tmp_time;                   // time used by decoding the last frame

    int structure         = (*dec_picture)->structure;
    int slice_type        = (*dec_picture)->slice_type;
    int frame_poc         = (*dec_picture)->frame_poc;  
    int refpic            = (*dec_picture)->used_for_reference;
    int qp                = (*dec_picture)->qp;
    int pic_num           = (*dec_picture)->pic_num;
    int is_idr            = (*dec_picture)->idr_flag;
    int chroma_format_idc = (*dec_picture)->chroma_format_idc;

    // report
    char cslice_type[9];  

    if (p_Inp->silent == FALSE) {
        if (structure == TOP_FIELD || structure == FRAME) {
            if (slice_type == I_SLICE && is_idr) // IDR picture
                strcpy(cslice_type,"IDR");
            else if (slice_type == I_SLICE) // I picture
                strcpy(cslice_type," I ");
            else if (slice_type == P_SLICE) // P pictures
                strcpy(cslice_type," P ");
            else if (slice_type == SP_SLICE) // SP pictures
                strcpy(cslice_type,"SP ");
            else if (slice_type == SI_SLICE)
                strcpy(cslice_type,"SI ");
            else if (refpic) // stored B pictures
                strcpy(cslice_type," B ");
            else // B pictures
                strcpy(cslice_type," b ");

            if (structure == FRAME)
                strncat(cslice_type,")       ",8-strlen(cslice_type));
        } else if (structure == BOTTOM_FIELD) {
            if (slice_type == I_SLICE && is_idr) // IDR picture
                strncat(cslice_type,"|IDR)",8-strlen(cslice_type));
            else if (slice_type == I_SLICE) // I picture
                strncat(cslice_type,"| I )",8-strlen(cslice_type));
            else if (slice_type == P_SLICE) // P pictures
                strncat(cslice_type,"| P )",8-strlen(cslice_type));
            else if (slice_type == SP_SLICE) // SP pictures
                strncat(cslice_type,"|SP )",8-strlen(cslice_type));
            else if (slice_type == SI_SLICE)
                strncat(cslice_type,"|SI )",8-strlen(cslice_type));
            else if (refpic) // stored B pictures
                strncat(cslice_type,"| B )",8-strlen(cslice_type));
            else // B pictures
                strncat(cslice_type,"| b )",8-strlen(cslice_type));   
        }
    }

    if (structure == FRAME || structure == BOTTOM_FIELD) {
        gettime (&(p_Vid->end_time));              // end time

        tmp_time  = timediff(&(p_Vid->start_time), &(p_Vid->end_time));
        p_Vid->tot_time += tmp_time;
        tmp_time  = timenorm(tmp_time);
        sprintf(yuvFormat,"%s", yuv_types[chroma_format_idc]);

        if (p_Inp->silent == FALSE) {
            SNRParameters   *snr = p_Vid->snr;
            if (p_Vid->p_ref != -1)
                fprintf(stdout,"%05d(%s%5d %5d %5d %8.4f %8.4f %8.4f  %s %7d\n",
                        p_Vid->frame_no, cslice_type, frame_poc, pic_num, qp, snr->snr[0], snr->snr[1], snr->snr[2], yuvFormat, (int) tmp_time);
            else
                fprintf(stdout,"%05d(%s%5d %5d %5d                             %s %7d\n",
                        p_Vid->frame_no, cslice_type, frame_poc, pic_num, qp, yuvFormat, (int)tmp_time);
        } else
            fprintf(stdout,"Completed Decoding frame %05d.\r",snr->frame_ctr);

        fflush(stdout);

        if (slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic) { // I or P pictures
#if (MVC_EXTENSION_ENABLE)
            if((p_Vid->ppSliceList[0])->view_id!=0)
#endif
                ++(p_Vid->number);
        } else
            ++(p_Vid->Bframe_ctr);    // B pictures
        ++(snr->frame_ctr);

#if (MVC_EXTENSION_ENABLE)
        if ((p_Vid->ppSliceList[0])->view_id != 0)
#endif
            ++(p_Vid->g_nFrame);   
    }
}

void exit_picture(VideoParameters *p_Vid, StorablePicture **dec_picture)
{
    // return if the last picture has already been finished
    if (*dec_picture == NULL ||
        (p_Vid->num_dec_mb != p_Vid->PicSizeInMbs &&
         (p_Vid->yuv_format != YUV444 || !p_Vid->separate_colour_plane_flag)))
        return;

#if (DISABLE_ERC == 0)
    erc_picture(p_Vid, dec_picture);
#endif

    pic_deblock(p_Vid, *dec_picture);

    if ((*dec_picture)->mb_aff_frame_flag)
        MbAffPostProc(p_Vid);

    if (p_Vid->structure != FRAME)
        p_Vid->number /= 2;
#if (MVC_EXTENSION_ENABLE)
    if ((*dec_picture)->used_for_reference || ((*dec_picture)->inter_view_flag == 1))
        pad_dec_picture(p_Vid, *dec_picture);
#endif
#if MVC_EXTENSION_ENABLE
    store_picture_in_dpb(p_Vid->p_Dpb_layer[(*dec_picture)->view_id], *dec_picture);
#endif

    if (p_Vid->last_has_mmco_5)
        p_Vid->pre_frame_num = 0;

    status_picture(p_Vid, dec_picture);

    *dec_picture = NULL;
}
