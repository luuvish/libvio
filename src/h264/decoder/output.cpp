#include "global.h"
#include "input_parameters.h"
#include "h264decoder.h"
#include "dpb.h"
#include "image.h"
#include "memalloc.h"
#include "sei.h"
#include "output.h"


static int testEndian(void)
{
    short s;
    byte *p;

    p = (byte *)&s;
    s = 1;

    return (*p == 0);
}

static void img2buf_normal(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride)
{
    int size = 0;

    // sizeof (imgpel) > sizeof(char)
    // little endian
    if (sizeof (imgpel) < symbol_size_in_bytes) {
        // this should not happen. we should not have smaller imgpel than our source material.
        size = sizeof(imgpel);
        // clear buffer
        int twidth  = size_x - crop_left - crop_right;
        int theight = size_y - crop_top - crop_bottom;
        for (int j = 0; j < theight; j++)
            memset(buf + j * iOutStride, 0, twidth * symbol_size_in_bytes);
    } else
        size = symbol_size_in_bytes;

    if ((crop_top || crop_bottom || crop_left || crop_right) || (size != 1)) {
        for (int i = crop_top; i < size_y - crop_bottom; i++) {
            int ipos = (i - crop_top) * iOutStride;
            for (int j = crop_left; j < size_x - crop_right; j++)
                memcpy(buf + (ipos + (j - crop_left) * symbol_size_in_bytes), &imgX[i][j], size);
        }
    } else {
        for (int j = 0; j < size_y; j++) {  
            imgpel* cur_pixel = imgX[j];
            unsigned char* pDst = buf + j * iOutStride;
            for (int i = 0; i < size_x; i++)
                *(pDst++) = (unsigned char)*(cur_pixel++);
        }
    }
}

static void img2buf_byte(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride)
{
    int twidth  = size_x - crop_left - crop_right;
    int theight = size_y - crop_top - crop_bottom;
    imgpel **img = &imgX[crop_top];
    for (int i = 0; i < theight; i++) {
        memcpy(buf, *img++ + crop_left, twidth);
        buf += iOutStride;
    }
}

static void img2buf_endian(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride)
{
    // big endian
    switch (symbol_size_in_bytes) {
    case 1:
        for (int i = crop_top; i < size_y - crop_bottom; i++) {
            for (int j = crop_left; j < size_x - crop_right; j++) {
                unsigned char ui8 = (unsigned char)imgX[i][j];
                buf[j - crop_left + (i - crop_top) * iOutStride] = ui8;
            }
        }
        break;
    case 2:
        for (int i = crop_top; i < size_y - crop_bottom; i++) {
            for (int j = crop_left; j < size_x - crop_right; j++) {
                uint16_t tmp16 = (uint16_t)imgX[i][j];
                uint16_t ui16  = (uint16_t)((tmp16 >> 8) | ((tmp16 & 0xFF) << 8));
                memcpy(buf + (j - crop_left + (i - crop_top) * iOutStride) * 2, &ui16, 2);
            }
        }
        break;
    case 4:
        for (int i = crop_top; i < size_y - crop_bottom; i++) {
            for (int j = crop_left; j < size_x - crop_right; j++) {
                unsigned long tmp32 = (unsigned long)imgX[i][j];
                unsigned long ui32  = (unsigned long)(((tmp32 & 0xFF00  ) << 8) | ((tmp32 &       0xFF) << 24) |
                                                      ((tmp32 & 0xFF0000) >> 8) | ((tmp32 & 0xFF000000) >> 24));
                memcpy(buf + ((j - crop_left + ((i-crop_top) * iOutStride)) * 4), &ui32, 4);
            }
        }
        break;
    default:
        error("writing only to formats of 8, 16 or 32 bit allowed on big endian architecture", 500);
        break;
    }  
}


static void allocate_p_dec_pic(VideoParameters *p_Vid, DecodedPicList *pDecPic, storable_picture *p, int iLumaSize, int iFrameSize, int iLumaSizeX, int iLumaSizeY, int iChromaSizeX, int iChromaSizeY)
{
    sps_t *sps = p_Vid->active_sps;
    int pic_unit_bitsize_on_disk;
    if (sps->BitDepthY > sps->BitDepthC || sps->chroma_format_idc == YUV400)
        pic_unit_bitsize_on_disk = (sps->BitDepthY > 8) ? 16 : 8;
    else
        pic_unit_bitsize_on_disk = (sps->BitDepthC > 8) ? 16 : 8;
    int symbol_size_in_bytes = ((pic_unit_bitsize_on_disk + 7) >> 3);
  
    if (pDecPic->pY)
        delete []pDecPic->pY;
    pDecPic->iBufSize = iFrameSize;
    pDecPic->pY = new uint8_t[pDecPic->iBufSize];

    pDecPic->pU = pDecPic->pY+iLumaSize;
    pDecPic->pV = pDecPic->pU + ((iFrameSize-iLumaSize)>>1);
    //init;
    pDecPic->iYUVFormat        = p->chroma_format_idc;
    pDecPic->iYUVStorageFormat = 0;
    pDecPic->iBitDepth         = pic_unit_bitsize_on_disk;
    pDecPic->iWidth            = iLumaSizeX;
    pDecPic->iHeight           = iLumaSizeY;
    pDecPic->iYBufStride       = iLumaSizeX*symbol_size_in_bytes;
    pDecPic->iUVBufStride      = iChromaSizeX*symbol_size_in_bytes;
}


static DecodedPicList* get_one_avail_dec_pic_from_list(DecodedPicList *pDecPicList, int b3D, int view_id)
{
    DecodedPicList *pPic = pDecPicList, *pPrior = NULL;
    if (b3D) {
        while (pPic && (pPic->bValid & (1 << view_id))) {
            pPrior = pPic;
            pPic = pPic->pNext;
        }
    } else {
        while (pPic && pPic->bValid) {
            pPrior = pPic;
            pPic = pPic->pNext;
        }
    }

    if (!pPic) {
        pPic = (DecodedPicList *)calloc(1, sizeof(*pPic));
        pPrior->pNext = pPic;
    }

    return pPic;
}

static void write_out_picture(VideoParameters *p_Vid, storable_picture *p, int p_out)
{
    InputParameters *p_Inp = p_Vid->p_Inp;
    DecodedPicList *pDecPic;
    sps_t *sps = p_Vid->active_sps;

    auto img2buf = img2buf_normal;
    if (sizeof(char) == sizeof(imgpel)) {
        int pic_unit_bitsize_on_disk = max(sps->BitDepthY, sps->BitDepthC) > 8 ? 16 : 8;
        int symbol_size_in_bytes = (pic_unit_bitsize_on_disk + 7) >> 3;
        if (sizeof(char) == symbol_size_in_bytes)
            img2buf = img2buf_byte;
        else
            img2buf = img2buf_normal;
    } else {
        if (testEndian())
            img2buf = img2buf_endian;
        else
            img2buf = img2buf_normal;
    }

    int pic_unit_bitsize_on_disk;
    if (sps->BitDepthY > sps->BitDepthC || sps->chroma_format_idc == YUV400)
        pic_unit_bitsize_on_disk = (sps->BitDepthY > 8) ? 16 : 8;
    else
        pic_unit_bitsize_on_disk = (sps->BitDepthC > 8) ? 16 : 8;

    static const int SubWidthC  [4]= { 1, 2, 2, 1};
    static const int SubHeightC [4]= { 1, 2, 1, 1};

    int crop_left, crop_right, crop_top, crop_bottom;
    int symbol_size_in_bytes = ((pic_unit_bitsize_on_disk + 7) >> 3);
    bool rgb_output = sps->vui_parameters.matrix_coefficients == 0;
    unsigned char *buf;
    int iLumaSize, iFrameSize;
    int iLumaSizeX, iLumaSizeY;
    int iChromaSizeX, iChromaSizeY;

    int ret;

    if (p->non_existing)
        return;

    // note: this tone-mapping is working for RGB format only. Sharp
    if (p->seiHasTone_mapping && rgb_output) {
        symbol_size_in_bytes = (p->tonemapped_bit_depth > 8) ? 2 : 1;
        tone_map(p->imgY, p->tone_mapping_lut, p->size_x, p->size_y);
        tone_map(p->imgUV[0], p->tone_mapping_lut, p->size_x_cr, p->size_y_cr);
        tone_map(p->imgUV[1], p->tone_mapping_lut, p->size_x_cr, p->size_y_cr);
    }

    // should this be done only once?
    if (p->frame_cropping_flag) {
        crop_left   = SubWidthC [p->chroma_format_idc] * p->frame_crop_left_offset;
        crop_right  = SubWidthC [p->chroma_format_idc] * p->frame_crop_right_offset;
        crop_top    = SubHeightC[p->chroma_format_idc] * (2 - p->frame_mbs_only_flag) * p->frame_crop_top_offset;
        crop_bottom = SubHeightC[p->chroma_format_idc] * (2 - p->frame_mbs_only_flag) * p->frame_crop_bottom_offset;
    } else
        crop_left = crop_right = crop_top = crop_bottom = 0;
    iChromaSizeX =  p->size_x_cr- p->frame_crop_left_offset -p->frame_crop_right_offset;
    iChromaSizeY = p->size_y_cr - ( 2 - p->frame_mbs_only_flag ) * p->frame_crop_top_offset -( 2 - p->frame_mbs_only_flag ) * p->frame_crop_bottom_offset;
    iLumaSizeX = p->size_x - crop_left-crop_right;
    iLumaSizeY = p->size_y - crop_top - crop_bottom;
    iLumaSize  = iLumaSizeX * iLumaSizeY * symbol_size_in_bytes;
    iFrameSize = (iLumaSizeX * iLumaSizeY + 2 * (iChromaSizeX * iChromaSizeY)) * symbol_size_in_bytes;

    // We need to further cleanup this function
    if (p_out == -1)
        return;

    // KS: this buffer should actually be allocated only once, but this is still much faster than the previous version
    pDecPic = get_one_avail_dec_pic_from_list(p_Vid->pDecOuputPic, 0, 0);
    if (pDecPic->pY == NULL || pDecPic->iBufSize < iFrameSize)
        allocate_p_dec_pic(p_Vid, pDecPic, p, iLumaSize, iFrameSize, iLumaSizeX, iLumaSizeY, iChromaSizeX, iChromaSizeY);

#if (MVC_EXTENSION_ENABLE)
    pDecPic->bValid = 1;
    pDecPic->iViewId = p->view_id >=0 ? p->view_id : -1;
#else
    pDecPic->bValid = 1;
#endif
    pDecPic->iPOC = p->frame_poc;
  
    if (NULL == pDecPic->pY)
        no_mem_exit("write_out_picture: buf");

    if (rgb_output) {
        buf = new unsigned char[p->size_x * p->size_y * symbol_size_in_bytes];
        crop_left   = p->frame_crop_left_offset;
        crop_right  = p->frame_crop_right_offset;
        crop_top    = ( 2 - p->frame_mbs_only_flag ) * p->frame_crop_top_offset;
        crop_bottom = ( 2 - p->frame_mbs_only_flag ) * p->frame_crop_bottom_offset;

        img2buf(p->imgUV[1], buf, p->size_x_cr, p->size_y_cr, symbol_size_in_bytes, crop_left, crop_right, crop_top, crop_bottom, pDecPic->iYBufStride);
        if (p_out >= 0) {
            ret = write(p_out, buf, (p->size_y_cr-crop_bottom-crop_top)*(p->size_x_cr-crop_right-crop_left)*symbol_size_in_bytes);
            if (ret != ((p->size_y_cr-crop_bottom-crop_top)*(p->size_x_cr-crop_right-crop_left)*symbol_size_in_bytes))
                error("write_out_picture: error writing to RGB file", 500);
        }

        if (p->frame_cropping_flag) {
            crop_left   = SubWidthC[p->chroma_format_idc] * p->frame_crop_left_offset;
            crop_right  = SubWidthC[p->chroma_format_idc] * p->frame_crop_right_offset;
            crop_top    = SubHeightC[p->chroma_format_idc]*( 2 - p->frame_mbs_only_flag ) * p->frame_crop_top_offset;
            crop_bottom = SubHeightC[p->chroma_format_idc]*( 2 - p->frame_mbs_only_flag ) * p->frame_crop_bottom_offset;
        } else
            crop_left = crop_right = crop_top = crop_bottom = 0;
        delete []buf;
    }

    buf = (pDecPic->bValid == 1) ? pDecPic->pY : pDecPic->pY + iLumaSizeX * symbol_size_in_bytes;

    img2buf(p->imgY, buf, p->size_x, p->size_y, symbol_size_in_bytes, crop_left, crop_right, crop_top, crop_bottom, pDecPic->iYBufStride);
    if (p_out >= 0) {
        ret = write(p_out, buf, (p->size_y-crop_bottom-crop_top)*(p->size_x-crop_right-crop_left)*symbol_size_in_bytes);
        if (ret != ((p->size_y-crop_bottom-crop_top)*(p->size_x-crop_right-crop_left)*symbol_size_in_bytes))
            error("write_out_picture: error writing to YUV file", 500);
    }

    if (p->chroma_format_idc != YUV400) {
        crop_left   = p->frame_crop_left_offset;
        crop_right  = p->frame_crop_right_offset;
        crop_top    = ( 2 - p->frame_mbs_only_flag ) * p->frame_crop_top_offset;
        crop_bottom = ( 2 - p->frame_mbs_only_flag ) * p->frame_crop_bottom_offset;
        buf = (pDecPic->bValid==1)? pDecPic->pU : pDecPic->pU + iChromaSizeX*symbol_size_in_bytes;
        img2buf(p->imgUV[0], buf, p->size_x_cr, p->size_y_cr, symbol_size_in_bytes, crop_left, crop_right, crop_top, crop_bottom, pDecPic->iUVBufStride);
        if (p_out >= 0) {
            ret = write(p_out, buf, (p->size_y_cr-crop_bottom-crop_top)*(p->size_x_cr-crop_right-crop_left)* symbol_size_in_bytes);
            if (ret != ((p->size_y_cr-crop_bottom-crop_top)*(p->size_x_cr-crop_right-crop_left)* symbol_size_in_bytes))
                error("write_out_picture: error writing to YUV file", 500);
        }

        if (!rgb_output) {
            buf = (pDecPic->bValid==1)? pDecPic->pV : pDecPic->pV + iChromaSizeX*symbol_size_in_bytes;
            img2buf(p->imgUV[1], buf, p->size_x_cr, p->size_y_cr, symbol_size_in_bytes, crop_left, crop_right, crop_top, crop_bottom, pDecPic->iUVBufStride);

            if (p_out >= 0) {
                ret = write(p_out, buf, (p->size_y_cr-crop_bottom-crop_top)*(p->size_x_cr-crop_right-crop_left)*symbol_size_in_bytes);
                if (ret != ((p->size_y_cr-crop_bottom-crop_top)*(p->size_x_cr-crop_right-crop_left)*symbol_size_in_bytes))
                    error("write_out_picture: error writing to YUV file", 500);
            }
        }
    } else {
        if (p_Inp->write_uv) {
            int i, j;
            imgpel cr_val = (imgpel) (1<<(p_Vid->active_sps->BitDepthY - 1));

            get_mem3Dpel(&p->imgUV, 1, p->size_y/2, p->size_x/2);
      
            for (j = 0; j < p->size_y / 2; j++) {
                for (i = 0; i < p->size_x / 2; i++)
                    p->imgUV[0][j][i] = cr_val;
            }

            // fake out U=V=128 to make a YUV 4:2:0 stream
            buf = new unsigned char[p->size_x * p->size_y * symbol_size_in_bytes];
            img2buf(p->imgUV[0], buf, p->size_x/2, p->size_y/2, symbol_size_in_bytes, crop_left/2, crop_right/2, crop_top/2, crop_bottom/2, pDecPic->iYBufStride/2);

            ret = write(p_out, buf, symbol_size_in_bytes * (p->size_y-crop_bottom-crop_top)/2 * (p->size_x-crop_right-crop_left)/2 );
            if (ret != (symbol_size_in_bytes * (p->size_y-crop_bottom-crop_top)/2 * (p->size_x-crop_right-crop_left)/2))
                error("write_out_picture: error writing to YUV file", 500);
            ret = write(p_out, buf, symbol_size_in_bytes * (p->size_y-crop_bottom-crop_top)/2 * (p->size_x-crop_right-crop_left)/2 );
            if (ret != (symbol_size_in_bytes * (p->size_y-crop_bottom-crop_top)/2 * (p->size_x-crop_right-crop_left)/2))
                error("write_out_picture: error writing to YUV file", 500);
            delete []buf;
            free_mem3Dpel(p->imgUV);
            p->imgUV = NULL;
        }
    }

    if (p_out >= 0)
        pDecPic->bValid = 0;
}

static void clear_picture(VideoParameters *p_Vid, storable_picture *p)
{
    sps_t *sps = p_Vid->active_sps;
    int i,j;

    for (i = 0; i < p->size_y; i++) {
        for (j = 0; j < p->size_x; j++)
            p->imgY[i][j] = (imgpel) (1 << (sps->BitDepthY - 1));
    }
    for (i = 0; i < p->size_y_cr; i++) {
        for (j = 0; j < p->size_x_cr; j++)
            p->imgUV[0][i][j] = (imgpel) (1 << (sps->BitDepthC - 1));
    }
    for (i = 0; i < p->size_y_cr; i++) {
        for (j = 0; j < p->size_x_cr; j++)
            p->imgUV[1][i][j] = (imgpel) (1 << (sps->BitDepthC - 1));
    }
}

static void write_unpaired_field(VideoParameters *p_Vid, frame_store* fs, int p_out)
{
    storable_picture *p;
    assert (fs->is_used < 3);

    if (fs->is_used & 0x01) {
        // we have a top field
        // construct an empty bottom field
        p = fs->top_field;
        fs->bottom_field = alloc_storable_picture(p_Vid, BOTTOM_FIELD, p->size_x, 2*p->size_y, p->size_x_cr, 2*p->size_y_cr, 1);
        fs->bottom_field->chroma_format_idc = p->chroma_format_idc;
        clear_picture(p_Vid, fs->bottom_field);
        dpb_combine_field_yuv(p_Vid, fs);
#if (MVC_EXTENSION_ENABLE)
        fs->frame->view_id = fs->view_id;
#endif
        write_out_picture(p_Vid, fs->frame, p_out);
    }

    if (fs->is_used & 0x02) {
        // we have a bottom field
        // construct an empty top field
        p = fs->bottom_field;
        fs->top_field = alloc_storable_picture(p_Vid, TOP_FIELD, p->size_x, 2*p->size_y, p->size_x_cr, 2*p->size_y_cr, 1);
        fs->top_field->chroma_format_idc = p->chroma_format_idc;
        clear_picture(p_Vid, fs->top_field);
        fs->top_field->frame_cropping_flag = fs->bottom_field->frame_cropping_flag;
        if (fs ->top_field->frame_cropping_flag) {
            fs->top_field->frame_crop_top_offset = fs->bottom_field->frame_crop_top_offset;
            fs->top_field->frame_crop_bottom_offset = fs->bottom_field->frame_crop_bottom_offset;
            fs->top_field->frame_crop_left_offset = fs->bottom_field->frame_crop_left_offset;
            fs->top_field->frame_crop_right_offset = fs->bottom_field->frame_crop_right_offset;
        }
        dpb_combine_field_yuv(p_Vid, fs);
#if (MVC_EXTENSION_ENABLE)
        fs->frame->view_id = fs->view_id;
#endif
        write_out_picture(p_Vid, fs->frame, p_out);
    }

    fs->is_used = 3;
}

static void flush_direct_output(VideoParameters *p_Vid, int p_out)
{
    write_unpaired_field(p_Vid, p_Vid->out_buffer, p_out);

    free_storable_picture(p_Vid->out_buffer->frame);
    p_Vid->out_buffer->frame = NULL;
    free_storable_picture(p_Vid->out_buffer->top_field);
    p_Vid->out_buffer->top_field = NULL;
    free_storable_picture(p_Vid->out_buffer->bottom_field);
    p_Vid->out_buffer->bottom_field = NULL;
    p_Vid->out_buffer->is_used = 0;
}


void write_stored_frame( VideoParameters *p_Vid, frame_store *fs, int p_out)
{
    // make sure no direct output field is pending
    flush_direct_output(p_Vid, p_out);

    if (fs->is_used < 3)
        write_unpaired_field(p_Vid, fs, p_out);
    else {
        if (fs->recovery_frame)
            p_Vid->recovery_flag = 1;
        if (!p_Vid->non_conforming_stream || p_Vid->recovery_flag)
            write_out_picture(p_Vid, fs->frame, p_out);
    }

    fs->is_output = 1;
}

void direct_output(VideoParameters *p_Vid, storable_picture *p, int p_out)
{
    if (p->structure == FRAME) {
        // we have a frame (or complementary field pair)
        // so output it directly
        flush_direct_output(p_Vid, p_out);
        write_out_picture (p_Vid, p, p_out);
        p_Vid->calculate_frame_no(p);
        free_storable_picture(p);
        return;
    }

    if (p->structure == TOP_FIELD) {
        if (p_Vid->out_buffer->is_used & 1)
            flush_direct_output(p_Vid, p_out);
        p_Vid->out_buffer->top_field = p;
        p_Vid->out_buffer->is_used |= 1;
    }

    if (p->structure == BOTTOM_FIELD) {
        if (p_Vid->out_buffer->is_used & 2)
            flush_direct_output(p_Vid, p_out);
        p_Vid->out_buffer->bottom_field = p;
        p_Vid->out_buffer->is_used |= 2;
    }

    if (p_Vid->out_buffer->is_used == 3) {
        // we have both fields, so output them
        dpb_combine_field_yuv(p_Vid, p_Vid->out_buffer);
#if (MVC_EXTENSION_ENABLE)
        p_Vid->out_buffer->frame->view_id = p_Vid->out_buffer->view_id;
#endif
        write_out_picture(p_Vid, p_Vid->out_buffer->frame, p_out);

        p_Vid->calculate_frame_no(p);
        free_storable_picture(p_Vid->out_buffer->frame);
        p_Vid->out_buffer->frame = NULL;
        free_storable_picture(p_Vid->out_buffer->top_field);
        p_Vid->out_buffer->top_field = NULL;
        free_storable_picture(p_Vid->out_buffer->bottom_field);
        p_Vid->out_buffer->bottom_field = NULL;
        p_Vid->out_buffer->is_used = 0;
    }
}
