#include "global.h"
#include "input_parameters.h"
#include "h264decoder.h"
#include "dpb.h"
#include "frame_buffer.h"
#include "memalloc.h"
#include "sei.h"
#include "output.h"

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>


static int testEndian(void)
{
    short s;
    byte *p;

    p = (byte *)&s;
    s = 1;

    return (*p == 0);
}

static void img2buf_byte(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride)
{
    int twidth  = size_x - crop_left - crop_right;
    int theight = size_y - crop_top - crop_bottom;
    imgpel** img = &imgX[crop_top];
    for (int i = 0; i < theight; i++) {
        memcpy(buf, *img++ + crop_left, twidth);
        buf += iOutStride;
    }
}

// little endian
static void img2buf_le(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride)
{
    if (sizeof (imgpel) < symbol_size_in_bytes) {
        int twidth  = size_x - crop_left - crop_right;
        int theight = size_y - crop_top - crop_bottom;
        for (int j = 0; j < theight; j++)
            memset(buf + j * iOutStride, 0, twidth * symbol_size_in_bytes);
    }

    int size = min<int>(sizeof(imgpel), symbol_size_in_bytes);

    if (size != 1) {
        for (int i = crop_top; i < size_y - crop_bottom; i++) {
            int ipos = (i - crop_top) * iOutStride;
            for (int j = crop_left; j < size_x - crop_right; j++)
                memcpy(buf + (ipos + (j - crop_left) * symbol_size_in_bytes), &imgX[i][j], size);
        }
    } else {
        for (int j = crop_top; j < size_y - crop_bottom; j++) {  
            uint8_t* pDst = buf + (j - crop_top) * iOutStride;
            imgpel* cur_pixel = &imgX[j][crop_left];
            for (int i = crop_left; i < size_x - crop_right; i++)
                *(pDst++) = (uint8_t)*(cur_pixel++);
        }
    }
}

// big endian
static void img2buf_be(imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes, int crop_left, int crop_right, int crop_top, int crop_bottom, int iOutStride)
{
    switch (symbol_size_in_bytes) {
    case 1:
        for (int i = crop_top; i < size_y - crop_bottom; i++) {
            for (int j = crop_left; j < size_x - crop_right; j++) {
                uint8_t ui8 = (uint8_t)imgX[i][j];
                buf[(j - crop_left) + (i - crop_top) * iOutStride] = ui8;
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
                uint32_t tmp32 = (uint32_t)imgX[i][j];
                uint32_t ui32  = (uint32_t)(((tmp32 & 0x000000FF) << 24) |
                                            ((tmp32 & 0x0000FF00) <<  8) |
                                            ((tmp32 & 0x00FF0000) >>  8) |
                                            ((tmp32 & 0xFF000000) >> 24));
                memcpy(buf + (j - crop_left + (i - crop_top) * iOutStride) * 4, &ui32, 4);
            }
        }
        break;
    default:
        error("writing only to formats of 8, 16 or 32 bit allowed on big endian architecture", 500);
        break;
    }  
}


static void write_out_picture(VideoParameters *p_Vid, storable_picture *p, int p_out)
{
    InputParameters* p_Inp = p_Vid->p_Inp;
    sps_t& sps = *p_Vid->active_sps;

    int pic_unit_bitsize_on_disk = max(sps.BitDepthY, sps.BitDepthC) > 8 ? 16 : 8;
    int symbol_size_in_bytes = (pic_unit_bitsize_on_disk + 7) >> 3;

    decltype(img2buf_le)* img2buf;
    if (sizeof(char) == sizeof(imgpel)) {
        if (sizeof(char) == symbol_size_in_bytes)
            img2buf = img2buf_byte;
        else
            img2buf = img2buf_le;
    } else {
        if (testEndian())
            img2buf = img2buf_be;
        else
            img2buf = img2buf_le;
    }

    bool rgb_output = sps.vui_parameters.matrix_coefficients == 0;

    if (p->non_existing)
        return;

    int size_x_l = sps.PicWidthInMbs    * 16;
    int size_y_l = sps.FrameHeightInMbs * 16;
    int size_x_c = sps.PicWidthInMbs    * sps.MbWidthC;
    int size_y_c = sps.FrameHeightInMbs * sps.MbHeightC;

    // note: this tone-mapping is working for RGB format only. Sharp
    if (p->seiHasTone_mapping && rgb_output) {
        symbol_size_in_bytes = (p->tonemapped_bit_depth > 8) ? 2 : 1;
        tone_map(p->imgY,     p->tone_mapping_lut, size_x_l, size_y_l);
        tone_map(p->imgUV[0], p->tone_mapping_lut, size_x_c, size_y_c);
        tone_map(p->imgUV[1], p->tone_mapping_lut, size_x_c, size_y_c);
    }

    // should this be done only once?
    int crop_left_l = 0, crop_right_l = 0, crop_top_l = 0, crop_bottom_l = 0;
    int crop_left_c = 0, crop_right_c = 0, crop_top_c = 0, crop_bottom_c = 0;
    if (sps.frame_cropping_flag) {
        crop_left_c   = sps.frame_crop_left_offset;
        crop_right_c  = sps.frame_crop_right_offset;
        crop_top_c    = sps.frame_crop_top_offset    * (2 - sps.frame_mbs_only_flag);
        crop_bottom_c = sps.frame_crop_bottom_offset * (2 - sps.frame_mbs_only_flag);
        crop_left_l   = sps.SubWidthC  * crop_left_c;
        crop_right_l  = sps.SubWidthC  * crop_right_c;
        crop_top_l    = sps.SubHeightC * crop_top_c;
        crop_bottom_l = sps.SubHeightC * crop_bottom_c;
    }
    int iLumaSizeX   = size_x_l - (crop_left_l + crop_right_l);
    int iLumaSizeY   = size_y_l - (crop_top_l + crop_bottom_l);
    int iLumaSize    = iLumaSizeX * iLumaSizeY * symbol_size_in_bytes;
    int iChromaSizeX = size_x_c - (crop_left_c + crop_right_c);
    int iChromaSizeY = size_y_c - (crop_top_c + crop_bottom_c);
    int iChromaSize  = iChromaSizeX * iChromaSizeY * symbol_size_in_bytes;
    int iFrameSize   = iLumaSize + 2 * iChromaSize;

    // We need to further cleanup this function
    if (p_out == -1)
        return;

    if (!p_Vid->pDecOuputPic.pY) {
        p_Vid->pDecOuputPic.pY = new uint8_t[iFrameSize];
        p_Vid->pDecOuputPic.pU = p_Vid->pDecOuputPic.pY + iLumaSize;
        p_Vid->pDecOuputPic.pV = p_Vid->pDecOuputPic.pU + iChromaSize;
    }

    if (rgb_output) {
        uint8_t* buf = new uint8_t[size_x_l * size_y_l * symbol_size_in_bytes];
        img2buf(p->imgUV[1], buf, size_x_c, size_y_c, symbol_size_in_bytes,
                crop_left_c, crop_right_c, crop_top_c, crop_bottom_c, iLumaSizeX * symbol_size_in_bytes);
        if (iChromaSize != write(p_out, buf, iChromaSize))
            error("write_out_picture: error writing to RGB file", 500);
        delete []buf;
    }

    img2buf(p->imgY, p_Vid->pDecOuputPic.pY, size_x_l, size_y_l, symbol_size_in_bytes,
            crop_left_l, crop_right_l, crop_top_l, crop_bottom_l, iLumaSizeX * symbol_size_in_bytes);
    if (iLumaSize != write(p_out, p_Vid->pDecOuputPic.pY, iLumaSize))
        error("write_out_picture: error writing to YUV file", 500);

    if (sps.chroma_format_idc != YUV400) {
        img2buf(p->imgUV[0], p_Vid->pDecOuputPic.pU, size_x_c, size_y_c, symbol_size_in_bytes,
                crop_left_c, crop_right_c, crop_top_c, crop_bottom_c, iChromaSizeX * symbol_size_in_bytes);
        if (iChromaSize != write(p_out, p_Vid->pDecOuputPic.pU, iChromaSize))
            error("write_out_picture: error writing to YUV file", 500);

        if (!rgb_output) {
            img2buf(p->imgUV[1], p_Vid->pDecOuputPic.pV, size_x_c, size_y_c, symbol_size_in_bytes,
                    crop_left_c, crop_right_c, crop_top_c, crop_bottom_c, iChromaSizeX * symbol_size_in_bytes);
            if (iChromaSize != write(p_out, p_Vid->pDecOuputPic.pV, iChromaSize))
                error("write_out_picture: error writing to YUV file", 500);
        }
    } else if (p_Inp->write_uv) {
        get_mem2Dpel(&p->imgUV[0], size_y_l / 2, size_x_l / 2);
  
        imgpel cr_val = (imgpel)(1 << (sps.BitDepthY - 1));
        for (int j = 0; j < size_y_l / 2; j++) {
            for (int i = 0; i < size_x_l / 2; i++)
                p->imgUV[0][j][i] = cr_val;
        }

        // fake out U=V=128 to make a YUV 4:2:0 stream
        uint8_t* buf = new uint8_t[size_x_l * size_y_l * symbol_size_in_bytes];
        img2buf(p->imgUV[0], buf, size_x_l/2, size_y_l/2, symbol_size_in_bytes,
                crop_left_l/2, crop_right_l/2, crop_top_l/2, crop_bottom_l/2, iLumaSizeX * symbol_size_in_bytes / 2);
        if (iLumaSize / 4 != write(p_out, buf, iLumaSize / 4))
            error("write_out_picture: error writing to YUV file", 500);
        if (iLumaSize / 4 != write(p_out, buf, iLumaSize / 4))
            error("write_out_picture: error writing to YUV file", 500);
        delete []buf;

        free_mem2Dpel(p->imgUV[0]);
        p->imgUV[0] = nullptr;
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
        fs->bottom_field = new storable_picture(p_Vid, BOTTOM_FIELD,
            p->size_x, p->size_y * 2, p->size_x_cr, p->size_y_cr * 2, 1);
        fs->bottom_field->clear_picture(p_Vid);
        fs->dpb_combine_field_yuv(p_Vid);
#if (MVC_EXTENSION_ENABLE)
        fs->frame->slice.view_id = fs->view_id;
#endif
        write_out_picture(p_Vid, fs->frame, p_out);
    }

    if (fs->is_used & 0x02) {
        // we have a bottom field
        // construct an empty top field
        p = fs->bottom_field;
        fs->top_field = new storable_picture(p_Vid, TOP_FIELD,
            p->size_x, p->size_y * 2, p->size_x_cr, p->size_y_cr * 2, 1);
        fs->top_field->clear_picture(p_Vid);
        fs->dpb_combine_field_yuv(p_Vid);
#if (MVC_EXTENSION_ENABLE)
        fs->frame->slice.view_id = fs->view_id;
#endif
        write_out_picture(p_Vid, fs->frame, p_out);
    }

    fs->is_used = 3;
}

static void flush_direct_output(VideoParameters *p_Vid, int p_out)
{
    write_unpaired_field(p_Vid, p_Vid->out_buffer, p_out);

    if (p_Vid->out_buffer->frame) {
        delete p_Vid->out_buffer->frame;
        p_Vid->out_buffer->frame = nullptr;
    }
    if (p_Vid->out_buffer->top_field) {
        delete p_Vid->out_buffer->top_field;
        p_Vid->out_buffer->top_field = nullptr;
    }
    if (p_Vid->out_buffer->bottom_field) {
        delete p_Vid->out_buffer->bottom_field;
        p_Vid->out_buffer->bottom_field = nullptr;
    }
    p_Vid->out_buffer->is_used = 0;
}


void write_stored_frame(VideoParameters *p_Vid, frame_store *fs, int p_out)
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
    if (p->slice.structure == FRAME) {
        // we have a frame (or complementary field pair)
        // so output it directly
        flush_direct_output(p_Vid, p_out);
        write_out_picture(p_Vid, p, p_out);
        p_Vid->calculate_frame_no(p);
        if (p) {
            delete p;
            p = nullptr;
        }
        return;
    }

    if (p->slice.structure == TOP_FIELD) {
        if (p_Vid->out_buffer->is_used & 1)
            flush_direct_output(p_Vid, p_out);
        p_Vid->out_buffer->top_field = p;
        p_Vid->out_buffer->is_used |= 1;
    }

    if (p->slice.structure == BOTTOM_FIELD) {
        if (p_Vid->out_buffer->is_used & 2)
            flush_direct_output(p_Vid, p_out);
        p_Vid->out_buffer->bottom_field = p;
        p_Vid->out_buffer->is_used |= 2;
    }

    if (p_Vid->out_buffer->is_used == 3) {
        // we have both fields, so output them
        p_Vid->out_buffer->dpb_combine_field_yuv(p_Vid);
#if (MVC_EXTENSION_ENABLE)
        p_Vid->out_buffer->frame->slice.view_id = p_Vid->out_buffer->view_id;
#endif
        write_out_picture(p_Vid, p_Vid->out_buffer->frame, p_out);

        p_Vid->calculate_frame_no(p);
        if (p_Vid->out_buffer->frame) {
            delete p_Vid->out_buffer->frame;
            p_Vid->out_buffer->frame = nullptr;
        }
        if (p_Vid->out_buffer->top_field) {
            delete p_Vid->out_buffer->top_field;
            p_Vid->out_buffer->top_field = nullptr;
        }
        if (p_Vid->out_buffer->bottom_field) {
            delete p_Vid->out_buffer->bottom_field;
            p_Vid->out_buffer->bottom_field = nullptr;
        }
        p_Vid->out_buffer->is_used = 0;
    }
}
