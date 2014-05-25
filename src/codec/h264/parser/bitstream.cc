/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 *
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 *
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : bitstream.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

#include "memalloc.h"
#include "bitstream.h"
#include "sets.h"


namespace vio  {
namespace h264 {


static void nal_unit_header_mvc_extension(nal_unit_t& nal)
{
    nal.non_idr_flag     = (nal.rbsp_byte[1] >> 6) & 1;
    nal.priority_id      = (nal.rbsp_byte[1]     ) & 63;
    nal.view_id          = (nal.rbsp_byte[2] << 2) | ((nal.rbsp_byte[3] >> 6) & 3);
    nal.temporal_id      = (nal.rbsp_byte[3] >> 3) & 7;
    nal.anchor_pic_flag  = (nal.rbsp_byte[3] >> 2) & 1;
    nal.inter_view_flag  = (nal.rbsp_byte[3] >> 1) & 1;
    nal.reserved_one_bit = (nal.rbsp_byte[3]     ) & 1;

    if (nal.reserved_one_bit != 1)
        printf("Nalu Header MVC Extension: reserved_one_bit is not 1!\n");

    if (nal.nal_unit_type == 20)
        nal.nal_unit_type = 1;
}

static void nal_unit_header_svc_extension(nal_unit_t& nal)
{
    //to be implemented for Annex G;
}

// 7.3.1 NAL unit syntax

static void nal_unit(nal_unit_t& nal)
{
    nal.forbidden_zero_bit = (nal.rbsp_byte[0] >> 7) & 1;
    nal.nal_ref_idc        = (nal.rbsp_byte[0] >> 5) & 3;
    nal.nal_unit_type      = (nal.rbsp_byte[0] & 0x1f);

    nal.mvc_extension_flag = 0;
    nal.svc_extension_flag = 0;

    if (nal.nal_unit_type == 14 || nal.nal_unit_type == 20 || nal.nal_unit_type == 21) {
        nal.svc_extension_flag = (nal.rbsp_byte[1] >> 7) & 1;
        nal.mvc_extension_flag = ~nal.svc_extension_flag;
        if (nal.svc_extension_flag)
            nal_unit_header_svc_extension(nal);
        else
            nal_unit_header_mvc_extension(nal);
    }
}

static int NALUtoRBSP(nal_unit_t& nal)
{
    int nalUnitHeaderBytes = 1;
    if (nal.nal_unit_type == 14 || nal.nal_unit_type == 20 || nal.nal_unit_type == 21)
        nalUnitHeaderBytes += 3;

    if (nal.num_bytes_in_nal_unit < nalUnitHeaderBytes)
        return nal.num_bytes_in_rbsp = nal.num_bytes_in_nal_unit;

    int count = 0;
    int j = nalUnitHeaderBytes;

    for (int i = nalUnitHeaderBytes; i < nal.num_bytes_in_nal_unit; ++i) {
        //starting from begin_bytepos to avoid header information
        //in NAL unit, 0x000000, 0x000001 or 0x000002 shall not occur at any byte-aligned position
        if (count == 2 && nal.rbsp_byte[i] < 0x03)
            return nal.num_bytes_in_rbsp = -1;
        if (count == 2 && nal.rbsp_byte[i] == 0x03) {
            //check the 4th byte after 0x000003, except when cabac_zero_word is used, in which case the last three bytes of this NAL unit must be 0x000003
            if (i < nal.num_bytes_in_nal_unit - 1 && nal.rbsp_byte[i + 1] > 0x03)
                return nal.num_bytes_in_rbsp = -1;
            //if cabac_zero_word is used, the final byte of this NAL unit(0x03) is discarded, and the last two bytes of RBSP must be 0x0000
            if (i == nal.num_bytes_in_nal_unit - 1)
                return nal.num_bytes_in_rbsp = j;

            ++i;
            count = 0;
        }
        nal.rbsp_byte[j] = nal.rbsp_byte[i];
        if (nal.rbsp_byte[i] == 0x00)
            ++count;
        else
            count = 0;
        ++j;
    }

    return nal.num_bytes_in_rbsp = j;
}


struct annex_b_t {
    static const int MAX_IOBUF_SIZE = 512 * 1024;

    int         BitStreamFile;
    bool        is_eof;
    size_t      iobuf_size;
    uint8_t*    iobuf_data;
    size_t      rdbuf_size;
    uint8_t*    rdbuf_data;

    int32_t     nextstartcodebytes;
    uint8_t*    Buf;

                annex_b_t(uint32_t max_size);
                ~annex_b_t();

    void        open (const char* fn);
    void        close();
    void        reset();

    annex_b_t& operator>>(nal_unit_t& nal);
    uint32_t    get_nalu(nal_unit_t& nal);

    inline uint32_t getChunk();
    inline uint8_t  getfbyte();
    inline bool     FindStartCode(uint8_t* Buf, uint32_t zeros_in_startcode);
};


annex_b_t::annex_b_t(uint32_t max_size)
{
    this->is_eof = false;
    this->iobuf_size = 0;
    this->iobuf_data = nullptr;
    this->Buf = new uint8_t[max_size];
}

annex_b_t::~annex_b_t()
{
    delete this->Buf;
}


void annex_b_t::open(const char* fn)
{
    if (this->iobuf_data)
        error(500, "open_annex_b: tried to open Annex B file twice");
    if ((this->BitStreamFile = ::open(fn, O_RDONLY)) == -1)
        error(500, "Cannot open Annex B ByteStream file '%s'", fn);

    this->is_eof = false;
    this->iobuf_size = annex_b_t::MAX_IOBUF_SIZE * sizeof(uint8_t);
    this->iobuf_data = new uint8_t[this->iobuf_size];
    this->getChunk();
}

void annex_b_t::close()
{
    if (this->BitStreamFile != -1) {
        ::close(this->BitStreamFile);
        this->BitStreamFile = -1;
    }
    delete this->iobuf_data;
    this->iobuf_data = NULL;
}

void annex_b_t::reset()
{
    this->is_eof     = false;
    this->rdbuf_size = 0;
    this->rdbuf_data = this->iobuf_data;
}


annex_b_t& annex_b_t::operator>>(nal_unit_t& nal)
{
    this->get_nalu(nal);
    return *this;
}

uint32_t annex_b_t::get_nalu(nal_unit_t& nal)
{
    uint32_t pos = 0;
    uint8_t* pBuf = this->Buf;
    //uint8_t* pBuf = nal.buf;

    if (this->nextstartcodebytes != 0) {
        for (int i = 0; i < this->nextstartcodebytes - 1; ++i) {
            *pBuf++ = 0;
            pos++;
        }
        *pBuf++ = 1;
        pos++;
    } else {
        while (!this->is_eof) {
            pos++;
            if ((*pBuf++ = this->getfbyte()) != 0)
                break;
        }
    }

    if (this->is_eof) {
        nal.num_bytes_in_nal_unit = pos == 0 ? 0 : -1;
        return nal.num_bytes_in_nal_unit;
    }

    if (*(pBuf - 1) != 1 || pos < 3) {
        nal.num_bytes_in_nal_unit = -1;
        return nal.num_bytes_in_nal_unit;
    }

    int LeadingZero8BitsCount = pos;
    bool info2 = false, info3 = false;

    //pBuf = nal.buf;

    bool StartCodeFound = false;
    while (!StartCodeFound) {
        if (this->is_eof) {
            pBuf -= 2;
            while (*pBuf-- == 0)
                pos--;

            nal.num_bytes_in_nal_unit = (pos - 1) - LeadingZero8BitsCount;
            memcpy(nal.rbsp_byte, this->Buf + LeadingZero8BitsCount, nal.num_bytes_in_nal_unit);
            this->nextstartcodebytes = 0;
            nal_unit(nal);
            return pos - 1;
        }

        pos++;
        *pBuf++ = this->getfbyte();
        info3 = this->FindStartCode(pBuf - 4, 3);
        if (!info3) {
            info2 = this->FindStartCode(pBuf - 3, 2);
            StartCodeFound = info2;
        } else
            StartCodeFound = true;
    }

    // Here, we have found another start code (and read length of startcode bytes more than we should
    // have.  Hence, go back in the file
    if (info3) { //if the detected start code is 00 00 01, trailing_zero_8bits is sure not to be present
        pBuf -= 5;
        while (*pBuf-- == 0)
            pos--;
        this->nextstartcodebytes = 4;
    } else if (info2)
        this->nextstartcodebytes = 3;
    else {
        printf(" Panic: Error in next start code search \n");
        return -1;
    }

    pos -= this->nextstartcodebytes;

    // Here the leading zeros(if any), Start code, the complete NALU, trailing zeros(if any)
    // and the next start code is in the Buf.
    // The size of Buf is pos - rewind, pos are the number of bytes excluding the next
    // start code, and (pos) - LeadingZero8BitsCount
    // is the size of the NALU.

    nal.num_bytes_in_nal_unit = pos - LeadingZero8BitsCount;
    memcpy(nal.rbsp_byte, this->Buf + LeadingZero8BitsCount, nal.num_bytes_in_nal_unit);
    nal.lost_packets = 0;
    nal_unit(nal);

    return pos;
}


inline uint32_t annex_b_t::getChunk()
{
    size_t reads = ::read(this->BitStreamFile, this->iobuf_data, this->iobuf_size);
    if (0 == reads) {
        this->is_eof = true;
        return 0;
    }

    this->rdbuf_size = reads;
    this->rdbuf_data = this->iobuf_data;
    return reads;
}

inline uint8_t annex_b_t::getfbyte()
{
    if (0 == this->rdbuf_size) {
        if (0 == this->getChunk())
            return 0;
    }

    --(this->rdbuf_size);
    return *this->rdbuf_data++;
}

inline bool annex_b_t::FindStartCode(uint8_t* Buf, uint32_t zeros_in_startcode)
{
    for (int i = 0; i < zeros_in_startcode; ++i) {
        if (*Buf++ != 0)
            return false;
    }
    return *Buf == 1;
}


void bitstream_t::open(const char* name, type format, uint32_t max_size)
{
    this->FileFormat = format;

    switch (format) {
    case type::RTP:
        open_rtp(name, &this->BitStreamFile);
        break;
    case type::ANNEX_B:
    default:
        this->annex_b = new annex_b_t(max_size);
        this->annex_b->open(name);
        break;
    }
}

void bitstream_t::close()
{
    switch (this->FileFormat) {
    case type::RTP:
        close_rtp(&this->BitStreamFile);
        break;
    case type::ANNEX_B:
    default:
        this->annex_b->close();
        delete this->annex_b;
        break;
    }
}


bitstream_t& bitstream_t::operator>>(nal_unit_t& nal)
{
    int ret;

    switch (this->FileFormat) {
    case type::RTP:
        ret = get_nalu_from_rtp(nal, this->BitStreamFile);
        break;
    case type::ANNEX_B:
    default:
        ret = this->annex_b->get_nalu(nal);
        break;
    }

    if (ret < 0) {
        error(601, "Error while getting the NALU in file format %s, exit\n",
                   this->FileFormat == type::ANNEX_B ? "Annex B" : "RTP");
    }
    if (ret == 0) {
        nal.num_bytes_in_rbsp = 0;
        return *this;
    }

    ret = NALUtoRBSP(nal);

    if (ret < 0)
        error(602, "Invalid startcode emulation prevention found.");

    // Got a NALU
    if (nal.forbidden_zero_bit)
        error(603, "Found NALU with forbidden_zero_bit set, bit error?");

    return *this;
}


}
}
