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
#include "data_partition.h"


namespace arrow {
namespace video {
namespace h264  {


struct annex_b_t {
    static const int IOBUFFERSIZE = 512 * 1024;

    int32_t     BitStreamFile;
    uint8_t*    iobuffer;
    uint8_t*    iobufferread;
    int32_t     bytesinbuffer;
    bool        is_eof;
    int32_t     iIOBufferSize;

    int32_t     IsFirstByteStreamNALU;
    int32_t     nextstartcodebytes;
    uint8_t*    Buf;  

                annex_b_t(uint32_t max_size);
                ~annex_b_t();

    void        open (const char* fn);
    void        close();
    void        reset();

    uint32_t    get_nalu(nalu_t* nalu);

    inline uint32_t getChunk();
    inline uint8_t  getfbyte();
    inline bool     FindStartCode(uint8_t* Buf, uint32_t zeros_in_startcode);
};


annex_b_t::annex_b_t(uint32_t max_size)
{
    this->Buf = new uint8_t[max_size];
    if (!this->Buf)
        error("malloc_annex_b: Buf", 101);
}

annex_b_t::~annex_b_t()
{
    delete this->Buf;
}


void annex_b_t::open(const char* fn)
{
    if (this->iobuffer)
        error("open_annex_b: tried to open Annex B file twice", 500);
    if ((this->BitStreamFile = ::open(fn, O_RDONLY)) == -1) {
        snprintf(errortext, ET_SIZE, "Cannot open Annex B ByteStream file '%s'", fn);
        error(errortext, 500);
    }

    this->iIOBufferSize = annex_b_t::IOBUFFERSIZE * sizeof(uint8_t);
    this->iobuffer = new uint8_t[this->iIOBufferSize];
    if (!this->iobuffer)
        error("open_annex_b: cannot allocate IO buffer", 500);
    this->is_eof = false;
    this->getChunk();
}

void annex_b_t::close()
{
    if (this->BitStreamFile != -1) {
        ::close(this->BitStreamFile);
        this->BitStreamFile = -1;
    }
    delete this->iobuffer;
    this->iobuffer = NULL;
}

void annex_b_t::reset()
{
    this->is_eof        = false;
    this->bytesinbuffer = 0;
    this->iobufferread  = this->iobuffer;
}


uint32_t annex_b_t::get_nalu(nalu_t* nalu)
{
    bool info2 = false, info3 = false;
    bool StartCodeFound = false;
    int  LeadingZero8BitsCount = 0;
    uint32_t pos = 0;
    uint8_t* pBuf = this->Buf;

    if (this->nextstartcodebytes != 0) {
        for (int i = 0; i < this->nextstartcodebytes - 1; i++) {
            (*pBuf++) = 0;
            pos++;
        }
        (*pBuf++) = 1;
        pos++;
    } else {
        while (!this->is_eof) {
            pos++;
            if ((*(pBuf++) = this->getfbyte()) != 0)
                break;
        }
    }

    if (this->is_eof) {
        if (pos == 0)
            return 0;
        else {
            printf( "get_annex_b_NALU can't read start code\n");
            return -1;
        }
    }

    if (*(pBuf - 1) != 1 || pos < 3) {
        printf("get_annex_b_NALU: no Start Code at the beginning of the NALU, return -1\n");
        return -1;
    }

    if (pos == 3)
        nalu->startcodeprefix_len = 3;
    else {
        LeadingZero8BitsCount = pos - 4;
        nalu->startcodeprefix_len = 4;
    }

    //the 1st byte stream NAL unit can has leading_zero_8bits, but subsequent ones are not
    //allowed to contain it since these zeros(if any) are considered trailing_zero_8bits
    //of the previous byte stream NAL unit.
    if (!this->IsFirstByteStreamNALU && LeadingZero8BitsCount > 0) {
        printf("get_annex_b_NALU: The leading_zero_8bits syntax can only be present in the first byte stream NAL unit, return -1\n");
        return -1;
    }

    LeadingZero8BitsCount = pos;
    this->IsFirstByteStreamNALU = 0;

    while (!StartCodeFound) {
        if (this->is_eof) {
            pBuf -= 2;
            while (*(pBuf--) == 0)
                pos--;

            nalu->len = (pos - 1) - LeadingZero8BitsCount;
            memcpy(nalu->buf, this->Buf + LeadingZero8BitsCount, nalu->len);
            nalu->forbidden_bit = (*(nalu->buf) >> 7) & 1;
            nalu->nal_ref_idc   = (*(nalu->buf) >> 5) & 3;
            nalu->nal_unit_type = (NaluType) ((*(nalu->buf)) & 0x1f);
            this->nextstartcodebytes = 0;

            return pos - 1;
        }

        pos++;
        *(pBuf++) = this->getfbyte();    
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
        while (*(pBuf--) == 0)
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

    nalu->len = pos - LeadingZero8BitsCount;
    memcpy(nalu->buf, this->Buf + LeadingZero8BitsCount, nalu->len);
    nalu->forbidden_bit = (*(nalu->buf) >> 7) & 1;
    nalu->nal_ref_idc   = (*(nalu->buf) >> 5) & 3;
    nalu->nal_unit_type = (NaluType) ((*(nalu->buf)) & 0x1f);
    nalu->lost_packets  = 0;

    return pos;
}


inline uint32_t annex_b_t::getChunk()
{
    uint32_t readbytes = ::read(this->BitStreamFile, this->iobuffer, this->iIOBufferSize); 
    if (0 == readbytes) {
        this->is_eof = true;
        return 0;
    }

    this->bytesinbuffer = readbytes;
    this->iobufferread = this->iobuffer;
    return readbytes;
}

inline uint8_t annex_b_t::getfbyte()
{
    if (0 == this->bytesinbuffer) {
        if (0 == this->getChunk())
            return 0;
    }

    this->bytesinbuffer--;
    return *this->iobufferread++;
}

inline bool annex_b_t::FindStartCode(uint8_t* Buf, uint32_t zeros_in_startcode)
{
    for (int i = 0; i < zeros_in_startcode; i++) {
        if (*(Buf++) != 0)
            return false;
    }

    if (*Buf != 1)
        return false;

    return true;
}


void bitstream_t::open(const char* name, type format, uint32_t max_size)
{
    this->FileFormat           = format;
    this->LastAccessUnitExists = 0;
    this->NALUCount            = 0;

    switch (format) {
    case type::RTP:
        open_rtp(name, &this->BitStreamFile);
        break;
    case type::ANNEX_B:
    default:
        this->annex_b = new annex_b_t(max_size);
        if (!this->annex_b) {
            snprintf(errortext, ET_SIZE, "Memory allocation for Annex_B file failed");
            error(errortext,100);
        }
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

void bitstream_t::reset()
{
    if (this->FileFormat == type::ANNEX_B)
        this->annex_b->reset(); 
}


static int EBSPtoRBSP(byte *streamBuffer, int end_bytepos, int begin_bytepos)
{
    //Start code and Emulation Prevention need this to be defined in identical manner at encoder and decoder
    #define ZEROBYTES_SHORTSTARTCODE 2 //indicates the number of zero bytes in the short start-code prefix

    int i, j, count;
    count = 0;

    if (end_bytepos < begin_bytepos)
        return end_bytepos;

    j = begin_bytepos;

    for (i = begin_bytepos; i < end_bytepos; ++i) {
        //starting from begin_bytepos to avoid header information
        //in NAL unit, 0x000000, 0x000001 or 0x000002 shall not occur at any byte-aligned position
        if (count == ZEROBYTES_SHORTSTARTCODE && streamBuffer[i] < 0x03) 
            return -1;
        if (count == ZEROBYTES_SHORTSTARTCODE && streamBuffer[i] == 0x03) {
            //check the 4th byte after 0x000003, except when cabac_zero_word is used, in which case the last three bytes of this NAL unit must be 0x000003
            if ((i < end_bytepos - 1) && (streamBuffer[i+1] > 0x03))
                return -1;
            //if cabac_zero_word is used, the final byte of this NAL unit(0x03) is discarded, and the last two bytes of RBSP must be 0x0000
            if (i == end_bytepos - 1)
                return j;

            ++i;
            count = 0;
        }
        streamBuffer[j] = streamBuffer[i];
        if (streamBuffer[i] == 0x00)
            ++count;
        else
            count = 0;
        ++j;
    }

    return j;
}

static int NALUtoRBSP(nalu_t *nalu)
{
    assert(nalu != NULL);

    nalu->len = EBSPtoRBSP(nalu->buf, nalu->len, 1);

    return nalu->len;
}

int bitstream_t::read_next_nalu(nalu_t* nalu)
{
    int ret;

    switch (this->FileFormat) {
    case type::RTP:
        ret = get_nalu_from_rtp(nalu, this->BitStreamFile);
        break;   
    case type::ANNEX_B:
    default:
        ret = this->annex_b->get_nalu(nalu);
        break;
    }

    if (ret < 0) {
        snprintf(errortext, ET_SIZE, "Error while getting the NALU in file format %s, exit\n",
                 this->FileFormat == type::ANNEX_B ? "Annex B" : "RTP");
        error(errortext, 601);
    }
    if (ret == 0)
        return 0;

    //In some cases, zero_byte shall be present. If current NALU is a VCL NALU, we can't tell
    //whether it is the first VCL NALU at this point, so only non-VCL NAL unit is checked here.
    this->CheckZeroByteNonVCL(nalu);

    ret = NALUtoRBSP(nalu);

    if (ret < 0)
        error("Invalid startcode emulation prevention found.", 602);

    // Got a NALU
    if (nalu->forbidden_bit)
        error("Found NALU with forbidden_bit set, bit error?", 603);

    return nalu->len;
}

void bitstream_t::CheckZeroByteNonVCL(nalu_t* nalu)
{
    int CheckZeroByte = 0;

    //This function deals only with non-VCL NAL units
    if (nalu->nal_unit_type >= 1 && nalu->nal_unit_type <= 5)
        return;

    //for SPS and PPS, zero_byte shall exist
    if (nalu->nal_unit_type == NALU_TYPE_SPS || nalu->nal_unit_type == NALU_TYPE_PPS)
        CheckZeroByte = 1;
    //check the possibility of the current NALU to be the start of a new access unit, according to 7.4.1.2.3
    if (nalu->nal_unit_type == NALU_TYPE_AUD || nalu->nal_unit_type == NALU_TYPE_SPS ||
        nalu->nal_unit_type == NALU_TYPE_PPS || nalu->nal_unit_type == NALU_TYPE_SEI ||
        (nalu->nal_unit_type >= 13 && nalu->nal_unit_type <= 18)) {
        if (this->LastAccessUnitExists) {
            this->LastAccessUnitExists = 0; //deliver the last access unit to decoder
            this->NALUCount = 0;
        }
    }
    this->NALUCount++;
    //for the first NAL unit in an access unit, zero_byte shall exists
    if (this->NALUCount == 1)
        CheckZeroByte = 1;
    if (CheckZeroByte && nalu->startcodeprefix_len == 3)
        printf("Warning: zero_byte shall exist\n");
        //because it is not a very serious problem, we do not exit here
}

void bitstream_t::CheckZeroByteVCL(nalu_t* nalu)
{
    int CheckZeroByte = 0;

    //This function deals only with VCL NAL units
    if (!(nalu->nal_unit_type >= NALU_TYPE_SLICE &&
          nalu->nal_unit_type <= NALU_TYPE_IDR))
        return;

    if (this->LastAccessUnitExists)
        this->NALUCount = 0;
    this->NALUCount++;
    //the first VCL NAL unit that is the first NAL unit after last VCL NAL unit indicates
    //the start of a new access unit and hence the first NAL unit of the new access unit.           (sounds like a tongue twister :-)
    if (this->NALUCount == 1)
        CheckZeroByte = 1;
    this->LastAccessUnitExists = 1;
    if (CheckZeroByte && nalu->startcodeprefix_len == 3)
        printf("warning: zero_byte shall exist\n");
        //because it is not a very serious problem, we do not exit here
}


}
}
}
