/*!
 ************************************************************************
 * \file vlc.c
 *
 * \brief
 *    VLC support functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Detlev Marpe
 *    - Gabi Blaettermann
 ************************************************************************
 */

#include "global.h"
#include "macroblock.h"
#include "bitstream.h"
#include "bitstream_elements.h"


bool bit_stream_dec::more_rbsp_data(void)
{
    byte *buffer       = this->streamBuffer;
    int   totbitoffset = this->frame_bitoffset;
    int   bytecount    = this->bitstream_length;

    long byteoffset = (totbitoffset >> 3);      // byte from start of buffer
    // there is more until we're in the last byte
    if (byteoffset < (bytecount - 1)) 
        return true;
    else {
        int bitoffset   = (7 - (totbitoffset & 0x07));      // bit from start of byte
        byte *cur_byte  = &(buffer[byteoffset]);
        // read one bit
        int ctr_bit     = ctr_bit = ((*cur_byte)>> (bitoffset--)) & 0x01;      // control bit for current bit posision

        //assert (byteoffset<bytecount);       

        // a stop bit has to be one
        if (ctr_bit == 0)
            return true;
        else {
            int cnt = 0;

            while (bitoffset >= 0 && !cnt)
                cnt |= ((*cur_byte) >> (bitoffset--)) & 0x01;   // set up control bit

            return cnt ? true : false;
        }
    }
}

uint32_t bit_stream_dec::next_bits(uint8_t n)
{
    byte *buffer       = this->streamBuffer;
    int   totbitoffset = this->frame_bitoffset;
    int   bitcount     = (this->bitstream_length << 3) + 7;

    if (totbitoffset + n > bitcount) 
        return 0;

    int bitoffset  = 7 - (totbitoffset & 0x07); // bit from start of byte
    int byteoffset = (totbitoffset >> 3); // byte from start of buffer
    byte *curbyte  = &(buffer[byteoffset]);
    uint32_t inf = 0;

    while (n--) {
        inf <<= 1;
        inf |= ((*curbyte) >> (bitoffset--)) & 0x01;    
        if (bitoffset == -1) { //Move onto next byte to get all of numbits
            curbyte++;
            bitoffset = 7;
        }
        // Above conditional could also be avoided using the following:
        // curbyte   -= (bitoffset >> 3);
        // bitoffset &= 0x07;
    }

    return inf;
}

uint32_t bit_stream_dec::read_bits(uint8_t n)
{
    uint32_t value = next_bits(n);
    this->frame_bitoffset += n;
    return value;
}

uint32_t bit_stream_dec::u(uint8_t n, const char *name)
{
    return read_bits(n);
}

int32_t bit_stream_dec::i(uint8_t n, const char *name)
{
    uint32_t value = read_bits(n);
    return -(value & (1 << (n - 1))) | value;
}

uint32_t bit_stream_dec::f(uint8_t n, const char *name)
{
    return read_bits(n);
}

uint32_t bit_stream_dec::b(uint8_t n, const char *name)
{
    return read_bits(n);
}

uint32_t bit_stream_dec::ue(const char *name)
{
    int leadingZeroBits = -1;
    uint32_t b;
    uint32_t codeNum;

    for (b = 0; !b; leadingZeroBits++)
        b = read_bits(1);

    codeNum = (1 << leadingZeroBits) - 1 + read_bits(leadingZeroBits);
    return codeNum;
}

int32_t bit_stream_dec::se(const char *name)
{
    uint32_t codeNum = ue();
    return (codeNum % 2 ? 1 : -1) * ((codeNum + 1) / 2);
}

uint32_t bit_stream_dec::ae(const char *name)
{
    return 0;
}

uint32_t bit_stream_dec::ce(const char *name)
{
    return 0;
}

uint32_t bit_stream_dec::me(const char *name)
{
    return 0;
}

uint32_t bit_stream_dec::te(const char *name)
{
    return 0;
}





static int GetVLCSymbol(int *info, Bitstream *currStream)
{
    byte *buffer       = currStream->streamBuffer;
    int   totbitoffset = currStream->frame_bitoffset;
    int   bytecount    = currStream->bitstream_length;

    long byteoffset = (totbitoffset >> 3);         // byte from start of buffer
    int  bitoffset  = (7 - (totbitoffset & 0x07)); // bit from start of byte
    int  bitcounter = 1;
    int  len        = 0;
    byte *cur_byte  = &(buffer[byteoffset]);
    int  ctr_bit    = ((*cur_byte) >> (bitoffset)) & 0x01;  // control bit for current bit posision

    while (ctr_bit == 0) {                 // find leading 1 bit
        len++;
        bitcounter++;
        bitoffset--;
        bitoffset &= 0x07;
        cur_byte  += (bitoffset == 7);
        byteoffset+= (bitoffset == 7);      
        ctr_bit    = ((*cur_byte) >> (bitoffset)) & 0x01;
    }

    if (byteoffset + ((len + 7) >> 3) > bytecount)
        return -1;
    else {
        // make infoword
        int inf = 0;  // shortest possible code is 1, then info is always 0    

        while (len--) {
            bitoffset--;
            bitoffset &= 0x07;
            cur_byte  += (bitoffset == 7);
            bitcounter++;
            inf <<= 1;
            inf |= ((*cur_byte) >> (bitoffset)) & 0x01;
        }

        *info = inf;
        currStream->frame_bitoffset += bitcounter;
        return bitcounter; // return absolute offset in bit from start of frame
    }
}



void linfo_ue(int len, int info, int *value1, int *dummy)
{
  //assert ((len >> 1) < 32);
  *value1 = (int) (((unsigned int) 1 << (len >> 1)) + (unsigned int) (info) - 1);
}

void linfo_se(int len,  int info, int *value1, int *dummy)
{
  //assert ((len >> 1) < 32);
  unsigned int n = ((unsigned int) 1 << (len >> 1)) + (unsigned int) info - 1;
  *value1 = (n + 1) >> 1;
  if((n & 0x01) == 0)                           // lsb is signed bit
    *value1 = -*value1;
}


int readSyntaxElement_UVLC(mb_t *currMB, SyntaxElement *sym, struct datapartition_dec *dP)
{
    sym->len = GetVLCSymbol(&sym->inf, dP->bitstream);
    if (sym->len == -1)
        return -1;

    sym->mapping(sym->len, sym->inf, &(sym->value1), &(sym->value2));
    return 1;
}

int readSyntaxElement_CABAC(mb_t *currMB, SyntaxElement *se, struct datapartition_dec *dP)
{
    cabac_engine_t *dep_dp = &dP->bitstream->de_cabac;
    se->reading(currMB, se, dep_dp);
    return 0;
}

/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer
 *
 * \param buffer
 *    containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param info
 *    returns value of the read bits
 * \param bitcount
 *    total bytes in bitstream
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */
int GetBits(byte buffer[], int totbitoffset, int *info, int bitcount, int numbits)
{
    if (totbitoffset + numbits > bitcount) 
        return -1;
    else {
        int bitoffset  = 7 - (totbitoffset & 0x07); // bit from start of byte
        int byteoffset = (totbitoffset >> 3); // byte from start of buffer
        int bitcounter = numbits;
        byte *curbyte  = &(buffer[byteoffset]);
        int inf = 0;

        while (numbits--) {
            inf <<= 1;
            inf |= ((*curbyte) >> (bitoffset--)) & 0x01;    
            if (bitoffset == -1) { //Move onto next byte to get all of numbits
                curbyte++;
                bitoffset = 7;
            }
            // Above conditional could also be avoided using the following:
            // curbyte   -= (bitoffset >> 3);
            // bitoffset &= 0x07;
        }
        *info = inf;

        return bitcounter; // return absolute offset in bit from start of frame
    }
}
