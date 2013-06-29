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


// Note that all NA values are filled with 0

/*!
 *************************************************************************************
 * \brief
 *    read_ue_v, reads an ue(v) syntax element, the length in bits is stored in
 *    the global p_Dec->UsedBits variable
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int read_ue_v(const char *tracestring, Bitstream *bitstream, int *used_bits)
{
    SyntaxElement symbol;
    symbol.type    = SE_HEADER;
    symbol.mapping = linfo_ue;   // Mapping rule

    readSyntaxElement_VLC(&symbol, bitstream);
    *used_bits += symbol.len;
    return symbol.value1;
}


/*!
 *************************************************************************************
 * \brief
 *    read_ue_v, reads an se(v) syntax element, the length in bits is stored in
 *    the global p_Dec->UsedBits variable
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int read_se_v(const char *tracestring, Bitstream *bitstream, int *used_bits)
{
    SyntaxElement symbol;
    symbol.type    = SE_HEADER;
    symbol.mapping = linfo_se;   // Mapping rule: signed integer

    readSyntaxElement_VLC(&symbol, bitstream);
    *used_bits += symbol.len;
    return symbol.value1;
}


/*!
 *************************************************************************************
 * \brief
 *    read_ue_v, reads an u(v) syntax element, the length in bits is stored in
 *    the global p_Dec->UsedBits variable
 *
 * \param LenInBits
 *    length of the syntax element
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int read_u_v(int LenInBits, const char *tracestring, Bitstream *bitstream, int *used_bits)
{
    SyntaxElement symbol;
    symbol.inf     = 0;
    symbol.type    = SE_HEADER;
    symbol.mapping = linfo_ue;   // Mapping rule
    symbol.len     = LenInBits;

    readSyntaxElement_FLC(&symbol, bitstream);
    *used_bits += symbol.len;
    return symbol.inf;
}

/*!
 *************************************************************************************
 * \brief
 *    read_i_v, reads an i(v) syntax element, the length in bits is stored in
 *    the global p_Dec->UsedBits variable
 *
 * \param LenInBits
 *    length of the syntax element
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int read_i_v(int LenInBits, const char *tracestring, Bitstream *bitstream, int *used_bits)
{
    SyntaxElement symbol;
    symbol.inf     = 0;
    symbol.type    = SE_HEADER;
    symbol.mapping = linfo_ue;   // Mapping rule
    symbol.len     = LenInBits;

    readSyntaxElement_FLC(&symbol, bitstream);
    *used_bits += symbol.len;
    // can be negative
    symbol.inf = -( symbol.inf & (1 << (LenInBits - 1)) ) | symbol.inf;
    return symbol.inf;
}


/*!
 *************************************************************************************
 * \brief
 *    read_ue_v, reads an u(1) syntax element, the length in bits is stored in
 *    the global p_Dec->UsedBits variable
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
Boolean read_u_1 (const char *tracestring, Bitstream *bitstream, int *used_bits)
{
  return (Boolean) read_u_v (1, tracestring, bitstream, used_bits);
}


/*!
 ************************************************************************
 * \brief
 *    mapping rule for ue(v) syntax elements
 * \par Input:
 *    lenght and info
 * \par Output:
 *    number in the code table
 ************************************************************************
 */
void linfo_ue(int len, int info, int *value1, int *dummy)
{
  //assert ((len >> 1) < 32);
  *value1 = (int) (((unsigned int) 1 << (len >> 1)) + (unsigned int) (info) - 1);
}

/*!
 ************************************************************************
 * \brief
 *    mapping rule for se(v) syntax elements
 * \par Input:
 *    lenght and info
 * \par Output:
 *    signed mvd
 ************************************************************************
 */
void linfo_se(int len,  int info, int *value1, int *dummy)
{
  //assert ((len >> 1) < 32);
  unsigned int n = ((unsigned int) 1 << (len >> 1)) + (unsigned int) info - 1;
  *value1 = (n + 1) >> 1;
  if((n & 0x01) == 0)                           // lsb is signed bit
    *value1 = -*value1;
}


/*!
 ************************************************************************
 * \brief
 *  read one exp-golomb VLC symbol
 *
 * \param buffer
 *    containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param  info
 *    returns the value of the symbol
 * \param bytecount
 *    buffer length
 * \return
 *    bits read
 ************************************************************************
 */
static int GetVLCSymbol(byte buffer[], int totbitoffset, int *info, int bytecount)
{
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
        return bitcounter; // return absolute offset in bit from start of frame
    }
}

/*!
 ************************************************************************
 * \brief
 *    read next UVLC codeword from UVLC-partition and
 *    map it to the corresponding syntax element
 ************************************************************************
 */
int readSyntaxElement_VLC(SyntaxElement *sym, Bitstream *currStream)
{
    sym->len =  GetVLCSymbol (currStream->streamBuffer, currStream->frame_bitoffset, &(sym->inf), currStream->bitstream_length);
    if (sym->len == -1)
        return -1;

    currStream->frame_bitoffset += sym->len;
    sym->mapping(sym->len, sym->inf, &(sym->value1), &(sym->value2));

    return 1;
}


/*!
 ************************************************************************
 * \brief
 *    read next UVLC codeword from UVLC-partition and
 *    map it to the corresponding syntax element
 ************************************************************************
 */
int readSyntaxElement_UVLC(Macroblock *currMB, SyntaxElement *sym, struct datapartition_dec *dP)
{
    return readSyntaxElement_VLC(sym, dP->bitstream);
}


/*!
 ************************************************************************
 * \brief
 *    test if bit buffer contains only stop bit
 *
 * \param buffer
 *    buffer containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param bytecount
 *    buffer length
 * \return
 *    true if more bits available
 ************************************************************************
 */
int more_rbsp_data (byte buffer[],int totbitoffset,int bytecount)
{
  long byteoffset = (totbitoffset >> 3);      // byte from start of buffer
  // there is more until we're in the last byte
  if (byteoffset < (bytecount - 1)) 
    return TRUE;
  else
  {
    int bitoffset   = (7 - (totbitoffset & 0x07));      // bit from start of byte
    byte *cur_byte  = &(buffer[byteoffset]);
    // read one bit
    int ctr_bit     = ctr_bit = ((*cur_byte)>> (bitoffset--)) & 0x01;      // control bit for current bit posision

    //assert (byteoffset<bytecount);       

    // a stop bit has to be one
    if (ctr_bit==0) 
      return TRUE;  
    else
    {
      int cnt = 0;

      while (bitoffset>=0 && !cnt)
      {
        cnt |= ((*cur_byte)>> (bitoffset--)) & 0x01;   // set up control bit
      }

      return (cnt);
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    read FLC codeword from UVLC-partition
 ************************************************************************
 */
int readSyntaxElement_FLC(SyntaxElement *sym, Bitstream *currStream)
{
  int BitstreamLengthInBits  = (currStream->bitstream_length << 3) + 7;
  
  if ((GetBits(currStream->streamBuffer, currStream->frame_bitoffset, &(sym->inf), BitstreamLengthInBits, sym->len)) < 0)
    return -1;

  sym->value1 = sym->inf;
  currStream->frame_bitoffset += sym->len; // move bitstream pointer

  return 1;
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

/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer
 *
 * \param buffer
 *    buffer containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param bitcount
 *    total bytes in bitstream
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */

int ShowBits(byte buffer[], int totbitoffset, int bitcount, int numbits)
{
    if (totbitoffset + numbits > bitcount)
        return -1;
    else {
        int bitoffset  = 7 - (totbitoffset & 0x07); // bit from start of byte
        int byteoffset = (totbitoffset >> 3); // byte from start of buffer
        byte *curbyte  = &(buffer[byteoffset]);
        int inf        = 0;

        while (numbits--) {
            inf <<= 1;
            inf |= ((*curbyte)>> (bitoffset--)) & 0x01;

            if (bitoffset == -1) { //Move onto next byte to get all of numbits
                curbyte++;
                bitoffset = 7;
            }
        }
        return inf; // return absolute offset in bit from start of frame
    }
}
