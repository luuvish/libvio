/*!
 ***********************************************************************
 * \file read_comp_cavlc.c
 *
 * \brief
 *     Read Coefficient Components (CAVLC version)
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

#include "global.h"
#include "slice.h"
#include "bitstream_elements.h"
#include "bitstream.h"
#include "macroblock.h"
#include "mb_read.h"
#include "transform.h"
#include "neighbour.h"

#define IS_I16MB(MB)    ((MB)->mb_type==I16MB  || (MB)->mb_type==IPCM)
#define IS_DIRECT(MB)   ((MB)->mb_type==0     && (currSlice->slice_type == B_SLICE ))

#define TOTRUN_NUM       15
#define RUNBEFORE_NUM_M1  6


/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer (Threshold based)
 *
 * \param inf
 *    bytes to extract numbits from with bitoffset already applied
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */

static inline int ShowBitsThres (int inf, int numbits)
{
  return ((inf) >> ((sizeof(byte) * 24) - (numbits)));
}


/*!
 ************************************************************************
 * \brief
 *    code from bitstream (2d tables)
 ************************************************************************
 */

static int code_from_bitstream_2d(SyntaxElement *sym,
                                  Bitstream *currStream,
                                  const byte *lentab,
                                  const byte *codtab,
                                  int tabwidth,
                                  int tabheight,
                                  int *code)
{
  int i, j;
  const byte *len = &lentab[0], *cod = &codtab[0];

  int *frame_bitoffset = &currStream->frame_bitoffset;
  byte *buf            = &currStream->streamBuffer[*frame_bitoffset >> 3];

  //Apply bitoffset to three bytes (maximum that may be traversed by ShowBitsThres)
  unsigned int inf = ((*buf) << 16) + (*(buf + 1) << 8) + *(buf + 2); //Even at the end of a stream we will still be pulling out of allocated memory as alloc is done by MAX_CODED_FRAME_SIZE
  inf <<= (*frame_bitoffset & 0x07);                                  //Offset is constant so apply before extracting different numbers of bits
  inf  &= 0xFFFFFF;                                                   //Arithmetic shift so wipe any sign which may be extended inside ShowBitsThres
  
  // this VLC decoding method is not optimized for speed
  for (j = 0; j < tabheight; j++) 
  {
    for (i = 0; i < tabwidth; i++)
    {
      if ((*len == 0) || (ShowBitsThres(inf, (int) *len) != *cod))
      {
        ++len;
        ++cod;
      }
      else
      {
        sym->len = *len;
        *frame_bitoffset += *len; // move bitstream pointer
        *code = *cod;             
        sym->value1 = i;
        sym->value2 = j;        
        return 0;                 // found code and return 
      }
    }
  }
  return -1;  // failed to find code
}


/*!
 ************************************************************************
 * \brief
 *    read NumCoeff/TrailingOnes codeword from UVLC-partition
 ************************************************************************
 */

static int readSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *sym,  
                                           Bitstream *currStream,
                                           char *type)
{
  int frame_bitoffset        = currStream->frame_bitoffset;
  int BitstreamLengthInBytes = currStream->bitstream_length;
  int BitstreamLengthInBits  = (BitstreamLengthInBytes << 3) + 7;
  byte *buf                  = currStream->streamBuffer;

  static const byte lentab[3][4][17] =
  {
    {   // 0702
      { 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
      { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
      { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
      { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16},
    },
    {
      { 2, 6, 6, 7, 8, 8, 9,11,11,12,12,12,13,13,13,14,14},
      { 0, 2, 5, 6, 6, 7, 8, 9,11,11,12,12,13,13,14,14,14},
      { 0, 0, 3, 6, 6, 7, 8, 9,11,11,12,12,13,13,13,14,14},
      { 0, 0, 0, 4, 4, 5, 6, 6, 7, 9,11,11,12,13,13,13,14},
    },
    {
      { 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9,10,10,10,10},
      { 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10},
      { 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,10},
      { 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9,10,10,10},
    },
  };

  static const byte codtab[3][4][17] =
  {
    {
      { 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7,4},
      { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10,6},
      { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9,5},
      { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12,8},
    },
    {
      { 3,11, 7, 7, 7, 4, 7,15,11,15,11, 8,15,11, 7, 9,7},
      { 0, 2, 7,10, 6, 6, 6, 6,14,10,14,10,14,10,11, 8,6},
      { 0, 0, 3, 9, 5, 5, 5, 5,13, 9,13, 9,13, 9, 6,10,5},
      { 0, 0, 0, 5, 4, 6, 8, 4, 4, 4,12, 8,12,12, 8, 1,4},
    },
    {
      {15,15,11, 8,15,11, 9, 8,15,11,15,11, 8,13, 9, 5,1},
      { 0,14,15,12,10, 8,14,10,14,14,10,14,10, 7,12, 8,4},
      { 0, 0,13,14,11, 9,13, 9,13,10,13, 9,13, 9,11, 7,3},
      { 0, 0, 0,12,11,10, 9, 8,13,12,12,12, 8,12,10, 6,2},
    },
  };

  int retval = 0, code;
  int vlcnum = sym->value1;
  // vlcnum is the index of Table used to code coeff_token
  // vlcnum==3 means (8<=nC) which uses 6bit FLC

  if (vlcnum == 3)
  {
    // read 6 bit FLC
    code = ShowBits(buf, frame_bitoffset, BitstreamLengthInBits, 6);
    currStream->frame_bitoffset += 6;
    sym->value2 = (code & 3);
    sym->value1 = (code >> 2);

    if (!sym->value1 && sym->value2 == 3)
    {
      // #c = 0, #t1 = 3 =>  #c = 0
      sym->value2 = 0;
    }
    else
      sym->value1++;

    sym->len = 6;
  }
  else
  {
    retval = code_from_bitstream_2d(sym, currStream, lentab[vlcnum][0], codtab[vlcnum][0], 17, 4, &code);
    if (retval)
    {
      printf("ERROR: failed to find NumCoeff/TrailingOnes\n");
      exit(-1);
    }
  }

  return retval;
}


/*!
 ************************************************************************
 * \brief
 *    read NumCoeff/TrailingOnes codeword from UVLC-partition ChromaDC
 ************************************************************************
 */
static int readSyntaxElement_NumCoeffTrailingOnesChromaDC(VideoParameters *p_Vid, SyntaxElement *sym,  Bitstream *currStream)
{
  static const byte lentab[3][4][17] =
  {
    //YUV420
    {{ 2, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 1, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 3, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV422
    {{ 1, 7, 7, 9, 9,10,11,12,13, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 2, 7, 7, 9,10,11,12,12, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 3, 7, 7, 9,10,11,12, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 5, 6, 7, 7,10,11, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV444
    {{ 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
    { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
    { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
    { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16}}
  };

  static const byte codtab[3][4][17] =
  {
    //YUV420
    {{ 1, 7, 4, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 1, 6, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV422
    {{ 1,15,14, 7, 6, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 1,13,12, 5, 6, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 1,11,10, 4, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 1, 1, 9, 8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV444
    {{ 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7, 4},
    { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10, 6},
    { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9, 5},
    { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12, 8}}

  };

  int code;
  int yuv = p_Vid->active_sps->chroma_format_idc - 1;
  int retval = code_from_bitstream_2d(sym, currStream, &lentab[yuv][0][0], &codtab[yuv][0][0], 17, 4, &code);

  if (retval)
  {
    printf("ERROR: failed to find NumCoeff/TrailingOnes ChromaDC\n");
    exit(-1);
  }

  return retval;
}

/*!
 ************************************************************************
 * \brief
 *    read Level VLC0 codeword from UVLC-partition
 ************************************************************************
 */
static int readSyntaxElement_Level_VLC0(SyntaxElement *sym, Bitstream *currStream)
{
  int frame_bitoffset        = currStream->frame_bitoffset;
  int BitstreamLengthInBytes = currStream->bitstream_length;
  int BitstreamLengthInBits  = (BitstreamLengthInBytes << 3) + 7;
  byte *buf                  = currStream->streamBuffer;
  int len = 1, sign = 0, level = 0, code = 1;

  while (!ShowBits(buf, frame_bitoffset++, BitstreamLengthInBits, 1))
    len++;

  if (len < 15)
  {
    sign  = (len - 1) & 1;
    level = ((len - 1) >> 1) + 1;
  }
  else if (len == 15)
  {
    // escape code
    code <<= 4;
    code |= ShowBits(buf, frame_bitoffset, BitstreamLengthInBits, 4);
    len  += 4;
    frame_bitoffset += 4;
    sign = (code & 0x01);
    level = ((code >> 1) & 0x07) + 8;
  }
  else if (len >= 16)
  {
    // escape code
    int addbit = (len - 16);
    int offset = (2048 << addbit) - 2032;
    len   -= 4;
    code   = ShowBits(buf, frame_bitoffset, BitstreamLengthInBits, len);
    sign   = (code & 0x01);
    frame_bitoffset += len;    
    level = (code >> 1) + offset;

    code |= (1 << (len)); // for display purpose only
    len += addbit + 16;
 }

  sym->inf = (sign) ? -level : level ;
  sym->len = len;

  currStream->frame_bitoffset = frame_bitoffset;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    read Level VLC codeword from UVLC-partition
 ************************************************************************
 */
static int readSyntaxElement_Level_VLCN(SyntaxElement *sym, int vlc, Bitstream *currStream)
{
  int frame_bitoffset        = currStream->frame_bitoffset;
  int BitstreamLengthInBytes = currStream->bitstream_length;
  int BitstreamLengthInBits  = (BitstreamLengthInBytes << 3) + 7;
  byte *buf                  = currStream->streamBuffer;

  int levabs, sign;
  int len = 1;
  int code = 1, sb;

  int shift = vlc - 1;

  // read pre zeros
  while (!ShowBits(buf, frame_bitoffset ++, BitstreamLengthInBits, 1))
    len++;

  frame_bitoffset -= len;

  if (len < 16)
  {
    levabs = ((len - 1) << shift) + 1;

    // read (vlc-1) bits -> suffix
    if (shift)
    {
      sb =  ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, shift);
      code = (code << (shift) )| sb;
      levabs += sb;
      len += (shift);
    }

    // read 1 bit -> sign
    sign = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, 1);
    code = (code << 1)| sign;
    len ++;
  }
  else // escape
  {
    int addbit = len - 5;
    int offset = (1 << addbit) + (15 << shift) - 2047;

    sb = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, addbit);
    code = (code << addbit ) | sb;
    len   += addbit;

    levabs = sb + offset;
    
    // read 1 bit -> sign
    sign = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, 1);

    code = (code << 1)| sign;

    len++;
  }

  sym->inf = (sign)? -levabs : levabs;
  sym->len = len;

  currStream->frame_bitoffset = frame_bitoffset + len;

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    read Total Zeros codeword from UVLC-partition
 ************************************************************************
 */
static int readSyntaxElement_TotalZeros(SyntaxElement *sym,  Bitstream *currStream)
{
  static const byte lentab[TOTRUN_NUM][16] =
  {

    { 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},
    { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},
    { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},
    { 5,3,4,4,3,3,3,4,3,4,5,5,5},
    { 4,4,4,3,3,3,3,3,4,5,4,5},
    { 6,5,3,3,3,3,3,3,4,3,6},
    { 6,5,3,3,3,2,3,4,3,6},
    { 6,4,5,3,2,2,3,3,6},
    { 6,6,4,2,2,3,2,5},
    { 5,5,3,2,2,2,4},
    { 4,4,3,3,1,3},
    { 4,4,2,1,3},
    { 3,3,1,2},
    { 2,2,1},
    { 1,1},
  };

  static const byte codtab[TOTRUN_NUM][16] =
  {
    {1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
    {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
    {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
    {3,7,5,4,6,5,4,3,3,2,2,1,0},
    {5,4,3,7,6,5,4,3,2,1,1,0},
    {1,1,7,6,5,4,3,2,1,1,0},
    {1,1,5,4,3,3,2,1,1,0},
    {1,1,1,3,3,2,2,1,0},
    {1,0,1,3,2,1,1,1,},
    {1,0,1,3,2,1,1,},
    {0,1,1,2,1,3},
    {0,1,1,1,1},
    {0,1,1,1},
    {0,1,1},
    {0,1},
  };

  int code;
  int vlcnum = sym->value1;
  int retval = code_from_bitstream_2d(sym, currStream, &lentab[vlcnum][0], &codtab[vlcnum][0], 16, 1, &code);

  if (retval)
  {
    printf("ERROR: failed to find Total Zeros !cdc\n");
    exit(-1);
  }

  return retval;
}

/*!
 ************************************************************************
 * \brief
 *    read Total Zeros Chroma DC codeword from UVLC-partition
 ************************************************************************
 */
static int readSyntaxElement_TotalZerosChromaDC(VideoParameters *p_Vid, SyntaxElement *sym,  Bitstream *currStream)
{
  static const byte lentab[3][TOTRUN_NUM][16] =
  {
    //YUV420
   {{ 1,2,3,3},
    { 1,2,2},
    { 1,1}},
    //YUV422
   {{ 1,3,3,4,4,4,5,5},
    { 3,2,3,3,3,3,3},
    { 3,3,2,2,3,3},
    { 3,2,2,2,3},
    { 2,2,2,2},
    { 2,2,1},
    { 1,1}},
    //YUV444
   {{ 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},
    { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},
    { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},
    { 5,3,4,4,3,3,3,4,3,4,5,5,5},
    { 4,4,4,3,3,3,3,3,4,5,4,5},
    { 6,5,3,3,3,3,3,3,4,3,6},
    { 6,5,3,3,3,2,3,4,3,6},
    { 6,4,5,3,2,2,3,3,6},
    { 6,6,4,2,2,3,2,5},
    { 5,5,3,2,2,2,4},
    { 4,4,3,3,1,3},
    { 4,4,2,1,3},
    { 3,3,1,2},
    { 2,2,1},
    { 1,1}}
  };

  static const byte codtab[3][TOTRUN_NUM][16] =
  {
    //YUV420
   {{ 1,1,1,0},
    { 1,1,0},
    { 1,0}},
    //YUV422
   {{ 1,2,3,2,3,1,1,0},
    { 0,1,1,4,5,6,7},
    { 0,1,1,2,6,7},
    { 6,0,1,2,7},
    { 0,1,2,3},
    { 0,1,1},
    { 0,1}},
    //YUV444
   {{1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
    {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
    {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
    {3,7,5,4,6,5,4,3,3,2,2,1,0},
    {5,4,3,7,6,5,4,3,2,1,1,0},
    {1,1,7,6,5,4,3,2,1,1,0},
    {1,1,5,4,3,3,2,1,1,0},
    {1,1,1,3,3,2,2,1,0},
    {1,0,1,3,2,1,1,1,},
    {1,0,1,3,2,1,1,},
    {0,1,1,2,1,3},
    {0,1,1,1,1},
    {0,1,1,1},
    {0,1,1},
    {0,1}}
  };

  int code;
  int yuv = p_Vid->active_sps->chroma_format_idc - 1;
  int vlcnum = sym->value1;
  int retval = code_from_bitstream_2d(sym, currStream, &lentab[yuv][vlcnum][0], &codtab[yuv][vlcnum][0], 16, 1, &code);

  if (retval)
  {
    printf("ERROR: failed to find Total Zeros\n");
    exit(-1);
  }

  return retval;
}


/*!
 ************************************************************************
 * \brief
 *    read  Run codeword from UVLC-partition
 ************************************************************************
 */
static int readSyntaxElement_Run(SyntaxElement *sym, Bitstream *currStream)
{
  static const byte lentab[TOTRUN_NUM][16] =
  {
    {1,1},
    {1,2,2},
    {2,2,2,2},
    {2,2,2,3,3},
    {2,2,3,3,3,3},
    {2,3,3,3,3,3,3},
    {3,3,3,3,3,3,3,4,5,6,7,8,9,10,11},
  };

  static const byte codtab[TOTRUN_NUM][16] =
  {
    {1,0},
    {1,1,0},
    {3,2,1,0},
    {3,2,1,1,0},
    {3,2,3,2,1,0},
    {3,0,1,3,2,5,4},
    {7,6,5,4,3,2,1,1,1,1,1,1,1,1,1},
  };
  int code;
  int vlcnum = sym->value1;
  int retval = code_from_bitstream_2d(sym, currStream, &lentab[vlcnum][0], &codtab[vlcnum][0], 16, 1, &code);

  if (retval)
  {
    printf("ERROR: failed to find Run\n");
    exit(-1);
  }

  return retval;
}


/*!
 ************************************************************************
 * \brief
 *    Reads coeff of an 4x4 block (CAVLC)
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */
static void read_coeff_4x4_CAVLC(Macroblock *currMB,
                                 int block_type,
                                 int i, int j, int levarr[16], int runarr[16],
                                 int *number_coefficients)
{
  Slice *currSlice = currMB->p_Slice;
  sps_t *sps = currSlice->active_sps;
  VideoParameters *p_Vid = currMB->p_Vid;
  int mb_nr = currMB->mbAddrX;
  SyntaxElement currSE;
  DataPartition *dP;
  const byte *partMap = assignSE2partition[currSlice->dp_mode];
  Bitstream *currStream;

  int k, code, vlcnum;
  int numcoeff = 0, numtrailingones;
  int level_two_or_higher;
  int numones, totzeros, abslevel, cdc=0, cac=0;
  int zerosleft, ntr, dptype = 0;
  int max_coeff_num = 0, nnz;
  char type[15];
  static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768};    // maximum vlc = 6

    int num_cdc_coeff;
    if (sps->chroma_format_idc != YUV400)
        num_cdc_coeff = (((1 << sps->chroma_format_idc) & (~0x1)) << 1);
    else
        num_cdc_coeff = 0;

  switch (block_type)
  {
  case LUMA:
    max_coeff_num = 16;
    dptype = (currMB->is_intra_block == TRUE) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case LUMA_INTRA16x16DC:
    max_coeff_num = 16;
    dptype = SE_LUM_DC_INTRA;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case LUMA_INTRA16x16AC:
    max_coeff_num = 15;
    dptype = SE_LUM_AC_INTRA;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case CHROMA_DC:
    max_coeff_num = num_cdc_coeff;
    cdc = 1;
    dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case CHROMA_AC:
    max_coeff_num = 15;
    cac = 1;
    dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  default:
    error ("read_coeff_4x4_CAVLC: invalid block type", 600);
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  }

  currSE.type = dptype;
  dP = &(currSlice->partArr[partMap[dptype]]);
  currStream = dP->bitstream;  

  if (!cdc)
  {    
    // luma or chroma AC    
    nnz = (!cac) ? predict_nnz(currMB, LUMA, i<<2, j<<2) : predict_nnz_chroma(currMB, i, ((j-4)<<2));

    currSE.value1 = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));

    readSyntaxElement_NumCoeffTrailingOnes(&currSE, currStream, type);

    numcoeff        =  currSE.value1;
    numtrailingones =  currSE.value2;

    p_Vid->nz_coeff[mb_nr][0][j][i] = (byte) numcoeff;
  }
  else
  {
    // chroma DC
    readSyntaxElement_NumCoeffTrailingOnesChromaDC(p_Vid, &currSE, currStream);

    numcoeff        =  currSE.value1;
    numtrailingones =  currSE.value2;
  }

  memset(levarr, 0, max_coeff_num * sizeof(int));
  memset(runarr, 0, max_coeff_num * sizeof(int));

  numones = numtrailingones;
  *number_coefficients = numcoeff;

  if (numcoeff)
  {
    if (numtrailingones)
    {      
      code = currStream->f(numtrailingones);
      ntr = numtrailingones;
      for (k = numcoeff - 1; k > numcoeff - 1 - numtrailingones; k--)
      {
        ntr --;
        levarr[k] = (code>>ntr)&1 ? -1 : 1;
      }
    }

    // decode levels
    level_two_or_higher = (numcoeff > 3 && numtrailingones == 3)? 0 : 1;
    vlcnum = (numcoeff > 10 && numtrailingones < 3) ? 1 : 0;

    for (k = numcoeff - 1 - numtrailingones; k >= 0; k--)
    {
      if (vlcnum == 0)
        readSyntaxElement_Level_VLC0(&currSE, currStream);
      else
        readSyntaxElement_Level_VLCN(&currSE, vlcnum, currStream);

      if (level_two_or_higher)
      {
        currSE.inf += (currSE.inf > 0) ? 1 : -1;
        level_two_or_higher = 0;
      }

      levarr[k] = currSE.inf;
      abslevel = iabs(levarr[k]);
      if (abslevel  == 1)
        ++numones;

      // update VLC table
      if (abslevel  > incVlc[vlcnum])
        ++vlcnum;

      if (k == numcoeff - 1 - numtrailingones && abslevel >3)
        vlcnum = 2;      
    }

    if (numcoeff < max_coeff_num)
    {
      // decode total run
      vlcnum = numcoeff - 1;
      currSE.value1 = vlcnum;

      if (cdc)
        readSyntaxElement_TotalZerosChromaDC(p_Vid, &currSE, currStream);
      else
        readSyntaxElement_TotalZeros(&currSE, currStream);

      totzeros = currSE.value1;
    }
    else
    {
      totzeros = 0;
    }

    // decode run before each coefficient
    zerosleft = totzeros;
    i = numcoeff - 1;

    if (zerosleft > 0 && i > 0)
    {
      do
      {
        // select VLC for runbefore
        vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);

        currSE.value1 = vlcnum;

        readSyntaxElement_Run(&currSE, currStream);
        runarr[i] = currSE.value1;

        zerosleft -= runarr[i];
        i --;
      } while (zerosleft != 0 && i != 0);
    }
    runarr[i] = zerosleft;    
  } // if numcoeff
}

/*!
 ************************************************************************
 * \brief
 *    Reads coeff of an 4x4 block (CAVLC)
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */
static void read_coeff_4x4_CAVLC_444(Macroblock *currMB,
                                     int block_type,
                                     int i, int j, int levarr[16], int runarr[16],
                                     int *number_coefficients)
{
  Slice *currSlice = currMB->p_Slice;
  sps_t *sps = currSlice->active_sps;
  VideoParameters *p_Vid = currMB->p_Vid;
  int mb_nr = currMB->mbAddrX;
  SyntaxElement currSE;
  DataPartition *dP;
  const byte *partMap = assignSE2partition[currSlice->dp_mode];
  Bitstream *currStream;

  int k, code, vlcnum;
  int numcoeff = 0, numtrailingones;
  int level_two_or_higher;
  int numones, totzeros, abslevel, cdc=0, cac=0;
  int zerosleft, ntr, dptype = 0;
  int max_coeff_num = 0, nnz;
  char type[15];
  static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768};    // maximum vlc = 6

    int num_cdc_coeff;
    if (sps->chroma_format_idc != YUV400)
        num_cdc_coeff = (((1 << sps->chroma_format_idc) & (~0x1)) << 1);
    else
        num_cdc_coeff = 0;

  switch (block_type)
  {
  case LUMA:
    max_coeff_num = 16;
    dptype = (currMB->is_intra_block == TRUE) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case LUMA_INTRA16x16DC:
    max_coeff_num = 16;
    dptype = SE_LUM_DC_INTRA;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case LUMA_INTRA16x16AC:
    max_coeff_num = 15;
    dptype = SE_LUM_AC_INTRA;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case CB:
    max_coeff_num = 16;
    dptype = ((currMB->is_intra_block == TRUE)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    p_Vid->nz_coeff[mb_nr][1][j][i] = 0; 
    break;
  case CB_INTRA16x16DC:
    max_coeff_num = 16;
    dptype = SE_LUM_DC_INTRA;
    p_Vid->nz_coeff[mb_nr][1][j][i] = 0; 
    break;
  case CB_INTRA16x16AC:
    max_coeff_num = 15;
    dptype = SE_LUM_AC_INTRA;
    p_Vid->nz_coeff[mb_nr][1][j][i] = 0; 
    break;
  case CR:
    max_coeff_num = 16;
    dptype = ((currMB->is_intra_block == TRUE)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    p_Vid->nz_coeff[mb_nr][2][j][i] = 0; 
    break;
  case CR_INTRA16x16DC:
    max_coeff_num = 16;
    dptype = SE_LUM_DC_INTRA;
    p_Vid->nz_coeff[mb_nr][2][j][i] = 0; 
    break;
  case CR_INTRA16x16AC:
    max_coeff_num = 15;
    dptype = SE_LUM_AC_INTRA;
    p_Vid->nz_coeff[mb_nr][2][j][i] = 0; 
    break;        
  case CHROMA_DC:
    max_coeff_num = num_cdc_coeff;
    cdc = 1;
    dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  case CHROMA_AC:
    max_coeff_num = 15;
    cac = 1;
    dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  default:
    error ("read_coeff_4x4_CAVLC: invalid block type", 600);
    p_Vid->nz_coeff[mb_nr][0][j][i] = 0; 
    break;
  }

  currSE.type = dptype;
  dP = &(currSlice->partArr[partMap[dptype]]);
  currStream = dP->bitstream;  

  if (!cdc)
  {    
    // luma or chroma AC    
    if(block_type==LUMA || block_type==LUMA_INTRA16x16DC || block_type==LUMA_INTRA16x16AC ||block_type==CHROMA_AC)
    {
      nnz = (!cac) ? predict_nnz(currMB, LUMA, i<<2, j<<2) : predict_nnz_chroma(currMB, i, ((j-4)<<2));
    }
    else if (block_type==CB || block_type==CB_INTRA16x16DC || block_type==CB_INTRA16x16AC)
    {   
      nnz = predict_nnz(currMB, CB, i<<2, j<<2);
    }
    else
    { 
      nnz = predict_nnz(currMB, CR, i<<2, j<<2);
    }

    currSE.value1 = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));

    readSyntaxElement_NumCoeffTrailingOnes(&currSE, currStream, type);

    numcoeff        =  currSE.value1;
    numtrailingones =  currSE.value2;

    if(block_type==LUMA || block_type==LUMA_INTRA16x16DC || block_type==LUMA_INTRA16x16AC ||block_type==CHROMA_AC)
      p_Vid->nz_coeff[mb_nr][0][j][i] = (byte) numcoeff;
    else if (block_type==CB || block_type==CB_INTRA16x16DC || block_type==CB_INTRA16x16AC)
      p_Vid->nz_coeff[mb_nr][1][j][i] = (byte) numcoeff;
    else
      p_Vid->nz_coeff[mb_nr][2][j][i] = (byte) numcoeff;        
  }
  else
  {
    // chroma DC
    readSyntaxElement_NumCoeffTrailingOnesChromaDC(p_Vid, &currSE, currStream);

    numcoeff        =  currSE.value1;
    numtrailingones =  currSE.value2;
  }

  memset(levarr, 0, max_coeff_num * sizeof(int));
  memset(runarr, 0, max_coeff_num * sizeof(int));

  numones = numtrailingones;
  *number_coefficients = numcoeff;

  if (numcoeff)
  {
    if (numtrailingones)
    {      
      code = currStream->f(numtrailingones);
      ntr = numtrailingones;
      for (k = numcoeff - 1; k > numcoeff - 1 - numtrailingones; k--)
      {
        ntr --;
        levarr[k] = (code>>ntr)&1 ? -1 : 1;
      }
    }

    // decode levels
    level_two_or_higher = (numcoeff > 3 && numtrailingones == 3)? 0 : 1;
    vlcnum = (numcoeff > 10 && numtrailingones < 3) ? 1 : 0;

    for (k = numcoeff - 1 - numtrailingones; k >= 0; k--)
    {

      if (vlcnum == 0)
        readSyntaxElement_Level_VLC0(&currSE, currStream);
      else
        readSyntaxElement_Level_VLCN(&currSE, vlcnum, currStream);

      if (level_two_or_higher)
      {
        currSE.inf += (currSE.inf > 0) ? 1 : -1;
        level_two_or_higher = 0;
      }

      levarr[k] = currSE.inf;
      abslevel = iabs(levarr[k]);
      if (abslevel  == 1)
        ++numones;

      // update VLC table
      if (abslevel  > incVlc[vlcnum])
        ++vlcnum;

      if (k == numcoeff - 1 - numtrailingones && abslevel >3)
        vlcnum = 2;      
    }

    if (numcoeff < max_coeff_num)
    {
      // decode total run
      vlcnum = numcoeff - 1;
      currSE.value1 = vlcnum;

      if (cdc)
        readSyntaxElement_TotalZerosChromaDC(p_Vid, &currSE, currStream);
      else
        readSyntaxElement_TotalZeros(&currSE, currStream);

      totzeros = currSE.value1;
    }
    else
    {
      totzeros = 0;
    }

    // decode run before each coefficient
    zerosleft = totzeros;
    i = numcoeff - 1;

    if (zerosleft > 0 && i > 0)
    {
      do
      {
        // select VLC for runbefore
        vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);

        currSE.value1 = vlcnum;

        readSyntaxElement_Run(&currSE, currStream);
        runarr[i] = currSE.value1;

        zerosleft -= runarr[i];
        i --;
      } while (zerosleft != 0 && i != 0);
    }
    runarr[i] = zerosleft;    
  } // if numcoeff
}

/*!
************************************************************************
* \brief
*    Get coefficients (run/level) of 4x4 blocks in a MB
*    from the NAL (CABAC Mode)
************************************************************************
*/
static void read_comp_coeff_4x4_CAVLC(Macroblock *currMB, ColorPlane pl, int (*InvLevelScale4x4)[4], int qp_per, int cbp, byte **nzcoeff)
{
  int block_y, block_x, b8;
  int i, j, k;
  int i0, j0;
  int levarr[16] = {0}, runarr[16] = {0}, numcoeff;
  Slice *currSlice = currMB->p_Slice;
  const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
  const byte *pos_scan_4x4 = pos_scan4x4[0];
  int start_scan = IS_I16MB(currMB) ? 1 : 0;
  int64 *cur_cbp = &currMB->s_cbp[pl].blk;
  int cur_context; 
  int block_y4, block_x4;

  if (IS_I16MB(currMB))
  {
    if (pl == PLANE_Y)
      cur_context = LUMA_INTRA16x16AC;
    else if (pl == PLANE_U)
      cur_context = CB_INTRA16x16AC;
    else
      cur_context = CR_INTRA16x16AC;
  }
  else
  {
    if (pl == PLANE_Y)
      cur_context = LUMA;
    else if (pl == PLANE_U)
      cur_context = CB;
    else
      cur_context = CR;
  }


  for (block_y = 0; block_y < 4; block_y += 2) /* all modes */
  {
    block_y4 = block_y << 2;
    for (block_x = 0; block_x < 4; block_x += 2)
    {
      block_x4 = block_x << 2;
      b8 = (block_y + (block_x >> 1));

      if (cbp & (1 << b8))  // test if the block contains any coefficients
      {
        for (j = block_y4; j < block_y4 + 8; j += BLOCK_SIZE)
        {
          for (i = block_x4; i < block_x4 + 8; i += BLOCK_SIZE)
          {
            currSlice->read_coeff_4x4_CAVLC(currMB, cur_context, i >> 2, j >> 2, levarr, runarr, &numcoeff);
            pos_scan_4x4 = pos_scan4x4[start_scan];

            for (k = 0; k < numcoeff; ++k)
            {
              if (levarr[k] != 0)
              {
                pos_scan_4x4 += (runarr[k] << 1);

                i0 = *pos_scan_4x4++;
                j0 = *pos_scan_4x4++;

                // inverse quant for 4x4 transform only
                *cur_cbp |= i64_power2(j + (i >> 2));

                currSlice->cof[pl][j + j0][i + i0]= rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0])<<qp_per, 4);
              }
            }
          }
        }
      }
      else
      {
        nzcoeff[block_y    ][block_x    ] = 0;
        nzcoeff[block_y    ][block_x + 1] = 0;
        nzcoeff[block_y + 1][block_x    ] = 0;
        nzcoeff[block_y + 1][block_x + 1] = 0;
      }
    }
  }      
}

/*!
************************************************************************
* \brief
*    Get coefficients (run/level) of 4x4 blocks in a MB
*    from the NAL (CAVLC Lossless Mode)
************************************************************************
*/
static void read_comp_coeff_4x4_CAVLC_ls (Macroblock *currMB, ColorPlane pl, int (*InvLevelScale4x4)[4], int qp_per, int cbp, byte **nzcoeff)
{
  int block_y, block_x, b8;
  int i, j, k;
  int i0, j0;
  int levarr[16] = {0}, runarr[16] = {0}, numcoeff;
  Slice *currSlice = currMB->p_Slice;
  const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
  int start_scan = IS_I16MB(currMB) ? 1 : 0;
  int64 *cur_cbp = &currMB->s_cbp[pl].blk;
  int coef_ctr, cur_context; 

  if (IS_I16MB(currMB))
  {
    if (pl == PLANE_Y)
      cur_context = LUMA_INTRA16x16AC;
    else if (pl == PLANE_U)
      cur_context = CB_INTRA16x16AC;
    else
      cur_context = CR_INTRA16x16AC;
  }
  else
  {
    if (pl == PLANE_Y)
      cur_context = LUMA;
    else if (pl == PLANE_U)
      cur_context = CB;
    else
      cur_context = CR;
  }

  for (block_y=0; block_y < 4; block_y += 2) /* all modes */
  {
    for (block_x=0; block_x < 4; block_x += 2)
    {
      b8 = 2*(block_y>>1) + (block_x>>1);

      if (cbp & (1<<b8))  /* are there any coeff in current block at all */
      {
        for (j=block_y; j < block_y+2; ++j)
        {
          for (i=block_x; i < block_x+2; ++i)
          {
            currSlice->read_coeff_4x4_CAVLC(currMB, cur_context, i, j, levarr, runarr, &numcoeff);

            coef_ctr = start_scan - 1;

            for (k = 0; k < numcoeff; ++k)
            {
              if (levarr[k] != 0)
              {
                coef_ctr += runarr[k]+1;

                i0=pos_scan4x4[coef_ctr][0];
                j0=pos_scan4x4[coef_ctr][1];

                *cur_cbp |= i64_power2((j<<2) + i);
                currSlice->cof[pl][(j<<2) + j0][(i<<2) + i0]= levarr[k];
              }
            }
          }
        }
      }
      else
      {
        nzcoeff[block_y    ][block_x    ] = 0;
        nzcoeff[block_y    ][block_x + 1] = 0;
        nzcoeff[block_y + 1][block_x    ] = 0;
        nzcoeff[block_y + 1][block_x + 1] = 0;
      }
    }
  }    
}

/*!
************************************************************************
* \brief
*    Get coefficients (run/level) of 4x4 blocks in a MB
*    from the NAL (CABAC Mode)
************************************************************************
*/
static void read_comp_coeff_8x8_CAVLC (Macroblock *currMB, ColorPlane pl, int (*InvLevelScale8x8)[8], int qp_per, int cbp, byte **nzcoeff)
{
  int block_y, block_x, b4, b8;
  int block_y4, block_x4;
  int i, j, k;
  int i0, j0;
  int levarr[16] = {0}, runarr[16] = {0}, numcoeff;
  Slice *currSlice = currMB->p_Slice;
  const byte (*pos_scan8x8)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN8x8 : FIELD_SCAN8x8;
  int start_scan = IS_I16MB(currMB) ? 1 : 0;
  int64 *cur_cbp = &currMB->s_cbp[pl].blk;
  int coef_ctr, cur_context; 

  if (IS_I16MB(currMB))
  {
    if (pl == PLANE_Y)
      cur_context = LUMA_INTRA16x16AC;
    else if (pl == PLANE_U)
      cur_context = CB_INTRA16x16AC;
    else
      cur_context = CR_INTRA16x16AC;
  }
  else
  {
    if (pl == PLANE_Y)
      cur_context = LUMA;
    else if (pl == PLANE_U)
      cur_context = CB;
    else
      cur_context = CR;
  }

  for (block_y = 0; block_y < 4; block_y += 2) /* all modes */
  {
    block_y4 = block_y << 2;

    for (block_x = 0; block_x < 4; block_x += 2)
    {
      block_x4 = block_x << 2;
      b8 = block_y + (block_x>>1);

      if (cbp & (1<<b8))  /* are there any coeff in current block at all */
      {
        for (j = block_y; j < block_y + 2; ++j)
        {
          for (i = block_x; i < block_x + 2; ++i)
          {
            currSlice->read_coeff_4x4_CAVLC(currMB, cur_context, i, j, levarr, runarr, &numcoeff);

            coef_ctr = start_scan - 1;

            for (k = 0; k < numcoeff; ++k)
            {
              if (levarr[k] != 0)
              {
                coef_ctr += runarr[k] + 1;

                // do same as CABAC for deblocking: any coeff in the 8x8 marks all the 4x4s
                //as containing coefficients
                *cur_cbp |= 51 << (block_y4 + block_x);

                b4 = (coef_ctr << 2) + 2*(j - block_y) + (i - block_x);

                i0 = pos_scan8x8[b4][0];
                j0 = pos_scan8x8[b4][1];

                currSlice->mb_rres[pl][block_y4 +j0][block_x4 +i0] = rshift_rnd_sf((levarr[k] * InvLevelScale8x8[j0][i0])<<qp_per, 6); // dequantization
              }
            }//else (!currMB->luma_transform_size_8x8_flag)
          }
        }
      }
      else
      {
        nzcoeff[block_y    ][block_x    ] = 0;
        nzcoeff[block_y    ][block_x + 1] = 0;
        nzcoeff[block_y + 1][block_x    ] = 0;
        nzcoeff[block_y + 1][block_x + 1] = 0;
      }
    }
  }   
}

/*!
************************************************************************
* \brief
*    Get coefficients (run/level) of 8x8 blocks in a MB
*    from the NAL (CAVLC Lossless Mode)
************************************************************************
*/
static void read_comp_coeff_8x8_CAVLC_ls (Macroblock *currMB, ColorPlane pl, int (*InvLevelScale8x8)[8], int qp_per, int cbp, byte **nzcoeff)
{
  int block_y, block_x, b4, b8;
  int i, j, k;
  int levarr[16] = {0}, runarr[16] = {0}, numcoeff;
  Slice *currSlice = currMB->p_Slice;
  const byte (*pos_scan8x8)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN8x8 : FIELD_SCAN8x8;
  int start_scan = IS_I16MB(currMB) ? 1 : 0;
  int64 *cur_cbp = &currMB->s_cbp[pl].blk;
  int coef_ctr, cur_context; 

  if (IS_I16MB(currMB))
  {
    if (pl == PLANE_Y)
      cur_context = LUMA_INTRA16x16AC;
    else if (pl == PLANE_U)
      cur_context = CB_INTRA16x16AC;
    else
      cur_context = CR_INTRA16x16AC;
  }
  else
  {
    if (pl == PLANE_Y)
      cur_context = LUMA;
    else if (pl == PLANE_U)
      cur_context = CB;
    else
      cur_context = CR;
  }

  for (block_y=0; block_y < 4; block_y += 2) /* all modes */
  {
    for (block_x=0; block_x < 4; block_x += 2)
    {
      b8 = 2*(block_y>>1) + (block_x>>1);

      if (cbp & (1<<b8))  /* are there any coeff in current block at all */
      {
        int iz, jz;

        for (j=block_y; j < block_y+2; ++j)
        {
          for (i=block_x; i < block_x+2; ++i)
          {

            currSlice->read_coeff_4x4_CAVLC(currMB, cur_context, i, j, levarr, runarr, &numcoeff);

            coef_ctr = start_scan - 1;

            for (k = 0; k < numcoeff; ++k)
            {
              if (levarr[k] != 0)
              {
                coef_ctr += runarr[k]+1;

                // do same as CABAC for deblocking: any coeff in the 8x8 marks all the 4x4s
                //as containing coefficients
                *cur_cbp  |= 51 << ((block_y<<2) + block_x);

                b4 = 2*(j-block_y)+(i-block_x);

                iz=pos_scan8x8[coef_ctr*4+b4][0];
                jz=pos_scan8x8[coef_ctr*4+b4][1];

                currSlice->mb_rres[pl][block_y*4 +jz][block_x*4 +iz] = levarr[k];
              }
            }
          }
        }
      }
      else
      {
        nzcoeff[block_y    ][block_x    ] = 0;
        nzcoeff[block_y    ][block_x + 1] = 0;
        nzcoeff[block_y + 1][block_x    ] = 0;
        nzcoeff[block_y + 1][block_x + 1] = 0;
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Get coded block pattern and coefficients (run/level)
 *    from the NAL
 ************************************************************************
 */
static void read_CBP_and_coeffs_from_NAL_CAVLC_400(Macroblock *currMB)
{
  int k;
  int mb_nr = currMB->mbAddrX;
  int cbp;
  SyntaxElement currSE;
  DataPartition *dP = NULL;
  Slice *currSlice = currMB->p_Slice;
  const byte *partMap = assignSE2partition[currSlice->dp_mode];
  int i0, j0;

  int levarr[16], runarr[16], numcoeff;

  int qp_per, qp_rem;
  VideoParameters *p_Vid = currMB->p_Vid;

  int intra = (currMB->is_intra_block == TRUE);

  int need_transform_size_flag;

  int (*InvLevelScale4x4)[4] = NULL;
  int (*InvLevelScale8x8)[8] = NULL;
  // select scan type
  const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
  const byte *pos_scan_4x4 = pos_scan4x4[0];


  // read CBP if not new intra mode
  if (!IS_I16MB (currMB))
  {
    //=====   C B P   =====
    //---------------------
    currSE.type = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB) 
      ? SE_CBP_INTRA
      : SE_CBP_INTER;

    dP = &(currSlice->partArr[partMap[currSE.type]]);

    currSE.mapping = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)
      ? currSlice->linfo_cbp_intra
      : currSlice->linfo_cbp_inter;

    dP->readSyntaxElement(currMB, &currSE, dP);
    currMB->cbp = cbp = currSE.value1;


    //============= Transform size flag for INTER MBs =============
    //-------------------------------------------------------------
    need_transform_size_flag = (((currMB->mb_type >= 1 && currMB->mb_type <= 3)||
      (IS_DIRECT(currMB) && p_Vid->active_sps->direct_8x8_inference_flag) ||
      (currMB->NoMbPartLessThan8x8Flag))
      && currMB->mb_type != I8MB && currMB->mb_type != I4MB
      && (currMB->cbp&15)
      && currSlice->active_pps->transform_8x8_mode_flag);

    if (need_transform_size_flag)
    {
      currSE.type   =  SE_HEADER;
      dP = &(currSlice->partArr[partMap[SE_HEADER]]);

      currMB->luma_transform_size_8x8_flag = dP->bitstream->f(1);
    }

    //=====   DQUANT   =====
    //----------------------
    // Delta quant only if nonzero coeffs
    if (cbp !=0)
    {
      read_delta_quant(&currSE, dP, currMB, partMap, ((currMB->is_intra_block == FALSE)) ? SE_DELTA_QUANT_INTER : SE_DELTA_QUANT_INTRA);

      if (currSlice->dp_mode)
      {
        if ((currMB->is_intra_block == FALSE) && currSlice->dpC_NotPresent ) 
          currMB->dpl_flag = 1;

        if( intra && currSlice->dpB_NotPresent )
        {
          currMB->ei_flag = 1;
          currMB->dpl_flag = 1;
        }

        // check for prediction from neighbours
        check_dp_neighbors (currMB);
        if (currMB->dpl_flag)
        {
          cbp = 0; 
          currMB->cbp = cbp;
        }
      }
    }
  }
  else  // read DC coeffs for new intra modes
  {
    cbp = currMB->cbp;
  
    read_delta_quant(&currSE, dP, currMB, partMap, SE_DELTA_QUANT_INTRA);

    if (currSlice->dp_mode)
    {  
      if (currSlice->dpB_NotPresent)
      {
        currMB->ei_flag  = 1;
        currMB->dpl_flag = 1;
      }
      check_dp_neighbors (currMB);
      if (currMB->dpl_flag)
      {
        currMB->cbp = cbp = 0; 
      }
    }

    if (!currMB->dpl_flag)
    {
      pos_scan_4x4 = pos_scan4x4[0];

      currSlice->read_coeff_4x4_CAVLC(currMB, LUMA_INTRA16x16DC, 0, 0, levarr, runarr, &numcoeff);

      for(k = 0; k < numcoeff; ++k)
      {
        if (levarr[k] != 0)                     // leave if level == 0
        {
          pos_scan_4x4 += 2 * runarr[k];

          i0 = ((*pos_scan_4x4++) << 2);
          j0 = ((*pos_scan_4x4++) << 2);

          currSlice->cof[0][j0][i0] = levarr[k];// add new intra DC coeff
        }
      }


      if(currMB->is_lossless == FALSE)
        itrans_2(currMB, (ColorPlane) currSlice->colour_plane_id);// transform new intra DC
    }
  }

  update_qp(currMB, currSlice->SliceQpY);

  qp_per = currMB->qp_scaled[PLANE_Y] / 6;
  qp_rem = currMB->qp_scaled[PLANE_Y] % 6;

  InvLevelScale4x4 = intra? currSlice->InvLevelScale4x4_Intra[currSlice->colour_plane_id][qp_rem] : currSlice->InvLevelScale4x4_Inter[currSlice->colour_plane_id][qp_rem];
  InvLevelScale8x8 = intra? currSlice->InvLevelScale8x8_Intra[currSlice->colour_plane_id][qp_rem] : currSlice->InvLevelScale8x8_Inter[currSlice->colour_plane_id][qp_rem];

  // luma coefficients
  if (cbp)
  {
    if (!currMB->luma_transform_size_8x8_flag) // 4x4 transform
    {
      currMB->read_comp_coeff_4x4_CAVLC (currMB, PLANE_Y, InvLevelScale4x4, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
    }
    else // 8x8 transform
    {
      currMB->read_comp_coeff_8x8_CAVLC (currMB, PLANE_Y, InvLevelScale8x8, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
    }
  }
  else
  {
    memset(p_Vid->nz_coeff[mb_nr][0][0], 0, BLOCK_PIXELS * sizeof(byte));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Get coded block pattern and coefficients (run/level)
 *    from the NAL
 ************************************************************************
 */
static void read_CBP_and_coeffs_from_NAL_CAVLC_422(Macroblock *currMB)
{
  int i,j,k;
  int mb_nr = currMB->mbAddrX;
  int cbp;
  SyntaxElement currSE;
  DataPartition *dP = NULL;
  Slice *currSlice = currMB->p_Slice;
  sps_t *sps = currSlice->active_sps;
  const byte *partMap = assignSE2partition[currSlice->dp_mode];
  int coef_ctr, i0, j0, b8;
  int ll;
  int levarr[16], runarr[16], numcoeff;

  int qp_per, qp_rem;
  VideoParameters *p_Vid = currMB->p_Vid;

  int uv; 
  int qp_per_uv[2];
  int qp_rem_uv[2];

  int intra = (currMB->is_intra_block == TRUE);

  int b4;
  //StorablePicture *dec_picture = currSlice->dec_picture;
  int m6[4];

  int need_transform_size_flag;

  int (*InvLevelScale4x4)[4] = NULL;
  int (*InvLevelScale8x8)[8] = NULL;
  // select scan type
  const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
  const byte *pos_scan_4x4 = pos_scan4x4[0];

    int num_blk8x8_uv;
    if (sps->chroma_format_idc != YUV400)
        num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    else
        num_blk8x8_uv = 0;
    int num_uv_blocks;
    if (sps->chroma_format_idc != YUV400)
        num_uv_blocks = (((1 << sps->chroma_format_idc) & (~0x1)) >> 1);
    else
        num_uv_blocks = 0;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;


  // read CBP if not new intra mode
  if (!IS_I16MB (currMB))
  {
    //=====   C B P   =====
    //---------------------
    currSE.type = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB) 
      ? SE_CBP_INTRA
      : SE_CBP_INTER;

    dP = &(currSlice->partArr[partMap[currSE.type]]);

    currSE.mapping = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)
      ? currSlice->linfo_cbp_intra
      : currSlice->linfo_cbp_inter;

    dP->readSyntaxElement(currMB, &currSE, dP);
    currMB->cbp = cbp = currSE.value1;


    //============= Transform size flag for INTER MBs =============
    //-------------------------------------------------------------
    need_transform_size_flag = (((currMB->mb_type >= 1 && currMB->mb_type <= 3)||
      (IS_DIRECT(currMB) && p_Vid->active_sps->direct_8x8_inference_flag) ||
      (currMB->NoMbPartLessThan8x8Flag))
      && currMB->mb_type != I8MB && currMB->mb_type != I4MB
      && (currMB->cbp&15)
      && currSlice->active_pps->transform_8x8_mode_flag);

    if (need_transform_size_flag)
    {
      currSE.type   =  SE_HEADER;
      dP = &(currSlice->partArr[partMap[SE_HEADER]]);

      currMB->luma_transform_size_8x8_flag = dP->bitstream->f(1);
    }

    //=====   DQUANT   =====
    //----------------------
    // Delta quant only if nonzero coeffs
    if (cbp !=0)
    {
      read_delta_quant(&currSE, dP, currMB, partMap, ((currMB->is_intra_block == FALSE)) ? SE_DELTA_QUANT_INTER : SE_DELTA_QUANT_INTRA);

      if (currSlice->dp_mode)
      {
        if ((currMB->is_intra_block == FALSE) && currSlice->dpC_NotPresent ) 
          currMB->dpl_flag = 1;

        if( intra && currSlice->dpB_NotPresent )
        {
          currMB->ei_flag = 1;
          currMB->dpl_flag = 1;
        }

        // check for prediction from neighbours
        check_dp_neighbors (currMB);
        if (currMB->dpl_flag)
        {
          cbp = 0; 
          currMB->cbp = cbp;
        }
      }
    }
  }
  else  // read DC coeffs for new intra modes
  {
    cbp = currMB->cbp;

    read_delta_quant(&currSE, dP, currMB, partMap, SE_DELTA_QUANT_INTRA);

    if (currSlice->dp_mode)
    {  
      if (currSlice->dpB_NotPresent)
      {
        currMB->ei_flag  = 1;
        currMB->dpl_flag = 1;
      }
      check_dp_neighbors (currMB);
      if (currMB->dpl_flag)
      {
        currMB->cbp = cbp = 0; 
      }
    }

    if (!currMB->dpl_flag)
    {
      pos_scan_4x4 = pos_scan4x4[0];

      currSlice->read_coeff_4x4_CAVLC(currMB, LUMA_INTRA16x16DC, 0, 0, levarr, runarr, &numcoeff);

      for(k = 0; k < numcoeff; ++k)
      {
        if (levarr[k] != 0)                     // leave if level == 0
        {
          pos_scan_4x4 += 2 * runarr[k];

          i0 = ((*pos_scan_4x4++) << 2);
          j0 = ((*pos_scan_4x4++) << 2);

          currSlice->cof[0][j0][i0] = levarr[k];// add new intra DC coeff
        }
      }


      if(currMB->is_lossless == FALSE)
        itrans_2(currMB, (ColorPlane) currSlice->colour_plane_id);// transform new intra DC
    }
  }

  update_qp(currMB, currSlice->SliceQpY);

  qp_per = currMB->qp_scaled[currSlice->colour_plane_id] / 6;
  qp_rem = currMB->qp_scaled[currSlice->colour_plane_id] % 6;

  //init quant parameters for chroma 
  for(i=0; i < 2; ++i)
  {
    qp_per_uv[i] = currMB->qp_scaled[i + 1] / 6;
    qp_rem_uv[i] = currMB->qp_scaled[i + 1] % 6;
  }

  InvLevelScale4x4 = intra? currSlice->InvLevelScale4x4_Intra[currSlice->colour_plane_id][qp_rem] : currSlice->InvLevelScale4x4_Inter[currSlice->colour_plane_id][qp_rem];
  InvLevelScale8x8 = intra? currSlice->InvLevelScale8x8_Intra[currSlice->colour_plane_id][qp_rem] : currSlice->InvLevelScale8x8_Inter[currSlice->colour_plane_id][qp_rem];

  // luma coefficients
  if (cbp)
  {
    if (!currMB->luma_transform_size_8x8_flag) // 4x4 transform
    {
      currMB->read_comp_coeff_4x4_CAVLC (currMB, PLANE_Y, InvLevelScale4x4, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
    }
    else // 8x8 transform
    {
      currMB->read_comp_coeff_8x8_CAVLC (currMB, PLANE_Y, InvLevelScale8x8, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
    }
  }
  else
  {
    memset(p_Vid->nz_coeff[mb_nr][0][0], 0, BLOCK_PIXELS * sizeof(byte));
  }

  //========================== CHROMA DC ============================
  //-----------------------------------------------------------------
  // chroma DC coeff
  if(cbp>15)
  {    
    for (ll=0;ll<3;ll+=2)
    {
      int (*InvLevelScale4x4)[4] = NULL;
      uv = ll>>1;
      {
        int **imgcof = currSlice->cof[PLANE_U + uv];
        int m3[2][4] = {{0,0,0,0},{0,0,0,0}};
        int m4[2][4] = {{0,0,0,0},{0,0,0,0}};
        int qp_per_uv_dc = (currMB->qpc[uv] + 3 + sps->QpBdOffsetC) / 6;       //for YUV422 only
        int qp_rem_uv_dc = (currMB->qpc[uv] + 3 + sps->QpBdOffsetC) % 6;       //for YUV422 only
        if (intra)
          InvLevelScale4x4 = currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv_dc];
        else 
          InvLevelScale4x4 = currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv_dc];


        //===================== CHROMA DC YUV422 ======================
        currSlice->read_coeff_4x4_CAVLC(currMB, CHROMA_DC, 0, 0, levarr, runarr, &numcoeff);
        coef_ctr=-1;
        for(k = 0; k < numcoeff; ++k)
        {
          if (levarr[k] != 0)
          {
            currMB->s_cbp[0].blk |= ((int64)0xff0000) << (ll<<2);
            coef_ctr += runarr[k]+1;
            i0 = SCAN_YUV422[coef_ctr][0];
            j0 = SCAN_YUV422[coef_ctr][1];

            m3[i0][j0]=levarr[k];
          }
        }

        // inverse CHROMA DC YUV422 transform
        // horizontal
        if(currMB->is_lossless == FALSE)
        {
          m4[0][0] = m3[0][0] + m3[1][0];
          m4[0][1] = m3[0][1] + m3[1][1];
          m4[0][2] = m3[0][2] + m3[1][2];
          m4[0][3] = m3[0][3] + m3[1][3];

          m4[1][0] = m3[0][0] - m3[1][0];
          m4[1][1] = m3[0][1] - m3[1][1];
          m4[1][2] = m3[0][2] - m3[1][2];
          m4[1][3] = m3[0][3] - m3[1][3];

          for (i = 0; i < 2; ++i)
          {
            m6[0] = m4[i][0] + m4[i][2];
            m6[1] = m4[i][0] - m4[i][2];
            m6[2] = m4[i][1] - m4[i][3];
            m6[3] = m4[i][1] + m4[i][3];

            imgcof[ 0][i<<2] = m6[0] + m6[3];
            imgcof[ 4][i<<2] = m6[1] + m6[2];
            imgcof[ 8][i<<2] = m6[1] - m6[2];
            imgcof[12][i<<2] = m6[0] - m6[3];
          }//for (i=0;i<2;++i)

          for(j = 0;j < mb_cr_size_y; j += BLOCK_SIZE)
          {
            for(i=0;i < mb_cr_size_x;i+=BLOCK_SIZE)
            {
              imgcof[j][i] = rshift_rnd_sf((imgcof[j][i] * InvLevelScale4x4[0][0]) << qp_per_uv_dc, 6);
            }
          }
        }
        else
        {
          for(j=0;j<4;++j)
          {
            currSlice->cof[PLANE_U + uv][j<<2][0] = m3[0][j];
            currSlice->cof[PLANE_U + uv][j<<2][4] = m3[1][j];
          }
        }

      }
    }//for (ll=0;ll<3;ll+=2)    
  }

  //========================== CHROMA AC ============================
  //-----------------------------------------------------------------
  // chroma AC coeff, all zero fram start_scan
  if (cbp<=31)
  {
    memset(p_Vid->nz_coeff [mb_nr ][1][0], 0, 2 * BLOCK_PIXELS * sizeof(byte));
  }
  else
  {
    if(currMB->is_lossless == FALSE)
    {
      for (b8=0; b8 < num_blk8x8_uv; ++b8)
      {
        currMB->is_v_block = uv = (b8 > (num_uv_blocks - 1 ));
        InvLevelScale4x4 = intra ? currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv[uv]] : currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv[uv]];

        for (b4=0; b4 < 4; ++b4)
        {
          i = cofuv_blk_x[1][b8][b4];
          j = cofuv_blk_y[1][b8][b4];

          currSlice->read_coeff_4x4_CAVLC(currMB, CHROMA_AC, i + 2*uv, j + 4, levarr, runarr, &numcoeff);
          coef_ctr = 0;

          for(k = 0; k < numcoeff;++k)
          {
            if (levarr[k] != 0)
            {
              currMB->s_cbp[0].blk |= i64_power2(cbp_blk_chroma[b8][b4]);
              coef_ctr += runarr[k] + 1;

              i0=pos_scan4x4[coef_ctr][0];
              j0=pos_scan4x4[coef_ctr][1];

              currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0])<<qp_per_uv[uv], 4);
            }
          }
        }
      }        
    }
    else
    {
      for (b8=0; b8 < num_blk8x8_uv; ++b8)
      {
        currMB->is_v_block = uv = (b8 > (num_uv_blocks - 1 ));

        for (b4=0; b4 < 4; ++b4)
        {
          i = cofuv_blk_x[1][b8][b4];
          j = cofuv_blk_y[1][b8][b4];

          currSlice->read_coeff_4x4_CAVLC(currMB, CHROMA_AC, i + 2*uv, j + 4, levarr, runarr, &numcoeff);
          coef_ctr = 0;

          for(k = 0; k < numcoeff;++k)
          {
            if (levarr[k] != 0)
            {
              currMB->s_cbp[0].blk |= i64_power2(cbp_blk_chroma[b8][b4]);
              coef_ctr += runarr[k] + 1;

              i0=pos_scan4x4[coef_ctr][0];
              j0=pos_scan4x4[coef_ctr][1];

              currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = levarr[k];
            }
          }
        }
      }        
    }
  } //if (dec_picture->chroma_format_idc != YUV400)
}

/*!
 ************************************************************************
 * \brief
 *    Get coded block pattern and coefficients (run/level)
 *    from the NAL
 ************************************************************************
 */
static void read_CBP_and_coeffs_from_NAL_CAVLC_444(Macroblock *currMB)
{
  int i,k;
  int mb_nr = currMB->mbAddrX;
  int cbp;
  SyntaxElement currSE;
  DataPartition *dP = NULL;
  Slice *currSlice = currMB->p_Slice;
  const byte *partMap = assignSE2partition[currSlice->dp_mode];
  int coef_ctr, i0, j0;
  int levarr[16], runarr[16], numcoeff;

  int qp_per, qp_rem;
  VideoParameters *p_Vid = currMB->p_Vid;

  int uv; 
  int qp_per_uv[3];
  int qp_rem_uv[3];

  int intra = (currMB->is_intra_block == TRUE);

  int need_transform_size_flag;

  int (*InvLevelScale4x4)[4] = NULL;
  int (*InvLevelScale8x8)[8] = NULL;
  // select scan type
  const byte (*pos_scan4x4)[2] = (!currSlice->field_pic_flag && (!currMB->mb_field_decoding_flag)) ? SNGL_SCAN : FIELD_SCAN;
  const byte *pos_scan_4x4 = pos_scan4x4[0];

  // read CBP if not new intra mode
  if (!IS_I16MB (currMB))
  {
    //=====   C B P   =====
    //---------------------
    currSE.type = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB) 
      ? SE_CBP_INTRA
      : SE_CBP_INTER;

    dP = &(currSlice->partArr[partMap[currSE.type]]);

    currSE.mapping = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)
      ? currSlice->linfo_cbp_intra
      : currSlice->linfo_cbp_inter;

    dP->readSyntaxElement(currMB, &currSE, dP);
    currMB->cbp = cbp = currSE.value1;


    //============= Transform size flag for INTER MBs =============
    //-------------------------------------------------------------
    need_transform_size_flag = (((currMB->mb_type >= 1 && currMB->mb_type <= 3)||
      (IS_DIRECT(currMB) && p_Vid->active_sps->direct_8x8_inference_flag) ||
      (currMB->NoMbPartLessThan8x8Flag))
      && currMB->mb_type != I8MB && currMB->mb_type != I4MB
      && (currMB->cbp&15)
      && currSlice->active_pps->transform_8x8_mode_flag);

    if (need_transform_size_flag)
    {
      currSE.type   =  SE_HEADER;
      dP = &(currSlice->partArr[partMap[SE_HEADER]]);

      currMB->luma_transform_size_8x8_flag = dP->bitstream->f(1);
    }

    //=====   DQUANT   =====
    //----------------------
    // Delta quant only if nonzero coeffs
    if (cbp !=0)
    {
      read_delta_quant(&currSE, dP, currMB, partMap, ((currMB->is_intra_block == FALSE)) ? SE_DELTA_QUANT_INTER : SE_DELTA_QUANT_INTRA);

      if (currSlice->dp_mode)
      {
        if ((currMB->is_intra_block == FALSE) && currSlice->dpC_NotPresent ) 
          currMB->dpl_flag = 1;

        if( intra && currSlice->dpB_NotPresent )
        {
          currMB->ei_flag = 1;
          currMB->dpl_flag = 1;
        }

        // check for prediction from neighbours
        check_dp_neighbors (currMB);
        if (currMB->dpl_flag)
        {
          cbp = 0; 
          currMB->cbp = cbp;
        }
      }
    }
  }
  else  // read DC coeffs for new intra modes
  {
    cbp = currMB->cbp;

    read_delta_quant(&currSE, dP, currMB, partMap, SE_DELTA_QUANT_INTRA);

    if (currSlice->dp_mode)
    {  
      if (currSlice->dpB_NotPresent)
      {
        currMB->ei_flag  = 1;
        currMB->dpl_flag = 1;
      }
      check_dp_neighbors (currMB);
      if (currMB->dpl_flag)
      {
        currMB->cbp = cbp = 0; 
      }
    }

    if (!currMB->dpl_flag)
    {
      pos_scan_4x4 = pos_scan4x4[0];

      currSlice->read_coeff_4x4_CAVLC(currMB, LUMA_INTRA16x16DC, 0, 0, levarr, runarr, &numcoeff);

      for(k = 0; k < numcoeff; ++k)
      {
        if (levarr[k] != 0)                     // leave if level == 0
        {
          pos_scan_4x4 += 2 * runarr[k];

          i0 = ((*pos_scan_4x4++) << 2);
          j0 = ((*pos_scan_4x4++) << 2);

          currSlice->cof[0][j0][i0] = levarr[k];// add new intra DC coeff
        }
      }


      if(currMB->is_lossless == FALSE)
        itrans_2(currMB, (ColorPlane) currSlice->colour_plane_id);// transform new intra DC
    }
  }

  update_qp(currMB, currSlice->SliceQpY);

  qp_per = currMB->qp_scaled[currSlice->colour_plane_id] / 6;
  qp_rem = currMB->qp_scaled[currSlice->colour_plane_id] % 6;

  //init quant parameters for chroma 
  for(i=PLANE_U; i <= PLANE_V; ++i)
  {
    qp_per_uv[i] = currMB->qp_scaled[i] / 6;
    qp_rem_uv[i] = currMB->qp_scaled[i] % 6;
  }

  InvLevelScale4x4 = intra? currSlice->InvLevelScale4x4_Intra[currSlice->colour_plane_id][qp_rem] : currSlice->InvLevelScale4x4_Inter[currSlice->colour_plane_id][qp_rem];
  InvLevelScale8x8 = intra? currSlice->InvLevelScale8x8_Intra[currSlice->colour_plane_id][qp_rem] : currSlice->InvLevelScale8x8_Inter[currSlice->colour_plane_id][qp_rem];

  // luma coefficients
  if (cbp)
  {
    if (!currMB->luma_transform_size_8x8_flag) // 4x4 transform
    {
      currMB->read_comp_coeff_4x4_CAVLC (currMB, PLANE_Y, InvLevelScale4x4, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
    }
    else // 8x8 transform
    {
      currMB->read_comp_coeff_8x8_CAVLC (currMB, PLANE_Y, InvLevelScale8x8, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
    }
  }
  else
  {
    memset(p_Vid->nz_coeff[mb_nr][0][0], 0, BLOCK_PIXELS * sizeof(byte));
  }

  for (uv = PLANE_U; uv <= PLANE_V; ++uv )
  {
    /*----------------------16x16DC Luma_Add----------------------*/
    if (IS_I16MB (currMB)) // read DC coeffs for new intra modes       
    {
      if (uv == PLANE_U)
        currSlice->read_coeff_4x4_CAVLC(currMB, CB_INTRA16x16DC, 0, 0, levarr, runarr, &numcoeff);
      else
        currSlice->read_coeff_4x4_CAVLC(currMB, CR_INTRA16x16DC, 0, 0, levarr, runarr, &numcoeff);

      coef_ctr=-1;

      for(k = 0; k < numcoeff; ++k)
      {
        if (levarr[k] != 0)                     // leave if level == 0
        {
          coef_ctr += runarr[k] + 1;

          i0 = pos_scan4x4[coef_ctr][0];
          j0 = pos_scan4x4[coef_ctr][1];
          currSlice->cof[uv][j0<<2][i0<<2] = levarr[k];// add new intra DC coeff
        } //if leavarr[k]
      } //k loop

      if(currMB->is_lossless == FALSE)
      {
        itrans_2(currMB, (ColorPlane) (uv)); // transform new intra DC
      }
    } //IS_I16MB

    update_qp(currMB, currSlice->SliceQpY);

    //init constants for every chroma qp offset
    qp_per_uv[uv] = currMB->qp_scaled[uv] / 6;
    qp_rem_uv[uv] = currMB->qp_scaled[uv] % 6;

    InvLevelScale4x4 = intra? currSlice->InvLevelScale4x4_Intra[uv][qp_rem_uv[uv]] : currSlice->InvLevelScale4x4_Inter[uv][qp_rem_uv[uv]];
    InvLevelScale8x8 = intra? currSlice->InvLevelScale8x8_Intra[uv][qp_rem_uv[uv]] : currSlice->InvLevelScale8x8_Inter[uv][qp_rem_uv[uv]];

    if (!currMB->luma_transform_size_8x8_flag) // 4x4 transform
    {
      currMB->read_comp_coeff_4x4_CAVLC (currMB, (ColorPlane) (uv), InvLevelScale4x4, qp_per_uv[uv], cbp, p_Vid->nz_coeff[mb_nr][uv]);
    }
    else // 8x8 transform
    {
      currMB->read_comp_coeff_8x8_CAVLC (currMB, (ColorPlane) (uv), InvLevelScale8x8, qp_per_uv[uv], cbp, p_Vid->nz_coeff[mb_nr][uv]);
    }   
  }   
}

/*!
 ************************************************************************
 * \brief
 *    Get coded block pattern and coefficients (run/level)
 *    from the NAL
 ************************************************************************
 */
static void read_CBP_and_coeffs_from_NAL_CAVLC_420(Macroblock *currMB)
{
    int mb_nr = currMB->mbAddrX;
    int cbp;
    SyntaxElement currSE;
    VideoParameters *p_Vid = currMB->p_Vid;
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    DataPartition *dP = NULL;
    const byte *partMap = assignSE2partition[currSlice->dp_mode];

    int qp_per, qp_rem;
    int qp_per_uv[2];
    int qp_rem_uv[2];

    // select scan type
    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    const byte *pos_scan_4x4 = pos_scan4x4[0];

    int num_blk8x8_uv;
    if (sps->chroma_format_idc != YUV400)
        num_blk8x8_uv = (1 << sps->chroma_format_idc) & (~(0x1));
    else
        num_blk8x8_uv = 0;
    int num_uv_blocks;
    if (sps->chroma_format_idc != YUV400)
        num_uv_blocks = (((1 << sps->chroma_format_idc) & (~0x1)) >> 1);
    else
        num_uv_blocks = 0;

    // read CBP if not new intra mode
    if (!IS_I16MB(currMB)) {
        //=====   C B P   =====
        //---------------------
        currSE.type = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB) 
                      ? SE_CBP_INTRA : SE_CBP_INTER;
        currSE.mapping = (currMB->mb_type == I4MB || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)
                         ? currSlice->linfo_cbp_intra : currSlice->linfo_cbp_inter;
        dP = &(currSlice->partArr[partMap[currSE.type]]);

        dP->readSyntaxElement(currMB, &currSE, dP);
        currMB->cbp = cbp = currSE.value1;

        //============= Transform size flag for INTER MBs =============
        //-------------------------------------------------------------
        int need_transform_size_flag =
            ((currMB->mb_type >= 1 && currMB->mb_type <= 3) ||
             (IS_DIRECT(currMB) && p_Vid->active_sps->direct_8x8_inference_flag) ||
             currMB->NoMbPartLessThan8x8Flag) &&
            (currMB->mb_type != I8MB) && (currMB->mb_type != I4MB) &&
            (currMB->cbp & 15) && currSlice->active_pps->transform_8x8_mode_flag;

        if (need_transform_size_flag) {
            currSE.type = SE_HEADER;
            dP = &(currSlice->partArr[partMap[SE_HEADER]]);

            currMB->luma_transform_size_8x8_flag = dP->bitstream->f(1);
        }

        //=====   DQUANT   =====
        //----------------------
        // Delta quant only if nonzero coeffs
        if (cbp != 0) {
            int intra = currMB->is_intra_block == TRUE;
            read_delta_quant(&currSE, dP, currMB, partMap,
                currMB->is_intra_block == FALSE ? SE_DELTA_QUANT_INTER : SE_DELTA_QUANT_INTRA);

            if (currSlice->dp_mode) {
                if (currMB->is_intra_block == FALSE && currSlice->dpC_NotPresent)
                    currMB->dpl_flag = 1;
                if (intra && currSlice->dpB_NotPresent) {
                    currMB->ei_flag  = 1;
                    currMB->dpl_flag = 1;
                }
                check_dp_neighbors(currMB);
                if (currMB->dpl_flag)
                    currMB->cbp = cbp = 0;
            }
        }
    } else {
        cbp = currMB->cbp;  
        read_delta_quant(&currSE, dP, currMB, partMap, SE_DELTA_QUANT_INTRA);

        if (currSlice->dp_mode) {
            if (currSlice->dpB_NotPresent) {
                currMB->ei_flag  = 1;
                currMB->dpl_flag = 1;
            }
            check_dp_neighbors(currMB);
            if (currMB->dpl_flag)
                currMB->cbp = cbp = 0; 
        }

        if (!currMB->dpl_flag) {
            int k, i0, j0;
            int levarr[16], runarr[16], numcoeff;
            pos_scan_4x4 = pos_scan4x4[0];

            currSlice->read_coeff_4x4_CAVLC(currMB, LUMA_INTRA16x16DC, 0, 0, levarr, runarr, &numcoeff);

            for (k = 0; k < numcoeff; ++k) {
                if (levarr[k] != 0) { // leave if level == 0
                    pos_scan_4x4 += 2 * runarr[k];

                    i0 = (*pos_scan_4x4++) << 2;
                    j0 = (*pos_scan_4x4++) << 2;

                    currSlice->cof[0][j0][i0] = levarr[k]; // add new intra DC coeff
                }
            }

            if (currMB->is_lossless == FALSE)
                itrans_2(currMB, (ColorPlane)currSlice->colour_plane_id); // transform new intra DC
        }
    }

    update_qp(currMB, currSlice->SliceQpY);

    qp_per = currMB->qp_scaled[currSlice->colour_plane_id] / 6;
    qp_rem = currMB->qp_scaled[currSlice->colour_plane_id] % 6;

    int i;
    //init quant parameters for chroma 
    for (i = 0; i < 2; ++i) {
        qp_per_uv[i] = currMB->qp_scaled[i + 1] / 6;
        qp_rem_uv[i] = currMB->qp_scaled[i + 1] % 6;
    }

    // luma coefficients
    if (cbp) {
        int intra = currMB->is_intra_block == TRUE;
        if (!currMB->luma_transform_size_8x8_flag) { // 4x4 transform
            int (*InvLevelScale4x4)[4] = NULL;
            InvLevelScale4x4 = intra ? currSlice->InvLevelScale4x4_Intra[currSlice->colour_plane_id][qp_rem]
                                     : currSlice->InvLevelScale4x4_Inter[currSlice->colour_plane_id][qp_rem];
            currMB->read_comp_coeff_4x4_CAVLC(currMB, PLANE_Y, InvLevelScale4x4, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
        } else { // 8x8 transform
            int (*InvLevelScale8x8)[8] = NULL;
            InvLevelScale8x8 = intra ? currSlice->InvLevelScale8x8_Intra[currSlice->colour_plane_id][qp_rem]
                                     : currSlice->InvLevelScale8x8_Inter[currSlice->colour_plane_id][qp_rem];
            currMB->read_comp_coeff_8x8_CAVLC(currMB, PLANE_Y, InvLevelScale8x8, qp_per, cbp, p_Vid->nz_coeff[mb_nr][PLANE_Y]);
        }
    } else
        memset(p_Vid->nz_coeff[mb_nr][0][0], 0, BLOCK_PIXELS * sizeof(byte));

    //========================== CHROMA DC ============================
    //-----------------------------------------------------------------
    // chroma DC coeff
    if (cbp > 15) {
        int ll, uv;
        int k;
        int coef_ctr;
        int (*InvLevelScale4x4)[4] = NULL;
        int intra = currMB->is_intra_block == TRUE;
        int levarr[16], runarr[16], numcoeff;
        int cofu[16];

        for (ll = 0; ll < 3; ll += 2) {
            uv = ll >> 1;

            InvLevelScale4x4 = intra ? currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv[uv]]
                                     : currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv[uv]];
            //===================== CHROMA DC YUV420 ======================
            memset(cofu, 0, 4 * sizeof(int));
            coef_ctr = -1;

            currSlice->read_coeff_4x4_CAVLC(currMB, CHROMA_DC, 0, 0, levarr, runarr, &numcoeff);

            for (k = 0; k < numcoeff; ++k) {
                if (levarr[k] != 0) {
                    currMB->s_cbp[0].blk |= 0xf0000 << (ll << 1);
                    coef_ctr += runarr[k] + 1;
                    cofu[coef_ctr] = levarr[k];
                }
            }

            int smb = (p_Vid->type == SP_SLICE && currMB->is_intra_block == FALSE) ||
                      (p_Vid->type == SI_SLICE && currMB->mb_type == SI4MB);
            if (smb || currMB->is_lossless == TRUE) { // check to see if MB type is SPred or SIntra4x4
                currSlice->cof[PLANE_U + uv][0][0] = cofu[0];
                currSlice->cof[PLANE_U + uv][0][4] = cofu[1];
                currSlice->cof[PLANE_U + uv][4][0] = cofu[2];
                currSlice->cof[PLANE_U + uv][4][4] = cofu[3];
            } else {
                int temp[4];
                ihadamard2x2(cofu, temp);

                currSlice->cof[PLANE_U + uv][0][0] = ((temp[0] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
                currSlice->cof[PLANE_U + uv][0][4] = ((temp[1] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
                currSlice->cof[PLANE_U + uv][4][0] = ((temp[2] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
                currSlice->cof[PLANE_U + uv][4][4] = ((temp[3] * InvLevelScale4x4[0][0]) << qp_per_uv[uv]) >> 5;
            }
        }
    }

    //========================== CHROMA AC ============================
    //-----------------------------------------------------------------
    // chroma AC coeff, all zero fram start_scan
    if (cbp <= 31)
        memset(p_Vid->nz_coeff[mb_nr][1][0], 0, 2 * BLOCK_PIXELS * sizeof(byte));
    else {
        int b8, b4;
        int uv;
        int k, i, j;
        int coef_ctr, i0, j0;
        int (*InvLevelScale4x4)[4] = NULL;
        int intra = currMB->is_intra_block == TRUE;
        int levarr[16], runarr[16], numcoeff;

        for (b8 = 0; b8 < num_blk8x8_uv; ++b8) {
            currMB->is_v_block = uv = b8 > (num_uv_blocks - 1);
            if (currMB->is_lossless == FALSE)
                InvLevelScale4x4 = intra ? currSlice->InvLevelScale4x4_Intra[PLANE_U + uv][qp_rem_uv[uv]]
                                         : currSlice->InvLevelScale4x4_Inter[PLANE_U + uv][qp_rem_uv[uv]];

            for (b4 = 0; b4 < 4; ++b4) {
                i = cofuv_blk_x[0][b8][b4];
                j = cofuv_blk_y[0][b8][b4];

                currSlice->read_coeff_4x4_CAVLC(currMB, CHROMA_AC, i + 2*uv, j + 4, levarr, runarr, &numcoeff);
                coef_ctr = 0;

                for (k = 0; k < numcoeff; ++k) {
                    if (levarr[k] != 0) {
                        currMB->s_cbp[0].blk |= i64_power2(cbp_blk_chroma[b8][b4]);
                        coef_ctr += runarr[k] + 1;

                        i0 = pos_scan4x4[coef_ctr][0];
                        j0 = pos_scan4x4[coef_ctr][1];

                        if (currMB->is_lossless == FALSE)
                            currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0])<<qp_per_uv[uv], 4);
                        else
                            currSlice->cof[PLANE_U + uv][(j<<2) + j0][(i<<2) + i0] = levarr[k];
                    }
                }
            }
        }
    }
}

/*!
************************************************************************
* \brief
*    setup coefficient reading functions for CAVLC
*
************************************************************************
*/

void set_read_CBP_and_coeffs_cavlc(Slice *currSlice)
{
    if (currSlice->p_Vid->active_sps->chroma_format_idc == YUV444 &&
        currSlice->p_Vid->active_sps->separate_colour_plane_flag == 0)
        currSlice->read_coeff_4x4_CAVLC = read_coeff_4x4_CAVLC_444;
    else
        currSlice->read_coeff_4x4_CAVLC = read_coeff_4x4_CAVLC;

    switch (currSlice->p_Vid->active_sps->chroma_format_idc) {
    case YUV444:
        if (currSlice->p_Vid->active_sps->separate_colour_plane_flag == 0)
            currSlice->read_CBP_and_coeffs_from_NAL = read_CBP_and_coeffs_from_NAL_CAVLC_444;
        else
            currSlice->read_CBP_and_coeffs_from_NAL = read_CBP_and_coeffs_from_NAL_CAVLC_400;
        break;
    case YUV422:
        currSlice->read_CBP_and_coeffs_from_NAL = read_CBP_and_coeffs_from_NAL_CAVLC_422;
        break;
    case YUV420:
        currSlice->read_CBP_and_coeffs_from_NAL = read_CBP_and_coeffs_from_NAL_CAVLC_420;
        break;
    case YUV400:
        currSlice->read_CBP_and_coeffs_from_NAL = read_CBP_and_coeffs_from_NAL_CAVLC_400;
        break;
    default:
        assert(1);
        currSlice->read_CBP_and_coeffs_from_NAL = NULL;
        break;
    }
}

void set_read_comp_coeff_cavlc(Macroblock *currMB)
{
    if (currMB->is_lossless == FALSE) {
        currMB->read_comp_coeff_4x4_CAVLC = read_comp_coeff_4x4_CAVLC;
        currMB->read_comp_coeff_8x8_CAVLC = read_comp_coeff_8x8_CAVLC;
    } else {
        currMB->read_comp_coeff_4x4_CAVLC = read_comp_coeff_4x4_CAVLC_ls;
        currMB->read_comp_coeff_8x8_CAVLC = read_comp_coeff_8x8_CAVLC_ls;
    }
}
