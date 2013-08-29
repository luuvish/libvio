/*!
 *************************************************************************************
 * \file cabac.c
 *
 * \brief
 *    CABAC entropy coding routines
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe
 **************************************************************************************
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "bitstream.h"
#include "bitstream_cabac.h"


// Table 9-44 Specification of rangeTabLPS depending on pStateIdx and qCodIRangeIdx
static const uint8_t rangeTabLPS[64][4] = {
    { 128, 176, 208, 240 }, { 128, 167, 197, 227 }, { 128, 158, 187, 216 }, { 123, 150, 178, 205 },
    { 116, 142, 169, 195 }, { 111, 135, 160, 185 }, { 105, 128, 152, 175 }, { 100, 122, 144, 166 },
    {  95, 116, 137, 158 }, {  90, 110, 130, 150 }, {  85, 104, 123, 142 }, {  81,  99, 117, 135 },
    {  77,  94, 111, 128 }, {  73,  89, 105, 122 }, {  69,  85, 100, 116 }, {  66,  80,  95, 110 },
    {  62,  76,  90, 104 }, {  59,  72,  86,  99 }, {  56,  69,  81,  94 }, {  53,  65,  77,  89 },
    {  51,  62,  73,  85 }, {  48,  59,  69,  80 }, {  46,  56,  66,  76 }, {  43,  53,  63,  72 },
    {  41,  50,  59,  69 }, {  39,  48,  56,  65 }, {  37,  45,  54,  62 }, {  35,  43,  51,  59 },
    {  33,  41,  48,  56 }, {  32,  39,  46,  53 }, {  30,  37,  43,  50 }, {  29,  35,  41,  48 },
    {  27,  33,  39,  45 }, {  26,  31,  37,  43 }, {  24,  30,  35,  41 }, {  23,  28,  33,  39 },
    {  22,  27,  32,  37 }, {  21,  26,  30,  35 }, {  20,  24,  29,  33 }, {  19,  23,  27,  31 },
    {  18,  22,  26,  30 }, {  17,  21,  25,  28 }, {  16,  20,  23,  27 }, {  15,  19,  22,  25 },
    {  14,  18,  21,  24 }, {  14,  17,  20,  23 }, {  13,  16,  19,  22 }, {  12,  15,  18,  21 },
    {  12,  14,  17,  20 }, {  11,  14,  16,  19 }, {  11,  13,  15,  18 }, {  10,  12,  15,  17 },
    {  10,  12,  14,  16 }, {   9,  11,  13,  15 }, {   9,  11,  12,  14 }, {   8,  10,  12,  14 },
    {   8,   9,  11,  13 }, {   7,   9,  11,  12 }, {   7,   9,  10,  12 }, {   7,   8,  10,  11 },
    {   6,   8,   9,  11 }, {   6,   7,   9,  10 }, {   6,   7,   8,   9 }, {   2,   2,   2,   2 }
};

// Table 9-45 State transition table
static const uint8_t transIdxLPS[64] = {
     0,  0,  1,  2,  2,  4,  4,  5,
     6,  7,  8,  9,  9, 11, 11, 12,
    13, 13, 15, 15, 16, 16, 18, 18,
    19, 19, 21, 21, 22, 22, 23, 24,
    24, 25, 26, 26, 27, 27, 28, 29,
    29, 30, 30, 30, 31, 32, 32, 33,
    33, 33, 34, 34, 35, 35, 35, 36,
    36, 36, 37, 37, 37, 38, 38, 63
};

static const uint8_t transIdxMPS[64] = {
     1,  2,  3,  4,  5,  6,  7,  8,
     9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24,
    25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56,
    57, 58, 59, 60, 61, 62, 62, 63
};

static const uint8_t renorm_table_32[32] = {
    6, 5, 4, 4, 3, 3, 3, 3,
    2, 2, 2, 2, 2, 2, 2, 2,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
};


static inline unsigned int getbyte(cabac_engine_t *dep)
{     
    return dep->Dcodestrm[(*dep->Dcodestrm_len)++];
}

static inline unsigned int getword(cabac_engine_t *dep)
{
    int *len = dep->Dcodestrm_len;
    byte *p_code_strm = &dep->Dcodestrm[*len];
    *len += 2;
    return ((*p_code_strm<<8) | *(p_code_strm + 1));
}


void cabac_engine_t::init(unsigned char *code_buffer, int firstbyte, int *code_len)
{
    this->Dcodestrm      = code_buffer;
    this->Dcodestrm_len  = code_len;
    *this->Dcodestrm_len = firstbyte;

    // lookahead of 2 bytes: always make sure that bitstream buffer
    // contains 2 more bytes than actual bitstream
    this->Dvalue = (getbyte(this) << 16) | getword(this);
    this->DbitsLeft = 15;
    this->codIRange = 510;
    //this->codIOffset = read_bits(9);
}

bool cabac_engine_t::decode_decision(cabac_context_t* ctx)
{  
    uint8_t qCodIRangeIdx = (this->codIRange >> 6) & 3;
    uint8_t codIRangeLPS  = rangeTabLPS[ctx->pStateIdx][qCodIRangeIdx];
    uint8_t binVal;

    this->codIRange -= codIRangeLPS;

    if (this->Dvalue < (this->codIRange << this->DbitsLeft)) {
        binVal = ctx->valMPS;
        ctx->pStateIdx = transIdxMPS[ctx->pStateIdx];
    } else {
        binVal = !ctx->valMPS;
        // this->codIOffset -= this->codIRange;
        this->Dvalue -= (this->codIRange << this->DbitsLeft);
        this->codIRange = codIRangeLPS;
        if (ctx->pStateIdx == 0)
            ctx->valMPS = 1 - ctx->valMPS; 
        ctx->pStateIdx = transIdxLPS[ctx->pStateIdx];
    }

    this->renormD();

    return binVal;
}

bool cabac_engine_t::decode_bypass()
{
    uint8_t binVal;

    //this->codIOffset = (this->codIOffset << 1) | read_bits(1);
    if (--(this->DbitsLeft) == 0) {
        // lookahead of 2 bytes: always make sure that bitstream buffer
        // contains 2 more bytes than actual bitstream
        this->Dvalue = (this->Dvalue << 16) | getword(this);
        this->DbitsLeft = 16;
    }

    if (this->Dvalue < (this->codIRange << this->DbitsLeft))
        binVal = 0;
    else {
        binVal = 1;
        //this->codIOffset -= this->codIRange;
        this->Dvalue -= (this->codIRange << this->DbitsLeft);
    }

    return binVal;
}

bool cabac_engine_t::decode_terminate()
{
    uint8_t binVal;

    this->codIRange -= 2;

    if (this->Dvalue < (this->codIRange << this->DbitsLeft)) {
        binVal = 0;
        this->renormD();
    } else
        binVal = 1;

    return binVal;
}

void cabac_engine_t::renormD()
{
    if (this->codIRange < 256) {
        int renorm = renorm_table_32[(this->codIRange >> 3) & 0x1F];
        this->codIRange <<= renorm;
        //this->CodIOffset = (this->CodIOffset << renorm) | read_bits(renorm);
        this->DbitsLeft -= renorm;

        if (this->DbitsLeft <= 0) {
            // lookahead of 2 bytes: always make sure that bitstream buffer
            // contains 2 more bytes than actual bitstream
            this->Dvalue = (this->Dvalue << 16) | getword(this);
            this->DbitsLeft += 16;
        }
    }
}

uint32_t cabac_engine_t::u(cabac_context_t* ctx, int ctx_offset)
{
    uint32_t bins = 0;
    //int bins = -1;

    if (!this->decode_decision(ctx))
        return bins;

    ctx += ctx_offset;
    for (int b = 1; b; bins++)
        b = this->decode_decision(ctx);
    return bins;
}

uint32_t cabac_engine_t::tu(cabac_context_t* ctx, int ctx_offset, unsigned int max_symbol)
{
    unsigned int bins = this->decode_decision(ctx);

    if (bins == 0 || max_symbol == 0)
        return bins;

    unsigned int l;
    ctx += ctx_offset;
    bins = 0;
    do {
        l = this->decode_decision(ctx);
        ++bins;
    } while (l != 0 && bins < max_symbol);

    if (l != 0 && bins == max_symbol)
        ++bins;
    return bins;
}


int cabac_startcode_follows(slice_t *currSlice, int eos_bit)
{
    unsigned int  bit;

    if (eos_bit) {
        Bitstream *bitstream = currSlice->partArr[0].bitstream;
        cabac_engine_t *dep_dp = &bitstream->de_cabac;

        bit = dep_dp->decode_terminate(); //GB
    } else
        bit = 0;

    return (bit == 1 ? 1 : 0);
}

int uvlc_startcode_follows(slice_t *currSlice, int dummy)
{
    Bitstream *bitstream = currSlice->partArr[0].bitstream;

    return !bitstream->more_rbsp_data();
}
