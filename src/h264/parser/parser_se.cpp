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
 *  File      : parser_se.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "neighbour.h"

#include "parser.h"

#include <functional>


namespace vio  {
namespace h264 {


Parser::SyntaxElement::SyntaxElement(mb_t& _mb) :
    sps { *_mb.p_Slice->active_sps },
    pps { *_mb.p_Slice->active_pps },
    slice { *_mb.p_Slice },
    mb { _mb },
    cavlc { _mb.p_Slice->parser.partArr[0] },
    cabac { _mb.p_Slice->parser.partArr[0].de_cabac },
    contexts { _mb.p_Slice->parser.mot_ctx }
{
}

Parser::SyntaxElement::~SyntaxElement()
{
}


static int ref_idx_ctxIdxInc(mb_t* mb, uint8_t list, uint8_t x0, uint8_t y0)
{
    slice_t* slice = mb->p_Slice;

    nb_t nbA = slice->neighbour.get_neighbour(slice, false, mb->mbAddrX, {x0 * 4 - 1, y0 * 4});
    nb_t nbB = slice->neighbour.get_neighbour(slice, false, mb->mbAddrX, {x0 * 4, y0 * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb->slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb->slice_nr ? nbB.mb : nullptr;

    int condTermFlagA = 0;
    int condTermFlagB = 0;
    int ctxIdxInc;

#define IS_DIRECT(MB) ((MB)->mb_type == 0 && (slice->slice_type == B_SLICE))
    if (nbA.mb) {
        int b8a = ((nbA.y / 4) & 2) + ((nbA.x / 8) & 1);
        auto mv_info = &slice->dec_picture->mv_info[nbA.y / 4][nbA.x / 4];
        if (!(nbA.mb->mb_type == IPCM || IS_DIRECT(nbA.mb) ||
             (nbA.mb->SubMbType[b8a] == 0 && nbA.mb->SubMbPredMode[b8a] == 2))) {
            if (slice->MbaffFrameFlag && !mb->mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag)
                condTermFlagA = (mv_info->ref_idx[list] > 1 ? 1 : 0);
            else
                condTermFlagA = (mv_info->ref_idx[list] > 0 ? 1 : 0);
        }
    }
    if (nbB.mb) {
        int b8b = ((nbB.y / 4) & 2) + ((nbB.x / 8) & 1);
        auto mv_info = &slice->dec_picture->mv_info[nbB.y / 4][nbB.x / 4];
        if (!(nbB.mb->mb_type == IPCM || IS_DIRECT(nbB.mb) ||
             (nbB.mb->SubMbType[b8b] == 0 && nbB.mb->SubMbPredMode[b8b] == 2))) {
            if (slice->MbaffFrameFlag && !mb->mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag)
                condTermFlagB = (mv_info->ref_idx[list] > 1 ? 1 : 0);
            else
                condTermFlagB = (mv_info->ref_idx[list] > 0 ? 1 : 0);
        }
    }
#undef IS_DIRECT

    ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    return ctxIdxInc;
}

static int mvd_ctxIdxInc(mb_t* mb, uint8_t list, uint8_t x0, uint8_t y0, bool comp)
{
    slice_t* slice = mb->p_Slice;

    nb_t nbA = slice->neighbour.get_neighbour(slice, false, mb->mbAddrX, {x0 * 4 - 1, y0 * 4});
    nb_t nbB = slice->neighbour.get_neighbour(slice, false, mb->mbAddrX, {x0 * 4, y0 * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb->slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb->slice_nr ? nbB.mb : nullptr;

    int absMvdCompA = 0;
    int absMvdCompB = 0;

    if (nbA.mb) {
        auto mvd = list == 0 ? nbA.mb->mvd_l0 : nbA.mb->mvd_l1;
        absMvdCompA = abs(mvd[(nbA.y & 15) / 4][(nbA.x & 15) / 4][comp]);
        if (slice->MbaffFrameFlag && comp) {
            if (!mb->mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag)
                absMvdCompA *= 2;
            else if (mb->mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag)
                absMvdCompA /= 2;
        }
    }
    if (nbB.mb) {
        auto mvd = list == 0 ? nbB.mb->mvd_l0 : nbB.mb->mvd_l1;
        absMvdCompB = abs(mvd[(nbB.y & 15) / 4][(nbB.x & 15) / 4][comp]);
        if (slice->MbaffFrameFlag && comp) {
            if (!mb->mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag)
                absMvdCompB *= 2;
            else if (mb->mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag)
                absMvdCompB /= 2;
        }
    }

    int absMvdSum = absMvdCompA + absMvdCompB;
    int ctxIdxInc = (absMvdSum < 3 ? 0 : absMvdSum <= 32 ? 1 : 2);

    return ctxIdxInc;
}

static int cbp_ctxIdxInc(mb_t* mb, uint8_t x0, uint8_t y0, uint8_t coded_block_pattern)
{
    slice_t* slice = mb->p_Slice;

    nb_t nbA = slice->neighbour.get_neighbour(slice, false, mb->mbAddrX, {x0 * 4 - 1, y0 * 4});
    nb_t nbB = slice->neighbour.get_neighbour(slice, false, mb->mbAddrX, {x0 * 4, y0 * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb->slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb->slice_nr ? nbB.mb : nullptr;

    int cbp_a = 0x3F, cbp_b = 0x3F;
    int cbp_a_idx = 0, cbp_b_idx = 0;
    if (x0 == 0) {
        if (nbA.mb && nbA.mb->mb_type != IPCM) {
            cbp_a = nbA.mb->CodedBlockPatternLuma;
            cbp_a_idx = (((nbA.y & 15) / 4) & ~1) + 1;
        }
    } else {
        cbp_a = coded_block_pattern;
        cbp_a_idx = y0;
    }
    if (y0 == 0) {
        if (nbB.mb && nbB.mb->mb_type != IPCM) {
            cbp_b = nbB.mb->CodedBlockPatternLuma;
            cbp_b_idx = (x0 / 2) + 2;
        }
    } else {
        cbp_b = coded_block_pattern;
        cbp_b_idx = (x0 / 2);
    }

    int condTermFlagA = (cbp_a & (1 << cbp_a_idx)) == 0 ? 1 : 0;
    int condTermFlagB = (cbp_b & (1 << cbp_b_idx)) == 0 ? 1 : 0;
    int ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    return ctxIdxInc;
}


uint32_t Parser::SyntaxElement::mb_skip_run()
{
    uint32_t mb_skip_run = 0;

    if (!pps.entropy_coding_mode_flag)
        mb_skip_run = cavlc.ue();

    return mb_skip_run;
}

bool Parser::SyntaxElement::mb_skip_flag()
{
    bool mb_skip_flag = 0;

    if (pps.entropy_coding_mode_flag) {
        cabac_context_t* ctx = contexts.skip_contexts;

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        int condTermFlagA = mbA && !mbA->mb_skip_flag ? 1 : 0;
        int condTermFlagB = mbB && !mbB->mb_skip_flag ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        mb_skip_flag = cabac.decode_decision(ctx + ctxIdxInc);
    }

    return mb_skip_flag;
}

bool Parser::SyntaxElement::mb_field_decoding_flag()
{
    bool mb_field_decoding_flag;

    if (!pps.entropy_coding_mode_flag)
        mb_field_decoding_flag = cavlc.f(1);
    else {
        cabac_context_t* ctx = contexts.mb_aff_contexts;

        int topMbAddr = slice.MbaffFrameFlag ? mb.mbAddrX & ~1 : mb.mbAddrX;

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, topMbAddr, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, topMbAddr, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        int condTermFlagA = mbA && mbA->mb_field_decoding_flag ? 1 : 0;
        int condTermFlagB = mbB && mbB->mb_field_decoding_flag ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        mb_field_decoding_flag = cabac.decode_decision(ctx + ctxIdxInc);
    }

    return mb_field_decoding_flag;
}

uint8_t Parser::SyntaxElement::mb_type()
{
    uint8_t mb_type;

    if (!pps.entropy_coding_mode_flag) {
        mb_type = cavlc.ue();
        mb_type += (slice.slice_type == P_slice || slice.slice_type == SP_slice) ? 1 : 0;
    } else {
        auto reading =
            slice.slice_type == I_slice || slice.slice_type == SI_slice ? std::mem_fn(&Parser::SyntaxElement::mb_type_i_slice) :
            slice.slice_type == P_slice || slice.slice_type == SP_slice ? std::mem_fn(&Parser::SyntaxElement::mb_type_p_slice) :
                                                                          std::mem_fn(&Parser::SyntaxElement::mb_type_b_slice);
        mb_type = reading(this);
    }

    return mb_type;
}

uint8_t Parser::SyntaxElement::mb_type_i_slice()
{
    uint8_t mb_type = 0;

    if (slice.slice_type == SI_slice) {
        cabac_context_t* ctx = contexts.mb_type_contexts; // ctxIdxOffset = 0

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        int condTermFlagA = mbA && mbA->mb_type != SI4MB ? 1 : 0;
        int condTermFlagB = mbB && mbB->mb_type != SI4MB ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        mb_type = cabac.decode_decision(ctx + ctxIdxInc);
    }

    if (slice.slice_type == I_slice || mb_type == 1) {
        cabac_context_t* ctx = contexts.mb_type_contexts + 3; // ctxIdxOffset = 3

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        int condTermFlagA = mbA && mbA->mb_type != I4MB && mbA->mb_type != I8MB ? 1 : 0;
        int condTermFlagB = mbB && mbB->mb_type != I4MB && mbB->mb_type != I8MB ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        if (cabac.decode_decision(ctx + ctxIdxInc)) {
            if (!cabac.decode_terminate()) {
                mb_type += 1;
                mb_type += cabac.decode_decision(ctx + 3) * 12;
                if (cabac.decode_decision(ctx + 4))
                    mb_type += cabac.decode_decision(ctx + 5) * 4 + 4;
                mb_type += cabac.decode_decision(ctx + 6) * 2;
                mb_type += cabac.decode_decision(ctx + 7);
            } else
                mb_type += 25;
        }
    }

    return mb_type;
}

uint8_t Parser::SyntaxElement::mb_type_p_slice()
{
    uint8_t mb_type = 1;

    cabac_context_t* ctx = contexts.mb_type_contexts; // ctxIdxOffset = 14
    if (!cabac.decode_decision(ctx + 0)) {
        if (!cabac.decode_decision(ctx + 1))
            mb_type += cabac.decode_decision(ctx + 2) * 3;
        else
            mb_type -= cabac.decode_decision(ctx + 3) - 2;
    } else
        mb_type += 5;

    if (mb_type == 6) {
        ctx = contexts.mb_type_contexts + 3; // ctxIdxOffset = 17

        if (cabac.decode_decision(ctx + 0)) {
            if (!cabac.decode_terminate()) {
                mb_type += 1;
                mb_type += cabac.decode_decision(ctx + 1) * 12;
                if (cabac.decode_decision(ctx + 2))
                    mb_type += cabac.decode_decision(ctx + 2) * 4 + 4;
                mb_type += cabac.decode_decision(ctx + 3) * 2;
                mb_type += cabac.decode_decision(ctx + 3);
            } else
                mb_type += 25;
        }
    }

    return mb_type;
}

uint8_t Parser::SyntaxElement::mb_type_b_slice()
{
    uint8_t mb_type = 0;

    mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->mb_type != 0 ? 1 : 0;
    int condTermFlagB = mbB && mbB->mb_type != 0 ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    cabac_context_t* ctx = contexts.mb_type_contexts; // ctxIdxOffset = 27
    if (cabac.decode_decision(ctx + ctxIdxInc)) {
        mb_type = 1;
        if (!cabac.decode_decision(ctx + 3))
            mb_type += cabac.decode_decision(ctx + 5);
        else {
            mb_type += 2;
            if (!cabac.decode_decision(ctx + 4)) {
                mb_type += cabac.decode_decision(ctx + 5) * 4;
                mb_type += cabac.decode_decision(ctx + 5) * 2;
                mb_type += cabac.decode_decision(ctx + 5);
            } else {
                mb_type += 9;
                mb_type += cabac.decode_decision(ctx + 5) * 8; 
                mb_type += cabac.decode_decision(ctx + 5) * 4;
                mb_type += cabac.decode_decision(ctx + 5) * 2;
                if (mb_type < 22)
                    mb_type += cabac.decode_decision(ctx + 5);

                if (mb_type == 22)
                    mb_type = 23;
                else if (mb_type == 24)  
                    mb_type = 11;
                else if (mb_type == 26)  
                    mb_type = 22;
            }
        }
    }

    if (mb_type == 23) {
        ctx = contexts.mb_type_contexts + 5; // ctxIdxOffset = 32
        if (cabac.decode_decision(ctx + 0)) {
            if (!cabac.decode_terminate()) {
                mb_type += 1;
                mb_type += cabac.decode_decision(ctx + 1) * 12;
                if (cabac.decode_decision(ctx + 2))
                    mb_type += cabac.decode_decision(ctx + 2) * 4 + 4;
                mb_type += cabac.decode_decision(ctx + 3) * 2;
                mb_type += cabac.decode_decision(ctx + 3);
            } else
                mb_type += 25;
        }
    }

    return mb_type;
}

uint8_t Parser::SyntaxElement::sub_mb_type()
{
    uint8_t sub_mb_type;

    if (!pps.entropy_coding_mode_flag) 
        sub_mb_type = cavlc.ue();
    else {
        if (slice.slice_type == P_slice || slice.slice_type == SP_slice)
            sub_mb_type = this->sub_mb_type_p_slice();
        else
            sub_mb_type = this->sub_mb_type_b_slice();
    }

    return sub_mb_type;
}

uint8_t Parser::SyntaxElement::sub_mb_type_p_slice()
{
    cabac_context_t* ctx = contexts.b8_type_contexts; // ctxIdxOffset = 21

    uint8_t sub_mb_type = 0;

    if (!cabac.decode_decision(ctx + 0)) {
        sub_mb_type += 1;
        if (cabac.decode_decision(ctx + 1))
            sub_mb_type += cabac.decode_decision(ctx + 2) ? 1 : 2;
    }

    return sub_mb_type;
}

uint8_t Parser::SyntaxElement::sub_mb_type_b_slice()
{
    cabac_context_t* ctx = contexts.b8_type_contexts; // ctxIdxOffset = 36

    uint8_t sub_mb_type = 0;

    if (cabac.decode_decision(ctx)) {
        sub_mb_type += 1;
        if (cabac.decode_decision(ctx + 1)) {
            sub_mb_type += 2;
            if (cabac.decode_decision(ctx + 2)) {
                sub_mb_type += 4;
                if (cabac.decode_decision(ctx + 3))
                    sub_mb_type += 4;
                else
                    sub_mb_type += cabac.decode_decision(ctx + 3) * 2;
            } else
                sub_mb_type += cabac.decode_decision(ctx + 3) * 2;
        }
        sub_mb_type += cabac.decode_decision(ctx + 3);
    }

    return sub_mb_type;
}

bool Parser::SyntaxElement::transform_size_8x8_flag()
{
    bool transform_size_8x8_flag;

    if (!pps.entropy_coding_mode_flag)
        transform_size_8x8_flag = cavlc.f(1);
    else {
        cabac_context_t* ctx = contexts.transform_size_contexts;

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        int condTermFlagA = mbA && mbA->transform_size_8x8_flag ? 1 : 0;
        int condTermFlagB = mbB && mbB->transform_size_8x8_flag ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        transform_size_8x8_flag = cabac.decode_decision(ctx + ctxIdxInc);
    }

    return transform_size_8x8_flag;
}

int8_t Parser::SyntaxElement::intra_pred_mode()
{
    int8_t intra_pred_mode = 0;

    if (!pps.entropy_coding_mode_flag) {
        if (cavlc.f(1))
            intra_pred_mode = -1;
        else
            intra_pred_mode = cavlc.f(3);
    } else {
        cabac_context_t* ctx = contexts.ipr_contexts;

        if (cabac.decode_decision(ctx))
            intra_pred_mode = -1;
        else {
            uint8_t ctxIdxIncs[] = { 0 };
            intra_pred_mode = cabac.fl(ctx + 1, ctxIdxIncs, 0, 7);
        }
    }

    return intra_pred_mode;
}

uint8_t Parser::SyntaxElement::intra_chroma_pred_mode()
{
    uint8_t intra_chroma_pred_mode;

    if (!pps.entropy_coding_mode_flag)
        intra_chroma_pred_mode = cavlc.ue();
    else {
        cabac_context_t* ctx = contexts.cipr_contexts;

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        int condTermFlagA = mbA && mbA->intra_chroma_pred_mode != 0 && mbA->mb_type != IPCM ? 1 : 0;
        int condTermFlagB = mbB && mbB->intra_chroma_pred_mode != 0 && mbB->mb_type != IPCM ? 1 : 0;
        uint8_t ctxIdxInc = condTermFlagA + condTermFlagB;
        uint8_t ctxIdxIncs[] = { ctxIdxInc, 3, 3 };

        intra_chroma_pred_mode = cabac.tu(ctx, ctxIdxIncs, 1, 3);
    }

    return intra_chroma_pred_mode;
}

uint8_t Parser::SyntaxElement::ref_idx_l(uint8_t list, uint8_t x0, uint8_t y0)
{
    bool refidx_present =
        slice.slice_type == B_slice || !mb.allrefzero || mb.mb_type != P8x8;
    int num_ref_idx_active = list == LIST_0 ?
        slice.num_ref_idx_l0_active_minus1 + 1 :
        slice.num_ref_idx_l1_active_minus1 + 1;

    if (!refidx_present || num_ref_idx_active <= 1)
        return 0;

    uint8_t ref_idx = 0;

    if (!pps.entropy_coding_mode_flag) {
        if (num_ref_idx_active == 2)
            ref_idx = 1 - cavlc.f(1);
        else
            ref_idx = cavlc.ue();
    } else {
        cabac_context_t* ctx = contexts.ref_no_contexts;
        uint8_t ctxIdxInc = ref_idx_ctxIdxInc(&mb, list, x0, y0);
        uint8_t ctxIdxIncs[] = { ctxIdxInc, 4, 5 };

        ref_idx = cabac.u(ctx, ctxIdxIncs, 2);
    }

    return ref_idx;
}

int16_t Parser::SyntaxElement::mvd_l(uint8_t list, uint8_t x0, uint8_t y0, uint8_t comp)
{
    int16_t mvd = 0;

    if (!pps.entropy_coding_mode_flag) 
        mvd = cavlc.se();
    else {
        cabac_context_t* ctx = (comp == 0) ? contexts.mvd_x_contexts : contexts.mvd_y_contexts;
        uint8_t ctxIdxInc = mvd_ctxIdxInc(&mb, list, x0, y0, comp);
        uint8_t ctxIdxIncs[] = { ctxIdxInc, 3, 4, 5, 6 };

        mvd = cabac.ueg(ctx, ctxIdxIncs, 4, 9, 3);
    }

    return mvd;
}

uint8_t Parser::SyntaxElement::coded_block_pattern()
{
    uint8_t coded_block_pattern = 0;

    if (!pps.entropy_coding_mode_flag) {
        //! gives CBP value from codeword number, both for intra and inter
        static const uint8_t NCBP[2][48][2] = {
            { { 15,  0 } , {  0,  1 } , {  7,  2 } , { 11,  4 } , { 13,  8 } , { 14,  3 },
              {  3,  5 } , {  5, 10 } , { 10, 12 } , { 12, 15 } , {  1,  7 } , {  2, 11 },
              {  4, 13 } , {  8, 14 } , {  6,  6 } , {  9,  9 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } },
            { { 47,  0 } , { 31, 16 } , { 15,  1 } , {  0,  2 } , { 23,  4 } , { 27,  8 },
              { 29, 32 } , { 30,  3 } , {  7,  5 } , { 11, 10 } , { 13, 12 } , { 14, 15 },
              { 39, 47 } , { 43,  7 } , { 45, 11 } , { 46, 13 } , { 16, 14 } , {  3,  6 },
              {  5,  9 } , { 10, 31 } , { 12, 35 } , { 19, 37 } , { 21, 42 } , { 26, 44 },
              { 28, 33 } , { 35, 34 } , { 37, 36 } , { 42, 40 } , { 44, 39 } , {  1, 43 },
              {  2, 45 } , {  4, 46 } , {  8, 17 } , { 17, 18 } , { 18, 20 } , { 20, 24 },
              { 24, 19 } , {  6, 21 } , {  9, 26 } , { 22, 28 } , { 25, 23 } , { 32, 27 },
              { 33, 29 } , { 34, 30 } , { 36, 22 } , { 40, 25 } , { 38, 38 } , { 41, 41 } }
        };

        bool normal  = (sps.chroma_format_idc == 0 || sps.chroma_format_idc == 3 ? 0 : 1);
        bool inter   = (mb.is_intra_block ? 0 : 1);
        int  cbp_idx = cavlc.ue();
        coded_block_pattern = NCBP[normal][cbp_idx][inter];
    } else {
        for (int mb_y = 0; mb_y < 4; mb_y += 2) {
            for (int mb_x = 0; mb_x < 4; mb_x += 2) {
                int ctxIdxInc = cbp_ctxIdxInc(&mb, mb_x, mb_y, coded_block_pattern);

                cabac_context_t* ctx = contexts.cbp_l_contexts;
                if (cabac.decode_decision(ctx + ctxIdxInc))
                    coded_block_pattern += (1 << (mb_y + (mb_x >> 1)));
            }
        }

        if (sps.chroma_format_idc != YUV400 && sps.chroma_format_idc != YUV444) {
            cabac_context_t* ctx = contexts.cbp_c_contexts;

            mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
            mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
            mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
            mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

            int condTermFlagA = mbA && (mbA->mb_type == IPCM || mbA->CodedBlockPatternChroma) ? 1 : 0;
            int condTermFlagB = mbB && (mbB->mb_type == IPCM || mbB->CodedBlockPatternChroma) ? 1 : 0;
            int ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

            if (cabac.decode_decision(ctx + ctxIdxInc)) {
                condTermFlagA = mbA && (mbA->mb_type == IPCM || mbA->CodedBlockPatternChroma == 2) ? 1 : 0;
                condTermFlagB = mbB && (mbB->mb_type == IPCM || mbB->CodedBlockPatternChroma == 2) ? 1 : 0;
                ctxIdxInc = condTermFlagA + 2 * condTermFlagB + 4;

                coded_block_pattern += cabac.decode_decision(ctx + ctxIdxInc) ? 32 : 16;
            }
        }
    }

    return coded_block_pattern;
}

int8_t Parser::SyntaxElement::mb_qp_delta()
{
    int8_t mb_qp_delta;

    if (!pps.entropy_coding_mode_flag)
        mb_qp_delta = cavlc.se();
    else {
        cabac_context_t* ctx = contexts.delta_qp_contexts;
        uint8_t ctxIdxInc = slice.parser.last_dquant != 0 ? 1 : 0;
        uint8_t ctxIdxIncs[] = { ctxIdxInc, 2, 3 };

        mb_qp_delta = cabac.u(ctx, ctxIdxIncs, 2);
        if (mb_qp_delta & 1)
            mb_qp_delta = ((mb_qp_delta + 1) >> 1);
        else
            mb_qp_delta = -((mb_qp_delta + 1) >> 1);
    }

    return mb_qp_delta;
}


}
}
