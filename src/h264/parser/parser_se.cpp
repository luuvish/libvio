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
    ctxidx { _mb },
    cavlc { _mb.p_Slice->parser.partArr[0] },
    cabac { _mb.p_Slice->parser.partArr[0].de_cabac },
    contexts { _mb.p_Slice->parser.mot_ctx }
{
}

Parser::SyntaxElement::~SyntaxElement()
{
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
        int ctxIdxInc = ctxidx.mb_skip_flag();
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
        int ctxIdxInc = ctxidx.mb_field_decoding_flag();
        mb_field_decoding_flag = cabac.decode_decision(ctx + ctxIdxInc);
    }

    return mb_field_decoding_flag;
}

uint8_t Parser::SyntaxElement::mb_type()
{
    uint8_t mb_type;

    if (!pps.entropy_coding_mode_flag) {
        mb_type = cavlc.ue();
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
        int ctxIdxInc = ctxidx.mb_type_si_slice();
        mb_type = cabac.decode_decision(ctx + ctxIdxInc);
    }

    if (slice.slice_type == I_slice || mb_type == 1) {
        cabac_context_t* ctx = contexts.mb_type_contexts + 3; // ctxIdxOffset = 3
        int ctxIdxInc = ctxidx.mb_type_i_slice();
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
    uint8_t mb_type = 0;

    cabac_context_t* ctx = contexts.mb_type_contexts; // ctxIdxOffset = 14
    if (!cabac.decode_decision(ctx + 0)) {
        if (!cabac.decode_decision(ctx + 1))
            mb_type += cabac.decode_decision(ctx + 2) * 3;
        else
            mb_type -= cabac.decode_decision(ctx + 3) - 2;
    } else
        mb_type += 5;

    if (mb_type == 5) {
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

    int ctxIdxInc = ctxidx.mb_type_b_slice();

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
        int ctxIdxInc = ctxidx.transform_size_8x8_flag();
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
        uint8_t ctxIdxInc = ctxidx.intra_chroma_pred_mode();
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
        uint8_t ctxIdxInc = ctxidx.ref_idx_l(list, x0, y0);
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
        uint8_t ctxIdxInc = ctxidx.mvd_l(list, x0, y0, comp);
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
                int ctxIdxInc = ctxidx.coded_block_pattern_luma(mb_x, mb_y, coded_block_pattern);

                cabac_context_t* ctx = contexts.cbp_l_contexts;
                if (cabac.decode_decision(ctx + ctxIdxInc))
                    coded_block_pattern += (1 << (mb_y + (mb_x >> 1)));
            }
        }

        if (sps.chroma_format_idc != YUV400 && sps.chroma_format_idc != YUV444) {
            cabac_context_t* ctx = contexts.cbp_c_contexts;
            int inc = ctxidx.coded_block_pattern_chroma();
            int ctxIdxInc[2] = { inc & 3, (inc >> 2) & 7 };
            if (cabac.decode_decision(ctx + ctxIdxInc[0]))
                coded_block_pattern += cabac.decode_decision(ctx + ctxIdxInc[1]) ? 32 : 16;
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

// Table 9-5 coeff_token mapping to TotalCoeff(coeff_token) and TrailingOnes(coeff_token)
static const uint8_t coeff_token_length[5][4][17] = {
    //  0 <= nC < 2
    {{  1,  6,  8,  9, 10, 11, 13, 13, 13, 14, 14, 15, 15, 16, 16, 16, 16 },
     {  0,  2,  6,  8,  9, 10, 11, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16 },
     {  0,  0,  3,  7,  8,  9, 10, 11, 13, 13, 14, 14, 15, 15, 16, 16, 16 },
     {  0,  0,  0,  5,  6,  7,  8,  9, 10, 11, 13, 14, 14, 15, 15, 16, 16 }},
    // 2 <= nC < 4
    {{  2,  6,  6,  7,  8,  8,  9, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14 },
     {  0,  2,  5,  6,  6,  7,  8,  9, 11, 11, 12, 12, 13, 13, 14, 14, 14 },
     {  0,  0,  3,  6,  6,  7,  8,  9, 11, 11, 12, 12, 13, 13, 13, 14, 14 },
     {  0,  0,  0,  4,  4,  5,  6,  6,  7,  9, 11, 11, 12, 13, 13, 13, 14 }},
    // 4 <= nC < 8
    {{  4,  6,  6,  6,  7,  7,  7,  7,  8,  8,  9,  9,  9, 10, 10, 10, 10 },
     {  0,  4,  5,  5,  5,  5,  6,  6,  7,  8,  8,  9,  9,  9, 10, 10, 10 },
     {  0,  0,  4,  5,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9, 10, 10, 10 },
     {  0,  0,  0,  4,  4,  4,  4,  4,  5,  6,  7,  8,  8,  9, 10, 10, 10 }},
    // nC == -1 
    {{  2,  6,  6,  6,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  1,  6,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  3,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }},
    // nC == -2 
    {{  1,  7,  7,  9,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  2,  7,  7,  9, 10, 11, 12, 12,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  3,  7,  7,  9, 10, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  5,  6,  7,  7, 10, 11,  0,  0,  0,  0,  0,  0,  0,  0 }}
};

static const uint8_t coeff_token_code[5][4][17] = {
    // 0 <= nC < 2
    {{  1,  5,  7,  7,  7,  7, 15, 11,  8, 15, 11, 15, 11, 15, 11,  7,  4 },
     {  0,  1,  4,  6,  6,  6,  6, 14, 10, 14, 10, 14, 10,  1, 14, 10,  6 },
     {  0,  0,  1,  5,  5,  5,  5,  5, 13,  9, 13,  9, 13,  9, 13,  9,  5 },
     {  0,  0,  0,  3,  3,  4,  4,  4,  4,  4, 12, 12,  8, 12,  8, 12,  8 }},
    // 2 <= nC < 4
    {{  3, 11,  7,  7,  7,  4,  7, 15, 11, 15, 11,  8, 15, 11,  7,  9,  7 },
     {  0,  2,  7, 10,  6,  6,  6,  6, 14, 10, 14, 10, 14, 10, 11,  8,  6 },
     {  0,  0,  3,  9,  5,  5,  5,  5, 13,  9, 13,  9, 13,  9,  6, 10,  5 },
     {  0,  0,  0,  5,  4,  6,  8,  4,  4,  4, 12,  8, 12, 12,  8,  1,  4 }},
    // 4 <= nC < 8
    {{ 15, 15, 11,  8, 15, 11,  9,  8, 15, 11, 15, 11,  8, 13,  9,  5,  1 },
     {  0, 14, 15, 12, 10,  8, 14, 10, 14, 14, 10, 14, 10,  7, 12,  8,  4 },
     {  0,  0, 13, 14, 11,  9, 13,  9, 13, 10, 13,  9, 13,  9, 11,  7,  3 },
     {  0,  0,  0, 12, 11, 10,  9,  8, 13, 12, 12, 12,  8, 12, 10,  6,  2 }},
    // nC == -1
    {{  1,  7,  4,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  1,  6,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  1,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }},
    // nC == -2
    {{  1, 15, 14,  7,  6,  7,  7,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  1, 13, 12,  5,  6,  6,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  1, 11, 10,  4,  5,  5,  4,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  1,  1,  9,  8,  4,  4,  0,  0,  0,  0,  0,  0,  0,  0 }}
};

// Table 9-7 total_zeros tables for 4x4 blocks with tzVlcIndex 1 to 7
// Table 9-8 total_zeros tables for 4x4 blocks with tzVlcIndex 8 to 15
// Table 9-9 total_zeros tables for chroma DC 2x2 and 2x4 blocks
static const uint8_t total_zeros_length[3][15][16] = {
    // YUV420
    {{ 1, 2, 3, 3 },
     { 1, 2, 2 },
     { 1, 1 }},
    // YUV422
    {{ 1, 3, 3, 4, 4, 4, 5, 5 },
     { 3, 2, 3, 3, 3, 3, 3 },
     { 3, 3, 2, 2, 3, 3 },
     { 3, 2, 2, 2, 3 },
     { 2, 2, 2, 2 },
     { 2, 2, 1 },
     { 1, 1 }},
    // YUV444
    {{ 1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9 },
     { 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6 },
     { 4, 3, 3, 3, 4, 4, 3, 3, 4, 5, 5, 6, 5, 6 },
     { 5, 3, 4, 4, 3, 3, 3, 4, 3, 4, 5, 5, 5 },
     { 4, 4, 4, 3, 3, 3, 3, 3, 4, 5, 4, 5 },
     { 6, 5, 3, 3, 3, 3, 3, 3, 4, 3, 6 },
     { 6, 5, 3, 3, 3, 2, 3, 4, 3, 6 },
     { 6, 4, 5, 3, 2, 2, 3, 3, 6 },
     { 6, 6, 4, 2, 2, 3, 2, 5 },
     { 5, 5, 3, 2, 2, 2, 4 },
     { 4, 4, 3, 3, 1, 3 },
     { 4, 4, 2, 1, 3 },
     { 3, 3, 1, 2 },
     { 2, 2, 1 },
     { 1, 1 }}
};

static const uint8_t total_zeros_code[3][15][16] = {
    // YUV420
    {{ 1, 1 , 1 , 0 },
     { 1, 1 , 0 },
     { 1, 0 }},
    // YUV422
    {{ 1, 2, 3, 2, 3, 1, 1, 0 },
     { 0, 1, 1, 4, 5, 6, 7 },
     { 0, 1, 1, 2, 6, 7 },
     { 6, 0, 1, 2, 7 },
     { 0, 1, 2, 3 },
     { 0, 1, 1 },
     { 0, 1 }},
    // YUV444
    {{ 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1 },
     { 7, 6, 5, 4, 3, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0 },
     { 5, 7, 6, 5, 4, 3, 4, 3, 2, 3, 2, 1, 1, 0 },
     { 3, 7, 5, 4, 6, 5, 4, 3, 3, 2, 2, 1, 0 },
     { 5, 4, 3, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
     { 1, 1, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
     { 1, 1, 5, 4, 3, 3, 2, 1, 1, 0 },
     { 1, 1, 1, 3, 3, 2, 2, 1, 0 },
     { 1, 0, 1, 3, 2, 1, 1, 1 },
     { 1, 0, 1, 3, 2, 1, 1 },
     { 0, 1, 1, 2, 1, 3 },
     { 0, 1, 1, 1, 1 },
     { 0, 1, 1, 1 },
     { 0, 1, 1 },
     { 0, 1 }}
};

// Table 9-10 Tables for run_before
static const uint8_t run_before_length[15][16] = {
    { 1, 1 },
    { 1, 2, 2 },
    { 2, 2, 2, 2 },
    { 2, 2, 2, 3, 3 },
    { 2, 2, 3, 3, 3, 3 },
    { 2, 3, 3, 3, 3, 3, 3 },
    { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
};

static const uint8_t run_before_code[15][16] = {
    { 1, 0 },
    { 1, 1, 0 },
    { 3, 2, 1, 0 },
    { 3, 2, 1, 1, 0 },
    { 3, 2, 3, 2, 1, 0 },
    { 3, 0, 1, 3, 2, 5, 4 },
    { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
};


uint8_t Parser::SyntaxElement::coeff_token(int nC)
{
    data_partition_t* dp = &slice.parser.partArr[slice.parser.dp_mode ? (mb.is_intra_block ? 1 : 2) : 0];

    if (nC >= 8) {
        int code = dp->read_bits(6);
        int TotalCoeff   = (code >> 2);
        int TrailingOnes = (code & 3);
        if (TotalCoeff == 0 && TrailingOnes == 3)
            TrailingOnes = 0;
        else
            TotalCoeff++;
        return (TotalCoeff << 2) | (TrailingOnes);
    }

    int tab = (nC == -2) ? 4 : (nC == -1) ? 3 : (nC < 2) ? 0 : (nC < 4) ? 1 : (nC < 8) ? 2 : 5;

    for (int TrailingOnes = 0; TrailingOnes < 4; TrailingOnes++) {
        for (int TotalCoeff = 0; TotalCoeff < 17; TotalCoeff++) {
            int length = coeff_token_length[tab][TrailingOnes][TotalCoeff];
            int code   = coeff_token_code  [tab][TrailingOnes][TotalCoeff];
            if (length > 0 && dp->next_bits(length) == code) {
                dp->read_bits(length);
                return (TotalCoeff << 2) | (TrailingOnes);
            }
        }
    }

    assert(false);
    return -1;
}

uint8_t Parser::SyntaxElement::total_zeros(int yuv, int tzVlcIndex)
{
    data_partition_t* dp = &slice.parser.partArr[slice.parser.dp_mode ? (mb.is_intra_block ? 1 : 2) : 0];

    int tab = tzVlcIndex - 1;

    for (int total_zeros = 0; total_zeros < 16; total_zeros++) {
        int length = total_zeros_length[yuv][tab][total_zeros];
        int code   = total_zeros_code  [yuv][tab][total_zeros];
        if (length > 0 && dp->next_bits(length) == code) {
            dp->read_bits(length);
            return total_zeros;
        }
    }

    assert(false);
    return -1;
}

uint8_t Parser::SyntaxElement::run_before(uint8_t zerosLeft)
{
    data_partition_t* dp = &slice.parser.partArr[slice.parser.dp_mode ? (mb.is_intra_block ? 1 : 2) : 0];

    int tab = min<int>(zerosLeft, 7) - 1;

    for (int run_before = 0; run_before < 16; run_before++) {
        int length = run_before_length[tab][run_before];
        int code   = run_before_code  [tab][run_before];
        if (length > 0 && dp->next_bits(length) == code) {
            dp->read_bits(length);
            return run_before;
        }
    }

    assert(false);
    return -1;
}


}
}
