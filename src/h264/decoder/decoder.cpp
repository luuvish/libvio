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
 *  File      : decoder.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "dpb.h"
#include "slice.h"
#include "macroblock.h"

#include "decoder.h"
#include "intra_prediction.h"
#include "inter_prediction.h"
#include "transform.h"
#include "deblock.h"


namespace vio  {
namespace h264 {


Decoder::Decoder() :
    intra_prediction { new IntraPrediction },
    inter_prediction { new InterPrediction },
    transform        { new Transform       },
    deblock          { new Deblock         }
{
}

Decoder::~Decoder()
{
    delete this->intra_prediction;
    delete this->inter_prediction;
    delete this->transform;
    delete this->deblock;
}


void Decoder::init(slice_t& slice)
{
    this->inter_prediction->init(slice);
    //this->transform->init(slice);
}

void Decoder::assign_quant_params(slice_t& slice)
{
    this->transform->init(slice);
}


void Decoder::decode(mb_t& mb)
{
    slice_t& slice = *mb.p_Slice;
    const sps_t& sps = *slice.active_sps;

    this->inter_prediction->set_chroma_vector(mb);

    this->decode_one_component(mb, PLANE_Y);

    if (sps.ChromaArrayType == 3) {
        this->decode_one_component(mb, PLANE_U);
        this->decode_one_component(mb, PLANE_V);

        slice.parser.is_reset_coeff    = false;
        slice.parser.is_reset_coeff_cr = false;
    }
}

void Decoder::coeff_luma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    this->transform->coeff_luma_dc(mb, pl, x0, y0, runarr, levarr);
}
void Decoder::coeff_luma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    this->transform->coeff_luma_ac(mb, pl, x0, y0, runarr, levarr);
}
void Decoder::coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    this->transform->coeff_chroma_dc(mb, pl, x0, y0, runarr, levarr);
}
void Decoder::coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    this->transform->coeff_chroma_ac(mb, pl, x0, y0, runarr, levarr);
}

void Decoder::transform_luma_dc(mb_t* mb, ColorPlane pl)
{
    this->transform->transform_luma_dc(mb, pl);
}
void Decoder::transform_chroma_dc(mb_t* mb, ColorPlane pl)
{
    this->transform->transform_chroma_dc(mb, pl);
}

void Decoder::deblock_filter(slice_t& slice)
{
    this->deblock->deblock(slice.p_Vid);
}

void Decoder::get_block_luma(storable_picture *curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y, imgpel block[16][16],
                             int shift_x,int maxold_x, int maxold_y, ColorPlane pl, mb_t* mb)
{
    this->inter_prediction->get_block_luma(curr_ref, x_pos, y_pos, block_size_x, block_size_y, block,
                                           shift_x, maxold_x, maxold_y, pl, mb);
}

void Decoder::decode_one_component(mb_t& mb, ColorPlane curr_plane)
{
    const slice_t& slice = *mb.p_Slice;
    const sps_t& sps = *slice.active_sps;

    VideoParameters* p_Vid = slice.p_Vid;

    if (sps.ChromaArrayType == 3) {
        if (!mb.is_intra_block) {
            storable_picture *vidref = p_Vid->no_reference_picture;
            int noref = slice.framepoc < p_Vid->recovery_poc;
            for (int j = 0; j < 6; j++) {
                for (int i = 0; i < slice.listXsize[j] ; i++) {
                    storable_picture *curr_ref = slice.listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = (curr_plane == PLANE_Y) ? curr_ref->imgY : curr_ref->imgUV[curr_plane-1];
                    }
                }
            }
        }
    }

    if (mb.mb_type == I_PCM)
        this->mb_pred_ipcm(mb, curr_plane);
    else if (mb.mb_type == I_16x16 || mb.mb_type == I_4x4 || mb.mb_type == I_8x8)
        this->mb_pred_intra(mb, curr_plane);
    else
        this->mb_pred_inter(mb, curr_plane);
}

void Decoder::mb_pred_ipcm(mb_t& mb, ColorPlane curr_plane)
{
    slice_t& slice = *mb.p_Slice;
    const sps_t& sps = *slice.active_sps;
    storable_picture* dec_picture = slice.dec_picture;

    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16 ; ++j)
            dec_picture->imgY[mb.mb.y * 16 + i][mb.mb.x * 16 + j] = (imgpel) this->transform->cof[0][i][j];
    }

    if (sps.ChromaArrayType != 0) {
        for (int k = 0; k < 2; ++k) {
            for (int i = 0; i < sps.MbHeightC; ++i) {
                for (int j = 0;j < sps.MbWidthC; ++j)
                    dec_picture->imgUV[k][mb.mb.y * sps.MbHeightC + i][mb.mb.x * sps.MbWidthC + j] = (imgpel) this->transform->cof[k + 1][i][j];
            }
        }
    }
}

void Decoder::mb_pred_intra(mb_t& mb, ColorPlane curr_plane)
{
    slice_t& slice = *mb.p_Slice;
    const sps_t& sps = *slice.active_sps;

    int blockoffset = mb.mb_type == I_4x4 ? 1 : mb.mb_type == I_8x8 ? 4 : 16;

    for (int block4x4 = 0; block4x4 < 16; block4x4 += blockoffset) {
        int ioff = ((block4x4 / 4) % 2) * 8 + ((block4x4 % 4) % 2) * 4;
        int joff = ((block4x4 / 4) / 2) * 8 + ((block4x4 % 4) / 2) * 4;

        if (mb.mb_type == I_4x4)
            this->intra_prediction->intra_pred_4x4(&mb, curr_plane, ioff, joff);
        else if (mb.mb_type == I_8x8)
            this->intra_prediction->intra_pred_8x8(&mb, curr_plane, ioff, joff);
        else if (mb.mb_type == I_16x16)
            this->intra_prediction->intra_pred_16x16(&mb, curr_plane, ioff, joff);

        if (mb.mb_type == I_4x4)
            this->transform->inverse_transform_4x4(&mb, curr_plane, ioff, joff);
        else if (mb.mb_type == I_8x8)
            this->transform->inverse_transform_8x8(&mb, curr_plane, ioff, joff);
        else if (mb.mb_type == I_16x16)
            this->transform->inverse_transform_16x16(&mb, curr_plane, ioff, joff);
    }

    if (mb.mb_type == I_16x16 || mb.CodedBlockPatternLuma != 0 || mb.CodedBlockPatternChroma != 0)
        slice.parser.is_reset_coeff = false;

    if (sps.chroma_format_idc != YUV400 && sps.chroma_format_idc != YUV444) {
        this->intra_prediction->intra_pred_chroma(&mb);

        for (int uv = 0; uv < 2; uv++)
            this->transform->inverse_transform_chroma(&mb, (ColorPlane)(uv + 1));
        if (mb.CodedBlockPatternChroma)
            slice.parser.is_reset_coeff_cr = false;
    }
}


//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

void Decoder::mb_pred_inter(mb_t& mb, ColorPlane curr_plane)
{
    slice_t& slice = *mb.p_Slice;
    const sps_t& sps = *slice.active_sps;

    bool b_inter_8x8 = (slice.slice_type == B_slice && mb.mb_type == B_8x8);

    int step_h0 = BLOCK_STEP[mb.mb_type][0];
    int step_v0 = BLOCK_STEP[mb.mb_type][1];
    if (mb.mb_type == PSKIP) {
        step_h0 = slice.slice_type == B_slice ? 2 : 4;
        step_v0 = slice.slice_type == B_slice ? 2 : 4;
    }

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = mb.SubMbType    [block8x8];
            int pred_dir = mb.SubMbPredMode[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];
            if (mv_mode == 0) {
                step_h4 = sps.direct_8x8_inference_flag ? 2 : 1;
                step_v4 = sps.direct_8x8_inference_flag ? 2 : 1;
            }

            if (b_inter_8x8 && slice.direct_spatial_mv_pred_flag) {
                auto mv_info = &slice.dec_picture->mv_info[mb.mb.y * 4 + j0][mb.mb.x * 4 + i0];
                pred_dir = (mv_info->ref_idx[LIST_1] < 0) ? 0 : (mv_info->ref_idx[LIST_0] < 0) ? 1 : 2;
            }

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    this->inter_prediction->perform_mc(&mb, curr_plane, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    if (slice.slice_type == SP_slice)
        this->transform->inverse_transform_sp(&mb, curr_plane);
    else
        this->transform->inverse_transform_inter(&mb, curr_plane);
    if (mb.CodedBlockPatternLuma != 0 || mb.CodedBlockPatternChroma != 0)
        slice.parser.is_reset_coeff = false;
}


}
}
