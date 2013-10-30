/*
 * =============================================================================
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
 * =============================================================================
 *
 *  File      : intra_prediction.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * =============================================================================
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "decoder.h"


namespace vio  {
namespace h264 {


void IntraPrediction::init(slice_t& slice)
{
    this->sets.pic   = slice.dec_picture;
    this->sets.sps   = slice.active_sps;
    this->sets.pps   = slice.active_pps;
    this->sets.slice = &slice;
}

void IntraPrediction::intra_pred_4x4(mb_t& mb, int comp, int xO, int yO)
{
    Intra4x4 samples {this->sets, mb, comp, xO, yO};

    px_t* pred = &this->sets.slice->mb_pred[comp][yO][xO];

    int i4x4 = ((yO / 4) / 2) * 8 + ((yO / 4) % 2) * 2 + ((xO / 4) / 2) * 4 + ((xO / 4) % 2);
    switch (mb.Intra4x4PredMode[i4x4]) {
    case Intra_4x4_Vertical:
        return samples.vertical(pred);
    case Intra_4x4_Horizontal:
        return samples.horizontal(pred);
    case Intra_4x4_DC:
        return samples.dc(pred);
    case Intra_4x4_Diagonal_Down_Left:
        return samples.diagonal_down_left(pred);
    case Intra_4x4_Diagonal_Down_Right:
        return samples.diagonal_down_right(pred);
    case Intra_4x4_Vertical_Right:
        return samples.vertical_right(pred);
    case Intra_4x4_Horizontal_Down:
        return samples.horizontal_down(pred);
    case Intra_4x4_Vertical_Left:
        return samples.vertical_left(pred);
    case Intra_4x4_Horizontal_Up:
        return samples.horizontal_up(pred);
    }
}

void IntraPrediction::intra_pred_8x8(mb_t& mb, int comp, int xO, int yO)
{
    Intra8x8 samples {this->sets, mb, comp, xO, yO};

    px_t* pred = &this->sets.slice->mb_pred[comp][yO][xO];

    int i8x8 = (yO / 8) * 2 + (xO / 8);
    switch (mb.Intra8x8PredMode[i8x8]) {
    case Intra_8x8_Vertical:
        return samples.vertical(pred);
    case Intra_8x8_Horizontal:
        return samples.horizontal(pred);
    case Intra_8x8_DC:
        return samples.dc(pred);
    case Intra_8x8_Diagonal_Down_Left:
        return samples.diagonal_down_left(pred);
    case Intra_8x8_Diagonal_Down_Right:
        return samples.diagonal_down_right(pred);
    case Intra_8x8_Vertical_Right:
        return samples.vertical_right(pred);
    case Intra_8x8_Horizontal_Down:
        return samples.horizontal_down(pred);
    case Intra_8x8_Vertical_Left:
        return samples.vertical_left(pred);
    case Intra_8x8_Horizontal_Up:
        return samples.horizontal_up(pred);
    }
}

void IntraPrediction::intra_pred_16x16(mb_t& mb, int comp)
{
    Intra16x16 samples {this->sets, mb, comp, 0, 0};

    px_t* pred = &this->sets.slice->mb_pred[comp][0][0];

    switch (mb.Intra16x16PredMode) {
    case Intra_16x16_Vertical:
        return samples.vertical(pred);
    case Intra_16x16_Horizontal:
        return samples.horizontal(pred);
    case Intra_16x16_DC:
        return samples.dc(pred);
    case Intra_16x16_Plane:
        return samples.plane(pred);
    }
}

void IntraPrediction::intra_pred_chroma(mb_t& mb, int comp)
{
    Chroma samples {this->sets, mb, comp, 0, 0};

    px_t* pred = &this->sets.slice->mb_pred[comp][0][0];

    switch (mb.intra_chroma_pred_mode) {
    case Intra_Chroma_DC:
        return samples.dc(pred);
    case Intra_Chroma_Horizontal:
        return samples.horizontal(pred);
    case Intra_Chroma_Vertical:
        return samples.vertical(pred);
    case Intra_Chroma_Plane:
        return samples.plane(pred);
    }
}


IntraPrediction::Intra4x4::Intra4x4(const sets_t& _sets, mb_t& mb, int comp, int xO, int yO) :
    sets {_sets}
{
    px_t** img = comp ? this->sets.pic->imgUV[comp - 1] : this->sets.pic->imgY;

    nb_t nbA[4];
    for (int i = 0; i < 4; i++) {
        nbA[i] = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO - 1, yO + i});
        nbA[i].mb = nbA[i].mb && nbA[i].mb->slice_nr == mb.slice_nr ? nbA[i].mb : nullptr;
    }
    nb_t nbB = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO    , yO - 1});
    nb_t nbC = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO + 4, yO - 1});
    nb_t nbD = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO - 1, yO - 1});
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbC.mb = nbC.mb && nbC.mb->slice_nr == mb.slice_nr ? nbC.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    nbC.mb = nbC.mb && !(xO == 4 && (yO == 4 || yO == 12)) ? nbC.mb : nullptr;

    if (this->sets.pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < 4; i++)
            available[0] &= nbA[i].mb && nbA[i].mb->is_intra_block ? 1 : 0;
        available[1] = nbB.mb && nbB.mb->is_intra_block ? 1 : 0;
        available[2] = nbC.mb && nbC.mb->is_intra_block ? 1 : 0;
        available[3] = nbD.mb && nbD.mb->is_intra_block ? 1 : 0;
    } else {
        available[0] = nbA[0].mb ? 1 : 0;
        available[1] = nbB.mb    ? 1 : 0;
        available[2] = nbC.mb    ? 1 : 0;
        available[3] = nbD.mb    ? 1 : 0;
    }

    if (available[3]) {
        px_t *pix = &img[nbD.y][nbD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < 4; y++)
            p(-1, y) = img[nbA[y].y][nbA[y].x];
    }
    if (available[1]) {
        px_t *pix = &img[nbB.y][nbB.x];
        for (int x = 0; x < 4; x++)
            p(x, -1) = pix[x];
        pix = &img[nbC.y][nbC.x - 4];
        for (int x = 4; x < 8; x++)
            p(x, -1) = available[2] ? pix[x] : p(3, -1);
        available[2] = available[1];
    }
}

void IntraPrediction::Intra4x4::vertical(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++)
            pred4x4L(x, y, pred) = p(x, -1);
    }
}

void IntraPrediction::Intra4x4::horizontal(px_t* pred)
{
    assert(this->available[0]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++)
            pred4x4L(x, y, pred) = p(-1, y);
    }
}

void IntraPrediction::Intra4x4::dc(px_t* pred)
{
    bool availA = this->available[0];
    bool availB = this->available[1];

    int sum = 0;
    if (availA || availB) {
        int round = 0 + (availA ? 2 : 0) + (availB ? 2 : 0);
        int shift = 1 + (availA ? 1 : 0) + (availB ? 1 : 0);
        if (availA) {
            for (int y = 0; y < 4; y++)
                sum += p(-1, y);
        }
        if (availB) {
            for (int x = 0; x < 4; x++)
                sum += p(x, -1);
        }
        sum = (sum + round) >> shift;
    } else
        sum = 1 << (this->sets.sps->BitDepthY - 1);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++)
            pred4x4L(x, y, pred) = sum;
    }
}

void IntraPrediction::Intra4x4::diagonal_down_left(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if (x == 3 && y == 3)
                pred4x4L(x, y, pred) = (p(x + y, -1) + 3 * p(x + y + 1, -1) + 2) >> 2;
            else
                pred4x4L(x, y, pred) =
                    (p(x + y, -1) + 2 * p(x + y + 1, -1) + p(x + y + 2, -1) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra4x4::diagonal_down_right(px_t* pred)
{
    assert(this->available[3]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if (x > y)
                pred4x4L(x, y, pred) =
                    (p(x - y - 2, -1) + 2 * p(x - y - 1, -1) + p(x - y, -1) + 2) >> 2;
            else if (x < y)
                pred4x4L(x, y, pred) =
                    (p(-1, y - x - 2) + 2 * p(-1, y - x - 1) + p(-1, y - x) + 2) >> 2;
            else
                pred4x4L(x, y, pred) = (p(0, -1) + 2 * p(-1, -1) + p(-1, 0) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra4x4::vertical_right(px_t* pred)
{
    assert(this->available[3]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            int zVR = 2 * x - y;
            if (zVR >= 0 && (zVR % 2) == 0)
                pred4x4L(x, y, pred) = (p(x - (y >> 1) - 1, -1) + p(x - (y >> 1), -1) + 1) >> 1;
            else if (zVR >= 0 && (zVR % 2) == 1)
                pred4x4L(x, y, pred) = (p(x - (y >> 1) - 2, -1) +
                    2 * p(x - (y >> 1) - 1, -1) + p(x - (y >> 1), -1) + 2) >> 2;
            else if (zVR == -1)
                pred4x4L(x, y, pred) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            else // zVR < -1
                pred4x4L(x, y, pred) = (p(-1, y - 2 * x - 1) +
                    2 * p(-1, y - 2 * x - 2) + p(-1, y - 2 * x - 3) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra4x4::horizontal_down(px_t* pred)
{
    assert(this->available[3]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            int zHD = 2 * y - x;
            if (zHD >= 0 && (zHD % 2) == 0)
                pred4x4L(x, y, pred) = (p(-1, y - (x >> 1) - 1) + p(-1, y - (x >> 1)) + 1) >> 1;
            else if (zHD >= 0 && (zHD % 2) == 1)
                pred4x4L(x, y, pred) = (p(-1, y - (x >> 1) - 2) +
                    2 * p(-1, y - (x >> 1) - 1) + p(-1, y - (x >> 1)) + 2) >> 2;
            else if (zHD == -1)
                pred4x4L(x, y, pred) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            else // zHD < -1
                pred4x4L(x, y, pred) = (p(x - 2 * y - 1, -1) +
                    2 * p(x - 2 * y - 2, -1) + p(x - 2 * y - 3, -1) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra4x4::vertical_left(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if ((y % 2) == 0)
                pred4x4L(x, y, pred) = (p(x + (y >> 1), -1) + p(x + (y >> 1) + 1, -1) + 1) >> 1;
            else // (y % 2) == 1
                pred4x4L(x, y, pred) = (p(x + (y >> 1), -1) +
                    2 * p(x + (y >> 1) + 1, -1) + p(x + (y >> 1) + 2, -1) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra4x4::horizontal_up(px_t* pred)
{
    assert(this->available[0]);

    int max = 4;
    int mHU = (max - 1) * 2 - 1;
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            int zHU = x + 2 * y;
            if (zHU < mHU && (zHU % 2) == 0)
                pred4x4L(x, y, pred) = (p(-1, y + (x >> 1)) + p(-1, y + (x >> 1) + 1) + 1) >> 1;
            else if (zHU < mHU && (zHU % 2) == 1)
                pred4x4L(x, y, pred) = (p(-1, y + (x >> 1)) +
                    2 * p(-1, y + (x >> 1) + 1) + p(-1, y + (x >> 1) + 2) + 2) >> 2;
            else if (zHU == mHU)
                pred4x4L(x, y, pred) = (p(-1, max - 2) + 3 * p(-1, max - 1) + 2) >> 2;
            else // zHU > mHU
                pred4x4L(x, y, pred) = p(-1, max - 1);
        }
    }
}

inline px_t& IntraPrediction::Intra4x4::pred4x4L(int x, int y, px_t* pred)
{
    return pred[y * 16 + x];
}

inline px_t& IntraPrediction::Intra4x4::p(int x, int y)
{
    return this->samples[(y + 1) * 9 + (x + 1)];
}


IntraPrediction::Intra8x8::Intra8x8(const sets_t& _sets, mb_t& mb, int comp, int xO, int yO) :
    sets {_sets}
{
    px_t** img = comp ? this->sets.pic->imgUV[comp - 1] : this->sets.pic->imgY;

    nb_t nbA[8];
    for (int i = 0; i < 8; i++) {
        nbA[i] = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO - 1, yO + i});
        nbA[i].mb = nbA[i].mb && nbA[i].mb->slice_nr == mb.slice_nr ? nbA[i].mb : nullptr;
    }
    nb_t nbB = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO    , yO - 1});
    nb_t nbC = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO + 8, yO - 1});
    nb_t nbD = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO - 1, yO - 1});
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbC.mb = nbC.mb && nbC.mb->slice_nr == mb.slice_nr ? nbC.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    nbC.mb = nbC.mb && !(xO == 8 && yO == 8) ? nbC.mb : nullptr;

    if (this->sets.pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < 8; i++)
            available[0] &= nbA[i].mb && nbA[i].mb->is_intra_block ? 1 : 0;
        available[1] = nbB.mb && nbB.mb->is_intra_block ? 1 : 0;
        available[2] = nbC.mb && nbC.mb->is_intra_block ? 1 : 0;
        available[3] = nbD.mb && nbD.mb->is_intra_block ? 1 : 0;
    } else {
        available[0] = nbA[0].mb ? 1 : 0;
        available[1] = nbB.mb    ? 1 : 0;
        available[2] = nbC.mb    ? 1 : 0;
        available[3] = nbD.mb    ? 1 : 0;
    }

    if (available[3]) {
        px_t *pix = &img[nbD.y][nbD.x];
        po(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < 8; y++)
            po(-1, y) = img[nbA[y].y][nbA[y].x];
    }
    if (available[1]) {
        px_t *pix = &img[nbB.y][nbB.x];
        for (int x = 0; x < 8; x++)
            po(x, -1) = pix[x];
        pix = &img[nbC.y][nbC.x - 8];
        for (int x = 8; x < 16; x++)
            po(x, -1) = available[2] ? pix[x] : po(7, -1);
        available[2] = available[1];
    }

    this->filtering();
}

void IntraPrediction::Intra8x8::filtering()
{
    bool availA = this->available[0];
    bool availB = this->available[1];
    bool availD = this->available[3];

    if (availB) {
        if (availD)
            p(0, -1) = (po(-1, -1) + 2 * po(0, -1) + po(1, -1) + 2) >> 2;
        else
            p(0, -1) = (3 * po(0, -1) + po(1, -1) + 2) >> 2;
        for (int x = 1; x < 15; x++)
            p(x, -1) = (po(x - 1, -1) + 2 * po(x, -1) + po(x + 1, -1) + 2) >> 2;
        p(15, -1) = (po(14, -1) + 3 * po(15, -1) + 2) >> 2;
    }
    if (availD) {
        if (availA && availB)
            p(-1, -1) = (po(0, -1) + 2 * po(-1, -1) + po(-1, 0) + 2) >> 2;
        else if (availB)
            p(-1, -1) = (3 * po(-1, -1) + po(0, -1) + 2) >> 2;
        else if (availA)
            p(-1, -1) = (3 * po(-1, -1) + po(-1, 0) + 2) >> 2;
        else
            p(-1, -1) = po(-1, -1);
    }
    if (availA) {
        if (availD)
            p(-1, 0) = (po(-1, -1) + 2 * po(-1, 0) + po(-1, 1) + 2) >> 2;
        else
            p(-1, 0) = (3 * po(-1, 0) + po(-1, 1) + 2) >> 2;
        for (int y = 1; y < 7; y++)
            p(-1, y) = (po(-1, y - 1) + 2 * po(-1, y) + po(-1, y + 1) + 2) >> 2;
        p(-1, 7) = (po(-1, 6) + 3 * po(-1, 7) + 2) >> 2;
    }
}

void IntraPrediction::Intra8x8::vertical(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++)
            pred8x8L(x, y, pred) = p(x, -1);
    }
}

void IntraPrediction::Intra8x8::horizontal(px_t* pred)
{
    assert(this->available[0]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++)
            pred8x8L(x, y, pred) = p(-1, y);
    }
}

void IntraPrediction::Intra8x8::dc(px_t* pred)
{
    bool availA = this->available[0];
    bool availB = this->available[1];

    int sum = 0;
    if (availA || availB) {
        int round = 0 + (availA ? 4 : 0) + (availB ? 4 : 0);
        int shift = 2 + (availA ? 1 : 0) + (availB ? 1 : 0);
        if (availA) {
            for (int y = 0; y < 8; y++)
                sum += p(-1, y);
        }
        if (availB) {
            for (int x = 0; x < 8; x++)
                sum += p(x, -1);
        }
        sum = (sum + round) >> shift;
    } else
        sum = 1 << (this->sets.sps->BitDepthY - 1);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++)
            pred8x8L(x, y, pred) = sum;
    }
}

void IntraPrediction::Intra8x8::diagonal_down_left(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            if (x == 7 && y == 7)
                pred8x8L(x, y, pred) = (p(x + y, -1) + 3 * p(x + y + 1, -1) + 2) >> 2;
            else
                pred8x8L(x, y, pred) =
                    (p(x + y, -1) + 2 * p(x + y + 1, -1) + p(x + y + 2, -1) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra8x8::diagonal_down_right(px_t* pred)
{
    assert(this->available[3]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            if (x > y)
                pred8x8L(x, y, pred) =
                    (p(x - y - 2, -1) + 2 * p(x - y - 1, -1) + p(x - y, -1) + 2) >> 2;
            else if (x < y)
                pred8x8L(x, y, pred) =
                    (p(-1, y - x - 2) + 2 * p(-1, y - x - 1) + p(-1, y - x) + 2) >> 2;
            else
                pred8x8L(x, y, pred) = (p(0, -1) + 2 * p(-1, -1) + p(-1, 0) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra8x8::vertical_right(px_t* pred)
{
    assert(this->available[3]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            int zVR = 2 * x - y;
            if (zVR >= 0 && (zVR % 2) == 0)
                pred8x8L(x, y, pred) = (p(x - (y >> 1) - 1, -1) + p(x - (y >> 1), -1) + 1) >> 1;
            else if (zVR >= 0 && (zVR % 2) == 1)
                pred8x8L(x, y, pred) = (p(x - (y >> 1) - 2, -1) +
                    2 * p(x - (y >> 1) - 1, -1) + p(x - (y >> 1), -1) + 2) >> 2;
            else if (zVR == -1)
                pred8x8L(x, y, pred) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            else // zVR < -1
                pred8x8L(x, y, pred) = (p(-1, y - 2 * x - 1) +
                    2 * p(-1, y - 2 * x - 2) + p(-1, y - 2 * x - 3) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra8x8::horizontal_down(px_t* pred)
{
    assert(this->available[3]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            int zHD = 2 * y - x;
            if (zHD >= 0 && (zHD % 2) == 0)
                pred8x8L(x, y, pred) = (p(-1, y - (x >> 1) - 1) + p(-1, y - (x >> 1)) + 1) >> 1;
            else if (zHD >= 0 && (zHD % 2) == 1)
                pred8x8L(x, y, pred) = (p(-1, y - (x >> 1) - 2) +
                    2 * p(-1, y - (x >> 1) - 1) + p(-1, y - (x >> 1)) + 2) >> 2;
            else if (zHD == -1)
                pred8x8L(x, y, pred) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            else // zHD < -1
                pred8x8L(x, y, pred) = (p(x - 2 * y - 1, -1) +
                    2 * p(x - 2 * y - 2, -1) + p(x - 2 * y - 3, -1) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra8x8::vertical_left(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            if ((y % 2) == 0)
                pred8x8L(x, y, pred) = (p(x + (y >> 1), -1) + p(x + (y >> 1) + 1, -1) + 1) >> 1;
            else // (y % 2) == 1
                pred8x8L(x, y, pred) = (p(x + (y >> 1), -1) +
                    2 * p(x + (y >> 1) + 1, -1) + p(x + (y >> 1) + 2, -1) + 2) >> 2;
        }
    }
}

void IntraPrediction::Intra8x8::horizontal_up(px_t* pred)
{
    assert(this->available[0]);

    int max = 8;
    int mHU = (max - 1) * 2 - 1;
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            int zHU = x + 2 * y;
            if (zHU < mHU && (zHU % 2) == 0)
                pred8x8L(x, y, pred) = (p(-1, y + (x >> 1)) + p(-1, y + (x >> 1) + 1) + 1) >> 1;
            else if (zHU < mHU && (zHU % 2) == 1)
                pred8x8L(x, y, pred) = (p(-1, y + (x >> 1)) +
                    2 * p(-1, y + (x >> 1) + 1) + p(-1, y + (x >> 1) + 2) + 2) >> 2;
            else if (zHU == mHU)
                pred8x8L(x, y, pred) = (p(-1, max - 2) + 3 * p(-1, max - 1) + 2) >> 2;
            else // zHU > mHU
                pred8x8L(x, y, pred) = p(-1, max - 1);
        }
    }
}

inline px_t& IntraPrediction::Intra8x8::pred8x8L(int x, int y, px_t* pred)
{
    return pred[y * 16 + x];
}

inline px_t& IntraPrediction::Intra8x8::po(int x, int y)
{
    return this->samples_lf[(y + 1) * 17 + (x + 1)];
}

inline px_t& IntraPrediction::Intra8x8::p(int x, int y)
{
    return this->samples[(y + 1) * 17 + (x + 1)];
}


IntraPrediction::Intra16x16::Intra16x16(const sets_t& _sets, mb_t& mb, int comp, int xO, int yO) :
    sets {_sets}
{
    px_t** img = comp ? this->sets.pic->imgUV[comp - 1] : this->sets.pic->imgY;

    nb_t nbA[16];
    for (int i = 0; i < 16; i++) {
        nbA[i] = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO - 1, yO + i});
        nbA[i].mb = nbA[i].mb && nbA[i].mb->slice_nr == mb.slice_nr ? nbA[i].mb : nullptr;
    }
    nb_t nbB = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO    , yO - 1});
    nb_t nbD = this->sets.slice->neighbour.get_neighbour(this->sets.slice, false, mb.mbAddrX, {xO - 1, yO - 1});
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    if (this->sets.pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < 16; i++)
            available[0] &= nbA[i].mb && nbA[i].mb->is_intra_block ? 1 : 0;
        available[1] = nbB.mb && nbB.mb->is_intra_block ? 1 : 0;
        available[2] = 0;
        available[3] = nbD.mb && nbD.mb->is_intra_block ? 1 : 0;
    } else {
        available[0] = nbA[0].mb ? 1 : 0;
        available[1] = nbB.mb    ? 1 : 0;
        available[2] = 0;
        available[3] = nbD.mb    ? 1 : 0;
    }

    if (available[3]) {
        px_t *pix = &img[nbD.y][nbD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < 16; y++)
            p(-1, y) = img[nbA[y].y][nbA[y].x];
    }
    if (available[1]) {
        px_t *pix = &img[nbB.y][nbB.x];
        for (int x = 0; x < 16; x++)
            p(x, -1) = pix[x];
    }
}

void IntraPrediction::Intra16x16::vertical(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++)
            predL(x, y, pred) = p(x, -1);
    }
}

void IntraPrediction::Intra16x16::horizontal(px_t* pred)
{
    assert(this->available[0]);

    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++)
            predL(x, y, pred) = p(-1, y);
    }
}

void IntraPrediction::Intra16x16::dc(px_t* pred)
{
    bool availA = this->available[0];
    bool availB = this->available[1];

    int sum = 0;
    if (availA || availB) {
        int round = 0 + (availA ? 8 : 0) + (availB ? 8 : 0);
        int shift = 3 + (availA ? 1 : 0) + (availB ? 1 : 0);
        if (availA) {
            for (int y = 0; y < 16; y++)
                sum += p(-1, y);
        }
        if (availB) {
            for (int x = 0; x < 16; x++)
                sum += p(x, -1);
        }
        sum = (sum + round) >> shift;
    } else
        sum = 1 << (this->sets.sps->BitDepthY - 1);

    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++)
            predL(x, y, pred) = sum;
    }
}

void IntraPrediction::Intra16x16::plane(px_t* pred)
{
    assert(this->available[3]);

    int H = 0;
    for (int x = 0; x < 8; x++)
        H += (x + 1) * (p(8 + x, -1) - p(6 - x, -1));
    int V = 0;
    for (int y = 0; y < 8; y++)
        V += (y + 1) * (p(-1, 8 + y) - p(-1, 6 - y));

    int a = 16 * (p(-1, 15) + p(15, -1));
    int b = (5 * H + 32) >> 6;
    int c = (5 * V + 32) >> 6;

    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++)
            predL(x, y, pred) = clip3(0, (1 << this->sets.sps->BitDepthY) - 1,
                (a + b * (x - 7) + c * (y - 7) + 16) >> 5);
    }
}

inline px_t& IntraPrediction::Intra16x16::predL(int x, int y, px_t* pred)
{
    return pred[y * 16 + x];
}

inline px_t& IntraPrediction::Intra16x16::p(int x, int y)
{
    return this->samples[(y + 1) * 17 + (x + 1)];
}


IntraPrediction::Chroma::Chroma(const sets_t& _sets, mb_t& mb, int comp, int xO, int yO) :
    sets {_sets}
{
    px_t** img = this->sets.pic->imgUV[comp - 1];

    nb_t nbA[16];
    for (int i = 0; i < this->sets.sps->MbHeightC; ++i) {
        nbA[i] = this->sets.slice->neighbour.get_neighbour(this->sets.slice, true, mb.mbAddrX, {xO - 1, yO + i});
        nbA[i].mb = nbA[i].mb && nbA[i].mb->slice_nr == mb.slice_nr ? nbA[i].mb : nullptr;
    }
    nb_t nbB = this->sets.slice->neighbour.get_neighbour(this->sets.slice, true, mb.mbAddrX, {xO    , yO - 1});
    nb_t nbD = this->sets.slice->neighbour.get_neighbour(this->sets.slice, true, mb.mbAddrX, {xO - 1, yO - 1});
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    if (this->sets.pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < this->sets.sps->MbHeightC / 2; ++i)
            available[0] &= nbA[i].mb && nbA[i].mb->is_intra_block ? 1 : 0;
        available[1] = nbB.mb && nbB.mb->is_intra_block ? 1 : 0;
        available[2] = 1;
        for (int i = this->sets.sps->MbHeightC / 2; i < this->sets.sps->MbHeightC; ++i)
            available[2] &= nbA[i].mb && nbA[i].mb->is_intra_block ? 1 : 0;
        available[3] = nbD.mb && nbD.mb->is_intra_block ? 1 : 0;
    } else {
        available[0] = nbA[0].mb ? 1 : 0;
        available[1] = nbB.mb    ? 1 : 0;
        available[2] = nbA[0].mb ? 1 : 0;
        available[3] = nbD.mb    ? 1 : 0;
    }

    if (available[3]) {
        px_t *pix = &img[nbD.y][nbD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < this->sets.sps->MbHeightC / 2; y++)
            p(-1, y) = img[nbA[y].y][nbA[y].x];
    }
    if (available[2]) {
        for (int y = this->sets.sps->MbHeightC / 2; y < this->sets.sps->MbHeightC; y++)
            p(-1, y) = img[nbA[y].y][nbA[y].x];
    }
    if (available[1]) {
        px_t *pix = &img[nbB.y][nbB.x];
        for (int x = 0; x < this->sets.sps->MbWidthC; x++)
            p(x, -1) = pix[x];
    }
}

void IntraPrediction::Chroma::dc4x4(px_t* pred, bool* available, int xO, int yO)
{
    bool availA = available[0];
    bool availB = available[1];

    int sum = 0;
    if (availA || availB) {
        int round = 0 + (availA ? 2 : 0) + (availB ? 2 : 0);
        int shift = 1 + (availA ? 1 : 0) + (availB ? 1 : 0);
        if (availA) {
            for (int y = 0; y < 4; y++)
                sum += p(-1, y + yO);
        }
        if (availB) {
            for (int x = 0; x < 4; x++)
                sum += p(x + xO, -1);
        }
        sum = (sum + round) >> shift;
    } else
        sum = 1 << (this->sets.sps->BitDepthC - 1);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++)
            predC(x + xO, y + yO, pred) = sum;
    }
}

void IntraPrediction::Chroma::dc(px_t* pred)
{
    for (int chroma4x4BlkIdx = 0;
         chroma4x4BlkIdx < (1 << (this->sets.sps->ChromaArrayType + 1));
         ++chroma4x4BlkIdx) {
        int xO = ((chroma4x4BlkIdx / 4) % (this->sets.sps->MbWidthC / 8)) * 8 +
                 ((chroma4x4BlkIdx % 4) % (this->sets.sps->MbWidthC / 4)) * 4;
        int yO = ((chroma4x4BlkIdx / 4) / (this->sets.sps->MbWidthC / 8)) * 8 +
                 ((chroma4x4BlkIdx % 4) / (this->sets.sps->MbWidthC / 4)) * 4;

        bool avail[4] = { 0 };
        if ((xO == 0 && yO == 0) || (xO > 0 && yO > 0)) {
            avail[0] = yO > 0 ? this->available[2] : this->available[0];
            avail[1] = this->available[1];
        } else if (xO > 0 && yO == 0) {
            avail[0] = this->available[1] ? 0 : this->available[0];
            avail[1] = this->available[1];
        } else if (xO == 0 && yO > 0) {
            avail[0] = this->available[2];
            avail[1] = this->available[2] ? 0 : this->available[1];
        }

        this->dc4x4(pred, avail, xO, yO);
    }
}

void IntraPrediction::Chroma::horizontal(px_t* pred)
{
    assert(this->available[0]);

    for (int y = 0; y < this->sets.sps->MbHeightC; y++) {
        for (int x = 0; x < this->sets.sps->MbWidthC; x++)
            predC(x, y, pred) = p(-1, y);
    }
}

void IntraPrediction::Chroma::vertical(px_t* pred)
{
    assert(this->available[1]);

    for (int y = 0; y < this->sets.sps->MbHeightC; y++) {
        for (int x = 0; x < this->sets.sps->MbWidthC; x++)
            predC(x, y, pred) = p(x, -1);
    }
}

void IntraPrediction::Chroma::plane(px_t* pred)
{
    assert(this->available[3]);

    int xCF = this->sets.sps->ChromaArrayType == 3 ? 4 : 0;
    int yCF = this->sets.sps->ChromaArrayType != 1 ? 4 : 0;

    int H = 0;
    for (int x = 0; x < 4 + xCF; x++)
        H += (x + 1) * (p(4 + xCF + x, -1) - p(2 + xCF - x, -1));
    int V = 0;
    for (int y = 0; y < 4 + yCF; y++)
        V += (y + 1) * (p(-1, 4 + yCF + y) - p(-1, 2 + yCF - y));

    int a = 16 * (p(-1, this->sets.sps->MbHeightC - 1) + p(this->sets.sps->MbWidthC - 1, -1));
    int b = ((34 - 29 * (this->sets.sps->ChromaArrayType == 3 ? 1 : 0)) * H + 32) >> 6;
    int c = ((34 - 29 * (this->sets.sps->ChromaArrayType != 1 ? 1 : 0)) * V + 32) >> 6;

    for (int y = 0; y < this->sets.sps->MbHeightC; y++) {
        for (int x = 0; x < this->sets.sps->MbWidthC; x++)
            predC(x, y, pred) = clip3(0, (1 << this->sets.sps->BitDepthC) - 1,
                (a + b * (x - 3 - xCF) + c * (y - 3 - yCF) + 16) >> 5);
    }
}

inline px_t& IntraPrediction::Chroma::predC(int x, int y, px_t* pred)
{
    return pred[y * 16 + x];
}

inline px_t& IntraPrediction::Chroma::p(int x, int y)
{
    return this->samples[(y + 1) * 17 + (x + 1)];
}


}
}
