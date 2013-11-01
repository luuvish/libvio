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
 *  File      : decoder.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * =============================================================================
 */

#ifndef _VIO_H264_DECODER_H_
#define _VIO_H264_DECODER_H_


namespace vio  {
namespace h264 {


class IntraPrediction {
public:
    enum {
        Intra_4x4_Vertical = 0,
        Intra_4x4_Horizontal,
        Intra_4x4_DC,
        Intra_4x4_Diagonal_Down_Left,
        Intra_4x4_Diagonal_Down_Right,
        Intra_4x4_Vertical_Right,
        Intra_4x4_Horizontal_Down,
        Intra_4x4_Vertical_Left,
        Intra_4x4_Horizontal_Up,

        Intra_8x8_Vertical = 0,
        Intra_8x8_Horizontal,
        Intra_8x8_DC,
        Intra_8x8_Diagonal_Down_Left,
        Intra_8x8_Diagonal_Down_Right,
        Intra_8x8_Vertical_Right,
        Intra_8x8_Horizontal_Down,
        Intra_8x8_Vertical_Left,
        Intra_8x8_Horizontal_Up,

        Intra_16x16_Vertical = 0,
        Intra_16x16_Horizontal,
        Intra_16x16_DC,
        Intra_16x16_Plane,

        Intra_Chroma_DC = 0,
        Intra_Chroma_Horizontal,
        Intra_Chroma_Vertical,
        Intra_Chroma_Plane
    };

    void init(slice_t& slice);

    void intra_pred_4x4   (mb_t& mb, int comp, int xO, int yO);
    void intra_pred_8x8   (mb_t& mb, int comp, int xO, int yO);
    void intra_pred_16x16 (mb_t& mb, int comp);
    void intra_pred_chroma(mb_t& mb, int comp);

protected:
    struct sets_t {
        storable_picture* pic;
        sps_t*            sps;
        pps_t*            pps;
        slice_t*          slice;
    };

    class Intra4x4 {
    public:
        Intra4x4(const sets_t& sets, mb_t& mb, int comp, int xO, int yO);

        void vertical           (px_t* pred);
        void horizontal         (px_t* pred);
        void dc                 (px_t* pred);
        void diagonal_down_left (px_t* pred);
        void diagonal_down_right(px_t* pred);
        void vertical_right     (px_t* pred);
        void horizontal_down    (px_t* pred);
        void vertical_left      (px_t* pred);
        void horizontal_up      (px_t* pred);

    protected:
        inline px_t& pred4x4L(int x, int y, px_t* pred);
        inline px_t& p       (int x, int y);

    private:
        bool available[4];
        px_t samples[9 * 9];

        const sets_t& sets;
    };

    class Intra8x8 {
    public:
        Intra8x8(const sets_t& sets, mb_t& mb, int comp, int xO, int yO);

        void vertical           (px_t* pred);
        void horizontal         (px_t* pred);
        void dc                 (px_t* pred);
        void diagonal_down_left (px_t* pred);
        void diagonal_down_right(px_t* pred);
        void vertical_right     (px_t* pred);
        void horizontal_down    (px_t* pred);
        void vertical_left      (px_t* pred);
        void horizontal_up      (px_t* pred);

    protected:
        inline px_t& pred8x8L(int x, int y, px_t* pred);
        inline px_t& po      (int x, int y);
        inline px_t& p       (int x, int y);
        void filtering();

    private:
        bool available[4];
        px_t samples_lf[17 * 17];
        px_t samples[17 * 17];

        const sets_t& sets;
    };

    class Intra16x16 {
    public:
        Intra16x16(const sets_t& sets, mb_t& mb, int comp, int xO, int yO);

        void vertical           (px_t* pred);
        void horizontal         (px_t* pred);
        void dc                 (px_t* pred);
        void plane              (px_t* pred);

    protected:
        inline px_t& predL(int x, int y, px_t* pred);
        inline px_t& p    (int x, int y);

    private:
        bool available[4];
        px_t samples[17 * 17];

        const sets_t& sets;
    };

    class Chroma {
    public:
        Chroma(const sets_t& sets, mb_t& mb, int comp, int xO, int yO);

        void dc4x4              (px_t* pred, bool* available, int xO, int yO);
        void dc                 (px_t* pred);
        void horizontal         (px_t* pred);
        void vertical           (px_t* pred);
        void plane              (px_t* pred);

    protected:
        inline px_t& predC(int x, int y, px_t* pred);
        inline px_t& p    (int x, int y);

    private:
        bool available[4];
        px_t samples[17 * 17];

        const sets_t& sets;
    };

private:
    sets_t sets;
};

class InterPrediction {
public:
    void        init(slice_t& slice);

    void        get_block_luma(storable_picture* curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y,
                    px_t block[16][16], int pl, mb_t& mb);

    void        inter_pred(mb_t& mb, int comp, int pred_dir, int i, int j, int block_size_x, int block_size_y);

protected:
    struct sets_t {
        storable_picture* pic;
        sps_t*            sps;
        pps_t*            pps;
        slice_t*          slice;
    };

    void        get_block_chroma(storable_picture* curr_ref, int mvLX[2],
                    int partWidthC, int partHeightC, px_t predPartLXC[16][16], int comp, mb_t& mb);

    void        mc_prediction(px_t* mb_pred,
                    px_t block[16][16], int block_size_y, int block_size_x,
                    mb_t& mb, int pl, short l0_refframe, int pred_dir);
    void        bi_prediction(px_t* mb_pred, 
                    px_t block_l0[16][16], px_t block_l1[16][16], int block_size_y, int block_size_x,
                    mb_t& mb, int pl, short l0_refframe, short l1_refframe);

    void        check_motion_vector_range(mb_t& mb, const mv_t *mv);
    int         CheckVertMV(mb_t *currMB, int vec_y, int block_size_y);

private:
    sets_t sets;
};

class Transform {
public:
    void        init(slice_t& slice);

    pos_t       inverse_scan_luma_dc  (mb_t* mb, int run);
    pos_t       inverse_scan_luma_ac  (mb_t* mb, int run);
    pos_t       inverse_scan_chroma_dc(mb_t* mb, int run);
    pos_t       inverse_scan_chroma_ac(mb_t* mb, int run);

    void        coeff_luma_dc  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_luma_ac  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);

    int         inverse_quantize(mb_t* mb, bool uv, ColorPlane pl, int i0, int j0, int levarr);

    void        transform_luma_dc       (mb_t* mb, ColorPlane pl);
    void        transform_chroma_dc     (mb_t* mb, ColorPlane pl);

    void        inverse_transform_4x4   (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void        inverse_transform_8x8   (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void        inverse_transform_16x16 (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void        inverse_transform_chroma(mb_t* mb, ColorPlane pl);

    void        inverse_transform_inter (mb_t* mb, ColorPlane pl);
    void        inverse_transform_sp    (mb_t* mb, ColorPlane pl);

    int         cof[3][16][16];

private:
    void        set_quant(slice_t& slice);

    void        ihadamard_2x2(int c[2][2], int f[2][2]);
    void        ihadamard_2x4(int c[4][2], int f[4][2]);
    void        ihadamard_4x4(int c[4][4], int f[4][4]);
    void        forward_4x4  (int p[16][16], int c[16][16], int pos_y, int pos_x);
    void        inverse_4x4  (int d[16][16], int r[16][16], int pos_y, int pos_x);
    void        inverse_8x8  (int d[16][16], int r[16][16], int pos_y, int pos_x);

    void        bypass_4x4   (int r[16][16], int f[16][16], int ioff, int joff, uint8_t pred_mode);
    void        bypass_8x8   (int r[16][16], int f[16][16], int ioff, int joff, uint8_t pred_mode);
    void        bypass_16x16 (int r[16][16], int f[16][16], int ioff, int joff, uint8_t pred_mode);
    void        bypass_chroma(int r[16][16], int f[16][16], int nW, int nH, uint8_t pred_mode);

    void        itrans_sp   (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void        itrans_sp_cr(mb_t* mb, ColorPlane pl);

    void        construction       (mb_t* mb, ColorPlane pl, int ioff, int joff, int nW, int nH);
    void        construction_16x16 (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void        construction_chroma(mb_t* mb, ColorPlane pl, int ioff, int joff);

    int         InvLevelScale4x4_Intra[3][6][4][4];
    int         InvLevelScale4x4_Inter[3][6][4][4];
    int         InvLevelScale8x8_Intra[3][6][8][8];
    int         InvLevelScale8x8_Inter[3][6][8][8];

    const int*  qmatrix[12];

    int         mb_rres[3][16][16];
    px_t        mb_rec [3][16][16];
};

class Deblock {
public:
    void init();
    void deblock(VideoParameters* p_Vid);

private:
    int  compare_mvs(const mv_t* mv0, const mv_t* mv1, int mvlimit);
    int  bs_compare_mvs(const pic_motion_params* mv_info_p, const pic_motion_params* mv_info_q, int mvlimit);

    void strength_vertical  (mb_t* MbQ, int edge);
    void strength_horizontal(mb_t* MbQ, int edge);
    void strength           (mb_t* mb);

    void filter_strong(px_t *pixQ, int width, int alpha, int beta, int bS, bool chromaStyleFilteringFlag);
    void filter_normal(px_t *pixQ, int width, int alpha, int beta, int bS, bool chromaStyleFilteringFlag, int tc0, int BitDepth);
    void filter_edge  (mb_t* MbQ, bool chromaEdgeFlag, ColorPlane pl, bool verticalEdgeFlag, bool fieldModeInFrameFilteringFlag, int edge);

    void filter_vertical  (mb_t* MbQ);
    void filter_horizontal(mb_t* MbQ);

    void init_neighbors       (VideoParameters *p_Vid);
    void make_frame_picture_JV(VideoParameters *p_Vid);
    void deblock_pic          (VideoParameters *p_Vid);
};


class Decoder {
public:
    Decoder();
    ~Decoder();

    void        init(slice_t& slice);

    void        assign_quant_params(slice_t& slice);

    void        decode(mb_t& mb);

    void        coeff_luma_dc  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_luma_ac  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);

    void        transform_luma_dc  (mb_t* mb, ColorPlane pl);
    void        transform_chroma_dc(mb_t* mb, ColorPlane pl);

    void        deblock_filter(slice_t& slice);

    // called in erc_do_p.cpp
    void        get_block_luma(storable_picture *curr_ref, int x_pos, int y_pos,
                               int block_size_x, int block_size_y, px_t block[16][16],
                               int pl, mb_t& mb);

protected:
    void        decode_one_component(mb_t& mb, ColorPlane curr_plane);
    void        mb_pred_ipcm        (mb_t& mb, ColorPlane curr_plane);
    void        mb_pred_intra       (mb_t& mb, ColorPlane curr_plane);
    void        mb_pred_inter       (mb_t& mb, ColorPlane curr_plane);

public:
    IntraPrediction* intra_prediction;
    InterPrediction* inter_prediction;
    Transform*       transform;
    Deblock*         deblock;
};


}
}


#endif // _VIO_H264_DECODER_H_
