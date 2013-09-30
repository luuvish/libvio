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
 *  File      : transform.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _VIO_H264_TRANSFORM_H_
#define _VIO_H264_TRANSFORM_H_


namespace vio  {
namespace h264 {


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

    int*        qmatrix[12];

    int         mb_rres[3][16][16];
    imgpel      mb_rec [3][16][16];
};


}
}


#endif // _VIO_H264_TRANSFORM_H_
