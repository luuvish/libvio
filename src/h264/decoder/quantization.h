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
 *  File      : quantization.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _QUANTIZATION_H_
#define _QUANTIZATION_H_


namespace vio  {
namespace h264 {


struct quantization_t {
    void assign_quant_params(slice_t* slice);

    void inverse_itrans_2  (mb_t* mb, ColorPlane pl, int** M4);
    void inverse_itrans_420(mb_t* mb, ColorPlane pl, int M4[4]);
    void inverse_itrans_422(mb_t* mb, ColorPlane pl, int** M4);

    void coeff_luma_dc  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void coeff_luma_ac  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);

    int  inverse_quantize(mb_t* mb, bool uv, ColorPlane pl, int i0, int j0, int levarr);

private:
    void set_quant(slice_t* slice);

    int         InvLevelScale4x4_Intra[3][6][4][4];
    int         InvLevelScale4x4_Inter[3][6][4][4];
    int         InvLevelScale8x8_Intra[3][6][8][8];
    int         InvLevelScale8x8_Inter[3][6][8][8];

    int*        qmatrix[12];
};


}
}


#endif /* _QUANTIZATION_H_ */
