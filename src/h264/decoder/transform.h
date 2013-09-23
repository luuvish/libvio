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

#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_


namespace vio  {
namespace h264 {


struct transform_t {
    void inverse_luma_dc         (mb_t* mb, ColorPlane pl);
    void inverse_chroma_dc       (mb_t* mb, ColorPlane pl);

    void inverse_transform_4x4   (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void inverse_transform_8x8   (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void inverse_transform_16x16 (mb_t* mb, ColorPlane pl, int ioff, int joff);
    void inverse_transform_chroma(mb_t* mb, ColorPlane pl);

    void inverse_transform_inter (mb_t* mb, ColorPlane pl);
    void inverse_transform_sp    (mb_t* mb, ColorPlane pl);

protected:
    void inverse_4x4(int **d, int **r, int pos_y, int pos_x);
    void inverse_8x8(int **d, int **r, int pos_y, int pos_x);

    void bypass_4x4   (int** r, int** f, int ioff, int joff, uint8_t pred_mode);
    void bypass_8x8   (int** r, int** f, int ioff, int joff, uint8_t pred_mode);
    void bypass_16x16 (int** r, int** f, int ioff, int joff, uint8_t pred_mode);
    void bypass_chroma(int** r, int** f, int nW, int nH, uint8_t pred_mode);

    void itrans_sp   (mb_t* mb, ColorPlane pl, int ioff, int joff);
};


}
}


#endif /* _TRANSFORM_H_ */
