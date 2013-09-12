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


static const byte QP_SCALE_CR[52] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35,
    35, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 39
};

struct transform_t {
	void inverse_luma_dc         (macroblock_t* mb, ColorPlane pl);
	void inverse_chroma_dc       (macroblock_t* mb, ColorPlane pl);

	void inverse_transform_4x4   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void inverse_transform_8x8   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void inverse_transform_16x16 (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void inverse_transform_chroma(macroblock_t* mb, ColorPlane pl);

	void inverse_transform_inter (macroblock_t* mb, ColorPlane pl, int smb);

protected:
	void Inv_Residual_trans_4x4   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void Inv_Residual_trans_8x8   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void Inv_Residual_trans_16x16 (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void Inv_Residual_trans_Chroma(macroblock_t* mb, ColorPlane pl, int ioff, int joff);

	void itrans4x4_ls(macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans4x4   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans8x8   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans16x16 (macroblock_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans_sp   (macroblock_t* mb, ColorPlane pl, int ioff, int joff);

	void iMBtrans4x4(macroblock_t* mb, ColorPlane pl, int smb);
	void iMBtrans8x8(macroblock_t* mb, ColorPlane pl, int smb);
};

extern transform_t transform;


}
}


#endif /* _TRANSFORM_H_ */
