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

	void inverse_transform_inter (mb_t* mb, ColorPlane pl, int smb);

protected:
	void Inv_Residual_trans_4x4   (mb_t* mb, ColorPlane pl, int ioff, int joff);
	void Inv_Residual_trans_8x8   (mb_t* mb, ColorPlane pl, int ioff, int joff);
	void Inv_Residual_trans_16x16 (mb_t* mb, ColorPlane pl, int ioff, int joff);
	void Inv_Residual_trans_Chroma(mb_t* mb, ColorPlane pl, int ioff, int joff);

	void itrans4x4_ls(mb_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans8x8_ls(mb_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans4x4   (mb_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans8x8   (mb_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans16x16 (mb_t* mb, ColorPlane pl, int ioff, int joff);
	void itrans_sp   (mb_t* mb, ColorPlane pl, int ioff, int joff);

	void iMBtrans4x4(mb_t* mb, ColorPlane pl, int smb);
	void iMBtrans8x8(mb_t* mb, ColorPlane pl, int smb);
};

extern transform_t transform;


}
}


#endif /* _TRANSFORM_H_ */
