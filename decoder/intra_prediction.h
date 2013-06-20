
/*!
 *************************************************************************************
 * \file mb_prediction.h
 *
 * \brief
 *    Functions for macroblock prediction
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>  
 *************************************************************************************
 */

#ifndef _INTRA_PREDICTION_H_
#define _INTRA_PREDICTION_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
pred_l[xO + x, yO + y] = pred4x4_l[x, y]
Intra4x4PredMode[luma4x4BlkIdx]
Intra8x8PredMode[luma8x8BlkIdx]
Intra16x16PredMode
intra_chroma_pred_mode
*/
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

    Intra_Chroma_Vertical = 0,
    Intra_Chroma_Horizontal,
    Intra_Chroma_DC,
    Intra_Chroma_Plane
};

typedef enum {
    VERT_PRED            = 0,
    HOR_PRED             = 1,
    DC_PRED              = 2,
    DIAG_DOWN_LEFT_PRED  = 3,
    DIAG_DOWN_RIGHT_PRED = 4,
    VERT_RIGHT_PRED      = 5,
    HOR_DOWN_PRED        = 6,
    VERT_LEFT_PRED       = 7,
    HOR_UP_PRED          = 8
} I4x4PredModes;

// 16x16 intra prediction modes
typedef enum {
    VERT_PRED_16   = 0,
    HOR_PRED_16    = 1,
    DC_PRED_16     = 2,
    PLANE_16       = 3
} I16x16PredModes;

// 8x8 chroma intra prediction modes
typedef enum {
    DC_PRED_8     =  0,
    HOR_PRED_8    =  1,
    VERT_PRED_8   =  2,
    PLANE_8       =  3
} I8x8PredModes;

void set_intra_prediction_modes(Slice *currSlice);
void update_direct_types(Slice *currSlice);

int mb_pred_intra4x4  (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int mb_pred_intra16x16(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture);
int mb_pred_intra8x8  (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int mb_pred_ipcm      (Macroblock *currMB);

#ifdef __cplusplus
}
#endif

#endif /* _INTRA_PREDICTION_H_ */
