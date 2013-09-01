#ifndef _INTRA_PREDICTION_H_
#define _INTRA_PREDICTION_H_


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

struct intra_prediction_t {
    void intra_pred_4x4(mb_t *currMB, ColorPlane pl, int ioff, int joff);
    void intra_pred_8x8(mb_t *currMB, ColorPlane pl, int ioff, int joff);
    void intra_pred_16x16(mb_t *currMB, ColorPlane pl, int ioff, int joff);
    void intra_pred_chroma(mb_t *currMB);
};


#endif /* _INTRA_PREDICTION_H_ */
