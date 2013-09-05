#ifndef _INTER_PREDICTION_H_
#define _INTER_PREDICTION_H_


#include "global.h"
#include "dpb.h"

struct macroblock_t;

void get_block_luma(StorablePicture *curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y, imgpel **block,
                    int shift_x,int maxold_x,int maxold_y, ColorPlane pl, macroblock_t* mb);

void prepare_direct_params(macroblock_t* mb, StorablePicture *dec_picture, MotionVector *pmvl0, MotionVector *pmvl1,char *l0_rFrame, char *l1_rFrame);
void perform_mc           (macroblock_t* mb, ColorPlane pl, StorablePicture *dec_picture, int pred_dir, int i, int j, int block_size_x, int block_size_y);

void get_direct8x8temporal(macroblock_t* mb, StorablePicture *dec_picture, int block8x8);
void get_direct4x4temporal(macroblock_t* mb, StorablePicture *dec_picture, int block8x8);
void get_direct8x8spatial (macroblock_t* mb, StorablePicture *dec_picture);
void get_direct4x4spatial (macroblock_t* mb, StorablePicture *dec_picture);
int  get_inter8x8         (macroblock_t* mb, StorablePicture *dec_picture, int block8x8);

void update_direct_mv_info(macroblock_t* mb);


struct inter_prediction_t {
	void motion_compensation(macroblock_t *mb);
};


#endif /* _INTER_PREDICTION_H_ */
