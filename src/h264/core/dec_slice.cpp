#include <math.h>
#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "fmo.h"
#include "data_partition.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"

#include "sei.h"
#include "output.h"
#include "neighbour.h"
#include "memalloc.h"
#include "macroblock.h"
#include "mb_read.h"

#include "intra_prediction.h"
#include "deblock.h"

#include "erc_api.h"
#include "dpb.h"



static void fill_wp_params(slice_t *currSlice)
{
    if (currSlice->slice_type == B_SLICE) {
        int i, j, k;
        int comp;
        int log_weight_denom;
        int tb, td;  
        int tx,DistScaleFactor;

        int max_l0_ref = currSlice->num_ref_idx_l0_active_minus1 + 1;
        int max_l1_ref = currSlice->num_ref_idx_l1_active_minus1 + 1;

        if (currSlice->active_pps->weighted_bipred_idc == 2) {
            currSlice->luma_log2_weight_denom = 5;
            currSlice->chroma_log2_weight_denom = 5;

            for (i = 0; i < MAX_REFERENCE_PICTURES; ++i) {
                for (comp = 0; comp < 3; ++comp) {
                    log_weight_denom = (comp == 0) ? currSlice->luma_log2_weight_denom : currSlice->chroma_log2_weight_denom;
                    currSlice->wp_weight[0][i][comp] = 1 << log_weight_denom;
                    currSlice->wp_weight[1][i][comp] = 1 << log_weight_denom;
                    currSlice->wp_offset[0][i][comp] = 0;
                    currSlice->wp_offset[1][i][comp] = 0;
                }
            }
        }

        for (i = 0; i < max_l0_ref; ++i) {
            for (j = 0; j < max_l1_ref; ++j) {
                for (comp = 0; comp < 3; ++comp) {
                    log_weight_denom = (comp == 0) ? currSlice->luma_log2_weight_denom : currSlice->chroma_log2_weight_denom;
                    if (currSlice->active_pps->weighted_bipred_idc == 1) {
                        currSlice->wbp_weight[0][i][j][comp] =  currSlice->wp_weight[0][i][comp];
                        currSlice->wbp_weight[1][i][j][comp] =  currSlice->wp_weight[1][j][comp];
                    } else if (currSlice->active_pps->weighted_bipred_idc == 2) {
                        td = iClip3(-128,127,currSlice->listX[LIST_1][j]->poc - currSlice->listX[LIST_0][i]->poc);
                        if (td == 0 || currSlice->listX[LIST_1][j]->is_long_term || currSlice->listX[LIST_0][i]->is_long_term) {
                            currSlice->wbp_weight[0][i][j][comp] = 32;
                            currSlice->wbp_weight[1][i][j][comp] = 32;
                        } else {
                            tb = iClip3(-128,127,currSlice->ThisPOC - currSlice->listX[LIST_0][i]->poc);

                            tx = (16384 + iabs(td/2))/td;
                            DistScaleFactor = iClip3(-1024, 1023, (tx*tb + 32 )>>6);

                            currSlice->wbp_weight[1][i][j][comp] = DistScaleFactor >> 2;
                            currSlice->wbp_weight[0][i][j][comp] = 64 - currSlice->wbp_weight[1][i][j][comp];
                            if (currSlice->wbp_weight[1][i][j][comp] < -64 || currSlice->wbp_weight[1][i][j][comp] > 128) {
                                currSlice->wbp_weight[0][i][j][comp] = 32;
                                currSlice->wbp_weight[1][i][j][comp] = 32;
                                currSlice->wp_offset[0][i][comp] = 0;
                                currSlice->wp_offset[1][j][comp] = 0;
                            }
                        }
                    }
                }
            }
        }

        if (currSlice->MbaffFrameFlag) {
            for (i = 0; i < 2 * max_l0_ref; ++i) {
                for (j = 0; j < 2 * max_l1_ref; ++j) {
                    for (comp = 0; comp < 3; ++comp) {
                        for (k = 2; k < 6; k += 2) {
                            currSlice->wp_offset[k+0][i][comp] = currSlice->wp_offset[0][i>>1][comp];
                            currSlice->wp_offset[k+1][j][comp] = currSlice->wp_offset[1][j>>1][comp];

                            log_weight_denom = (comp == 0) ? currSlice->luma_log2_weight_denom : currSlice->chroma_log2_weight_denom;
                            if (currSlice->active_pps->weighted_bipred_idc == 1) {
                                currSlice->wbp_weight[k+0][i][j][comp] =  currSlice->wp_weight[0][i>>1][comp];
                                currSlice->wbp_weight[k+1][i][j][comp] =  currSlice->wp_weight[1][j>>1][comp];
                            } else if (currSlice->active_pps->weighted_bipred_idc == 2) {
                                td = iClip3(-128, 127, currSlice->listX[k+LIST_1][j]->poc - currSlice->listX[k+LIST_0][i]->poc);
                                if (td == 0 || currSlice->listX[k+LIST_1][j]->is_long_term || currSlice->listX[k+LIST_0][i]->is_long_term) {
                                    currSlice->wbp_weight[k+0][i][j][comp] =   32;
                                    currSlice->wbp_weight[k+1][i][j][comp] =   32;
                                } else {
                                    tb = iClip3(-128,127,((k==2)?currSlice->TopFieldOrderCnt:currSlice->BottomFieldOrderCnt) - currSlice->listX[k+LIST_0][i]->poc);

                                    tx = (16384 + iabs(td/2))/td;
                                    DistScaleFactor = iClip3(-1024, 1023, (tx*tb + 32 )>>6);

                                    currSlice->wbp_weight[k+1][i][j][comp] = DistScaleFactor >> 2;
                                    currSlice->wbp_weight[k+0][i][j][comp] = 64 - currSlice->wbp_weight[k+1][i][j][comp];
                                    if (currSlice->wbp_weight[k+1][i][j][comp] < -64 || currSlice->wbp_weight[k+1][i][j][comp] > 128) {
                                        currSlice->wbp_weight[k+1][i][j][comp] = 32;
                                        currSlice->wbp_weight[k+0][i][j][comp] = 32;
                                        currSlice->wp_offset[k+0][i][comp] = 0;
                                        currSlice->wp_offset[k+1][j][comp] = 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void compute_colocated(slice_t *currSlice, StorablePicture **listX[6])
{
    int i, j;

    VideoParameters *p_Vid = currSlice->p_Vid;

    if (currSlice->direct_spatial_mv_pred_flag == 0) {
        for (j = 0; j < 2 + (currSlice->MbaffFrameFlag * 4); j += 2) {
            for (i = 0; i < currSlice->listXsize[j]; i++) {
                int prescale, iTRb, iTRp;

                if (j == 0)
                    iTRb = iClip3( -128, 127, p_Vid->dec_picture->poc - listX[LIST_0 + j][i]->poc );
                else if (j == 2)
                    iTRb = iClip3( -128, 127, p_Vid->dec_picture->top_poc - listX[LIST_0 + j][i]->poc );
                else
                    iTRb = iClip3( -128, 127, p_Vid->dec_picture->bottom_poc - listX[LIST_0 + j][i]->poc );

                iTRp = iClip3( -128, 127,  listX[LIST_1 + j][0]->poc - listX[LIST_0 + j][i]->poc);

                if (iTRp != 0) {
                    prescale = ( 16384 + iabs( iTRp / 2 ) ) / iTRp;
                    currSlice->mvscale[j][i] = iClip3( -1024, 1023, ( iTRb * prescale + 32 ) >> 6 ) ;
                } else
                    currSlice->mvscale[j][i] = 9999;
            }
        }
    }
}

// this is intended to make get_block_luma faster by doing this at a more appropriate level
// i.e. per slice rather than per MB
static void init_cur_imgy(slice_t *currSlice, VideoParameters *p_Vid)
{
    int i, j;
    if (currSlice->active_sps->separate_colour_plane_flag != 0) {
        StorablePicture *vidref = p_Vid->no_reference_picture;
        int noref = (currSlice->framepoc < p_Vid->recovery_poc);
        if (currSlice->colour_plane_id == 0) {
            for (j = 0; j < 6; j++) {
                for (i = 0; i < MAX_LIST_SIZE; i++) {
                    StorablePicture *curr_ref = currSlice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = curr_ref->imgY;
                    }
                }
            }
        }
    } else {
        StorablePicture *vidref = p_Vid->no_reference_picture;
        int noref = (currSlice->framepoc < p_Vid->recovery_poc);
        int total_lists = currSlice->MbaffFrameFlag ? 6 :
                          currSlice->slice_type == B_SLICE ? 2 : 1;
        for (j = 0; j < total_lists; j++) {
            // note that if we always set this to MAX_LIST_SIZE, we avoid crashes with invalid ref_idx being set
            // since currently this is done at the slice level, it seems safe to do so.
            // Note for some reason I get now a mismatch between version 12 and this one in cabac. I wonder why.
            for (i = 0; i < MAX_LIST_SIZE; i++) {
                StorablePicture *curr_ref = currSlice->listX[j][i];
                if (curr_ref) {
                    curr_ref->no_ref = noref && (curr_ref == vidref);
                    curr_ref->cur_imgY = curr_ref->imgY;
                }
            }
        }
    }
}


bool slice_t::init()
{
    int i;

    VideoParameters *p_Vid = this->p_Vid;
    p_Vid->active_sps = this->active_sps;
    p_Vid->active_pps = this->active_pps;
    int current_header = this->current_header;

    this->mb_skip_run = -1;

    this->prescan_skip_read = false;
    this->prescan_mb_field_decoding_read = false;

    this->PrevQpY = this->SliceQpY;

    init_lists(this);

#if (MVC_EXTENSION_ENABLE)
    if (this->svc_extension_flag == 0 || this->svc_extension_flag == 1)
        reorder_lists_mvc(this, this->ThisPOC);
    else
        reorder_lists(this);

    if (this->fs_listinterview0) {
        free(this->fs_listinterview0);
        this->fs_listinterview0 = NULL;
    }
    if (this->fs_listinterview1) {
        free(this->fs_listinterview1);
        this->fs_listinterview1 = NULL;
    }
#endif

    if (!this->field_pic_flag)
        init_mbaff_lists(p_Vid, this);

    // update reference flags and set current p_Vid->ref_flag
    if (!(this->redundant_pic_cnt != 0 && p_Vid->previous_frame_num == this->frame_num)) {
        for (i = 16; i > 0; i--)
            this->ref_flag[i] = this->ref_flag[i-1];
    }
    this->ref_flag[0] = this->redundant_pic_cnt == 0 ? p_Vid->Is_primary_correct
                                                               : p_Vid->Is_redundant_correct;

    if (this->active_pps->entropy_coding_mode_flag) {
        this->mot_ctx->init(this->slice_type, this->cabac_init_idc, this->SliceQpY);
        this->last_dquant = 0;
    }

    if ((this->active_pps->weighted_bipred_idc > 0 && this->slice_type == B_SLICE) ||
        (this->active_pps->weighted_pred_flag && this->slice_type != I_SLICE))
        fill_wp_params(this);


    if (current_header != SOP && current_header != SOS)
        return false;

    if (this->active_sps->separate_colour_plane_flag != 0)
        change_plane_JV(p_Vid, this->colour_plane_id, this);
    else {
        this->mb_data     = p_Vid->mb_data;
        this->dec_picture = p_Vid->dec_picture;
        this->intra_block = p_Vid->intra_block;
    }

    if (this->slice_type == B_SLICE)
        compute_colocated(this, this->listX);

    if (this->slice_type != I_SLICE && this->slice_type != SI_SLICE)
        init_cur_imgy(this, p_Vid);

    return true;
}

void slice_t::decode()
{
    bool end_of_slice = 0;

    while (!end_of_slice) { // loop over macroblocks
        mb_t *currMB = &this->mb_data[this->current_mb_nr]; 

        // Initializes the current macroblock
        currMB->init(this);
        // Get the syntax elements from the NAL
        currMB->parse();
        currMB->decode();

        if (this->MbaffFrameFlag && currMB->mb_field_decoding_flag) {
            this->num_ref_idx_l0_active_minus1 = ((this->num_ref_idx_l0_active_minus1 + 1) >> 1) - 1;
            this->num_ref_idx_l1_active_minus1 = ((this->num_ref_idx_l1_active_minus1 + 1) >> 1) - 1;
        }

#if (DISABLE_ERC == 0)
        ercWriteMBMODEandMV(currMB);
#endif

        end_of_slice = currMB->close(this);
    }
}
