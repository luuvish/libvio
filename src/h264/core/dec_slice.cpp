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
#include "memalloc.h"
#include "macroblock.h"
#include "neighbour.h"

using vio::h264::mb_t;

#include "erc_api.h"
#include "dpb.h"
#include "ref_list.h"



static void compute_colocated(slice_t* currSlice, storable_picture* listX[6][33])
{
    int i, j;

    VideoParameters *p_Vid = currSlice->p_Vid;

    if (currSlice->direct_spatial_mv_pred_flag == 0) {
        for (j = 0; j < 2 + (currSlice->MbaffFrameFlag * 4); j += 2) {
            for (i = 0; i < currSlice->listXsize[j]; i++) {
                int curPoc = (j == 0) ? p_Vid->dec_picture->poc :
                             (j == 2) ? p_Vid->dec_picture->top_poc :
                                        p_Vid->dec_picture->bottom_poc;
                int pic0 = listX[LIST_0 + j][i]->poc;
                int pic1 = listX[LIST_1 + j][0]->poc;
                int tb = clip3(-128, 127, curPoc - pic0);
                int td = clip3(-128, 127, pic1 - pic0);

                if (td != 0) {
                    int tx = (16384 + abs(td / 2)) / td;
                    int DistScaleFactor = clip3(-1024, 1023, (tb * tx + 32) >> 6);
                    currSlice->mvscale[j][i] = DistScaleFactor;
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
        storable_picture *vidref = p_Vid->no_reference_picture;
        int noref = (currSlice->PicOrderCnt < p_Vid->recovery_poc);
        if (currSlice->colour_plane_id == 0) {
            for (j = 0; j < 6; j++) {
                for (i = 0; i < MAX_LIST_SIZE; i++) {
                    storable_picture *curr_ref = currSlice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                    }
                }
            }
        }
    } else {
        storable_picture *vidref = p_Vid->no_reference_picture;
        int noref = (currSlice->PicOrderCnt < p_Vid->recovery_poc);
        int total_lists = currSlice->MbaffFrameFlag ? 6 :
                          currSlice->slice_type == B_slice ? 2 : 1;
        for (j = 0; j < total_lists; j++) {
            // note that if we always set this to MAX_LIST_SIZE, we avoid crashes with invalid ref_idx being set
            // since currently this is done at the slice level, it seems safe to do so.
            // Note for some reason I get now a mismatch between version 12 and this one in cabac. I wonder why.
            for (i = 0; i < MAX_LIST_SIZE; i++) {
                storable_picture *curr_ref = currSlice->listX[j][i];
                if (curr_ref) {
                    curr_ref->no_ref = noref && (curr_ref == vidref);
                }
            }
        }
    }
}


bool slice_t::init()
{
    VideoParameters *p_Vid = this->p_Vid;
    p_Vid->active_sps = this->active_sps;
    p_Vid->active_pps = this->active_pps;
    int current_header = this->current_header;

    init_ref_lists(this);

    this->parser.init(*this);
    this->decoder.init(*this);


    if (current_header != SOP && current_header != SOS)
        return false;

    if (this->active_sps->separate_colour_plane_flag) {
        p_Vid->mb_data     = p_Vid->mb_data_JV    [this->colour_plane_id];
        p_Vid->dec_picture = p_Vid->dec_picture_JV[this->colour_plane_id];
    }
    this->neighbour.mb_data = p_Vid->mb_data;
    this->dec_picture = p_Vid->dec_picture;

    if (this->slice_type == B_slice)
        compute_colocated(this, this->listX);

    if (this->slice_type != I_slice && this->slice_type != SI_slice)
        init_cur_imgy(this, p_Vid);

    return true;
}

void slice_t::decode()
{
    bool end_of_slice = 0;

    while (!end_of_slice) { // loop over macroblocks
        mb_t& mb = this->neighbour.mb_data[this->parser.current_mb_nr]; 

        // Initializes the current macroblock
        mb.init(*this);
        // Get the syntax elements from the NAL
        this->parser.parse(mb);
        this->decoder.decode(mb);

        if (this->MbaffFrameFlag && mb.mb_field_decoding_flag) {
            this->num_ref_idx_l0_active_minus1 = ((this->num_ref_idx_l0_active_minus1 + 1) >> 1) - 1;
            this->num_ref_idx_l1_active_minus1 = ((this->num_ref_idx_l1_active_minus1 + 1) >> 1) - 1;
        }

#if (DISABLE_ERC == 0)
        ercWriteMBMODEandMV(&mb);
#endif

        end_of_slice = mb.close(*this);

        ++this->num_dec_mb;
    }
}
