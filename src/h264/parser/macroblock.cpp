#include "global.h"
#include "slice.h"
#include "dpb.h"
#include "data_partition.h"
#include "macroblock.h"
#include "fmo.h"
#include "image.h"
#include "neighbour.h"
#include "intra_prediction.h"
#include "transform.h"


using namespace vio::h264;


void macroblock_t::init(slice_t* slice)
{
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    this->p_Slice = slice;
    this->p_Vid   = slice->p_Vid;
    this->mbAddrX = slice->current_mb_nr;
    // Save the slice number of this macroblock. When the macroblock below
    // is coded it will use this to decide if prediction for above is possible
    this->slice_nr = (short) slice->current_slice_nr;

    /* Update coordinates of the current macroblock */
    if (slice->MbaffFrameFlag) {
        this->mb.x = (this->mbAddrX / 2) % sps->PicWidthInMbs;
        this->mb.y = (this->mbAddrX / 2) / sps->PicWidthInMbs * 2 + (this->mbAddrX % 2);
    } else {
        this->mb.x = this->mbAddrX % sps->PicWidthInMbs;
        this->mb.y = this->mbAddrX / sps->PicWidthInMbs;
    }

    this->is_intra_block          = 0;
    this->mb_skip_flag            = 0;
    this->mb_type                 = 0;
    this->mb_qp_delta             = 0;
    this->intra_chroma_pred_mode  = IntraPrediction::Intra_Chroma_DC;
    this->CodedBlockPatternLuma   = 0;
    this->CodedBlockPatternChroma = 0;

    // Reset syntax element entries in MB struct
    if (slice->slice_type != I_slice) {
        memset(this->mvd_l0, 0, sizeof(this->mvd_l0));
        if (slice->slice_type == B_slice)
            memset(this->mvd_l1, 0, sizeof(this->mvd_l1));
    }

    memset(this->cbp_blks, 0, sizeof(this->cbp_blks));
    memset(this->cbp_bits, 0, sizeof(this->cbp_bits));

    if (!slice->is_reset_coeff) {
        memset(slice->decoder.transform->cof[0][0], 0, 16 * 16 * sizeof(int));
        if (!slice->is_reset_coeff_cr) {
            memset(slice->decoder.transform->cof[1][0], 0, 2 * 16 * 16 * sizeof(int));
            slice->is_reset_coeff_cr = 1;
        }
        slice->is_reset_coeff = 1;
    }

    this->mixedModeEdgeFlag = 0;

    this->update_qp(slice->SliceQpY);

    CheckAvailabilityOfNeighbors(this);

    this->mb_field_decoding_flag = 0;
    if (slice->MbaffFrameFlag) {
        bool prevMbSkipped = (this->mbAddrX % 2 == 1) ?
            slice->mb_data[this->mbAddrX - 1].mb_skip_flag : 0;
        if (this->mbAddrX % 2 == 0 || prevMbSkipped) {
            if (this->mbAvailA)
                this->mb_field_decoding_flag = slice->mb_data[this->mbAddrA].mb_field_decoding_flag;
            else if (this->mbAvailB)
                this->mb_field_decoding_flag = slice->mb_data[this->mbAddrB].mb_field_decoding_flag;
        } else
            this->mb_field_decoding_flag = slice->mb_data[this->mbAddrX - 1].mb_field_decoding_flag;
    }

    if (pps->entropy_coding_mode_flag)
        CheckAvailabilityOfNeighborsCABAC(this);
}

bool macroblock_t::close(slice_t* slice)
{
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    bool eos_bit = (!slice->MbaffFrameFlag || slice->current_mb_nr % 2);
    bool startcode_follows;

    //! The if() statement below resembles the original code, which tested
    //! p_Vid->current_mb_nr == p_Vid->PicSizeInMbs.  Both is, of course, nonsense
    //! In an error prone environment, one can only be sure to have a new
    //! picture by checking the tr of the next slice header!
    ++(slice->num_dec_mb);

    int PicSizeInMbs = sps->PicWidthInMbs * (sps->FrameHeightInMbs / (1 + slice->field_pic_flag));
    if (slice->current_mb_nr == PicSizeInMbs - 1)
        return true;

    slice->current_mb_nr = FmoGetNextMBNr(slice->p_Vid, slice->current_mb_nr);

    if (pps->entropy_coding_mode_flag)
        startcode_follows = eos_bit && slice->partArr[0].de_cabac.decode_terminate();
    else
        startcode_follows = !slice->partArr[0].more_rbsp_data();

    if (slice->current_mb_nr == -1) { // End of slice_t group, MUST be end of slice
        assert(startcode_follows);
        return true;
    }

    if (!startcode_follows)
        return false;

    if (slice->slice_type == I_slice || slice->slice_type == SI_slice)
        return true;
    if (pps->entropy_coding_mode_flag)
        return true;
    if (slice->mb_skip_run <= 0)
        return true;

    return false;
}
