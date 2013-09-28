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


namespace vio  {
namespace h264 {


void macroblock_t::init(slice_t& slice)
{
    sps_t& sps = *slice.active_sps;
    mb_t& mb = *this;

    mb.p_Slice = &slice;
    mb.mbAddrX = slice.parser.current_mb_nr;
    // Save the slice number of this macroblock. When the macroblock below
    // is coded it will use this to decide if prediction for above is possible
    mb.slice_nr = (short) slice.current_slice_nr;

    /* Update coordinates of the current macroblock */
    if (slice.MbaffFrameFlag) {
        mb.mb.x = (mb.mbAddrX / 2) % sps.PicWidthInMbs;
        mb.mb.y = (mb.mbAddrX / 2) / sps.PicWidthInMbs * 2 + (mb.mbAddrX % 2);
    } else {
        mb.mb.x = mb.mbAddrX % sps.PicWidthInMbs;
        mb.mb.y = mb.mbAddrX / sps.PicWidthInMbs;
    }

    mb.is_intra_block          = 0;
    mb.mb_skip_flag            = 0;
    mb.mb_type                 = 0;
    mb.mb_qp_delta             = 0;
    mb.intra_chroma_pred_mode  = 0;
    mb.CodedBlockPatternLuma   = 0;
    mb.CodedBlockPatternChroma = 0;

    // Reset syntax element entries in MB struct
    if (slice.slice_type != I_slice) {
        memset(mb.mvd_l0, 0, sizeof(mb.mvd_l0));
        if (slice.slice_type == B_slice)
            memset(mb.mvd_l1, 0, sizeof(mb.mvd_l1));
    }

    memset(mb.cbp_blks, 0, sizeof(mb.cbp_blks));
    memset(mb.cbp_bits, 0, sizeof(mb.cbp_bits));

    if (!slice.parser.is_reset_coeff) {
        memset(slice.decoder.transform->cof[0][0], 0, 16 * 16 * sizeof(int));
        if (!slice.parser.is_reset_coeff_cr) {
            memset(slice.decoder.transform->cof[1][0], 0, 2 * 16 * 16 * sizeof(int));
            slice.parser.is_reset_coeff_cr = 1;
        }
        slice.parser.is_reset_coeff = 1;
    }

    slice.parser.update_qp(mb, slice.SliceQpY);

    mb.mb_field_decoding_flag = 0;
    if (slice.MbaffFrameFlag) {
        bool prevMbSkipped = (mb.mbAddrX % 2 == 1) ?
            slice.mb_data[mb.mbAddrX - 1].mb_skip_flag : 0;
        if (mb.mbAddrX % 2 == 0 || prevMbSkipped) {
            int topMbAddr = mb.mbAddrX & ~1;

            mb_t* mbA = slice.neighbour.get_mb(&slice, false, topMbAddr, {-1, 0});
            mb_t* mbB = slice.neighbour.get_mb(&slice, false, topMbAddr, {0, -1});
            mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
            mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

            if (mbA)
                mb.mb_field_decoding_flag = mbA->mb_field_decoding_flag;
            else if (mbB)
                mb.mb_field_decoding_flag = mbB->mb_field_decoding_flag;
        } else
            mb.mb_field_decoding_flag = slice.mb_data[mb.mbAddrX - 1].mb_field_decoding_flag;
    }
}

bool macroblock_t::close(slice_t& slice)
{
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;

    bool eos_bit = (!slice.MbaffFrameFlag || this->mbAddrX % 2);
    bool startcode_follows;

    //! The if() statement below resembles the original code, which tested
    //! mbAddrX == p_Vid->PicSizeInMbs.  Both is, of course, nonsense
    //! In an error prone environment, one can only be sure to have a new
    //! picture by checking the tr of the next slice header!

    int PicSizeInMbs = sps.PicWidthInMbs * (sps.FrameHeightInMbs / (1 + slice.field_pic_flag));
    if (this->mbAddrX == PicSizeInMbs - 1)
        return true;

    slice.parser.current_mb_nr = FmoGetNextMBNr(slice.p_Vid, slice.parser.current_mb_nr);

    if (pps.entropy_coding_mode_flag)
        startcode_follows = eos_bit && slice.parser.partArr[0].de_cabac.decode_terminate();
    else
        startcode_follows = !slice.parser.partArr[0].more_rbsp_data();

    if (slice.parser.current_mb_nr == -1) { // End of slice_t group, MUST be end of slice
        assert(startcode_follows);
        return true;
    }

    if (!startcode_follows)
        return false;

    if (slice.slice_type == I_slice || slice.slice_type == SI_slice)
        return true;
    if (pps.entropy_coding_mode_flag)
        return true;
    if (slice.parser.mb_skip_run <= 0)
        return true;

    return false;
}

    
}
}
