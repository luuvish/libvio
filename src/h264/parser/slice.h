/*!
 *************************************************************************************
 * \file header.h
 *
 * \brief
 *    Prototypes for header.c
 *************************************************************************************
 */

#ifndef _SLICE_H_
#define _SLICE_H_

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_NUM_REF_IDX 32

struct cabac_contexts_t;
struct macroblock_t;


enum {
    P_slice  = 0, // or 5
    B_slice  = 1, // or 6
    I_slice  = 2, // or 7
    SP_slice = 3, // or 8
    SI_slice = 4  // or 9
};


/*! Buffer structure for decoded referenc picture marking commands */
typedef struct DecRefPicMarking_s
{
  int                        memory_management_control_operation;
  int                        difference_of_pic_nums_minus1;
  int                        long_term_pic_num;
  int                        long_term_frame_idx;
  int                        max_long_term_frame_idx_plus1;
  struct DecRefPicMarking_s *Next;
} DecRefPicMarking_t;

//! slice_t
typedef struct slice_t {
    struct video_par         *p_Vid;
    struct inp_par           *p_Inp;
    pps_t *active_pps;
    sps_t *active_sps;
    int                       svc_extension_flag;

    // dpb pointer
    struct decoded_picture_buffer_t *p_Dpb;


    //slice property;
    int         idr_flag;
    int         nal_ref_idc;
    uint8_t     nal_unit_type;


    uint32_t    first_mb_in_slice;                                    // ue(v)
    uint8_t     slice_type;                                           // ue(v)
    uint8_t     pic_parameter_set_id;                                 // ue(v)
    uint8_t     colour_plane_id;                                      // u(2)
    uint32_t    frame_num;                                            // u(v)
    bool        field_pic_flag;                                       // u(1)
    bool        bottom_field_flag;                                    // u(1)
    uint16_t    idr_pic_id;                                           // ue(v)
    uint32_t    pic_order_cnt_lsb;                                    // u(v)
    int32_t     delta_pic_order_cnt_bottom;                           // se(v)
    int32_t     delta_pic_order_cnt[2];                               // se(v)
    uint8_t     redundant_pic_cnt;                                    // ue(v)
    bool        direct_spatial_mv_pred_flag;                          // u(1)
    bool        num_ref_idx_active_override_flag;                     // u(1)
    uint8_t     num_ref_idx_l0_active_minus1;                         // ue(v)
    uint8_t     num_ref_idx_l1_active_minus1;                         // ue(v)

    bool        ref_pic_list_modification_flag_l0;                    // u(1)
    bool        ref_pic_list_modification_flag_l1;                    // u(1)
    uint8_t     modification_of_pic_nums_idc[2][32+1];                // ue(v)
    uint32_t    abs_diff_pic_num_minus1     [2][32+1];                // ue(v)
    uint32_t    long_term_pic_num           [2][32+1];                // ue(v)
#if (MVC_EXTENSION_ENABLE)
    uint32_t    abs_diff_view_idx_minus1    [2][32+1];                // ue(v)
#endif

    uint8_t     luma_log2_weight_denom;                               // ue(v)
    uint8_t     chroma_log2_weight_denom;                             // ue(v)
    bool        luma_weight_l0_flag  [MAX_NUM_REF_IDX];               // u(1)
    int8_t      luma_weight_l0       [MAX_NUM_REF_IDX];               // se(v)
    int8_t      luma_offset_l0       [MAX_NUM_REF_IDX];               // se(v)
    bool        chroma_weight_l0_flag[MAX_NUM_REF_IDX];               // u(1)
    int8_t      chroma_weight_l0     [MAX_NUM_REF_IDX][2];            // se(v)
    int8_t      chroma_offset_l0     [MAX_NUM_REF_IDX][2];            // se(v)
    bool        luma_weight_l1_flag  [MAX_NUM_REF_IDX];               // u(1)
    int8_t      luma_weight_l1       [MAX_NUM_REF_IDX];               // se(v)
    int8_t      luma_offset_l1       [MAX_NUM_REF_IDX];               // se(v)
    bool        chroma_weight_l1_flag[MAX_NUM_REF_IDX];               // u(1)
    int8_t      chroma_weight_l1     [MAX_NUM_REF_IDX][2];            // se(v)
    int8_t      chroma_offset_l1     [MAX_NUM_REF_IDX][2];            // se(v)

    bool        no_output_of_prior_pics_flag;                         // u(1)
    bool        long_term_reference_flag;                             // u(1)
    bool        adaptive_ref_pic_marking_mode_flag;                   // u(1)
  //uint32_t    memory_management_control_operation;                  // ue(v)
  //uint32_t    difference_of_pic_nums_minus1;                        // ue(v)
  //uint32_t    long_term_pic_num;                                    // ue(v)
  //uint32_t    long_term_frame_idx;                                  // ue(v)
  //uint32_t    max_long_term_frame_idx_plus1;                        // ue(v)
    DecRefPicMarking_t       *dec_ref_pic_marking_buffer; //!< stores the memory management control operations

    uint8_t     cabac_init_idc;                                       // ue(v)
    int8_t      slice_qp_delta;                                       // se(v)
    bool        sp_for_switch_flag;                                   // u(1)
    int8_t      slice_qs_delta;                                       // se(v)
    uint8_t     disable_deblocking_filter_idc;                        // ue(v)
    int8_t      slice_alpha_c0_offset_div2;                           // se(v)
    int8_t      slice_beta_offset_div2;                               // se(v)
    uint32_t    slice_group_change_cycle;                             // u(v)

    PictureStructure    structure;     //!< Identify picture structure type
    bool        MbaffFrameFlag;
    uint32_t    PicHeightInMbs;
    uint32_t    PicHeightInSampleL;
    uint32_t    PicHeightInSampleC;
    uint32_t    PicSizeInMbs;
    uint32_t    MaxPicNum;
    uint32_t    CurrPicNum;
    int8_t      SliceQpY;
    int8_t      QsY;
    int8_t      FilterOffsetA;
    int8_t      FilterOffsetB;
    uint32_t    MapUnitsInSliceGroup0;

    int32_t     PicOrderCntMsb;
    uint32_t    FrameNumOffset;
    int32_t     TopFieldOrderCnt;
    int32_t     BottomFieldOrderCnt;
    // for POC mode 1:
    int         ThisPOC;
    int         framepoc;

    int8_t      PrevQpY;

    //weighted prediction
    unsigned short            weighted_pred_flag;
    unsigned short            weighted_bipred_idc;
    int                    ***wp_weight;  // weight in [list][index][component] order
    int                    ***wp_offset;  // offset in [list][index][component] order
    int                   ****wbp_weight; //weight in [list][fw_index][bw_index][component] order





    //information need to move to slice;
    unsigned int              current_mb_nr; // bitstream order
    unsigned int              num_dec_mb;
    short                     current_slice_nr;
    int                       mb_skip_run;                   //!< Current count of number of skipped macroblocks in a row
    int                       allrefzero;
    //end;


    int                       max_part_nr;
    int                       dp_mode;       //!< data partitioning mode
    int                       current_header;
    int                       last_dquant;

#if (MVC_EXTENSION_ENABLE)
    int                       view_id;
    int                       inter_view_flag;
    int                       anchor_pic_flag;

    NALUnitHeaderMVCExt_t     NaluHeaderMVCExt;
#endif

    //slice header information;
    int                       ref_flag[17]; //!< 0: i-th previous frame is incorrect
    char                      listXsize[6];
    struct storable_picture **listX[6];

    DataPartition            *partArr;      //!< array of partitions
    cabac_contexts_t         *mot_ctx;      //!< pointer to struct of context models for use in CABAC

    int                       mvscale[6][MAX_REFERENCE_PICTURES];


    int                       layer_id;


    int                       dpB_NotPresent;    //!< non-zero, if data partition B is lost
    int                       dpC_NotPresent;    //!< non-zero, if data partition C is lost

    Boolean                   is_reset_coeff;
    Boolean                   is_reset_coeff_cr;
    imgpel                 ***mb_pred; // IntraPrediction()
    imgpel                 ***mb_rec;
    int                    ***mb_rres;
    int                    ***cof;

    imgpel                  **tmp_block_l0; // InterPrediction()
    imgpel                  **tmp_block_l1; // InterPrediction()
    int                     **tmp_res;
    imgpel                  **tmp_block_l2; // InterPrediction()
    imgpel                  **tmp_block_l3; // InterPrediction()

    // Scaling matrix info
    int                       InvLevelScale4x4_Intra[3][6][4][4];
    int                       InvLevelScale4x4_Inter[3][6][4][4];
    int                       InvLevelScale8x8_Intra[3][6][8][8];
    int                       InvLevelScale8x8_Inter[3][6][8][8];

    int                      *qmatrix[12];

    // Cabac
    int                       coeff[64]; // one more for EOB
    int                       coeff_ctr;
    int                       pos;  


#if (MVC_EXTENSION_ENABLE)
    int                       listinterviewidx0;
    int                       listinterviewidx1;
    struct frame_store      **fs_listinterview0;
    struct frame_store      **fs_listinterview1;
#endif

    // for signalling to the neighbour logic that this is a deblocker call
    int                       max_mb_vmv_r; //!< maximum vertical motion vector range in luma quarter pixel units for the current level_idc

    int                       erc_mvperMB;
    struct macroblock_t     *mb_data;
    struct storable_picture  *dec_picture;
    char                     *intra_block;
    char                      chroma_vector_adjustment[6][32];

    bool        init();
    void        decode();
} slice_t;


void slice_header(slice_t *currSlice);
void ref_pic_list_modification(slice_t *currSlice);
void ref_pic_list_mvc_modification(slice_t *currSlice);
void pred_weight_table(slice_t *currSlice);
//void dec_ref_pic_marking(slice_t *currSlice);

void dec_ref_pic_marking(VideoParameters *p_Vid, Bitstream *currStream, slice_t *pSlice);

void decode_poc(VideoParameters *p_Vid, slice_t *pSlice);

#ifdef __cplusplus
}
#endif

#endif /* _SLICE_H_ */
