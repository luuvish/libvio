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

struct motion_info_context_t;
struct texture_info_context_t;
struct macroblock_dec;

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

//! Slice
typedef struct slice_t {
    struct video_par         *p_Vid;
    struct inp_par           *p_Inp;
    pic_parameter_set_rbsp_t *active_pps;
    seq_parameter_set_rbsp_t *active_sps;
    int                       svc_extension_flag;

    // dpb pointer
    struct decoded_picture_buffer *p_Dpb;


    //slice property;
    int                       idr_flag;
    int                       nal_reference_idc; //!< nal_reference_idc from NAL unit


    int                       start_mb_nr;   //!< MUST be set by NAL even in case of ei_flag == 1
    int                       slice_type;    //!< slice type
    int                       pic_parameter_set_id;   //!<the ID of the picture parameter set the slice is reffering to
    int                       colour_plane_id;          //!< colour_plane_id of the current coded slice
    unsigned int              frame_num;   //frame_num for this frame
    unsigned int              field_pic_flag;
    byte                      bottom_field_flag;
    PictureStructure          structure;     //!< Identify picture structure type
    int                       idr_pic_id;
    //the following is for slice header syntax elements of poc
    // for poc mode 0.
    unsigned int              pic_order_cnt_lsb;
    int                       delta_pic_order_cnt_bottom;
    // for poc mode 1.
    int                       delta_pic_order_cnt[2];
    int                       redundant_pic_cnt;
    int                       direct_spatial_mv_pred_flag;       //!< Indicator for direct mode type (1 for Spatial, 0 for Temporal)
    int                       num_ref_idx_active[2];             //!< number of available list references

    int                       ref_pic_list_reordering_flag[2];
    int                      *modification_of_pic_nums_idc[2];
    int                      *abs_diff_pic_num_minus1[2];
    int                      *long_term_pic_idx[2];
#if (MVC_EXTENSION_ENABLE)
    int                      *abs_diff_view_idx_minus1[2];
#endif
    int                       redundant_slice_ref_idx;  //!< reference index of redundant slice

    //weighted prediction
    unsigned short            weighted_pred_flag;
    unsigned short            weighted_bipred_idc;

    unsigned short            luma_log2_weight_denom;
    unsigned short            chroma_log2_weight_denom;
    int                    ***wp_weight;  // weight in [list][index][component] order
    int                    ***wp_offset;  // offset in [list][index][component] order
    int                   ****wbp_weight; //weight in [list][fw_index][bw_index][component] order
    short                     wp_round_luma;
    short                     wp_round_chroma;

    int                       no_output_of_prior_pics_flag;
    int                       long_term_reference_flag;
    int                       adaptive_ref_pic_buffering_flag;
    DecRefPicMarking_t       *dec_ref_pic_marking_buffer; //!< stores the memory management control operations

    int                       model_number;  //!< cabac model number
    int                       slice_qp_delta;
    int                       sp_switch;                //!< 1 for switching sp, 0 for normal sp  
    int                       slice_qs_delta;
    int                       qp;
    int                       qs;

    short                     DFDisableIdc;     //!< Disable deblocking filter on slice
    short                     DFAlphaC0Offset;  //!< Alpha and C0 offset for filtering slice
    short                     DFBetaOffset;     //!< Beta offset for filtering slice
    int                       slice_group_change_cycle;




    int                       Transform8x8Mode;
    Boolean                   chroma444_not_separate; //!< indicates chroma 4:4:4 coding with separate_colour_plane_flag equal to zero

    int                       toppoc;    //poc for this top field
    int                       bottompoc; //poc of bottom field of frame
    int                       framepoc;  //poc of this frame


    // ////////////////////////
    // for POC mode 0:
    signed   int              PicOrderCntMsb;
    // for POC mode 1:
    unsigned int              AbsFrameNum;
    int                       ThisPOC;

    //information need to move to slice;
    unsigned int              current_mb_nr; // bitstream order
    unsigned int              num_dec_mb;
    short                     current_slice_nr;
    int                       cod_counter;                   //!< Current count of number of skipped macroblocks in a row
    int                       allrefzero;
    //end;

    int                       mb_aff_frame_flag;

    int                       ei_flag;       //!< 0 if the partArr[0] contains valid information
    int                       end_mb_nr_plus1;
    int                       max_part_nr;
    int                       dp_mode;       //!< data partitioning mode
    int                       current_header;
    int                       next_header;
    int                       last_dquant;

#if (MVC_EXTENSION_ENABLE)
    int                       view_id;
    int                       inter_view_flag;
    int                       anchor_pic_flag;

    NALUnitHeaderMVCExt_t     NaluHeaderMVCExt;
#endif

    //slice header information;

    char                      listXsize[6];
    struct storable_picture **listX[6];

    DataPartition            *partArr;      //!< array of partitions
    struct motion_info_context_t  *mot_ctx;      //!< pointer to struct of context models for use in CABAC
    struct texture_info_context_t *tex_ctx;      //!< pointer to struct of context models for use in CABAC

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
    int                    ***fcf;

    int                       cofu[16];

    imgpel                  **tmp_block_l0;
    imgpel                  **tmp_block_l1;  
    int                     **tmp_res;
    imgpel                  **tmp_block_l2;
    imgpel                  **tmp_block_l3;  

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


    WPParams                **wp_params; // wp parameters in [list][index]


#if (MVC_EXTENSION_ENABLE)
    int                       listinterviewidx0;
    int                       listinterviewidx1;
    struct frame_store      **fs_listinterview0;
    struct frame_store      **fs_listinterview1;
#endif

    // for signalling to the neighbour logic that this is a deblocker call
    int                       max_mb_vmv_r; //!< maximum vertical motion vector range in luma quarter pixel units for the current level_idc
    int                       ref_flag[17]; //!< 0: i-th previous frame is incorrect

    int                       erc_mvperMB;
    struct macroblock_dec     *mb_data;
    struct storable_picture  *dec_picture;
    int                     **siblock;
    byte                    **ipredmode;
    char                     *intra_block;
    char                      chroma_vector_adjustment[6][32];

    void (*read_CBP_and_coeffs_from_NAL)(struct macroblock_dec *currMB);
    int  (*decode_one_component     )(struct macroblock_dec *currMB, ColorPlane curr_plane, imgpel **currImg, struct storable_picture *dec_picture);
    int  (*readSlice                )(struct video_par *, struct inp_par *);  
    int  (*nal_startcode_follows    )(struct slice_t *, int );
    void (*read_motion_info_from_NAL)(struct macroblock_dec *currMB);
    void (*read_one_macroblock      )(struct macroblock_dec *currMB);
    void (*interpret_mb_mode        )(struct macroblock_dec *currMB);
    void (*init_lists               )(struct slice_t *currSlice);

    void (*intra_pred_chroma)(struct macroblock_dec *currMB);
    int  (*intra_pred_4x4)   (struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff,int i4,int j4);
    int  (*intra_pred_8x8)   (struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
    int  (*intra_pred_16x16) (struct macroblock_dec *currMB, ColorPlane pl, int predmode);

    void (*linfo_cbp_intra      )(int len, int info, int *cbp, int *dummy);
    void (*linfo_cbp_inter      )(int len, int info, int *cbp, int *dummy);    
    void (*update_direct_mv_info)(struct macroblock_dec *currMB);
    void (*read_coeff_4x4_CAVLC )(struct macroblock_dec *currMB, int block_type, int i, int j, int levarr[16], int runarr[16], int *number_coefficients);
} Slice;


int FirstPartOfSliceHeader(Slice *currSlice);
int RestOfSliceHeader     (Slice *currSlice);

void dec_ref_pic_marking(VideoParameters *p_Vid, Bitstream *currStream, Slice *pSlice);

void decode_poc(VideoParameters *p_Vid, Slice *pSlice);

#ifdef __cplusplus
}
#endif

#endif /* _SLICE_H_ */
