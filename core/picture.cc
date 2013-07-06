
#include <math.h>
#include <limits.h>

#include "global.h"
#include "slice.h"
#include "image.h"
#include "fmo.h"
#include "bitstream_nal.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "parset.h"

#include "sei.h"
#include "output.h"
#include "neighbour.h"
#include "memalloc.h"
#include "macroblock.h"
#include "mb.h"
#include "mb_read.h"

#include "intra_prediction.h"
#include "deblock.h"

#include "biaridecod.h"

#include "erc_api.h"
#include "dpb_common.h"
#include "dpb_mvc.h"


#include "dec_slice.h"

#define MAX_QP          51

static inline int is_BL_profile(unsigned int profile_idc) 
{
  return ( profile_idc == FREXT_CAVLC444 || profile_idc == BASELINE || profile_idc == MAIN || profile_idc == EXTENDED ||
           profile_idc == FREXT_HP || profile_idc == FREXT_Hi10P || profile_idc == FREXT_Hi422 || profile_idc == FREXT_Hi444);
}
static inline int is_EL_profile(unsigned int profile_idc) 
{
  return ( (profile_idc == MVC_HIGH) || (profile_idc == STEREO_HIGH)
           );
}

static inline int is_MVC_profile(unsigned int profile_idc)
{
  return ( (0)
#if (MVC_EXTENSION_ENABLE)
  || (profile_idc == MVC_HIGH) || (profile_idc == STEREO_HIGH)
#endif
  );
}

static void init_qp_process(CodingParameters *cps, int bitdepth_qp_scale)
{
  int i;

  // We should allocate memory outside of this process since maybe we will have a change of SPS 
  // and we may need to recreate these. Currently should only support same bitdepth
  if (cps->qp_per_matrix == NULL)
    if ((cps->qp_per_matrix = (int*)malloc((MAX_QP + 1 +  bitdepth_qp_scale)*sizeof(int))) == NULL)
      no_mem_exit("init_qp_process: cps->qp_per_matrix");

  if (cps->qp_rem_matrix == NULL)
    if ((cps->qp_rem_matrix = (int*)malloc((MAX_QP + 1 +  bitdepth_qp_scale)*sizeof(int))) == NULL)
      no_mem_exit("init_qp_process: cps->qp_rem_matrix");

  for (i = 0; i < MAX_QP + bitdepth_qp_scale + 1; i++)
  {
    cps->qp_per_matrix[i] = i / 6;
    cps->qp_rem_matrix[i] = i % 6;
  }
}

static int init_global_buffers(VideoParameters *p_Vid, int layer_id)
{
  int memory_size=0;
  int i;
  CodingParameters *cps = p_Vid->p_EncodePar[layer_id];
  sps_t *sps = p_Vid->active_sps;
  BlockPos* PicPos;
  int FrameSizeInMbs = sps->PicWidthInMbs * sps->FrameHeightInMbs;

  if (p_Vid->global_init_done[layer_id])
  {
    free_layer_buffers(p_Vid, layer_id);
  }

  // allocate memory for reference frame in find_snr
  memory_size += get_mem2Dpel(&cps->imgY_ref, cps->height, cps->width);
  if (sps->chroma_format_idc != YUV400)
  {
    memory_size += get_mem3Dpel(&cps->imgUV_ref, 2, cps->height_cr, cps->width_cr);
  }
  else
    cps->imgUV_ref = NULL;

  // allocate memory in structure p_Vid
  if( (sps->separate_colour_plane_flag != 0) )
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      if(((cps->mb_data_JV[i]) = (Macroblock *) calloc(FrameSizeInMbs, sizeof(Macroblock))) == NULL)
        no_mem_exit("init_global_buffers: cps->mb_data_JV");
    }
    cps->mb_data = NULL;
  }
  else
  {
    if(((cps->mb_data) = (Macroblock *) calloc(FrameSizeInMbs, sizeof(Macroblock))) == NULL)
      no_mem_exit("init_global_buffers: cps->mb_data");
  }
  if(sps->separate_colour_plane_flag)
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      if(((cps->intra_block_JV[i]) = (char*) calloc(FrameSizeInMbs, sizeof(char))) == NULL)
        no_mem_exit("init_global_buffers: cps->intra_block_JV");
    }
    cps->intra_block = NULL;
  }
  else
  {
    if(((cps->intra_block) = (char*) calloc(FrameSizeInMbs, sizeof(char))) == NULL)
      no_mem_exit("init_global_buffers: cps->intra_block");
  }


  if(((cps->PicPos) = (BlockPos*) calloc(FrameSizeInMbs + 1, sizeof(BlockPos))) == NULL)
    no_mem_exit("init_global_buffers: PicPos");

  PicPos = cps->PicPos;
  for (i = 0; i < (int) FrameSizeInMbs + 1;++i)
  {
    PicPos[i].x = (short) (i % sps->PicWidthInMbs);
    PicPos[i].y = (short) (i / sps->PicWidthInMbs);
  }

  if(sps->separate_colour_plane_flag)
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      get_mem2D(&(cps->ipredmode_JV[i]), 4*sps->FrameHeightInMbs, 4*p_Vid->active_sps->PicWidthInMbs);
    }
    cps->ipredmode = NULL;
  }
  else
   memory_size += get_mem2D(&(cps->ipredmode), 4*sps->FrameHeightInMbs, 4*p_Vid->active_sps->PicWidthInMbs);

  // CAVLC mem
  memory_size += get_mem4D(&(cps->nz_coeff), FrameSizeInMbs, 3, BLOCK_SIZE, BLOCK_SIZE);
  if(sps->separate_colour_plane_flag)
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      get_mem2Dint(&(cps->siblock_JV[i]), sps->FrameHeightInMbs, sps->PicWidthInMbs);
      if(cps->siblock_JV[i]== NULL)
        no_mem_exit("init_global_buffers: p_Vid->siblock_JV");
    }
    cps->siblock = NULL;
  }
  else
  {
    memory_size += get_mem2Dint(&(cps->siblock), sps->FrameHeightInMbs, sps->PicWidthInMbs);
  }
  init_qp_process(cps, imax(sps->QpBdOffsetY, sps->QpBdOffsetC));

  if(layer_id == 0 )
    init_output(cps, ((cps->pic_unit_bitsize_on_disk+7) >> 3));
  else
    cps->img2buf = p_Vid->p_EncodePar[0]->img2buf;
  p_Vid->global_init_done[layer_id] = 1;

  return (memory_size);
}

static void updateMaxValue(FrameFormat *format)
{
  format->max_value[0] = (1 << format->bit_depth[0]) - 1;
  format->max_value_sq[0] = format->max_value[0] * format->max_value[0];
  format->max_value[1] = (1 << format->bit_depth[1]) - 1;
  format->max_value_sq[1] = format->max_value[1] * format->max_value[1];
  format->max_value[2] = (1 << format->bit_depth[2]) - 1;
  format->max_value_sq[2] = format->max_value[2] * format->max_value[2];
}

static void reset_format_info(sps_t *sps, VideoParameters *p_Vid, FrameFormat *source, FrameFormat *output)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  static const int SubWidthC  [4]= { 1, 2, 2, 1};
  static const int SubHeightC [4]= { 1, 2, 1, 1};

  int crop_left, crop_right;
  int crop_top, crop_bottom;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

  // cropping for luma
  if (sps->frame_cropping_flag)
  {
    crop_left   = SubWidthC [sps->chroma_format_idc] * sps->frame_crop_left_offset;
    crop_right  = SubWidthC [sps->chroma_format_idc] * sps->frame_crop_right_offset;
    crop_top    = SubHeightC[sps->chroma_format_idc] * ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_top_offset;
    crop_bottom = SubHeightC[sps->chroma_format_idc] * ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  source->width[0] = p_Vid->width - crop_left - crop_right;
  source->height[0] = p_Vid->height - crop_top - crop_bottom;

  // cropping for chroma
  if (sps->frame_cropping_flag)
  {
    crop_left   = sps->frame_crop_left_offset;
    crop_right  = sps->frame_crop_right_offset;
    crop_top    = ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_top_offset;
    crop_bottom = ( 2 - sps->frame_mbs_only_flag ) *  sps->frame_crop_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  if ((sps->chroma_format_idc==YUV400) && p_Inp->write_uv)
  {
    source->width[1]  = (source->width[0] >> 1);
    source->width[2]  = source->width[1];
    source->height[1] = (source->height[0] >> 1);
    source->height[2] = source->height[1];
  }
  else
  {
    source->width[1]  = p_Vid->width_cr - crop_left - crop_right;
    source->width[2]  = source->width[1];
    source->height[1] = p_Vid->height_cr - crop_top - crop_bottom;
    source->height[2] = source->height[1];
  }

  output->width[0]  = p_Vid->width;
  source->width[1]  = p_Vid->width_cr;
  source->width[2]  = p_Vid->width_cr;
  output->height[0] = p_Vid->height;
  output->height[1] = p_Vid->height_cr;
  output->height[2] = p_Vid->height_cr;

  source->size_cmp[0] = source->width[0] * source->height[0];
  source->size_cmp[1] = source->width[1] * source->height[1];
  source->size_cmp[2] = source->size_cmp[1];
  source->size        = source->size_cmp[0] + source->size_cmp[1] + source->size_cmp[2];
  source->mb_width    = source->width[0]  / MB_BLOCK_SIZE;
  source->mb_height   = source->height[0] / MB_BLOCK_SIZE;

  // output size (excluding padding)
  output->size_cmp[0] = output->width[0] * output->height[0];
  output->size_cmp[1] = output->width[1] * output->height[1];
  output->size_cmp[2] = output->size_cmp[1];
  output->size        = output->size_cmp[0] + output->size_cmp[1] + output->size_cmp[2];
  output->mb_width    = output->width[0]  / MB_BLOCK_SIZE;
  output->mb_height   = output->height[0] / MB_BLOCK_SIZE;


  output->bit_depth[0] = source->bit_depth[0] = sps->bit_depth_luma_minus8 + 8;
  output->bit_depth[1] = source->bit_depth[1] = sps->bit_depth_chroma_minus8 + 8;
  output->bit_depth[2] = source->bit_depth[2] = sps->bit_depth_chroma_minus8 + 8;  
  output->pic_unit_size_on_disk = (imax(output->bit_depth[0], output->bit_depth[1]) > 8) ? 16 : 8;
  output->pic_unit_size_shift3 = output->pic_unit_size_on_disk >> 3;

  output->frame_rate  = source->frame_rate;
  output->color_model = source->color_model;
  output->yuv_format  = source->yuv_format = (ColorFormat) sps->chroma_format_idc;

  output->auto_crop_bottom    = crop_bottom;
  output->auto_crop_right     = crop_right;
  output->auto_crop_bottom_cr = (crop_bottom * mb_cr_size_y) / MB_BLOCK_SIZE;
  output->auto_crop_right_cr  = (crop_right * mb_cr_size_x) / MB_BLOCK_SIZE;

  source->auto_crop_bottom    = output->auto_crop_bottom;
  source->auto_crop_right     = output->auto_crop_right;
  source->auto_crop_bottom_cr = output->auto_crop_bottom_cr;
  source->auto_crop_right_cr  = output->auto_crop_right_cr;

  updateMaxValue(source);
  updateMaxValue(output);

  if (p_Vid->first_sps) {
    p_Vid->first_sps = 0;
    if(!p_Inp->bDisplayDecParams) {
      fprintf(stdout,"Profile IDC  : %d\n", sps->profile_idc);
      fprintf(stdout,"Image Format : %dx%d (%dx%d)\n", source->width[0], source->height[0], p_Vid->width, p_Vid->height);
      if (sps->chroma_format_idc == YUV400)
        fprintf(stdout,"Color Format : 4:0:0 ");
      else if (sps->chroma_format_idc == YUV420)
        fprintf(stdout,"Color Format : 4:2:0 ");
      else if (sps->chroma_format_idc == YUV422)
        fprintf(stdout,"Color Format : 4:2:2 ");
      else
        fprintf(stdout,"Color Format : 4:4:4 ");

      fprintf(stdout,"(%d:%d:%d)\n", source->bit_depth[0], source->bit_depth[1], source->bit_depth[2]);
      fprintf(stdout,"--------------------------------------------------------------------------\n");
    }
    if (!p_Inp->silent)
    {
      fprintf(stdout,"POC must = frame# or field# for SNRs to be correct\n");
      fprintf(stdout,"--------------------------------------------------------------------------\n");
      fprintf(stdout,"  Frame          POC  Pic#   QP    SnrY     SnrU     SnrV   Y:U:V Time(ms)\n");
      fprintf(stdout,"--------------------------------------------------------------------------\n");
    }
  }
}

static void setup_layer_info(VideoParameters *p_Vid, sps_t *sps, LayerParameters *p_Lps)
{
  int layer_id = p_Lps->layer_id;
  p_Lps->p_Vid = p_Vid;
  p_Lps->p_Cps = p_Vid->p_EncodePar[layer_id];
  p_Lps->p_SPS = sps;
  p_Lps->p_Dpb = p_Vid->p_Dpb_layer[layer_id];
}

static void set_coding_par(sps_t *sps, CodingParameters *cps)
{
  // maximum vertical motion vector range in luma quarter pixel units
  cps->profile_idc = sps->profile_idc;
  if (sps->level_idc <= 10)
  {
    cps->max_vmv_r = 64 * 4;
  }
  else if (sps->level_idc <= 20)
  {
    cps->max_vmv_r = 128 * 4;
  }
  else if (sps->level_idc <= 30)
  {
    cps->max_vmv_r = 256 * 4;
  }
  else
  {
    cps->max_vmv_r = 512 * 4; // 512 pixels in quarter pixels
  }

  // Fidelity Range Extensions stuff (part 1)
  cps->width_cr        = 0;
  cps->height_cr       = 0;

  cps->width = sps->PicWidthInMbs * MB_BLOCK_SIZE;
  cps->height = sps->FrameHeightInMbs * MB_BLOCK_SIZE;  

  cps->iLumaPadX = MCBUF_LUMA_PAD_X;
  cps->iLumaPadY = MCBUF_LUMA_PAD_Y;
  cps->iChromaPadX = MCBUF_CHROMA_PAD_X;
  cps->iChromaPadY = MCBUF_CHROMA_PAD_Y;
  if (sps->chroma_format_idc == YUV420)
  {
    cps->width_cr  = (cps->width  >> 1);
    cps->height_cr = (cps->height >> 1);
  }
  else if (sps->chroma_format_idc == YUV422)
  {
    cps->width_cr  = (cps->width >> 1);
    cps->height_cr = cps->height;
    cps->iChromaPadY = MCBUF_CHROMA_PAD_Y*2;
  }
  else if (sps->chroma_format_idc == YUV444)
  {
    //YUV444
    cps->width_cr = cps->width;
    cps->height_cr = cps->height;
    cps->iChromaPadX = cps->iLumaPadX;
    cps->iChromaPadY = cps->iLumaPadY;
  }
  //pel bitdepth init

  if(sps->BitDepthY > sps->BitDepthC || sps->chroma_format_idc == YUV400)
    cps->pic_unit_bitsize_on_disk = (sps->BitDepthY > 8) ? 16 : 8;
  else
    cps->pic_unit_bitsize_on_disk = (sps->BitDepthC > 8) ? 16 : 8;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

  if (sps->chroma_format_idc != YUV400)
  {
    //for chrominance part
    cps->shiftpel_x  = mb_cr_size_x == 8 ? 3 : 2;
    cps->shiftpel_y  = mb_cr_size_y == 8 ? 3 : 2;
    cps->total_scale = cps->shiftpel_x + cps->shiftpel_y;
  }
  else
  {
    cps->shiftpel_x    = 0;
    cps->shiftpel_y    = 0;
    cps->total_scale   = 0;
  }

  cps->rgb_output =  (sps->vui_parameters.matrix_coefficients==0);
}

static void init_frext(VideoParameters *p_Vid)  //!< video parameters
{
    sps_t *sps = p_Vid->active_sps;
  //pel bitdepth init

  if(sps->BitDepthY > sps->BitDepthC || sps->chroma_format_idc == YUV400)
    p_Vid->pic_unit_bitsize_on_disk = (sps->BitDepthY > 8)? 16:8;
  else
    p_Vid->pic_unit_bitsize_on_disk = (sps->BitDepthC > 8)? 16:8;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

  if (sps->chroma_format_idc != YUV400)
  {
    //for chrominance part
    p_Vid->subpel_x    = mb_cr_size_x == 8 ? 7 : 3;
    p_Vid->subpel_y    = mb_cr_size_y == 8 ? 7 : 3;
    p_Vid->shiftpel_x  = mb_cr_size_x == 8 ? 3 : 2;
    p_Vid->shiftpel_y  = mb_cr_size_y == 8 ? 3 : 2;
    p_Vid->total_scale = p_Vid->shiftpel_x + p_Vid->shiftpel_y;
  }
  else
  {
    p_Vid->subpel_x      = 0;
    p_Vid->subpel_y      = 0;
    p_Vid->shiftpel_x    = 0;
    p_Vid->shiftpel_y    = 0;
    p_Vid->total_scale   = 0;
  }
}

static void set_global_coding_par(VideoParameters *p_Vid, CodingParameters *cps)
{
    p_Vid->width_cr        = 0;
    p_Vid->height_cr       = 0;
    p_Vid->max_vmv_r = cps->max_vmv_r;

    // Fidelity Range Extensions stuff (part 1)
    p_Vid->width = cps->width;
    p_Vid->height = cps->height;
    p_Vid->iLumaPadX = MCBUF_LUMA_PAD_X;
    p_Vid->iLumaPadY = MCBUF_LUMA_PAD_Y;
    p_Vid->iChromaPadX = MCBUF_CHROMA_PAD_X;
    p_Vid->iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (p_Vid->active_sps->chroma_format_idc == YUV420) {
        p_Vid->width_cr  = (p_Vid->width  >> 1);
        p_Vid->height_cr = (p_Vid->height >> 1);
    } else if (p_Vid->active_sps->chroma_format_idc == YUV422) {
        p_Vid->width_cr  = (p_Vid->width >> 1);
        p_Vid->height_cr = p_Vid->height;
        p_Vid->iChromaPadY = MCBUF_CHROMA_PAD_Y*2;
    } else if (p_Vid->active_sps->chroma_format_idc == YUV444) {
        //YUV444
        p_Vid->width_cr = p_Vid->width;
        p_Vid->height_cr = p_Vid->height;
        p_Vid->iChromaPadX = p_Vid->iLumaPadX;
        p_Vid->iChromaPadY = p_Vid->iLumaPadY;
    }

    init_frext(p_Vid);
}

void activate_sps (VideoParameters *p_Vid, sps_t *sps)
{
    InputParameters *p_Inp = p_Vid->p_Inp;  

    if (p_Vid->active_sps != sps) {
        if (p_Vid->dec_picture)
            // this may only happen on slice loss
            exit_picture(p_Vid, &p_Vid->dec_picture);
        p_Vid->active_sps = sps;

        if (p_Vid->dpb_layer_id == 0 && is_BL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[0]->init_done) {
            set_coding_par(sps, p_Vid->p_EncodePar[0]);
            setup_layer_info( p_Vid, sps, p_Vid->p_LayerPar[0]);
        } else if (p_Vid->dpb_layer_id == 1 && is_EL_profile(sps->profile_idc) && !p_Vid->p_Dpb_layer[1]->init_done) {
            set_coding_par(sps, p_Vid->p_EncodePar[1]);
            setup_layer_info(p_Vid, sps, p_Vid->p_LayerPar[1]);
        }

        //to be removed in future;
        set_global_coding_par(p_Vid, p_Vid->p_EncodePar[p_Vid->dpb_layer_id]);
        //end;

#if (MVC_EXTENSION_ENABLE)
        if (p_Vid->last_profile_idc != p_Vid->active_sps->profile_idc &&
            is_BL_profile(p_Vid->active_sps->profile_idc) && !p_Vid->p_Dpb_layer[0]->init_done) {
            init_global_buffers(p_Vid, 0);

            if (!p_Vid->no_output_of_prior_pics_flag) {
                flush_dpb(p_Vid->p_Dpb_layer[0]);
                flush_dpb(p_Vid->p_Dpb_layer[1]);
            }
            init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);
        } else if (p_Vid->last_profile_idc != p_Vid->active_sps->profile_idc &&
                   (is_MVC_profile(p_Vid->last_profile_idc) || is_MVC_profile(p_Vid->active_sps->profile_idc)) &&
                   (!p_Vid->p_Dpb_layer[1]->init_done)) {
            assert(p_Vid->p_Dpb_layer[0]->init_done);
            if (p_Vid->p_Dpb_layer[0]->init_done) {
                free_dpb(p_Vid->p_Dpb_layer[0]);
                init_dpb(p_Vid, p_Vid->p_Dpb_layer[0], 1);
            }
            init_global_buffers(p_Vid, 1);
          // for now lets re_init both buffers. Later, we should only re_init appropriate one
          // Note that we seem to be doing this for every frame which seems not good.
            init_dpb(p_Vid, p_Vid->p_Dpb_layer[1], 2);
        }
        p_Vid->last_pic_width_in_mbs_minus1 = p_Vid->active_sps->pic_width_in_mbs_minus1;  
        p_Vid->last_pic_height_in_map_units_minus1 = p_Vid->active_sps->pic_height_in_map_units_minus1;
        p_Vid->last_max_dec_frame_buffering = GetMaxDecFrameBuffering(p_Vid);
        p_Vid->last_profile_idc = p_Vid->active_sps->profile_idc;
#endif

#if (DISABLE_ERC == 0)
        ercInit(p_Vid, p_Vid->width, p_Vid->height, 1);
        if (p_Vid->dec_picture) {
            sps_t *sps = p_Vid->active_sps;
            Slice *currSlice = p_Vid->ppSliceList[0];
            int PicSizeInMbs = sps->PicWidthInMbs * (sps->FrameHeightInMbs / (1 + currSlice->field_pic_flag));
            ercReset(p_Vid->erc_errorVar, PicSizeInMbs, PicSizeInMbs, p_Vid->dec_picture->size_x);
            p_Vid->erc_mvperMB = 0;
        }
#endif
    }
  
    reset_format_info(sps, p_Vid, &p_Inp->source, &p_Inp->output);
}

void activate_pps(VideoParameters *p_Vid, pps_t *pps)
{  
    if (p_Vid->active_pps != pps) {
        if (p_Vid->dec_picture) {
            // this may only happen on slice loss
            exit_picture(p_Vid, &p_Vid->dec_picture);
        }

        p_Vid->active_pps = pps;
    }
}



void UseParameterSet (Slice *currSlice)
{
    VideoParameters *p_Vid = currSlice->p_Vid;
    int PicParsetId = currSlice->pic_parameter_set_id;  
    pps_t *pps = &p_Vid->PicParSet[PicParsetId];
    sps_t *sps = &p_Vid->SeqParSet[pps->seq_parameter_set_id];
    int i;

    if (pps->Valid != TRUE)
        printf ("Trying to use an invalid (uninitialized) Picture Parameter Set with ID %d, expect the unexpected...\n", PicParsetId);
#if (MVC_EXTENSION_ENABLE)
    if (currSlice->svc_extension_flag == -1) {
        if (sps->Valid != TRUE)
            printf ("PicParset %d references an invalid (uninitialized) Sequence Parameter Set with ID %d, expect the unexpected...\n", 
        PicParsetId, (int) pps->seq_parameter_set_id);
    } else {
        // Set SPS to the subset SPS parameters
        p_Vid->active_subset_sps = p_Vid->SubsetSeqParSet + pps->seq_parameter_set_id;
        sps = &(p_Vid->active_subset_sps->sps);
        if (p_Vid->SubsetSeqParSet[pps->seq_parameter_set_id].Valid != TRUE)
            printf ("PicParset %d references an invalid (uninitialized) Subset Sequence Parameter Set with ID %d, expect the unexpected...\n", 
                    PicParsetId, (int) pps->seq_parameter_set_id);
    }
#endif

    // In theory, and with a well-designed software, the lines above
    // are everything necessary.  In practice, we need to patch many values
    // in p_Vid-> (but no more in p_Inp-> -- these have been taken care of)

    // Set Sequence Parameter Stuff first
    if ((int) sps->pic_order_cnt_type < 0 || sps->pic_order_cnt_type > 2) {
        printf("invalid sps->pic_order_cnt_type = %d\n", (int) sps->pic_order_cnt_type);
        error("pic_order_cnt_type != 1", -1000);
    }

    if (sps->pic_order_cnt_type == 1) {
        if (sps->num_ref_frames_in_pic_order_cnt_cycle >= MAX_NUM_REF_FRAMES)
            error("num_ref_frames_in_pic_order_cnt_cycle too large",-1011);
    }
    p_Vid->dpb_layer_id = currSlice->layer_id;
    activate_sps(p_Vid, sps);
    activate_pps(p_Vid, pps);

    // currSlice->dp_mode is set by read_new_slice (NALU first byte available there)
    if (!pps->entropy_coding_mode_flag) {
        for (i = 0; i < 3; i++)
            currSlice->partArr[i].readSyntaxElement = readSyntaxElement_UVLC;      
    } else {
        for (i = 0; i < 3; i++)
            currSlice->partArr[i].readSyntaxElement = readSyntaxElement_CABAC;
    }
    p_Vid->type = currSlice->slice_type;
}

static void init_picture_decoding(VideoParameters *p_Vid)
{
    Slice *pSlice = p_Vid->ppSliceList[0];
    int j, i, iDeblockMode = 1;

    if (p_Vid->iSliceNumOfCurrPic >= MAX_NUM_SLICES)
        error ("Maximum number of supported slices exceeded. \nPlease recompile with increased value for MAX_NUM_SLICES", 200);

    if (p_Vid->pNextPPS->Valid && p_Vid->pNextPPS->pic_parameter_set_id == pSlice->pic_parameter_set_id) {
        pps_t tmpPPS;
        memcpy(&tmpPPS, &(p_Vid->PicParSet[pSlice->pic_parameter_set_id]), sizeof (pps_t));
        (p_Vid->PicParSet[pSlice->pic_parameter_set_id]).slice_group_id = NULL;
        MakePPSavailable(p_Vid, p_Vid->pNextPPS->pic_parameter_set_id, p_Vid->pNextPPS);
        memcpy(p_Vid->pNextPPS, &tmpPPS, sizeof (pps_t));
        tmpPPS.slice_group_id = NULL;
    }

    UseParameterSet(pSlice);
    if (pSlice->idr_flag)
        p_Vid->number = 0;

    p_Vid->structure = pSlice->structure;

    fmo_init(p_Vid, pSlice);

#if (MVC_EXTENSION_ENABLE)
    if (pSlice->layer_id > 0 && pSlice->svc_extension_flag == 0 && pSlice->NaluHeaderMVCExt.non_idr_flag == 0)
        idr_memory_management(p_Vid->p_Dpb_layer[pSlice->layer_id], p_Vid->dec_picture);
    update_ref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_ltref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
    update_pic_num(pSlice);
    i = pSlice->view_id;
#endif
    init_Deblock(p_Vid, pSlice->MbaffFrameFlag);
    //init mb_data;
    for (j = 0; j < p_Vid->iSliceNumOfCurrPic; j++) {
        if (p_Vid->ppSliceList[j]->disable_deblocking_filter_idc != 1)
            iDeblockMode = 0;
#if (MVC_EXTENSION_ENABLE)
        assert(p_Vid->ppSliceList[j]->view_id == i);
#endif
    }
    p_Vid->iDeblockMode = iDeblockMode;
}


void decode_picture(VideoParameters *p_Vid)
{
    Slice **ppSliceList = p_Vid->ppSliceList;
    int iSliceNo;

    p_Vid->num_dec_mb = 0;

    init_picture_decoding(p_Vid);

    for (iSliceNo = 0; iSliceNo < p_Vid->iSliceNumOfCurrPic; iSliceNo++) {
        Slice *currSlice = ppSliceList[iSliceNo];

        assert(currSlice->current_header != EOS);
        assert(currSlice->current_slice_nr == iSliceNo);

        if (!init_slice(currSlice))
            continue;
        decode_one_slice(currSlice);

        p_Vid->num_dec_mb  += currSlice->num_dec_mb;
        p_Vid->erc_mvperMB += currSlice->erc_mvperMB;
    }

#if MVC_EXTENSION_ENABLE
    p_Vid->last_dec_view_id = p_Vid->dec_picture->view_id;
#endif
    if (p_Vid->dec_picture->structure == FRAME)
        p_Vid->last_dec_poc = p_Vid->dec_picture->frame_poc;
    else if (p_Vid->dec_picture->structure == TOP_FIELD)
        p_Vid->last_dec_poc = p_Vid->dec_picture->top_poc;
    else if (p_Vid->dec_picture->structure == BOTTOM_FIELD)
        p_Vid->last_dec_poc = p_Vid->dec_picture->bottom_poc;

    exit_picture(p_Vid, &p_Vid->dec_picture);

    p_Vid->previous_frame_num = ppSliceList[0]->frame_num;
}


static void update_mbaff_macroblock_data(imgpel **cur_img, imgpel (*temp)[16], int x0, int width, int height)
{
    imgpel (*temp_evn)[16] = temp;
    imgpel (*temp_odd)[16] = temp + height; 
    imgpel **temp_img = cur_img;
    int y;

    for (y = 0; y < 2 * height; ++y)
        memcpy(*temp++, (*temp_img++ + x0), width * sizeof(imgpel));

    for (y = 0; y < height; ++y) {
        memcpy((*cur_img++ + x0), *temp_evn++, width * sizeof(imgpel));
        memcpy((*cur_img++ + x0), *temp_odd++, width * sizeof(imgpel));
    }
}

static void MbAffPostProc(VideoParameters *p_Vid)
{
    imgpel temp_buffer[32][16];

    StorablePicture *dec_picture = p_Vid->dec_picture;
    sps_t *sps = p_Vid->active_sps;
    imgpel ** imgY  = dec_picture->imgY;
    imgpel ***imgUV = dec_picture->imgUV;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    short i, x0, y0;

    for (i = 0; i < (int)dec_picture->PicSizeInMbs; i += 2) {
        if (dec_picture->motion.mb_field_decoding_flag[i]) {
            get_mb_pos(p_Vid, i, mb_size, &x0, &y0);
            update_mbaff_macroblock_data(imgY + y0, temp_buffer, x0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

            if (dec_picture->chroma_format_idc != YUV400) {
                x0 = (short) ((x0 * mb_cr_size_x) >> 4);
                y0 = (short) ((y0 * mb_cr_size_y) >> 4);

                update_mbaff_macroblock_data(imgUV[0] + y0, temp_buffer, x0, mb_cr_size_x, mb_cr_size_y);
                update_mbaff_macroblock_data(imgUV[1] + y0, temp_buffer, x0, mb_cr_size_x, mb_cr_size_y);
            }
        }
    }
}

void pad_buf(imgpel *pImgBuf, int iWidth, int iHeight, int iStride, int iPadX, int iPadY)
{
    int j;
    imgpel *pLine0 = pImgBuf - iPadX, *pLine;
    int i;
    for (i = -iPadX; i < 0; i++)
        pImgBuf[i] = *pImgBuf;
    for (i = 0; i < iPadX; i++)
        pImgBuf[i+iWidth] = *(pImgBuf+iWidth-1);

    for (j = -iPadY; j < 0; j++)
        memcpy(pLine0+j*iStride, pLine0, iStride*sizeof(imgpel));
    for (j = 1; j < iHeight; j++) {
        pLine = pLine0 + j*iStride;
        for (i = 0; i < iPadX; i++)
            pLine[i] = pLine[iPadX];
        pLine += iPadX+iWidth-1;
        for (i = 1; i < iPadX + 1; i++)
            pLine[i] = *pLine;
    }
    pLine = pLine0 + (iHeight - 1) * iStride;
    for (j = iHeight; j < iHeight + iPadY; j++)
        memcpy(pLine0+j*iStride,  pLine, iStride*sizeof(imgpel));
}

void pad_dec_picture(VideoParameters *p_Vid, StorablePicture *dec_picture)
{
    int iPadX = p_Vid->iLumaPadX;
    int iPadY = p_Vid->iLumaPadY;
    int iWidth = dec_picture->size_x;
    int iHeight = dec_picture->size_y;
    int iStride = dec_picture->iLumaStride;

    pad_buf(*dec_picture->imgY, iWidth, iHeight, iStride, iPadX, iPadY);

    if (dec_picture->chroma_format_idc != YUV400) {
        iPadX = p_Vid->iChromaPadX;
        iPadY = p_Vid->iChromaPadY;
        iWidth = dec_picture->size_x_cr;
        iHeight = dec_picture->size_y_cr;
        iStride = dec_picture->iChromaStride;
        pad_buf(*dec_picture->imgUV[0], iWidth, iHeight, iStride, iPadX, iPadY);
        pad_buf(*dec_picture->imgUV[1], iWidth, iHeight, iStride, iPadX, iPadY);
    }
}

static void status_picture(VideoParameters *p_Vid, StorablePicture **dec_picture)
{
    InputParameters *p_Inp = p_Vid->p_Inp;
    SNRParameters   *snr   = p_Vid->snr;

    char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};
    char yuvFormat[10];
    int64 tmp_time;                   // time used by decoding the last frame

    int structure         = (*dec_picture)->structure;
    int slice_type        = (*dec_picture)->slice_type;
    int frame_poc         = (*dec_picture)->frame_poc;  
    int refpic            = (*dec_picture)->used_for_reference;
    int qp                = (*dec_picture)->qp;
    int pic_num           = (*dec_picture)->pic_num;
    int is_idr            = (*dec_picture)->idr_flag;
    int chroma_format_idc = (*dec_picture)->chroma_format_idc;

    // report
    char cslice_type[9];  

    if (p_Inp->silent == FALSE) {
        if (structure == TOP_FIELD || structure == FRAME) {
            if (slice_type == I_SLICE && is_idr) // IDR picture
                strcpy(cslice_type,"IDR");
            else if (slice_type == I_SLICE) // I picture
                strcpy(cslice_type," I ");
            else if (slice_type == P_SLICE) // P pictures
                strcpy(cslice_type," P ");
            else if (slice_type == SP_SLICE) // SP pictures
                strcpy(cslice_type,"SP ");
            else if (slice_type == SI_SLICE)
                strcpy(cslice_type,"SI ");
            else if (refpic) // stored B pictures
                strcpy(cslice_type," B ");
            else // B pictures
                strcpy(cslice_type," b ");

            if (structure == FRAME)
                strncat(cslice_type,")       ",8-strlen(cslice_type));
        } else if (structure == BOTTOM_FIELD) {
            if (slice_type == I_SLICE && is_idr) // IDR picture
                strncat(cslice_type,"|IDR)",8-strlen(cslice_type));
            else if (slice_type == I_SLICE) // I picture
                strncat(cslice_type,"| I )",8-strlen(cslice_type));
            else if (slice_type == P_SLICE) // P pictures
                strncat(cslice_type,"| P )",8-strlen(cslice_type));
            else if (slice_type == SP_SLICE) // SP pictures
                strncat(cslice_type,"|SP )",8-strlen(cslice_type));
            else if (slice_type == SI_SLICE)
                strncat(cslice_type,"|SI )",8-strlen(cslice_type));
            else if (refpic) // stored B pictures
                strncat(cslice_type,"| B )",8-strlen(cslice_type));
            else // B pictures
                strncat(cslice_type,"| b )",8-strlen(cslice_type));   
        }
    }

    if (structure == FRAME || structure == BOTTOM_FIELD) {
        gettime (&(p_Vid->end_time));              // end time

        tmp_time  = timediff(&(p_Vid->start_time), &(p_Vid->end_time));
        p_Vid->tot_time += tmp_time;
        tmp_time  = timenorm(tmp_time);
        sprintf(yuvFormat,"%s", yuv_types[chroma_format_idc]);

        if (p_Inp->silent == FALSE) {
            SNRParameters   *snr = p_Vid->snr;
            if (p_Vid->p_ref != -1)
                fprintf(stdout,"%05d(%s%5d %5d %5d %8.4f %8.4f %8.4f  %s %7d\n",
                        p_Vid->frame_no, cslice_type, frame_poc, pic_num, qp, snr->snr[0], snr->snr[1], snr->snr[2], yuvFormat, (int) tmp_time);
            else
                fprintf(stdout,"%05d(%s%5d %5d %5d                             %s %7d\n",
                        p_Vid->frame_no, cslice_type, frame_poc, pic_num, qp, yuvFormat, (int)tmp_time);
        } else
            fprintf(stdout,"Completed Decoding frame %05d.\r",snr->frame_ctr);

        fflush(stdout);

        if (slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic) { // I or P pictures
#if (MVC_EXTENSION_ENABLE)
            if((p_Vid->ppSliceList[0])->view_id!=0)
#endif
                ++(p_Vid->number);
        } else
            ++(p_Vid->Bframe_ctr);    // B pictures
        ++(snr->frame_ctr);

#if (MVC_EXTENSION_ENABLE)
        if ((p_Vid->ppSliceList[0])->view_id != 0)
#endif
            ++(p_Vid->g_nFrame);   
    }
}

void exit_picture(VideoParameters *p_Vid, StorablePicture **dec_picture)
{
    sps_t *sps = p_Vid->active_sps;
    Slice *currSlice = p_Vid->ppSliceList[0];
    int PicSizeInMbs = sps->PicWidthInMbs * (sps->FrameHeightInMbs / (1 + currSlice->field_pic_flag));

    // return if the last picture has already been finished
    if (*dec_picture == NULL ||
        (p_Vid->num_dec_mb != PicSizeInMbs &&
         (sps->chroma_format_idc != YUV444 || !sps->separate_colour_plane_flag)))
        return;

#if (DISABLE_ERC == 0)
    erc_picture(p_Vid, dec_picture);
#endif

    pic_deblock(p_Vid, *dec_picture);

    if ((*dec_picture)->mb_aff_frame_flag)
        MbAffPostProc(p_Vid);

    if (p_Vid->structure != FRAME)
        p_Vid->number /= 2;
#if (MVC_EXTENSION_ENABLE)
    if ((*dec_picture)->used_for_reference || ((*dec_picture)->inter_view_flag == 1))
        pad_dec_picture(p_Vid, *dec_picture);
#endif
#if MVC_EXTENSION_ENABLE
    store_picture_in_dpb(p_Vid->p_Dpb_layer[(*dec_picture)->view_id], *dec_picture);
#endif

    if (p_Vid->last_has_mmco_5)
        p_Vid->pre_frame_num = 0;

    status_picture(p_Vid, dec_picture);

    *dec_picture = NULL;
}
