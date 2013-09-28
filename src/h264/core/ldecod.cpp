#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "data_partition.h"
#include "bitstream_cabac.h"
#include "image.h"
#include "memalloc.h"
#include "dpb.h"
#include "fmo.h"
#include "output.h"
#include "parset.h"
#include "sei.h"

using vio::h264::mb_t;

#include "erc_api.h"
#include "output.h"
#include "h264decoder.h"

#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"
#define TRACEFILE   "trace_dec.txt"

// Decoder definition. This should be the only global variable in the entire
// software. Global variables should be avoided.
DecoderParams  *p_Dec;
char errortext[ET_SIZE];


/*!
 ************************************************************************
 * \brief
 *    Error handling procedure. Print error message to stderr and exit
 *    with supplied code.
 * \param text
 *    Error message
 * \param code
 *    Exit code
 ************************************************************************
 */
void error(const char *text, int code)
{
  fprintf(stderr, "%s\n", text);
  if (p_Dec)
  {
    flush_dpb(p_Dec->p_Vid->p_Dpb_layer[0]);
#if (MVC_EXTENSION_ENABLE)
    flush_dpb(p_Dec->p_Vid->p_Dpb_layer[1]);
#endif
  }

  exit(code);
}



/*!
 ************************************************************************
 * \brief
 *    Reports the gathered information to appropriate outputs
 *
 * \par Input:
 *    InputParameters *p_Inp,
 *    VideoParameters *p_Vid,
 *    struct snr_par *stat
 *
 * \par Output:
 *    None
 ************************************************************************
 */
static void Report(VideoParameters *p_Vid)
{
  static const char yuv_formats[4][4]= { {"400"}, {"420"}, {"422"}, {"444"} };
  pps_t *active_pps = p_Vid->active_pps;
  InputParameters *p_Inp = p_Vid->p_Inp;
  SNRParameters   *snr   = p_Vid->snr;
#define OUTSTRING_SIZE 255
  char string[OUTSTRING_SIZE];
  FILE *p_log;

#ifndef WIN32
  time_t  now;
  struct tm *l_time;
#else
  char timebuf[128];
#endif

  // normalize time
  p_Vid->tot_time /= 1000;

  if (!p_Inp->silent)
  {
    fprintf(stdout,"-------------------- Average SNR all frames ------------------------------\n");
    fprintf(stdout," SNR Y(dB)           : %5.2f\n",snr->snra[0]);
    fprintf(stdout," SNR U(dB)           : %5.2f\n",snr->snra[1]);
    fprintf(stdout," SNR V(dB)           : %5.2f\n",snr->snra[2]);
    fprintf(stdout," Total decoding time : %.3f sec (%.3f fps)[%d frm/%lld ms]\n",p_Vid->tot_time*0.001,(snr->frame_ctr ) * 1000.0 / p_Vid->tot_time, snr->frame_ctr, p_Vid->tot_time);
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout," Exit JM %s decoder, ver %s ",JM, VERSION);
    fprintf(stdout,"\n");
  }
  else
  {
    fprintf(stdout,"\n----------------------- Decoding Completed -------------------------------\n");
    fprintf(stdout," Total decoding time : %.3f sec (%.3f fps)[%d frm/%lld  ms]\n",p_Vid->tot_time*0.001, (snr->frame_ctr) * 1000.0 / p_Vid->tot_time, snr->frame_ctr, p_Vid->tot_time);
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout," Exit JM %s decoder, ver %s ",JM, VERSION);
    fprintf(stdout,"\n");
  }

  // write to log file
  fprintf(stdout," Output status file                     : %s \n",LOGFILE);
  snprintf(string, OUTSTRING_SIZE, "%s", LOGFILE);

  if ((p_log=fopen(string,"r"))==0)                    // check if file exist
  {
    if ((p_log=fopen(string,"a"))==0)
    {
      snprintf(errortext, ET_SIZE, "Error open file %s for appending",string);
      error(errortext, 500);
    }
    else                                              // Create header to new file
    {
      fprintf(p_log," -------------------------------------------------------------------------------------------------------------------\n");
      fprintf(p_log,"|  Decoder statistics. This file is made first time, later runs are appended               |\n");
      fprintf(p_log," ------------------------------------------------------------------------------------------------------------------- \n");
      fprintf(p_log,"|   ver  | Date  | Time  |    Sequence        |#Img| Format  | YUV |Coding|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|\n");
      fprintf(p_log," -------------------------------------------------------------------------------------------------------------------\n");
    }
  }
  else
  {
    fclose(p_log);
    p_log=fopen(string,"a");                    // File exist,just open for appending
  }

  fprintf(p_log,"|%s/%-4s", VERSION, EXT_VERSION);

  now = time ((time_t *) NULL); // Get the system time and put it into 'now' as 'calender time'
  time (&now);
  l_time = localtime (&now);
  strftime (string, sizeof string, "%d-%b-%Y", l_time);
  fprintf(p_log,"| %1.5s |",string );

  strftime (string, sizeof string, "%H:%M:%S", l_time);
  fprintf(p_log,"| %1.5s |",string );

  fprintf(p_log,"%20.20s|",p_Inp->infile);

  fprintf(p_log,"%3d |",p_Vid->number);
  fprintf(p_log,"%4dx%-4d|", p_Vid->active_sps->PicWidthInMbs*16, p_Vid->active_sps->FrameHeightInMbs*16);
  fprintf(p_log," %s |", &(yuv_formats[p_Vid->active_sps->chroma_format_idc][0]));

  if (active_pps)
  {
    if (active_pps->entropy_coding_mode_flag)
      fprintf(p_log," CAVLC|");
    else
      fprintf(p_log," CABAC|");
  }

  fprintf(p_log,"%6.3f|",snr->snr1[0]);
  fprintf(p_log,"%6.3f|",snr->snr1[1]);
  fprintf(p_log,"%6.3f|",snr->snr1[2]);
  fprintf(p_log,"%6.3f|",snr->snra[0]);
  fprintf(p_log,"%6.3f|",snr->snra[1]);
  fprintf(p_log,"%6.3f|",snr->snra[2]);
  fprintf(p_log,"\n");
  fclose(p_log);

  snprintf(string, OUTSTRING_SIZE,"%s", DATADECFILE);
  p_log=fopen(string,"a");

  if(p_Vid->Bframe_ctr != 0) // B picture used
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      p_Vid->number, 0, p_Vid->ppSliceList[0]->SliceQpY,
      snr->snr1[0],
      snr->snr1[1],
      snr->snr1[2],
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snra[0],
      snr->snra[1],
      snr->snra[2],
      0,
      (double)0.001*p_Vid->tot_time/(p_Vid->number + p_Vid->Bframe_ctr - 1));
  }
  else
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      p_Vid->number, 0, p_Vid->ppSliceList[0] ? p_Vid->ppSliceList[0]->SliceQpY : 0,
      snr->snr1[0],
      snr->snr1[1],
      snr->snr1[2],
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snra[0],
      snr->snra[1],
      snr->snra[2],
      0,
      p_Vid->number ? ((double)0.001*p_Vid->tot_time/p_Vid->number) : 0.0);
  }
  fclose(p_log);
}


/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of frame size related global buffers
 *    buffers are defined in global.h, allocated memory is allocated in
 *    int init_global_buffers()
 *
 * \par Input:
 *    Input Parameters VideoParameters *p_Vid
 *
 * \par Output:
 *    none
 *
 ************************************************************************
 */
void free_layer_buffers(VideoParameters *p_Vid, int layer_id)
{
  CodingParameters *cps = p_Vid->p_EncodePar[layer_id];
  
  if(!p_Vid->global_init_done[layer_id])
    return;

  // free mem, allocated for structure p_Vid
  if( (p_Vid->active_sps->separate_colour_plane_flag != 0) )
  {
    int i;
    for(i=0; i<MAX_PLANE; i++)
    {
      free(cps->mb_data_JV[i]);
      cps->mb_data_JV[i] = NULL;
    }   
  }
  else
  {
    if (cps->mb_data != NULL)
    {
      free(cps->mb_data);
      cps->mb_data = NULL;
    }
  }

  p_Vid->global_init_done[layer_id] = 0;
}






static void alloc_video_params(VideoParameters **p_Vid)
{
    if ((*p_Vid = (VideoParameters *)calloc(1, sizeof(VideoParameters))) == NULL) 
        no_mem_exit("alloc_video_params: p_Vid");

    if (((*p_Vid)->old_slice = (OldSliceParams *)calloc(1, sizeof(OldSliceParams))) == NULL)
        no_mem_exit("alloc_video_params: p_Vid->old_slice");

    (*p_Vid)->old_slice->field_pic_flag            = 0;
    (*p_Vid)->old_slice->pps_id                    = INT_MAX;
    (*p_Vid)->old_slice->frame_num                 = INT_MAX;
    (*p_Vid)->old_slice->nal_ref_idc               = INT_MAX;
    (*p_Vid)->old_slice->idr_flag                  = 0;

    (*p_Vid)->old_slice->pic_oder_cnt_lsb          = UINT_MAX;
    (*p_Vid)->old_slice->delta_pic_oder_cnt_bottom = INT_MAX;

    (*p_Vid)->old_slice->delta_pic_order_cnt[0]    = INT_MAX;
    (*p_Vid)->old_slice->delta_pic_order_cnt[1]    = INT_MAX;

    if (((*p_Vid)->snr = (SNRParameters *)calloc(1, sizeof(SNRParameters))) == NULL)
        no_mem_exit("alloc_video_params: p_Vid->snr");  

    // Allocate new dpb buffer
    for (int i = 0; i < MAX_NUM_DPB_LAYERS; i++) {
        if (((*p_Vid)->p_Dpb_layer[i] = (dpb_t *)calloc(1, sizeof(dpb_t))) == NULL) 
            no_mem_exit("alloc_video_params: p_Vid->p_Dpb_layer[i]");
        (*p_Vid)->p_Dpb_layer[i]->layer_id = i;
        (*p_Vid)->p_Dpb_layer[i]->p_Vid = (*p_Vid);
        (*p_Vid)->p_Dpb_layer[i]->init_done = 0;

        if (((*p_Vid)->p_EncodePar[i] = (CodingParameters *)calloc(1, sizeof(CodingParameters))) == NULL)
            no_mem_exit("alloc_video_params:p_Vid->p_EncodePar[i]");
        ((*p_Vid)->p_EncodePar[i])->layer_id = i;

        if (((*p_Vid)->p_LayerPar[i] = (LayerParameters *)calloc(1, sizeof(LayerParameters))) == NULL)
            no_mem_exit("alloc_video_params:p_Vid->p_LayerPar[i]");
        ((*p_Vid)->p_LayerPar[i])->layer_id = i;
    }
    (*p_Vid)->global_init_done[0] = (*p_Vid)->global_init_done[1] = 0;

#if (ENABLE_OUTPUT_TONEMAPPING)  
    if (((*p_Vid)->seiToneMapping = (ToneMappingSEI*)calloc(1, sizeof(ToneMappingSEI))) == NULL)
        no_mem_exit("alloc_video_params: (*p_Vid)->seiToneMapping");  
#endif

    if (((*p_Vid)->ppSliceList = (slice_t **)calloc(MAX_NUM_DECSLICES, sizeof(slice_t *))) == NULL)
        no_mem_exit("alloc_video_params: p_Vid->ppSliceList");

    (*p_Vid)->iNumOfSlicesAllocated = MAX_NUM_DECSLICES;
    (*p_Vid)->pNextSlice = NULL;
    (*p_Vid)->nalu = new nalu_t(MAX_CODED_FRAME_SIZE);
    if (!(*p_Vid)->nalu)
        no_mem_exit("AllocNALU: n");
    (*p_Vid)->pDecOuputPic = (DecodedPicList *)calloc(1, sizeof(DecodedPicList));
    (*p_Vid)->pNextPPS = new pps_t;
    if (!(*p_Vid)->pNextPPS)
         no_mem_exit ("AllocPPS: PPS");
    (*p_Vid)->first_sps = 1;
}

static int alloc_decoder(DecoderParams **p_Dec)
{
    if ((*p_Dec = (DecoderParams *)calloc(1, sizeof(DecoderParams))) == NULL)
        no_mem_exit("alloc_decoder: p_Dec");

    alloc_video_params(&((*p_Dec)->p_Vid));

    if (((*p_Dec)->p_Inp = (InputParameters *)calloc(1, sizeof(InputParameters))) == NULL) 
        no_mem_exit("alloc_params: p_Inp");

    (*p_Dec)->p_Vid->p_Inp = (*p_Dec)->p_Inp;

    return 0;
}

static void init_video_params(VideoParameters *p_Vid)
{
    InputParameters *p_Inp = p_Vid->p_Inp;

    p_Vid->recovery_point = 0;
    p_Vid->recovery_point_found = 0;
    p_Vid->recovery_poc = 0x7fffffff; /* set to a max value */

    p_Vid->idr_psnr_number = p_Inp->ref_offset;
    p_Vid->psnr_number=0;

    p_Vid->number = 0;
    p_Vid->type = I_SLICE;

    p_Vid->g_nFrame = 0;
    // B pictures
    p_Vid->Bframe_ctr = p_Vid->snr->frame_ctr = 0;

    // time for total decoding session
    p_Vid->tot_time = 0;

    p_Vid->dec_picture = NULL;

    p_Vid->MbToSliceGroupMap = NULL;
    p_Vid->MapUnitToSliceGroupMap = NULL;

    p_Vid->out_buffer = NULL;
    p_Vid->pending_output = NULL;
    p_Vid->recovery_flag = 0;

#if (ENABLE_OUTPUT_TONEMAPPING)
    init_tone_mapping_sei(p_Vid->seiToneMapping);
#endif

    p_Vid->newframe = 0;
    p_Vid->previous_frame_num = 0;

    p_Vid->last_dec_layer_id = -1;
}

int OpenDecoder(InputParameters *p_Inp)
{
    int iRet;
    DecoderParams *pDecoder;

    iRet = alloc_decoder(&p_Dec);
    if (iRet)
        return (iRet | DEC_ERRMASK);

    pDecoder = p_Dec;
    memcpy(pDecoder->p_Inp, p_Inp, sizeof(InputParameters));
    pDecoder->p_Vid->conceal_mode = p_Inp->conceal_mode;
    pDecoder->p_Vid->ref_poc_gap  = p_Inp->ref_poc_gap;
    pDecoder->p_Vid->poc_gap      = p_Inp->poc_gap;

    int i;
    VideoParameters *p_Vid = pDecoder->p_Vid;
    // Set defaults
    p_Vid->p_out = -1;
    for (i = 0; i < MAX_VIEW_NUM; i++)
        p_Vid->p_out_mvc[i] = -1;

    if (p_Inp->DecodeAllLayers == 1)
        OpenOutputFiles(p_Vid, 0, 1);
    else { //Normal AVC      
        if (strcasecmp(p_Inp->outfile, "\"\"") != 0 && strlen(p_Inp->outfile) > 0) {
            if ((p_Vid->p_out_mvc[0] = open(p_Inp->outfile, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1) {
                snprintf(errortext, ET_SIZE, "Error open file %s ", p_Inp->outfile);
                error(errortext, 500);
            }
        }
        p_Vid->p_out = p_Vid->p_out_mvc[0];
    }

    if (strlen(pDecoder->p_Inp->reffile) > 0 && strcmp(pDecoder->p_Inp->reffile, "\"\"")) {
        if ((pDecoder->p_Vid->p_ref = open(pDecoder->p_Inp->reffile, O_RDONLY)) == -1) {
            fprintf(stdout, " Input reference file                   : %s does not exist \n", pDecoder->p_Inp->reffile);
            fprintf(stdout, "                                          SNR values are not available\n");
        }
    } else
        pDecoder->p_Vid->p_ref = -1;

    pDecoder->p_Vid->bitstream.open(
        pDecoder->p_Inp->infile,
        pDecoder->p_Inp->FileFormat ? bitstream_t::type::RTP : bitstream_t::type::ANNEX_B,
        pDecoder->p_Vid->nalu->max_size);

    init_video_params(pDecoder->p_Vid);
 
    init_out_buffer(pDecoder->p_Vid);

    pDecoder->p_Vid->active_sps = NULL;
    pDecoder->p_Vid->active_subset_sps = NULL;
    init_subset_sps_list(pDecoder->p_Vid->SubsetSeqParSet, MAXSPS);

    return DEC_OPEN_NOERR;
}


static void ClearDecPicList(VideoParameters *p_Vid)
{
    DecodedPicList *pPic = p_Vid->pDecOuputPic, *pPrior = NULL;
    //find the head first;
    while (pPic && !pPic->bValid) {
        pPrior = pPic;
        pPic = pPic->pNext;
    }

    if (pPic && (pPic != p_Vid->pDecOuputPic)) {
        //move all nodes before pPic to the end;
        DecodedPicList *pPicTail = pPic;
        while (pPicTail->pNext)
            pPicTail = pPicTail->pNext;

        pPicTail->pNext = p_Vid->pDecOuputPic;
        p_Vid->pDecOuputPic = pPic;
        pPrior->pNext = NULL;
    }
}

int DecodeOneFrame(DecodedPicList **ppDecPicList)
{
    int iRet;

    DecoderParams *pDecoder = p_Dec;
    ClearDecPicList(pDecoder->p_Vid);

    iRet = decode_one_frame(pDecoder);
    if (iRet == SOP)
        iRet = DEC_SUCCEED;
    else if (iRet == EOS)
        iRet = DEC_EOS;
    else
        iRet |= DEC_ERRMASK;

    *ppDecPicList = pDecoder->p_Vid->pDecOuputPic;
    return iRet;
}

int FinitDecoder(DecodedPicList **ppDecPicList)
{
    DecoderParams *pDecoder = p_Dec;
    if (!pDecoder)
        return DEC_GEN_NOERR;
    ClearDecPicList(pDecoder->p_Vid);

#if (MVC_EXTENSION_ENABLE)
    flush_dpb(pDecoder->p_Vid->p_Dpb_layer[0]);
    flush_dpb(pDecoder->p_Vid->p_Dpb_layer[1]);
#endif

    pDecoder->p_Vid->bitstream.reset();

    pDecoder->p_Vid->newframe = 0;
    pDecoder->p_Vid->previous_frame_num = 0;
    *ppDecPicList = pDecoder->p_Vid->pDecOuputPic;

    return DEC_GEN_NOERR;
}


static void free_global_buffers(VideoParameters *p_Vid)
{
    if (p_Vid->dec_picture) {
        free_storable_picture(p_Vid->dec_picture);
        p_Vid->dec_picture = NULL;
    }
#if MVC_EXTENSION_ENABLE
    if (p_Vid->active_subset_sps && p_Vid->active_subset_sps->sps.Valid &&
        (p_Vid->active_subset_sps->sps.profile_idc == MVC_HIGH || p_Vid->active_subset_sps->sps.profile_idc == STEREO_HIGH))
        free_img_data( p_Vid, &(p_Vid->tempData3) );
#endif
}

static void FreeDecPicList(DecodedPicList *pDecPicList)
{
    while (pDecPicList) {
        DecodedPicList *pPicNext = pDecPicList->pNext;
        if (pDecPicList->pY) {
            free(pDecPicList->pY);
            pDecPicList->pY = NULL;
            pDecPicList->pU = NULL;
            pDecPicList->pV = NULL;
        }
        free(pDecPicList);
        pDecPicList = pPicNext;
    }
}

static void free_slice(slice_t *currSlice)
{
    free_mem2Dint(currSlice->tmp_res);
    free_mem2Dpel(currSlice->tmp_block_l0);
    free_mem2Dpel(currSlice->tmp_block_l1);
    free_mem2Dpel(currSlice->tmp_block_l2);
    free_mem2Dpel(currSlice->tmp_block_l3);

    free_mem3Dpel(currSlice->mb_pred);

    free_mem3Dint(currSlice->wp_weight );
    free_mem3Dint(currSlice->wp_offset );
    free_mem4Dint(currSlice->wbp_weight);

    for (int i = 0; i < 6; i++) {
        if (currSlice->listX[i]) {
            free(currSlice->listX[i]);
            currSlice->listX[i] = NULL;
        }
    }
    while (currSlice->dec_ref_pic_marking_buffer) {
        DecRefPicMarking_t *tmp_drpm=currSlice->dec_ref_pic_marking_buffer;
        currSlice->dec_ref_pic_marking_buffer=tmp_drpm->Next;
        free(tmp_drpm);
    }

    delete currSlice;
}

static void free_img( VideoParameters *p_Vid)
{
    int i;
    if (p_Vid != NULL) {
#if (ENABLE_OUTPUT_TONEMAPPING)  
        if (p_Vid->seiToneMapping != NULL) {
            free(p_Vid->seiToneMapping);
            p_Vid->seiToneMapping = NULL;
        }
#endif

        // Free new dpb layers
        for (i = 0; i < MAX_NUM_DPB_LAYERS; i++) {
            if (p_Vid->p_Dpb_layer[i] != NULL) {
                free(p_Vid->p_Dpb_layer[i]);
                p_Vid->p_Dpb_layer[i] = NULL;
            }
            if (p_Vid->p_EncodePar[i]) {
                free(p_Vid->p_EncodePar[i]);
                p_Vid->p_EncodePar[i] = NULL;
            }
            if (p_Vid->p_LayerPar[i]) {
                free(p_Vid->p_LayerPar[i]);
                p_Vid->p_LayerPar[i] = NULL;
            }
        }
        if (p_Vid->snr != NULL) {
            free(p_Vid->snr);
            p_Vid->snr = NULL;
        }
        if (p_Vid->old_slice != NULL) {
            free(p_Vid->old_slice);
            p_Vid->old_slice = NULL;
        }
        if (p_Vid->pNextSlice) {
            free_slice(p_Vid->pNextSlice);
            p_Vid->pNextSlice = NULL;
        }
        if (p_Vid->ppSliceList) {
            for (int i = 0; i < p_Vid->iNumOfSlicesAllocated; i++) {
                if (p_Vid->ppSliceList[i])
                    free_slice(p_Vid->ppSliceList[i]);
            }
            free(p_Vid->ppSliceList);
        }
        if (p_Vid->nalu) {
            delete p_Vid->nalu;
            p_Vid->nalu = NULL;
        }
        //free memory;
        FreeDecPicList(p_Vid->pDecOuputPic);
        if (p_Vid->pNextPPS) {
            delete p_Vid->pNextPPS;
            p_Vid->pNextPPS = NULL;
        }

        free(p_Vid);
        p_Vid = NULL;
    }
}

int CloseDecoder()
{
    int i;

    DecoderParams *pDecoder = p_Dec;
    if (!pDecoder)
        return DEC_CLOSE_NOERR;
  
    Report(pDecoder->p_Vid);
    FmoFinit(pDecoder->p_Vid);

    free_layer_buffers(pDecoder->p_Vid, 0);
    free_layer_buffers(pDecoder->p_Vid, 1);
    free_global_buffers(pDecoder->p_Vid);

    pDecoder->p_Vid->bitstream.close();

#if (MVC_EXTENSION_ENABLE)
    for (i = 0; i < MAX_VIEW_NUM; i++) {
        if (pDecoder->p_Vid->p_out_mvc[i] != -1)
            close(pDecoder->p_Vid->p_out_mvc[i]);
    }
#endif

    if (pDecoder->p_Vid->p_ref != -1)
        close(pDecoder->p_Vid->p_ref);

#if (DISABLE_ERC == 0)
    ercClose(pDecoder->p_Vid, pDecoder->p_Vid->erc_errorVar);
#endif

    CleanUpPPS(pDecoder->p_Vid);
#if (MVC_EXTENSION_ENABLE)
    for (i = 0; i < MAXSPS; i++)
        reset_subset_sps(pDecoder->p_Vid->SubsetSeqParSet+i);
#endif

    for (i = 0; i < MAX_NUM_DPB_LAYERS; i++)
        free_dpb(pDecoder->p_Vid->p_Dpb_layer[i]);

    uninit_out_buffer(pDecoder->p_Vid);

    free_img(pDecoder->p_Vid);
    free(pDecoder->p_Inp);
    free(pDecoder);

    p_Dec = NULL;
    return DEC_CLOSE_NOERR;
}

#if (MVC_EXTENSION_ENABLE)
void OpenOutputFiles(VideoParameters *p_Vid, int view0_id, int view1_id)
{
    InputParameters *p_Inp = p_Vid->p_Inp;
    char out_ViewFileName[2][FILE_NAME_SIZE], chBuf[FILE_NAME_SIZE], *pch;  
    if (strcasecmp(p_Inp->outfile, "\"\"") != 0 && strlen(p_Inp->outfile) > 0) {
        strcpy(chBuf, p_Inp->outfile);
        pch = strrchr(chBuf, '.');
        if (pch)
            *pch = '\0';
        if (strcmp("nul", chBuf)) {
            sprintf(out_ViewFileName[0], "%s_ViewId%04d.yuv", chBuf, view0_id);
            sprintf(out_ViewFileName[1], "%s_ViewId%04d.yuv", chBuf, view1_id);
            if (p_Vid->p_out_mvc[0] >= 0) {
                close(p_Vid->p_out_mvc[0]);
                p_Vid->p_out_mvc[0] = -1;
            }
            if ((p_Vid->p_out_mvc[0] = open(out_ViewFileName[0], O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1) {
                snprintf(errortext, ET_SIZE, "Error open file %s ", out_ViewFileName[0]);
                fprintf(stderr, "%s\n", errortext);
                exit(500);
            }
      
            if (p_Vid->p_out_mvc[1] >= 0) {
                close(p_Vid->p_out_mvc[1]);
                p_Vid->p_out_mvc[1] = -1;
            }
            if ((p_Vid->p_out_mvc[1] = open(out_ViewFileName[1], O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1) {
                snprintf(errortext, ET_SIZE, "Error open file %s ", out_ViewFileName[1]);
                fprintf(stderr, "%s\n", errortext);
                exit(500);
            }
        }
    }
}
#endif
