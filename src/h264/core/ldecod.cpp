#include "global.h"
#include "input_parameters.h"
#include "h264decoder.h"
#include "report.h"

#include "slice.h"
#include "macroblock.h"
#include "interpret.h"
#include "bitstream_cabac.h"
#include "memalloc.h"
#include "dpb.h"
#include "output.h"
#include "sets.h"
#include "sei.h"

#include "erc_api.h"
#include "output.h"

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

#include <stdarg.h>


// Decoder definition. This should be the only global variable in the entire
// software. Global variables should be avoided.
DecoderParams  *p_Dec;


//void error(const char *text, int code)
void error(int code, const char* format, ...)
{
    va_list vl;
    va_start(vl, format);

    vfprintf(stderr, format, vl);

    if (p_Dec) {
        p_Dec->p_Vid->p_Dpb_layer[0]->flush();
#if (MVC_EXTENSION_ENABLE)
        p_Dec->p_Vid->p_Dpb_layer[1]->flush();
#endif
    }

    exit(code);
}



VideoParameters::VideoParameters()
{
    this->out_buffer = new pic_t {};
    this->snr        = new SNRParameters;

    // Allocate new dpb buffer
    for (int i = 0; i < MAX_NUM_DPB_LAYERS; i++) {
        this->p_Dpb_layer[i] = new dpb_t;
        this->p_Dpb_layer[i]->layer_id = i;
        this->p_Dpb_layer[i]->p_Vid = this;
        this->p_Dpb_layer[i]->init_done = 0;

        this->p_EncodePar[i] = new CodingParameters;
        this->p_EncodePar[i]->layer_id = i;

        this->p_LayerPar[i] = new LayerParameters;
        this->p_LayerPar[i]->layer_id = i;
    }
    this->global_init_done[0] = 0;
    this->global_init_done[1] = 0;

    this->seiToneMapping = new ToneMappingSEI;

    this->pNextSlice            = new slice_t;
    this->nalu                  = new nal_unit_t();
    this->pDecOuputPic.pY       = nullptr;
    this->pDecOuputPic.pU       = nullptr;
    this->pDecOuputPic.pV       = nullptr;
    this->pNextPPS              = new pps_t;

    this->recovery_flag         = 0;

    this->recovery_point        = false;
    this->recovery_point_found  = false;
    this->recovery_poc          = 0x7fffffff; /* set to a max value */

    this->number                = 0;
    this->type                  = I_slice;

    // B pictures
    this->snr->Bframe_ctr       = 0;
    this->snr->g_nFrame         = 0;
    this->snr->frame_ctr        = 0;
    this->snr->tot_time         = 0;

    this->dec_picture           = nullptr;

    init_tone_mapping_sei(this->seiToneMapping);

    this->newframe              = 0;
    this->previous_frame_num    = 0;

    this->last_dec_layer_id     = -1;
}

VideoParameters::~VideoParameters()
{
    delete this->out_buffer;
    delete this->snr;

    // Free new dpb layers
    for (int i = 0; i < MAX_NUM_DPB_LAYERS; i++) {
        delete this->p_Dpb_layer[i];
        delete this->p_EncodePar[i];
        delete this->p_LayerPar [i];
    }

    delete this->seiToneMapping;

    for (slice_t* slice : this->ppSliceList)
        delete slice;

    if (this->pNextSlice)
        delete this->pNextSlice;
    delete this->nalu;
    if (this->pDecOuputPic.pY)
        delete []this->pDecOuputPic.pY;
    delete this->pNextPPS;
}

#if (MVC_EXTENSION_ENABLE)
void VideoParameters::OpenOutputFiles(int view0_id, int view1_id)
{
    InputParameters *p_Inp = this->p_Inp;
    char out_ViewFileName[2][FILE_NAME_SIZE], chBuf[FILE_NAME_SIZE], *pch;  
    if (strcasecmp(p_Inp->outfile, "\"\"") != 0 && strlen(p_Inp->outfile) > 0) {
        strcpy(chBuf, p_Inp->outfile);
        pch = strrchr(chBuf, '.');
        if (pch)
            *pch = '\0';
        if (strcmp("nul", chBuf)) {
            sprintf(out_ViewFileName[0], "%s_ViewId%04d.yuv", chBuf, view0_id);
            sprintf(out_ViewFileName[1], "%s_ViewId%04d.yuv", chBuf, view1_id);
            if (this->p_out_mvc[0] >= 0) {
                close(this->p_out_mvc[0]);
                this->p_out_mvc[0] = -1;
            }
            if ((this->p_out_mvc[0] = open(out_ViewFileName[0], O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1) {
                error(500, "Error open file %s ", out_ViewFileName[0]);
            }
      
            if (this->p_out_mvc[1] >= 0) {
                close(this->p_out_mvc[1]);
                this->p_out_mvc[1] = -1;
            }
            if ((this->p_out_mvc[1] = open(out_ViewFileName[1], O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1) {
                error(500, "Error open file %s ", out_ViewFileName[1]);
            }
        }
    }
}
#endif


void DecoderParams::OpenDecoder(InputParameters *p_Inp)
{
    this->p_Vid = new VideoParameters;
    this->p_Vid->p_Inp = this->p_Inp = new InputParameters {};

    p_Dec = this;
    memcpy(this->p_Inp, p_Inp, sizeof(InputParameters));
    this->p_Vid->conceal_mode         = p_Inp->conceal_mode;
    this->p_Vid->snr->idr_psnr_number = p_Inp->ref_offset;

    // Set defaults
    this->p_Vid->p_out = -1;
    for (int i = 0; i < MAX_VIEW_NUM; i++)
        this->p_Vid->p_out_mvc[i] = -1;

    if (p_Inp->DecodeAllLayers == 1)
        this->p_Vid->OpenOutputFiles(0, 1);
    else { //Normal AVC      
        if (strcasecmp(p_Inp->outfile, "\"\"") != 0 && strlen(p_Inp->outfile) > 0) {
            if ((this->p_Vid->p_out_mvc[0] = open(p_Inp->outfile, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR)) == -1) {
                error(500, "Error open file %s ", p_Inp->outfile);
            }
        }
        this->p_Vid->p_out = this->p_Vid->p_out_mvc[0];
    }

    if (strlen(this->p_Inp->reffile) > 0 && strcmp(this->p_Inp->reffile, "\"\"")) {
        if ((this->p_Vid->p_ref = open(this->p_Inp->reffile, O_RDONLY)) == -1) {
            fprintf(stdout, " Input reference file                   : %s does not exist \n", this->p_Inp->reffile);
            fprintf(stdout, "                                          SNR values are not available\n");
        }
    } else
        this->p_Vid->p_ref = -1;

    this->p_Vid->bitstream.open(
        this->p_Inp->infile,
        this->p_Inp->FileFormat ? bitstream_t::type::RTP : bitstream_t::type::ANNEX_B,
        this->p_Vid->nalu->max_size);

    this->p_Vid->active_sps = NULL;
    this->p_Vid->active_subset_sps = NULL;
}


int DecoderParams::DecodeOneFrame()
{
    int iRet = this->decode_slice_headers();

    //init_picture(p_Vid->ppSliceList[0]);
    init_picture_decoding(this->p_Vid);

    this->decode_slice_datas();

    exit_picture(this->p_Vid);
    this->p_Vid->previous_frame_num = this->p_Vid->ppSliceList[0]->header.frame_num;

    if (iRet == SOP)
        iRet = DEC_SUCCEED;
    else if (iRet == EOS)
        iRet = DEC_EOS;
    else
        iRet |= DEC_ERRMASK;
    return iRet;
}

void DecoderParams::FinitDecoder()
{
#if (MVC_EXTENSION_ENABLE)
    this->p_Vid->p_Dpb_layer[0]->flush();
    this->p_Vid->p_Dpb_layer[1]->flush();
#endif

    this->p_Vid->newframe = 0;
    this->p_Vid->previous_frame_num = 0;
}


void free_layer_buffers(VideoParameters *p_Vid, int layer_id)
{
    CodingParameters *cps = p_Vid->p_EncodePar[layer_id];

    if (!p_Vid->global_init_done[layer_id])
        return;

    // free mem, allocated for structure p_Vid
    if (p_Vid->active_sps->separate_colour_plane_flag) {
        for (int i = 0; i < 3; i++) {
            delete []cps->mb_data_JV[i];
            cps->mb_data_JV[i] = nullptr;
        }
    } else {
        if (cps->mb_data) {
            delete []cps->mb_data;
            cps->mb_data = nullptr;
        }
    }

    p_Vid->global_init_done[layer_id] = 0;
}

static void free_global_buffers(VideoParameters *p_Vid)
{
    if (p_Vid->dec_picture) {
        delete p_Vid->dec_picture;
        p_Vid->dec_picture = nullptr;
    }
#if MVC_EXTENSION_ENABLE
    if (p_Vid->active_subset_sps && p_Vid->active_subset_sps->sps.Valid &&
        (p_Vid->active_subset_sps->sps.profile_idc == MVC_HIGH || p_Vid->active_subset_sps->sps.profile_idc == STEREO_HIGH))
        free_img_data(p_Vid, &p_Vid->tempData3);
#endif
}

void DecoderParams::CloseDecoder()
{
    this->p_Vid->report();

    free_layer_buffers(this->p_Vid, 0);
    free_layer_buffers(this->p_Vid, 1);
    free_global_buffers(this->p_Vid);

    this->p_Vid->bitstream.close();

#if (MVC_EXTENSION_ENABLE)
    for (int i = 0; i < MAX_VIEW_NUM; i++) {
        if (this->p_Vid->p_out_mvc[i] != -1)
            close(this->p_Vid->p_out_mvc[i]);
    }
#endif

    if (this->p_Vid->p_ref != -1)
        close(this->p_Vid->p_ref);

#if (DISABLE_ERC == 0)
    ercClose(this->p_Vid, this->p_Vid->erc_errorVar);
#endif

    for (int i = 0; i < MAX_NUM_DPB_LAYERS; i++)
        this->p_Vid->p_Dpb_layer[i]->free();

    delete this->p_Inp;
    delete this->p_Vid;
}
