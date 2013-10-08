#include <cstdint>

#include <stdio.h>

#include "global.h"
#include "input_parameters.h"
#include "report.h"
#include "slice.h"


#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"


void VideoParameters::calculate_frame_no(storable_picture *p)
{
    InputParameters *p_Inp = this->p_Inp;
    // calculate frame number
    int psnrPOC = p->poc / p_Inp->poc_scale;
  
    if (psnrPOC == 0)
        this->snr->idr_psnr_number = this->snr->g_nFrame * p_Inp->ref_poc_gap / p_Inp->poc_scale;

    this->snr->frame_no = this->snr->idr_psnr_number + psnrPOC;
}

void VideoParameters::status(storable_picture** dec_picture)
{
    InputParameters *p_Inp = this->p_Inp;
    SNRParameters   *snr   = this->snr;

    char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};
    char yuvFormat[10];

    int structure         = (*dec_picture)->structure;
    int slice_type        = (*dec_picture)->slice_type;
    int frame_poc         = (*dec_picture)->frame_poc;  
    int refpic            = (*dec_picture)->used_for_reference;
    int qp                = 0;
    int pic_num           = (*dec_picture)->PicNum;
    int is_idr            = (*dec_picture)->idr_flag;
    int chroma_format_idc = this->active_sps->chroma_format_idc;

    // report
    static char cslice_type[9] = { 0 };  

    if (!p_Inp->silent) {
        if (structure == TOP_FIELD || structure == FRAME) {
            if (slice_type == I_slice && is_idr) // IDR picture
                strcpy(cslice_type,"IDR");
            else if (slice_type == I_slice) // I picture
                strcpy(cslice_type," I ");
            else if (slice_type == P_slice) // P pictures
                strcpy(cslice_type," P ");
            else if (slice_type == SP_slice) // SP pictures
                strcpy(cslice_type,"SP ");
            else if (slice_type == SI_slice)
                strcpy(cslice_type,"SI ");
            else if (refpic) // stored B pictures
                strcpy(cslice_type," B ");
            else // B pictures
                strcpy(cslice_type," b ");

            if (structure == FRAME)
                strncat(cslice_type,")       ",8-strlen(cslice_type));
        } else if (structure == BOTTOM_FIELD) {
            if (slice_type == I_slice && is_idr) // IDR picture
                strncat(cslice_type,"|IDR)",8-strlen(cslice_type));
            else if (slice_type == I_slice) // I picture
                strncat(cslice_type,"| I )",8-strlen(cslice_type));
            else if (slice_type == P_slice) // P pictures
                strncat(cslice_type,"| P )",8-strlen(cslice_type));
            else if (slice_type == SP_slice) // SP pictures
                strncat(cslice_type,"|SP )",8-strlen(cslice_type));
            else if (slice_type == SI_slice)
                strncat(cslice_type,"|SI )",8-strlen(cslice_type));
            else if (refpic) // stored B pictures
                strncat(cslice_type,"| B )",8-strlen(cslice_type));
            else // B pictures
                strncat(cslice_type,"| b )",8-strlen(cslice_type));   
        }
    }

    if (structure == FRAME || structure == BOTTOM_FIELD) {
        snr->end_time = std::chrono::system_clock::now();
        int64_t tmp_time = std::chrono::duration_cast<std::chrono::microseconds>(snr->end_time - snr->start_time).count();
        snr->tot_time += tmp_time;

        sprintf(yuvFormat,"%s", yuv_types[chroma_format_idc]);

        if (!p_Inp->silent) {
            if (this->p_ref != -1)
                fprintf(stdout,"%05d(%s%5d %5d %5d %8.4f %8.4f %8.4f  %s %7d\n",
                        snr->frame_no, cslice_type, frame_poc, pic_num, qp, snr->snr[0], snr->snr[1], snr->snr[2], yuvFormat, (int) tmp_time);
            else
                fprintf(stdout,"%05d(%s%5d %5d %5d                             %s %7d\n",
                        snr->frame_no, cslice_type, frame_poc, pic_num, qp, yuvFormat, (int)(tmp_time/1000));
        } else
            fprintf(stdout,"Completed Decoding frame %05d.\r",snr->frame_ctr);

        fflush(stdout);

        if (slice_type == I_slice || slice_type == SI_slice || slice_type == P_slice || refpic) { // I or P pictures
#if (MVC_EXTENSION_ENABLE)
            if (this->ppSliceList[0]->view_id != 0)
#endif
                ++(this->number);
        } else
            ++(snr->Bframe_ctr);
        ++(snr->frame_ctr);

#if (MVC_EXTENSION_ENABLE)
        if ((this->ppSliceList[0])->view_id != 0)
#endif
            ++(snr->g_nFrame);   
    }
}

void VideoParameters::report()
{
    static const char yuv_formats[4][4]= { {"400"}, {"420"}, {"422"}, {"444"} };
    pps_t* active_pps = this->active_pps;
    InputParameters* p_Inp = this->p_Inp;
#define OUTSTRING_SIZE 255
    char string[OUTSTRING_SIZE];

    // normalize time
    snr->tot_time /= 1000;

    if (!p_Inp->silent) {
        fprintf(stdout, "-------------------- Average SNR all frames ------------------------------\n");
        fprintf(stdout, " SNR Y(dB)           : %5.2f\n", snr->snra[0]);
        fprintf(stdout, " SNR U(dB)           : %5.2f\n", snr->snra[1]);
        fprintf(stdout, " SNR V(dB)           : %5.2f\n", snr->snra[2]);
        fprintf(stdout, " Total decoding time : %.3f sec (%.3f fps)[%d frm/%lld ms]\n",
                snr->tot_time * 0.001, (snr->frame_ctr ) * 1000.0 / snr->tot_time, snr->frame_ctr, snr->tot_time);
        fprintf(stdout, "--------------------------------------------------------------------------\n");
        fprintf(stdout, " Exit JM %s decoder, ver %s ", JM, VERSION);
        fprintf(stdout, "\n");
    } else {
        fprintf(stdout, "\n----------------------- Decoding Completed -------------------------------\n");
        fprintf(stdout, " Total decoding time : %.3f sec (%.3f fps)[%d frm/%lld  ms]\n",
                snr->tot_time * 0.001, (snr->frame_ctr) * 1000.0 / snr->tot_time, snr->frame_ctr, snr->tot_time);
        fprintf(stdout, "--------------------------------------------------------------------------\n");
        fprintf(stdout, " Exit JM %s decoder, ver %s ", JM, VERSION);
        fprintf(stdout, "\n");
    }

    // write to log file
    fprintf(stdout, " Output status file                     : %s \n", LOGFILE);
    snprintf(string, OUTSTRING_SIZE, "%s", LOGFILE);

    FILE *p_log;
    if ((p_log = fopen(string,"r")) == 0) {
        if ((p_log = fopen(string,"a")) == 0) {
            snprintf(errortext, ET_SIZE, "Error open file %s for appending",string);
            error(errortext, 500);
        } else {
            fprintf(p_log, " --------------------------------------------------------------------------------------------------------------------\n");
            fprintf(p_log, "|  Decoder statistics. This file is made first time, later runs are appended                                        |\n");
            fprintf(p_log, " --------------------------------------------------------------------------------------------------------------------\n");
            fprintf(p_log, "|   ver  | Date  | Time  |    Sequence        |#Img| Format  | YUV |Coding|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|\n");
            fprintf(p_log, " --------------------------------------------------------------------------------------------------------------------\n");
        }
    } else {
        fclose(p_log);
        p_log = fopen(string, "a");
    }

    fprintf(p_log, "|%s/%-4s", VERSION, EXT_VERSION);

    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    time_t now_c = std::chrono::system_clock::to_time_t(now);
    struct tm* l_time = std::localtime(&now_c);
    strftime(string, sizeof string, "%d-%b-%Y", l_time);
    fprintf(p_log, "| %1.5s |", string);

    strftime(string, sizeof string, "%H:%M:%S", l_time);
    fprintf(p_log, "| %1.5s |", string);

    fprintf(p_log, "%20.20s|", p_Inp->infile);

    fprintf(p_log, "%3d |", this->number);
    fprintf(p_log, "%4dx%-4d|", this->active_sps->PicWidthInMbs * 16, this->active_sps->FrameHeightInMbs * 16);
    fprintf(p_log, " %s |", &yuv_formats[this->active_sps->chroma_format_idc][0]);

    if (active_pps) {
        if (active_pps->entropy_coding_mode_flag)
            fprintf(p_log, " CAVLC|");
        else
            fprintf(p_log, " CABAC|");
    }

    fprintf(p_log, "%6.3f|", snr->snr1[0]);
    fprintf(p_log, "%6.3f|", snr->snr1[1]);
    fprintf(p_log, "%6.3f|", snr->snr1[2]);
    fprintf(p_log, "%6.3f|", snr->snra[0]);
    fprintf(p_log, "%6.3f|", snr->snra[1]);
    fprintf(p_log, "%6.3f|", snr->snra[2]);
    fprintf(p_log, "\n");
    fclose(p_log);

    snprintf(string, OUTSTRING_SIZE, "%s", DATADECFILE);
    p_log = fopen(string, "a");

    if (snr->Bframe_ctr != 0) { // B picture used
        fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
                       "%2.2f %2.2f %2.2f %5d "
                       "%2.2f %2.2f %2.2f %5d %.3f\n",
                this->number, 0, this->ppSliceList[0]->SliceQpY,
                snr->snr1[0], snr->snr1[1], snr->snr1[2], 0,
                0.0, 0.0, 0.0, 0,
                snr->snra[0], snr->snra[1], snr->snra[2], 0,
                (double)0.001 * snr->tot_time / (this->number + snr->Bframe_ctr - 1));
    } else {
        fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
                       "%2.2f %2.2f %2.2f %5d "
                       "%2.2f %2.2f %2.2f %5d %.3f\n",
                this->number, 0, this->ppSliceList[0] ? this->ppSliceList[0]->SliceQpY : 0,
                snr->snr1[0], snr->snr1[1], snr->snr1[2], 0,
                0.0, 0.0, 0.0, 0,
                snr->snra[0], snr->snra[1], snr->snra[2], 0,
                this->number ? ((double)0.001 * snr->tot_time / this->number) : 0.0);
    }
    fclose(p_log);
}
