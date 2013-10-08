/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : Makefile
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 15, 2013    first release
 *
 * ===========================================================================
 */

#include <stdio.h>
#include <string.h>

#include "input_parameters.h"
#include "h264decoder.h"


static void Configure(InputParameters *p_Inp, int ac, char *av[])
{
    memset(p_Inp, 0, sizeof(InputParameters));
    strcpy(p_Inp->infile, "test.264");
    strcpy(p_Inp->outfile, "test_dec.yuv");
    strcpy(p_Inp->reffile, "test_rec.yuv");

    p_Inp->ParseCommand(ac, av);

    fprintf(stdout, "----------------------------- JM %s %s -----------------------------\n", VERSION, EXT_VERSION);
    if (!p_Inp->bDisplayDecParams) {
        fprintf(stdout, "--------------------------------------------------------------------------\n");
        fprintf(stdout, " Input H.264 bitstream                  : %s \n", p_Inp->infile);
        fprintf(stdout, " Output decoded YUV                     : %s \n", p_Inp->outfile);
        fprintf(stdout, " Input reference file                   : %s \n", p_Inp->reffile);
        fprintf(stdout, "--------------------------------------------------------------------------\n");
    }
}


int main(int argc, char** argv)
{
    int iRet;
    int iFramesDecoded = 0;
    InputParameters InputParams;
    DecoderParams Decoder;

    //get input parameters;
    Configure(&InputParams, argc, argv);
    //open decoder;
    Decoder.OpenDecoder(&InputParams);

    //decoding;
    do {
        iRet = Decoder.DecodeOneFrame();
        if (iRet == DEC_EOS || iRet == DEC_SUCCEED)
            iFramesDecoded++;
        else
            //error handling;
            fprintf(stderr, "Error in decoding process: 0x%x\n", iRet);
    } while ((iRet == DEC_SUCCEED) &&
             (Decoder.p_Inp->iDecFrmNum == 0 || iFramesDecoded < Decoder.p_Inp->iDecFrmNum));

    Decoder.FinitDecoder();
    Decoder.CloseDecoder();

    printf("%d frames are decoded.\n", iFramesDecoded);
    return 0;
}
