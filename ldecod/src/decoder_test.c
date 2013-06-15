
/*!
 ***********************************************************************
 *  \file
 *     decoder_test.c
 *  \brief
 *     H.264/AVC decoder test 
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Yuwen He       <yhe@dolby.com>
 ***********************************************************************
 */

#include <sys/stat.h>

#include "win32.h"
#include "h264decoder.h"
#include "configfile.h"


static void Configure(InputParameters *p_Inp, int ac, char *av[])
{
    memset(p_Inp, 0, sizeof(InputParameters));
    strcpy(p_Inp->infile, "test.264"); //! set default bitstream name
    strcpy(p_Inp->outfile, "test_dec.yuv"); //! set default output file name
    strcpy(p_Inp->reffile, "test_rec.yuv"); //! set default reference file name

    ParseCommand(p_Inp, ac, av);

    fprintf(stdout,"----------------------------- JM %s %s -----------------------------\n", VERSION, EXT_VERSION);
    if (!p_Inp->bDisplayDecParams) {
        fprintf(stdout,"--------------------------------------------------------------------------\n");
        fprintf(stdout," Input H.264 bitstream                  : %s \n",p_Inp->infile);
        fprintf(stdout," Output decoded YUV                     : %s \n",p_Inp->outfile);
        fprintf(stdout," Input reference file                   : %s \n",p_Inp->reffile);
        fprintf(stdout,"--------------------------------------------------------------------------\n");
    }
}

/*********************************************************
if bOutputAllFrames is 1, then output all valid frames to file onetime; 
else output the first valid frame and move the buffer to the end of list;
*********************************************************/
static int WriteOneFrame(DecodedPicList *pDecPic, int bOutputAllFrames)
{
    DecodedPicList *pPic = pDecPic;

    if (pPic &&
        ((pPic->iYUVStorageFormat == 2 && pPic->bValid == 3) ||
         (pPic->iYUVStorageFormat != 2 && pPic->bValid == 1))) {
        do {
            pPic->bValid = 0;
            pPic = pPic->pNext;
        } while (pPic != NULL && pPic->bValid && bOutputAllFrames);
    }

    return 0;
}

/*!
 ***********************************************************************
 * \brief
 *    main function for JM decoder
 ***********************************************************************
 */
int main(int argc, char **argv)
{
    int iRet;
    DecodedPicList *pDecPicList;
    int iFramesDecoded = 0;
    InputParameters InputParams;

    init_time();

    //get input parameters;
    Configure(&InputParams, argc, argv);
    //open decoder;
    iRet = OpenDecoder(&InputParams);
    if (iRet != DEC_OPEN_NOERR) {
        fprintf(stderr, "Open encoder failed: 0x%x!\n", iRet);
        return -1; //failed;
    }

    //decoding;
    do {
        iRet = DecodeOneFrame(&pDecPicList);
        if (iRet == DEC_EOS || iRet == DEC_SUCCEED) {
            //process the decoded picture, output or display;
            WriteOneFrame(pDecPicList, 0);
            iFramesDecoded++;
        } else {
            //error handling;
            fprintf(stderr, "Error in decoding process: 0x%x\n", iRet);
        }
    } while ((iRet == DEC_SUCCEED) &&
             (p_Dec->p_Inp->iDecFrmNum == 0 || iFramesDecoded < p_Dec->p_Inp->iDecFrmNum));

    iRet = FinitDecoder(&pDecPicList);
    WriteOneFrame(pDecPicList, 1);
    iRet = CloseDecoder();

    printf("%d frames are decoded.\n", iFramesDecoded);
    return 0;
}
