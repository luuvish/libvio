
/*!
 ************************************************************************
 * \file vlc.h
 *
 * \brief
 *    header for (CA)VLC coding functions
 *
 * \author
 *    Karsten Suehring
 *
 ************************************************************************
 */

#ifndef _VLC_H_
#define _VLC_H_

#ifdef __cplusplus
extern "C" {
#endif


extern int read_se_v (char *tracestring, Bitstream *bitstream, int *used_bits);
extern int read_ue_v (char *tracestring, Bitstream *bitstream, int *used_bits);
extern Boolean read_u_1 (char *tracestring, Bitstream *bitstream, int *used_bits);
extern int read_u_v (int LenInBits, char *tracestring, Bitstream *bitstream, int *used_bits);
extern int read_i_v (int LenInBits, char *tracestring, Bitstream *bitstream, int *used_bits);

// CAVLC mapping
extern void linfo_ue(int len, int info, int *value1, int *dummy);
extern void linfo_se(int len, int info, int *value1, int *dummy);

extern void linfo_cbp_intra_normal(int len,int info,int *cbp, int *dummy);
extern void linfo_cbp_inter_normal(int len,int info,int *cbp, int *dummy);
extern void linfo_cbp_intra_other(int len,int info,int *cbp, int *dummy);
extern void linfo_cbp_inter_other(int len,int info,int *cbp, int *dummy);

extern void linfo_levrun_inter(int len,int info,int *level,int *irun);
extern void linfo_levrun_c2x2(int len,int info,int *level,int *irun);

extern int  uvlc_startcode_follows(Slice *currSlice, int dummy);

extern int  readSyntaxElement_VLC (SyntaxElement *sym, Bitstream *currStream);
extern int  readSyntaxElement_UVLC(Macroblock *currMB, SyntaxElement *sym, struct datapartition_dec *dp);
extern int  readSyntaxElement_Intra4x4PredictionMode(SyntaxElement *sym, Bitstream   *currStream);

extern int  GetVLCSymbol (byte buffer[],int totbitoffset,int *info, int bytecount);
extern int  GetVLCSymbol_IntraMode (byte buffer[],int totbitoffset,int *info, int bytecount);

extern int readSyntaxElement_FLC                         (SyntaxElement *sym, Bitstream *currStream);
extern int readSyntaxElement_NumCoeffTrailingOnes        (SyntaxElement *sym,  Bitstream *currStream, char *type);
extern int readSyntaxElement_NumCoeffTrailingOnesChromaDC(VideoParameters *p_Vid, SyntaxElement *sym, Bitstream *currStream);
extern int readSyntaxElement_Level_VLC0                  (SyntaxElement *sym, Bitstream *currStream);
extern int readSyntaxElement_Level_VLCN                  (SyntaxElement *sym, int vlc, Bitstream *currStream);
extern int readSyntaxElement_TotalZeros                  (SyntaxElement *sym, Bitstream *currStream);
extern int readSyntaxElement_TotalZerosChromaDC          (VideoParameters *p_Vid, SyntaxElement *sym, Bitstream *currStream);
extern int readSyntaxElement_Run                         (SyntaxElement *sym, Bitstream *currStream);
extern int GetBits  (byte buffer[],int totbitoffset,int *info, int bitcount, int numbits);
extern int ShowBits (byte buffer[],int totbitoffset,int bitcount, int numbits);

extern int more_rbsp_data (byte buffer[],int totbitoffset,int bytecount);


#ifdef __cplusplus
}
#endif

#endif /* _VLC_H_ */
