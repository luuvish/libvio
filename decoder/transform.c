
/*!
 ***********************************************************************
 *  \file
 *      block.c
 *
 *  \brief
 *      Block functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Inge Lille-Langoy          <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg            <rickard.sjoberg@era.ericsson.se>
 ***********************************************************************
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "transform.h"
#include "image.h"
#include "neighbour.h"
#include "transform.h"
#include "quantization.h"
#include "memalloc.h"


#include "bitstream_elements.h"


// SP decoding parameter (EQ. 8-425)
static const int A[4][4] = {
  { 16, 20, 16, 20},
  { 20, 25, 20, 25},
  { 16, 20, 16, 20},
  { 20, 25, 20, 25}
};

static void sample_reconstruct (imgpel **curImg, imgpel **mpr, int **mb_rres, int mb_x, int opix_x, int width, int height, int max_imgpel_value, int dq_bits)
{
  imgpel *imgOrg, *imgPred;
  int    *m7;
  int i, j;

  for (j = 0; j < height; j++)
  {
    imgOrg = &curImg[j][opix_x];
    imgPred = &mpr[j][mb_x];
    m7 = &mb_rres[j][mb_x]; 
    for (i=0;i<width;i++)
      *imgOrg++ = (imgpel) iClip1( max_imgpel_value, rshift_rnd_sf(*m7++, dq_bits) + *imgPred++);
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Inverse 4x4 transformation, transforms cof to mb_rres
 ***********************************************************************
 */
void itrans4x4(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  Slice *currSlice = currMB->p_Slice;
  int    **mb_rres = currSlice->mb_rres[pl];

  inverse4x4(currSlice->cof[pl],mb_rres,joff,ioff);

  sample_reconstruct (&currSlice->mb_rec[pl][joff], &currSlice->mb_pred[pl][joff], &mb_rres[joff], ioff, ioff, BLOCK_SIZE, BLOCK_SIZE, currMB->p_Vid->max_pel_value_comp[pl], DQ_BITS);
}

/*!
 ****************************************************************************
 * \brief
 *    Inverse 4x4 lossless_qpprime transformation, transforms cof to mb_rres
 ****************************************************************************
 */
void itrans4x4_ls(Macroblock *currMB,   //!< current macroblock
                  ColorPlane pl,        //!< Color plane (for 4:4:4)                  
                  int ioff,             //!< index to 4x4 block
                  int joff)             //!< index to 4x4 block
{
  int i,j;

  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres [pl];

  for (j = joff; j < joff + BLOCK_SIZE; ++j)
  {
    for (i = ioff; i < ioff + BLOCK_SIZE; ++i)
    {      
      mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value, mb_pred[j][i] + mb_rres[j][i]);
    }
  }
}

/*!
************************************************************************
* \brief
*    Inverse residual DPCM for Intra lossless coding
*
************************************************************************
*/
void Inv_Residual_trans_4x4(Macroblock *currMB,   //!< current macroblock
                            ColorPlane pl,        //!< used color plane
                            int ioff,             //!< index to 4x4 block
                            int joff)             //!< index to 4x4 block
{
  int i,j;
  int temp[4][4];
  Slice *currSlice = currMB->p_Slice;
  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];
  int    **cof     = currSlice->cof[pl];

  if(currMB->ipmode_DPCM == VERT_PRED)
  {
    for(i=0; i<4; ++i)
    {
      temp[0][i] = cof[joff + 0][ioff + i];
      temp[1][i] = cof[joff + 1][ioff + i] + temp[0][i];
      temp[2][i] = cof[joff + 2][ioff + i] + temp[1][i];
      temp[3][i] = cof[joff + 3][ioff + i] + temp[2][i];
    }

    for(i=0; i<4; ++i)
    {
      mb_rres[joff    ][ioff + i]=temp[0][i];
      mb_rres[joff + 1][ioff + i]=temp[1][i];
      mb_rres[joff + 2][ioff + i]=temp[2][i];
      mb_rres[joff + 3][ioff + i]=temp[3][i];
    }
  }
  else if(currMB->ipmode_DPCM == HOR_PRED)
  {
    for(j=0; j<4; ++j)
    {
      temp[j][0] = cof[joff + j][ioff    ];
      temp[j][1] = cof[joff + j][ioff + 1] + temp[j][0];
      temp[j][2] = cof[joff + j][ioff + 2] + temp[j][1];
      temp[j][3] = cof[joff + j][ioff + 3] + temp[j][2];
    }

    for(j=0; j<4; ++j)
    {
      mb_rres[joff + j][ioff    ]=temp[j][0];
      mb_rres[joff + j][ioff + 1]=temp[j][1];
      mb_rres[joff + j][ioff + 2]=temp[j][2];
      mb_rres[joff + j][ioff + 3]=temp[j][3];
    }
  }
  else
  {
    for (j = joff; j < joff + BLOCK_SIZE; ++j)
      for (i = ioff; i < ioff + BLOCK_SIZE; ++i)
        mb_rres[j][i] = cof[j][i];
  }

  for (j = joff; j < joff + BLOCK_SIZE; ++j)
  {
    for (i = ioff; i < ioff + BLOCK_SIZE; ++i)
    {
      mb_rec[j][i] = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
  }
}

/*!
************************************************************************
* \brief
*    Inverse residual DPCM for Intra lossless coding
*
* \par Input:
*    ioff_x,joff_y: Block position inside a macro block (0,8).
************************************************************************
*/
//For residual DPCM
void Inv_Residual_trans_8x8(Macroblock *currMB, ColorPlane pl, int ioff,int joff)
{
  Slice *currSlice = currMB->p_Slice;
  int i, j;
  int temp[8][8];
  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];

  if(currMB->ipmode_DPCM == VERT_PRED)
  {
    for(i=0; i<8; ++i)
    {
      temp[0][i] = mb_rres[joff + 0][ioff + i];
      temp[1][i] = mb_rres[joff + 1][ioff + i] + temp[0][i];
      temp[2][i] = mb_rres[joff + 2][ioff + i] + temp[1][i];
      temp[3][i] = mb_rres[joff + 3][ioff + i] + temp[2][i];
      temp[4][i] = mb_rres[joff + 4][ioff + i] + temp[3][i];
      temp[5][i] = mb_rres[joff + 5][ioff + i] + temp[4][i];
      temp[6][i] = mb_rres[joff + 6][ioff + i] + temp[5][i];
      temp[7][i] = mb_rres[joff + 7][ioff + i] + temp[6][i];
    }
    for(i=0; i<8; ++i)
    {
      mb_rres[joff  ][ioff+i]=temp[0][i];
      mb_rres[joff+1][ioff+i]=temp[1][i];
      mb_rres[joff+2][ioff+i]=temp[2][i];
      mb_rres[joff+3][ioff+i]=temp[3][i];
      mb_rres[joff+4][ioff+i]=temp[4][i];
      mb_rres[joff+5][ioff+i]=temp[5][i];
      mb_rres[joff+6][ioff+i]=temp[6][i];
      mb_rres[joff+7][ioff+i]=temp[7][i];
    }
  }
  else if(currMB->ipmode_DPCM == HOR_PRED)//HOR_PRED
  {
    for(i=0; i<8; ++i)
    {
      temp[i][0] = mb_rres[joff + i][ioff + 0];
      temp[i][1] = mb_rres[joff + i][ioff + 1] + temp[i][0];
      temp[i][2] = mb_rres[joff + i][ioff + 2] + temp[i][1];
      temp[i][3] = mb_rres[joff + i][ioff + 3] + temp[i][2];
      temp[i][4] = mb_rres[joff + i][ioff + 4] + temp[i][3];
      temp[i][5] = mb_rres[joff + i][ioff + 5] + temp[i][4];
      temp[i][6] = mb_rres[joff + i][ioff + 6] + temp[i][5];
      temp[i][7] = mb_rres[joff + i][ioff + 7] + temp[i][6];
    }
    for(i=0; i<8; ++i)
    {
      mb_rres[joff+i][ioff+0]=temp[i][0];
      mb_rres[joff+i][ioff+1]=temp[i][1];
      mb_rres[joff+i][ioff+2]=temp[i][2];
      mb_rres[joff+i][ioff+3]=temp[i][3];
      mb_rres[joff+i][ioff+4]=temp[i][4];
      mb_rres[joff+i][ioff+5]=temp[i][5];
      mb_rres[joff+i][ioff+6]=temp[i][6];
      mb_rres[joff+i][ioff+7]=temp[i][7];
    }
  }

  for (j = joff; j < joff + BLOCK_SIZE*2; ++j)
  {
    for (i = ioff; i < ioff + BLOCK_SIZE*2; ++i)
    {
      mb_rec [j][i]  = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
  }
}



/*!
************************************************************************
* \brief
*    Inverse residual DPCM for Intra lossless coding
*
************************************************************************
*/
void Inv_Residual_trans_16x16(Macroblock *currMB,   //!< current macroblock
                              ColorPlane pl)        //!< used color plane
{
  int i,j;
  int temp[16][16];
  Slice *currSlice = currMB->p_Slice;
  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];
  int    **cof     = currSlice->cof[pl];

  if(currMB->ipmode_DPCM == VERT_PRED_16)
  {
    for(i=0; i<MB_BLOCK_SIZE; ++i)
    {
      temp[0][i] = cof[0][i];
      for(j = 1; j < MB_BLOCK_SIZE; j++)
        temp[j][i] = cof[j][i] + temp[j-1][i];
    }

    for(i=0; i<MB_BLOCK_SIZE; ++i)
    {
      for(j = 0; j < MB_BLOCK_SIZE; j++)
        mb_rres[j][i]=temp[j][i];
    }
  }
  else if(currMB->ipmode_DPCM == HOR_PRED_16)
  {
    for(j=0; j<MB_BLOCK_SIZE; ++j)
    {
      temp[j][ 0] = cof[j][ 0  ];
      for(i = 1; i < MB_BLOCK_SIZE; i++)
        temp[j][i] = cof[j][i] + temp[j][i-1];
    }

    for(j=0; j<MB_BLOCK_SIZE; ++j)
    {
      for(i = 0; i < MB_BLOCK_SIZE; ++i)
        mb_rres[j][i]=temp[j][i];
    }
  }
  else
  {
    for (j = 0; j < MB_BLOCK_SIZE; ++j)
      for (i = 0; i < MB_BLOCK_SIZE; ++i)
        mb_rres[j][i] = cof[j][i];
  }

  for (j = 0; j < MB_BLOCK_SIZE; ++j)
  {
    for (i = 0; i < MB_BLOCK_SIZE; ++i)
    {
      mb_rec[j][i] = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
  }
}


/*!
************************************************************************
* \brief
*    Inverse residual DPCM for Intra lossless coding
*
************************************************************************
*/
void Inv_Residual_trans_Chroma(Macroblock *currMB, int uv)  
{
  int i, j;
  int temp[16][16];
  Slice *currSlice = currMB->p_Slice;
  //imgpel **mb_pred = currSlice->mb_pred[uv+1];
  //imgpel **mb_rec  = currSlice->mb_rec[uv+1];
  int    **mb_rres = currSlice->mb_rres[uv+1];
  int    **cof     = currSlice->cof[uv+1];
  int width, height; 

  width = currMB->p_Vid->mb_cr_size_x;
  height = currMB->p_Vid->mb_cr_size_y;

  if(currMB->c_ipred_mode == VERT_PRED_8)
  {
    for(i=0; i<width; i++)
    {
      temp[0][i] = cof[0][i];
      for(j = 1; j < height; j++)
        temp[j][i] = temp[j-1][i] + cof[j][i];
    }
    for(i=0; i<width; i++)
    {
      for(j = 0; j < height; j++)
        mb_rres[j][i] = temp[j][i];
    }
  }
  else //HOR_PRED_8
  {
    for(i=0; i<height; i++)
    {
      temp[i][0] = cof[i][0];
      for(j = 1; j < width; j++)
        temp[i][j] = temp[i][j-1] + cof[i][j];
    }
    for(i=0; i<height; i++)
    {
      for(j = 0; j < width; j++)
        mb_rres[i][j] = temp[i][j];
    }
  }
}


/*!
 ***********************************************************************
 * \brief
 *    Luma DC inverse transform
 ***********************************************************************
 */ 
void itrans_2(Macroblock *currMB,    //!< current macroblock
              ColorPlane pl)         //!< used color plane
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int j;

  int transform_pl = (p_Vid->separate_colour_plane_flag != 0) ? PLANE_Y : pl;
  int **cof = currSlice->cof[transform_pl];
  int qp_scaled = currMB->qp_scaled[transform_pl];

  int qp_per = p_Vid->qp_per_matrix[ qp_scaled ];
  int qp_rem = p_Vid->qp_rem_matrix[ qp_scaled ];      

  int invLevelScale = currSlice->InvLevelScale4x4_Intra[pl][qp_rem][0][0];
  int **M4;
  get_mem2Dint(&M4, BLOCK_SIZE, BLOCK_SIZE);
  
  // horizontal
  for (j=0; j < 4;++j) 
  {
    M4[j][0]=cof[j<<2][0];
    M4[j][1]=cof[j<<2][4];
    M4[j][2]=cof[j<<2][8];
    M4[j][3]=cof[j<<2][12];
  }

  ihadamard4x4(M4, M4);

  // vertical
  for (j=0; j < 4;++j) 
  {
    cof[j<<2][0]  = rshift_rnd((( M4[j][0] * invLevelScale) << qp_per), 6);
    cof[j<<2][4]  = rshift_rnd((( M4[j][1] * invLevelScale) << qp_per), 6);
    cof[j<<2][8]  = rshift_rnd((( M4[j][2] * invLevelScale) << qp_per), 6);
    cof[j<<2][12] = rshift_rnd((( M4[j][3] * invLevelScale) << qp_per), 6);
  }

  free_mem2Dint(M4);
}


void itrans_sp(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;
  int i,j;  
  int ilev, icof;

  int qp = (currSlice->slice_type == SI_SLICE) ? currSlice->qs : currSlice->qp;
  int qp_per = p_Vid->qp_per_matrix[ qp ];
  int qp_rem = p_Vid->qp_rem_matrix[ qp ];

  int qp_per_sp = p_Vid->qp_per_matrix[ currSlice->qs ];
  int qp_rem_sp = p_Vid->qp_rem_matrix[ currSlice->qs ];
  int q_bits_sp = Q_BITS + qp_per_sp;

  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];
  int    **cof     = currSlice->cof[pl];
  int max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  const int (*InvLevelScale4x4)  [4] = dequant_coef[qp_rem];
  const int (*InvLevelScale4x4SP)[4] = dequant_coef[qp_rem_sp];  
  int **PBlock;  

  get_mem2Dint(&PBlock, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  for (j=0; j< BLOCK_SIZE; ++j)
  {
    PBlock[j][0] = mb_pred[j+joff][ioff    ];
    PBlock[j][1] = mb_pred[j+joff][ioff + 1];
    PBlock[j][2] = mb_pred[j+joff][ioff + 2];
    PBlock[j][3] = mb_pred[j+joff][ioff + 3];
  }

  forward4x4(PBlock, PBlock, 0, 0);

  if(currSlice->sp_switch || currSlice->slice_type==SI_SLICE)
  {    
    for (j=0;j<BLOCK_SIZE;++j)
    {
      for (i=0;i<BLOCK_SIZE;++i)
      {
        // recovering coefficient since they are already dequantized earlier
        icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
        //icof = ((cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i])>> (qp_per + 15)) ;
        // icof  = rshift_rnd_sf(cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i], qp_per + 15);
        ilev  = rshift_rnd_sf(iabs(PBlock[j][i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
        ilev  = isignab(ilev, PBlock[j][i]) + icof;
        cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
      }
    }
  }
  else
  {
    for (j=0;j<BLOCK_SIZE;++j)
    {
      for (i=0;i<BLOCK_SIZE;++i)
      {
        // recovering coefficient since they are already dequantized earlier
        icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
        //icof = cof[joff + j][ioff + i];
        //icof  = rshift_rnd_sf(cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i], qp_per + 15);
        ilev = PBlock[j][i] + ((icof * InvLevelScale4x4[j][i] * A[j][i] <<  qp_per) >> 6);
        ilev  = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
        //cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
        cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
      }
    }
  }

  inverse4x4(cof, mb_rres, joff, ioff);

  for (j=joff; j<joff +BLOCK_SIZE;++j)
  {
    mb_rec[j][ioff   ] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff   ], DQ_BITS));
    mb_rec[j][ioff+ 1] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff+ 1], DQ_BITS));
    mb_rec[j][ioff+ 2] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff+ 2], DQ_BITS));
    mb_rec[j][ioff+ 3] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff+ 3], DQ_BITS));
  }  

  free_mem2Dint(PBlock);
}


void itrans_sp_cr(Macroblock *currMB, int uv)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int i,j,ilev, icof, n2,n1;
  int mp1[BLOCK_SIZE];
  int qp_per,qp_rem;
  int qp_per_sp,qp_rem_sp,q_bits_sp;
  imgpel **mb_pred = currSlice->mb_pred[uv + 1];
  int    **cof = currSlice->cof[uv + 1];
  int **PBlock = new_mem2Dint(MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  qp_per    = p_Vid->qp_per_matrix[ ((currSlice->qp < 0 ? currSlice->qp : QP_SCALE_CR[currSlice->qp]))];
  qp_rem    = p_Vid->qp_rem_matrix[ ((currSlice->qp < 0 ? currSlice->qp : QP_SCALE_CR[currSlice->qp]))];

  qp_per_sp = p_Vid->qp_per_matrix[ ((currSlice->qs < 0 ? currSlice->qs : QP_SCALE_CR[currSlice->qs]))];
  qp_rem_sp = p_Vid->qp_rem_matrix[ ((currSlice->qs < 0 ? currSlice->qs : QP_SCALE_CR[currSlice->qs]))];
  q_bits_sp = Q_BITS + qp_per_sp;  

  if (currSlice->slice_type == SI_SLICE)
  {
    qp_per = qp_per_sp;
    qp_rem = qp_rem_sp;
  }

  for (j=0; j < p_Vid->mb_cr_size_y; ++j)
  {
    for (i=0; i < p_Vid->mb_cr_size_x; ++i)
    {
      PBlock[j][i] = mb_pred[j][i];
      mb_pred[j][i] = 0;
    }
  }

  for (n2=0; n2 < p_Vid->mb_cr_size_y; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 < p_Vid->mb_cr_size_x; n1 += BLOCK_SIZE)
    {
      forward4x4(PBlock, PBlock, n2, n1);
    }
  }

  //     2X2 transform of DC coeffs.
  mp1[0] = (PBlock[0][0] + PBlock[4][0] + PBlock[0][4] + PBlock[4][4]);
  mp1[1] = (PBlock[0][0] - PBlock[4][0] + PBlock[0][4] - PBlock[4][4]);
  mp1[2] = (PBlock[0][0] + PBlock[4][0] - PBlock[0][4] - PBlock[4][4]);
  mp1[3] = (PBlock[0][0] - PBlock[4][0] - PBlock[0][4] + PBlock[4][4]);

  if (currSlice->sp_switch || currSlice->slice_type == SI_SLICE)  
  {        
    for (n2=0; n2 < 2; ++n2 )
    {
      for (n1=0; n1 < 2; ++n1 )
      {
        //quantization fo predicted block
        ilev = rshift_rnd_sf(iabs (mp1[n1+n2*2]) * quant_coef[qp_rem_sp][0][0], q_bits_sp + 1);
        //addition
        ilev = isignab(ilev, mp1[n1+n2*2]) + cof[n2<<2][n1<<2];
        //dequantization
        mp1[n1+n2*2] =ilev * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
      }
    }

    for (n2 = 0; n2 < p_Vid->mb_cr_size_y; n2 += BLOCK_SIZE)
    {
      for (n1 = 0; n1 < p_Vid->mb_cr_size_x; n1 += BLOCK_SIZE)
      {
        for (j = 0; j < BLOCK_SIZE; ++j)
        {
          for (i = 0; i < BLOCK_SIZE; ++i)
          {
            // recovering coefficient since they are already dequantized earlier
            cof[n2 + j][n1 + i] = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];

            //quantization of the predicted block
            ilev = rshift_rnd_sf(iabs(PBlock[n2 + j][n1 + i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
            //addition of the residual
            ilev = isignab(ilev,PBlock[n2 + j][n1 + i]) + cof[n2 + j][n1 + i];
            // Inverse quantization
            cof[n2 + j][n1 + i] = ilev * dequant_coef[qp_rem_sp][j][i] << qp_per_sp;
          }
        }
      }
    }
  }
  else
  {
    for (n2=0; n2 < 2; ++n2 )
    {
      for (n1=0; n1 < 2; ++n1 )
      {
        ilev = mp1[n1+n2*2] + (((cof[n2<<2][n1<<2] * dequant_coef[qp_rem][0][0] * A[0][0]) << qp_per) >> 5);
        ilev = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][0][0], q_bits_sp + 1);
        //ilev = isignab(rshift_rnd_sf(iabs(ilev)* quant_coef[qp_rem_sp][0][0], q_bits_sp + 1), ilev);
        mp1[n1+n2*2] = ilev * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
      }
    }

    for (n2 = 0; n2 < p_Vid->mb_cr_size_y; n2 += BLOCK_SIZE)
    {
      for (n1 = 0; n1 < p_Vid->mb_cr_size_x; n1 += BLOCK_SIZE)
      {
        for (j = 0; j< BLOCK_SIZE; ++j)
        {
          for (i = 0; i< BLOCK_SIZE; ++i)
          {
            // recovering coefficient since they are already dequantized earlier
            //icof = ((((cof[n2 + j][n1 + i] << 4) + qp_per/2)>> qp_per) + dequant_coef[qp_rem][j][i]/2) / dequant_coef[qp_rem][j][i];
            icof = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];
            //dequantization and addition of the predicted block      
            ilev = PBlock[n2 + j][n1 + i] + ((icof * dequant_coef[qp_rem][j][i] * A[j][i] << qp_per) >> 6);
            //quantization and dequantization
            ilev = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
            cof[n2 + j][n1 + i] = ilev * dequant_coef[qp_rem_sp][j][i] << qp_per_sp;
            //printf( " %d %d %d\n", j, i, quant_coef[qp_rem_sp][j][i]);
          }
        }
      }
    }
  }

  cof[0][0] = (mp1[0] + mp1[1] + mp1[2] + mp1[3]) >> 1;
  cof[0][4] = (mp1[0] + mp1[1] - mp1[2] - mp1[3]) >> 1;
  cof[4][0] = (mp1[0] - mp1[1] + mp1[2] - mp1[3]) >> 1;
  cof[4][4] = (mp1[0] - mp1[1] - mp1[2] + mp1[3]) >> 1;

  free_mem2Dint(PBlock);
}

void iMBtrans4x4(Macroblock *currMB, ColorPlane pl, int smb)
{
  Slice *currSlice = currMB->p_Slice;
  //VideoParameters *p_Vid = currMB->p_Vid;

  StorablePicture *dec_picture = currMB->p_Slice->dec_picture;
  int jj, ii;
  int block8x8;
  int k;  

  imgpel **curr_img = pl ? dec_picture->imgUV[pl - 1]: dec_picture->imgY;

  // =============== 4x4 itrans ================
  // -------------------------------------------
  if (currMB->is_lossless && currMB->mb_type == I16MB)
  {
    Inv_Residual_trans_16x16(currMB, pl) ;
  }
  else if (smb || currMB->is_lossless == TRUE)
  {
    currMB->itrans_4x4 = (smb) ? itrans_sp : ((currMB->is_lossless == FALSE) ? itrans4x4 : Inv_Residual_trans_4x4);
    for (block8x8=0; block8x8 < MB_BLOCK_SIZE; block8x8 += 4)
    { 
      for (k = block8x8; k < block8x8 + 4; ++k )
      {
        jj = ((decode_block_scan[k] >> 2) & 3) << BLOCK_SHIFT;
        ii = (decode_block_scan[k] & 3) << BLOCK_SHIFT;

        currMB->itrans_4x4(currMB, pl, ii, jj);   // use integer transform and make 4x4 block mb_rres from prediction block mb_pred
      }
    }
  }
  else
  {
    int **cof = currSlice->cof[pl];
    int **mb_rres = currSlice->mb_rres[pl];

    if (currMB->is_intra_block == FALSE)
    {
      if (currMB->cbp & 0x01)
      {
        inverse4x4(cof, mb_rres, 0, 0);
        inverse4x4(cof, mb_rres, 0, 4);
        inverse4x4(cof, mb_rres, 4, 0);
        inverse4x4(cof, mb_rres, 4, 4);
      }
      if (currMB->cbp & 0x02)
      {
        inverse4x4(cof, mb_rres, 0, 8);
        inverse4x4(cof, mb_rres, 0, 12);
        inverse4x4(cof, mb_rres, 4, 8);
        inverse4x4(cof, mb_rres, 4, 12);
      }
      if (currMB->cbp & 0x04)
      {
        inverse4x4(cof, mb_rres, 8, 0);
        inverse4x4(cof, mb_rres, 8, 4);
        inverse4x4(cof, mb_rres, 12, 0);
        inverse4x4(cof, mb_rres, 12, 4);
      }
      if (currMB->cbp & 0x08)
      {
        inverse4x4(cof, mb_rres, 8, 8);
        inverse4x4(cof, mb_rres, 8, 12);
        inverse4x4(cof, mb_rres, 12, 8);
        inverse4x4(cof, mb_rres, 12, 12);
      }
    }
    else
    {
      for (jj = 0; jj < MB_BLOCK_SIZE; jj += BLOCK_SIZE)
      {
        inverse4x4(cof, mb_rres, jj, 0);
        inverse4x4(cof, mb_rres, jj, 4);
        inverse4x4(cof, mb_rres, jj, 8);
        inverse4x4(cof, mb_rres, jj, 12);
      }
    }
    sample_reconstruct (currSlice->mb_rec[pl], currSlice->mb_pred[pl], mb_rres, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE, currMB->p_Vid->max_pel_value_comp[pl], DQ_BITS);
  }

  // construct picture from 4x4 blocks
  copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_rec[pl], currMB->pix_x, 0);
}

void iMBtrans8x8(Macroblock *currMB, ColorPlane pl)
{
  //VideoParameters *p_Vid = currMB->p_Vid;
  StorablePicture *dec_picture = currMB->p_Slice->dec_picture;
  imgpel **curr_img = pl ? dec_picture->imgUV[pl - 1]: dec_picture->imgY;

  // Perform 8x8 idct
  if (currMB->cbp & 0x01) 
    itrans8x8(currMB, pl, 0, 0);
  else
    icopy8x8(currMB, pl, 0, 0);

  if (currMB->cbp & 0x02) 
    itrans8x8(currMB, pl, 8, 0);
  else
    icopy8x8(currMB, pl, 8, 0);

  if (currMB->cbp & 0x04) 
    itrans8x8(currMB, pl, 0, 8);
  else
    icopy8x8(currMB, pl, 0, 8);

  if (currMB->cbp & 0x08) 
    itrans8x8(currMB, pl, 8, 8);
  else
    icopy8x8(currMB, pl, 8, 8);

  copy_image_data_16x16(&curr_img[currMB->pix_y], currMB->p_Slice->mb_rec[pl], currMB->pix_x, 0);
}

void iTransform(Macroblock *currMB, ColorPlane pl, int smb)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  StorablePicture *dec_picture = currSlice->dec_picture;
  imgpel **curr_img;
  int uv = pl-1; 

  if ((currMB->cbp & 15) != 0 || smb)
  {
    if(currMB->luma_transform_size_8x8_flag == 0) // 4x4 inverse transform
    {
      iMBtrans4x4(currMB, pl, smb); 
    }
    else // 8x8 inverse transform
    {  
      iMBtrans8x8(currMB, pl);    
    }
  }
  else
  {
    curr_img = pl ? dec_picture->imgUV[uv] : dec_picture->imgY;
    copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_pred[pl], currMB->pix_x, 0);
  }
  if(smb)
    currSlice->is_reset_coeff = FALSE;

  if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
  {
    imgpel **curUV;
    int b8;
    int ioff, joff;
    imgpel **mb_rec;

    for(uv = PLANE_U; uv <= PLANE_V; ++uv)
    {
      // =============== 4x4 itrans ================
      // -------------------------------------------
      curUV = &dec_picture->imgUV[uv - 1][currMB->pix_c_y]; 
      mb_rec = currSlice->mb_rec[uv];

      if (!smb && (currMB->cbp >> 4))
      {
        if (currMB->is_lossless == FALSE)
        {
          const unsigned char *x_pos, *y_pos;

          for (b8 = 0; b8 < (p_Vid->num_uv_blocks); ++b8)
          {
            x_pos = subblk_offset_x[1][b8];
            y_pos = subblk_offset_y[1][b8];

            itrans4x4(currMB, (ColorPlane)uv, *x_pos++, *y_pos++);
            itrans4x4(currMB, (ColorPlane)uv, *x_pos++, *y_pos++);
            itrans4x4(currMB, (ColorPlane)uv, *x_pos++, *y_pos++);
            itrans4x4(currMB, (ColorPlane)uv, *x_pos  , *y_pos  );
          }
          sample_reconstruct (mb_rec, currSlice->mb_pred[uv], currSlice->mb_rres[uv], 0, 0, 
            p_Vid->mb_size[1][0], p_Vid->mb_size[1][1], currMB->p_Vid->max_pel_value_comp[uv], DQ_BITS);
        }
        else
        {
          const unsigned char *x_pos, *y_pos;
          for (b8 = 0; b8 < (p_Vid->num_uv_blocks); ++b8)
          {
            int i,j;
            x_pos = subblk_offset_x[1][b8];
            y_pos = subblk_offset_y[1][b8];

            for (i = 0 ; i < p_Vid->mb_cr_size_y ; i ++)
            {
              for (j = 0 ; j < p_Vid->mb_cr_size_x ; j ++)
              {
                currSlice->mb_rres[uv][i][j] = currSlice->cof[uv][i][j] ;
              }
            }

            itrans4x4_ls(currMB, (ColorPlane)uv, *x_pos++, *y_pos++);
            itrans4x4_ls(currMB, (ColorPlane)uv, *x_pos++, *y_pos++);
            itrans4x4_ls(currMB, (ColorPlane)uv, *x_pos++, *y_pos++);
            itrans4x4_ls(currMB, (ColorPlane)uv, *x_pos  , *y_pos  );
          }
        }
        copy_image_data(curUV, mb_rec, currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);

        currSlice->is_reset_coeff_cr = FALSE;
      }
      else if (smb)
      {
        currMB->itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : itrans4x4_ls;
        itrans_sp_cr(currMB, uv - 1);

        for (joff = 0; joff < p_Vid->mb_cr_size_y; joff += BLOCK_SIZE)
        {
          for(ioff = 0; ioff < p_Vid->mb_cr_size_x ;ioff += BLOCK_SIZE)
          {
            currMB->itrans_4x4(currMB, (ColorPlane)uv, ioff, joff);
          }
        }

        copy_image_data(curUV, mb_rec, currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
        currSlice->is_reset_coeff_cr = FALSE;
      }
      else 
      {
        copy_image_data(curUV, currSlice->mb_pred[uv], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
      }
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (16x16)
 *************************************************************************************
 */
void copy_image_data_16x16(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{
  int j;
  for(j = 0; j < MB_BLOCK_SIZE; j += 4)
  { 
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), MB_BLOCK_SIZE * sizeof (imgpel));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (8x8)
 *************************************************************************************
 */
void copy_image_data_8x8(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{  
  int j;
  for(j = 0; j < BLOCK_SIZE_8x8; j+=4)
  {
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
  }
}


/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (4x4)
 *************************************************************************************
 */
void copy_image_data_4x4(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{
  memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
  memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
  memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE * sizeof (imgpel));
  memcpy((*imgBuf1   + off1), (*imgBuf2   + off2), BLOCK_SIZE * sizeof (imgpel));
}


/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (8x8)
 *************************************************************************************
 */
void copy_image_data(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2, int width, int height)
{
  int j;
  for(j = 0; j < height; ++j)
  {
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), width * sizeof (imgpel));
  }
}




void forward4x4(int **block, int **tblock, int pos_y, int pos_x)
{
  int i, ii;  
  int tmp[16];
  int *pTmp = tmp, *pblock;
  int p0,p1,p2,p3;
  int t0,t1,t2,t3;

  // Horizontal
  for (i=pos_y; i < pos_y + BLOCK_SIZE; i++)
  {
    pblock = &block[i][pos_x];
    p0 = *(pblock++);
    p1 = *(pblock++);
    p2 = *(pblock++);
    p3 = *(pblock  );

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    *(pTmp++) =  t0 + t1;
    *(pTmp++) = (t3 << 1) + t2;
    *(pTmp++) =  t0 - t1;    
    *(pTmp++) =  t3 - (t2 << 1);
  }

  // Vertical 
  for (i=0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    p0 = *pTmp;
    p1 = *(pTmp += BLOCK_SIZE);
    p2 = *(pTmp += BLOCK_SIZE);
    p3 = *(pTmp += BLOCK_SIZE);

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    ii = pos_x + i;
    tblock[pos_y    ][ii] = t0 +  t1;
    tblock[pos_y + 1][ii] = t2 + (t3 << 1);
    tblock[pos_y + 2][ii] = t0 -  t1;
    tblock[pos_y + 3][ii] = t3 - (t2 << 1);
  }
}

void inverse4x4(int **tblock, int **block, int pos_y, int pos_x)
{
  int i, ii;  
  int tmp[16];
  int *pTmp = tmp, *pblock;
  int p0,p1,p2,p3;
  int t0,t1,t2,t3;

  // Horizontal
  for (i = pos_y; i < pos_y + BLOCK_SIZE; i++)
  {
    pblock = &tblock[i][pos_x];
    t0 = *(pblock++);
    t1 = *(pblock++);
    t2 = *(pblock++);
    t3 = *(pblock  );

    p0 =  t0 + t2;
    p1 =  t0 - t2;
    p2 = (t1 >> 1) - t3;
    p3 =  t1 + (t3 >> 1);

    *(pTmp++) = p0 + p3;
    *(pTmp++) = p1 + p2;
    *(pTmp++) = p1 - p2;
    *(pTmp++) = p0 - p3;
  }

  //  Vertical 
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    t0 = *pTmp;
    t1 = *(pTmp += BLOCK_SIZE);
    t2 = *(pTmp += BLOCK_SIZE);
    t3 = *(pTmp += BLOCK_SIZE);

    p0 = t0 + t2;
    p1 = t0 - t2;
    p2 =(t1 >> 1) - t3;
    p3 = t1 + (t3 >> 1);

    ii = i + pos_x;
    block[pos_y    ][ii] = p0 + p3;
    block[pos_y + 1][ii] = p1 + p2;
    block[pos_y + 2][ii] = p1 - p2;
    block[pos_y + 3][ii] = p0 - p3;
  }
}



void ihadamard4x4(int **tblock, int **block)
{
  int i;  
  int tmp[16];
  int *pTmp = tmp, *pblock;
  int p0,p1,p2,p3;
  int t0,t1,t2,t3;

  // Horizontal
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pblock = tblock[i];
    t0 = *(pblock++);
    t1 = *(pblock++);
    t2 = *(pblock++);
    t3 = *(pblock  );

    p0 = t0 + t2;
    p1 = t0 - t2;
    p2 = t1 - t3;
    p3 = t1 + t3;

    *(pTmp++) = p0 + p3;
    *(pTmp++) = p1 + p2;
    *(pTmp++) = p1 - p2;
    *(pTmp++) = p0 - p3;
  }

  //  Vertical 
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    t0 = *pTmp;
    t1 = *(pTmp += BLOCK_SIZE);
    t2 = *(pTmp += BLOCK_SIZE);
    t3 = *(pTmp += BLOCK_SIZE);

    p0 = t0 + t2;
    p1 = t0 - t2;
    p2 = t1 - t3;
    p3 = t1 + t3;
    
    block[0][i] = p0 + p3;
    block[1][i] = p1 + p2;
    block[2][i] = p1 - p2;
    block[3][i] = p0 - p3;
  }
}

void ihadamard4x2(int **tblock, int **block)
{
  int i;  
  int tmp[8];
  int *pTmp = tmp;
  int p0,p1,p2,p3;
  int t0,t1,t2,t3;

  // Horizontal
  *(pTmp++) = tblock[0][0] + tblock[1][0];
  *(pTmp++) = tblock[0][1] + tblock[1][1];
  *(pTmp++) = tblock[0][2] + tblock[1][2];
  *(pTmp++) = tblock[0][3] + tblock[1][3];

  *(pTmp++) = tblock[0][0] - tblock[1][0];
  *(pTmp++) = tblock[0][1] - tblock[1][1];
  *(pTmp++) = tblock[0][2] - tblock[1][2];
  *(pTmp  ) = tblock[0][3] - tblock[1][3];

  // Vertical
  pTmp = tmp;
  for (i = 0; i < 2; i++)
  {
    p0 = *(pTmp++);
    p1 = *(pTmp++);
    p2 = *(pTmp++);
    p3 = *(pTmp++);

    t0 = p0 + p2;
    t1 = p0 - p2;
    t2 = p1 - p3;
    t3 = p1 + p3;

    // coefficients (transposed)
    block[0][i] = t0 + t3;
    block[1][i] = t1 + t2;
    block[2][i] = t1 - t2;
    block[3][i] = t0 - t3;
  }
}

void ihadamard2x2(int tblock[4], int block[4])
{
  int t0,t1,t2,t3;

  t0 = tblock[0] + tblock[1];
  t1 = tblock[0] - tblock[1];
  t2 = tblock[2] + tblock[3];
  t3 = tblock[2] - tblock[3];

  block[0] = (t0 + t2);
  block[1] = (t1 + t3);
  block[2] = (t0 - t2);
  block[3] = (t1 - t3);
}

void inverse8x8(int **tblock, int **block, int pos_x)
{
  int i, ii;
  int tmp[64];
  int *pTmp = tmp, *pblock;
  int a0, a1, a2, a3;
  int p0, p1, p2, p3, p4, p5 ,p6, p7;  
  int b0, b1, b2, b3, b4, b5, b6, b7;

  // Horizontal  
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    pblock = &tblock[i][pos_x];
    p0 = *(pblock++);
    p1 = *(pblock++);
    p2 = *(pblock++);
    p3 = *(pblock++);
    p4 = *(pblock++);
    p5 = *(pblock++);
    p6 = *(pblock++);
    p7 = *(pblock  );

    a0 = p0 + p4;
    a1 = p0 - p4;
    a2 = p6 - (p2 >> 1);
    a3 = p2 + (p6 >> 1);

    b0 =  a0 + a3;
    b2 =  a1 - a2;
    b4 =  a1 + a2;
    b6 =  a0 - a3;

    a0 = -p3 + p5 - p7 - (p7 >> 1);    
    a1 =  p1 + p7 - p3 - (p3 >> 1);    
    a2 = -p1 + p7 + p5 + (p5 >> 1);    
    a3 =  p3 + p5 + p1 + (p1 >> 1);

    
    b1 =  a0 + (a3>>2);    
    b3 =  a1 + (a2>>2);    
    b5 =  a2 - (a1>>2);
    b7 =  a3 - (a0>>2);                

    *(pTmp++) = b0 + b7;
    *(pTmp++) = b2 - b5;
    *(pTmp++) = b4 + b3;
    *(pTmp++) = b6 + b1;
    *(pTmp++) = b6 - b1;
    *(pTmp++) = b4 - b3;
    *(pTmp++) = b2 + b5;
    *(pTmp++) = b0 - b7;
  }

  //  Vertical 
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    pTmp = tmp + i;
    p0 = *pTmp;
    p1 = *(pTmp += BLOCK_SIZE_8x8);
    p2 = *(pTmp += BLOCK_SIZE_8x8);
    p3 = *(pTmp += BLOCK_SIZE_8x8);
    p4 = *(pTmp += BLOCK_SIZE_8x8);
    p5 = *(pTmp += BLOCK_SIZE_8x8);
    p6 = *(pTmp += BLOCK_SIZE_8x8);
    p7 = *(pTmp += BLOCK_SIZE_8x8);

    a0 =  p0 + p4;
    a1 =  p0 - p4;
    a2 =  p6 - (p2>>1);
    a3 =  p2 + (p6>>1);

    b0 = a0 + a3;
    b2 = a1 - a2;
    b4 = a1 + a2;
    b6 = a0 - a3;

    a0 = -p3 + p5 - p7 - (p7 >> 1);
    a1 =  p1 + p7 - p3 - (p3 >> 1);
    a2 = -p1 + p7 + p5 + (p5 >> 1);
    a3 =  p3 + p5 + p1 + (p1 >> 1);


    b1 =  a0 + (a3 >> 2);
    b7 =  a3 - (a0 >> 2);
    b3 =  a1 + (a2 >> 2);
    b5 =  a2 - (a1 >> 2);

    ii = i + pos_x;
    block[0][ii] = b0 + b7;
    block[1][ii] = b2 - b5;
    block[2][ii] = b4 + b3;
    block[3][ii] = b6 + b1;
    block[4][ii] = b6 - b1;
    block[5][ii] = b4 - b3;
    block[6][ii] = b2 + b5;
    block[7][ii] = b0 - b7;
  }
}


static void recon8x8(int **m7, imgpel **mb_rec, imgpel **mpr, int max_imgpel_value, int ioff)
{
  int j;
  int    *m_tr  = NULL;
  imgpel *m_rec = NULL;
  imgpel *m_prd = NULL;

  for( j = 0; j < 8; j++)
  {
    m_tr = (*m7++) + ioff;
    m_rec = (*mb_rec++) + ioff;
    m_prd = (*mpr++) + ioff;

    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec++ = (imgpel) iClip1(max_imgpel_value, (*m_prd++) + rshift_rnd_sf(*m_tr++, DQ_BITS_8)); 
    *m_rec   = (imgpel) iClip1(max_imgpel_value, (*m_prd  ) + rshift_rnd_sf(*m_tr  , DQ_BITS_8)); 
  }
}

static void copy8x8(imgpel **mb_rec, imgpel **mpr, int ioff)
{
  int j;

  for( j = 0; j < 8; j++)
  {
    memcpy((*mb_rec++) + ioff, (*mpr++) + ioff, 8 * sizeof(imgpel));
  }
}

static void recon8x8_lossless(int **m7, imgpel **mb_rec, imgpel **mpr, int max_imgpel_value, int ioff)
{
  int i, j;
  for( j = 0; j < 8; j++)
  {
    for( i = ioff; i < ioff + 8; i++)
      (*mb_rec)[i] = (imgpel) iClip1(max_imgpel_value, ((*m7)[i] + (long)(*mpr)[i])); 
    mb_rec++;
    m7++;
    mpr++;
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Inverse 8x8 transformation
 ***********************************************************************
 */ 
void itrans8x8(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane       
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  Slice *currSlice = currMB->p_Slice;

  int    **m7     = currSlice->mb_rres[pl];

  if (currMB->is_lossless == TRUE)
  {
    recon8x8_lossless(&m7[joff], &currSlice->mb_rec[pl][joff], &currSlice->mb_pred[pl][joff], currMB->p_Vid->max_pel_value_comp[pl], ioff);
  }
  else
  {
    inverse8x8(&m7[joff], &m7[joff], ioff);
    recon8x8  (&m7[joff], &currSlice->mb_rec[pl][joff], &currSlice->mb_pred[pl][joff], currMB->p_Vid->max_pel_value_comp[pl], ioff);
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Inverse 8x8 transformation
 ***********************************************************************
 */ 
void icopy8x8(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane       
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  Slice *currSlice = currMB->p_Slice;

  copy8x8  (&currSlice->mb_rec[pl][joff], &currSlice->mb_pred[pl][joff], ioff);
}
