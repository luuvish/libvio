#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "biaridecod.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "bitstream_elements.h"
#include "bitstream_ctx.h"


/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the macroblock
 *    type info of a given MB.
 ************************************************************************
 */
static void readMB_typeInfo_CABAC_i_slice(Macroblock *currMB,  
                                          SyntaxElement *se,
                                          DecodingEnvironmentPtr dep_dp)
{
    Slice *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;

    int a = 0, b = 0;
    int act_ctx;
    int act_sym;
    int mode_sym;
    int curr_mb_type = 0;

    if (currSlice->slice_type == I_SLICE) { // INTRA-frame
        if (currMB->mb_up != NULL)
            b = (currMB->mb_up->mb_type != I4MB && currMB->mb_up->mb_type != I8MB ? 1 : 0 );
        if (currMB->mb_left != NULL)
            a = (currMB->mb_left->mb_type != I4MB && currMB->mb_left->mb_type != I8MB ? 1 : 0);

        act_ctx = a + b;
        act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
        se->context = act_ctx; // store context

        if (act_sym == 0) // 4x4 Intra
            curr_mb_type = act_sym;
        else { // 16x16 Intra
            mode_sym = biari_decode_final(dep_dp);
            if (mode_sym == 1)
                curr_mb_type = 25;
            else {
                act_sym = 1;
                act_ctx = 4;
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx); // decoding of AC/no AC
                act_sym += mode_sym * 12;
                act_ctx = 5;
                // decoding of cbp: 0,1,2
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                if (mode_sym != 0) {
                    act_ctx = 6;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += 4;
                    if (mode_sym != 0)
                        act_sym += 4;
                }
                // decoding of I pred-mode: 0,1,2,3
                act_ctx = 7;
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                act_sym += mode_sym * 2;
                act_ctx = 8;
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                act_sym += mode_sym;
                curr_mb_type = act_sym;
            }
        }
    } else if (currSlice->slice_type == SI_SLICE) { // SI-frame
        // special ctx's for SI4MB
        if (currMB->mb_up != NULL)
            b = (( (currMB->mb_up)->mb_type != SI4MB) ? 1 : 0 );
        if (currMB->mb_left != NULL)
            a = (( (currMB->mb_left)->mb_type != SI4MB) ? 1 : 0 );

        act_ctx = a + b;
        act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx);
        se->context = act_ctx; // store context

        if (act_sym == 0) //  SI 4x4 Intra
            curr_mb_type = 0;
        else { // analog INTRA_IMG
            if (currMB->mb_up != NULL)
                b = (( (currMB->mb_up)->mb_type != I4MB) ? 1 : 0 );
            if (currMB->mb_left != NULL)
                a = (( (currMB->mb_left)->mb_type != I4MB) ? 1 : 0 );

            act_ctx = a + b;
            act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
            se->context = act_ctx; // store context

            if (act_sym == 0) // 4x4 Intra
                curr_mb_type = 1;
            else { // 16x16 Intra
                mode_sym = biari_decode_final(dep_dp);
                if (mode_sym == 1)
                    curr_mb_type = 26;
                else {
                    act_sym = 2;
                    act_ctx = 4;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx); // decoding of AC/no AC
                    act_sym += mode_sym * 12;
                    act_ctx = 5;
                    // decoding of cbp: 0,1,2
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    if (mode_sym != 0) {
                        act_ctx = 6;
                        mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                        act_sym += 4;
                        if (mode_sym != 0)
                          act_sym += 4;
                    }
                    // decoding of I pred-mode: 0,1,2,3
                    act_ctx = 7;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += mode_sym * 2;
                    act_ctx = 8;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += mode_sym;
                    curr_mb_type = act_sym;
                }
            }
        }
    }

    se->value1 = curr_mb_type;
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the macroblock
 *    type info of a given MB.
 ************************************************************************
 */
static void readMB_typeInfo_CABAC_p_slice(Macroblock *currMB,  
                           SyntaxElement *se,
                           DecodingEnvironmentPtr dep_dp)
{
  Slice *currSlice = currMB->p_Slice;
  MotionInfoContexts *ctx = currSlice->mot_ctx;

  int act_ctx;
  int act_sym;
  int mode_sym;
  int curr_mb_type;
  BiContextType *mb_type_contexts = ctx->mb_type_contexts[1];

  if (biari_decode_symbol(dep_dp, &mb_type_contexts[4] ))
  {
    if (biari_decode_symbol(dep_dp, &mb_type_contexts[7] ))   
      act_sym = 7;
    else                                                              
      act_sym = 6;
  }
  else
  {
    if (biari_decode_symbol(dep_dp, &mb_type_contexts[5] ))
    {
      if (biari_decode_symbol(dep_dp, &mb_type_contexts[7] )) 
        act_sym = 2;
      else
        act_sym = 3;
    }
    else
    {
      if (biari_decode_symbol(dep_dp, &mb_type_contexts[6] ))
        act_sym = 4;
      else                                                            
        act_sym = 1;
    }
  }

  if (act_sym <= 6)
  {
    curr_mb_type = act_sym;
  }
  else  // additional info for 16x16 Intra-mode
  {
    mode_sym = biari_decode_final(dep_dp);
    if( mode_sym==1 )
    {
      curr_mb_type = 31;
    }
    else
    {
      act_ctx = 8;
      mode_sym =  biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx ); // decoding of AC/no AC
      act_sym += mode_sym*12;

      // decoding of cbp: 0,1,2
      act_ctx = 9;
      mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
      if (mode_sym != 0)
      {
        act_sym+=4;
        mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
        if (mode_sym != 0)
          act_sym+=4;
      }

      // decoding of I pred-mode: 0,1,2,3
      act_ctx = 10;
      mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
      act_sym += mode_sym*2;
      mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
      act_sym += mode_sym;
      curr_mb_type = act_sym;
    }
  }

  se->value1 = curr_mb_type;
}


/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the macroblock
 *    type info of a given MB.
 ************************************************************************
 */
static void readMB_typeInfo_CABAC_b_slice(Macroblock *currMB,  
                           SyntaxElement *se,
                           DecodingEnvironmentPtr dep_dp)
{
  Slice *currSlice = currMB->p_Slice;
  MotionInfoContexts *ctx = currSlice->mot_ctx;

  int a = 0, b = 0;
  int act_ctx;
  int act_sym;
  int mode_sym;
  int curr_mb_type;
  BiContextType *mb_type_contexts = ctx->mb_type_contexts[2];

  if (currMB->mb_up != NULL)
    b = (( (currMB->mb_up)->mb_type != 0) ? 1 : 0 );

  if (currMB->mb_left != NULL)
    a = (( (currMB->mb_left)->mb_type != 0) ? 1 : 0 );

  act_ctx = a + b;

  if (biari_decode_symbol (dep_dp, &mb_type_contexts[act_ctx]))
  {
    if (biari_decode_symbol (dep_dp, &mb_type_contexts[4]))
    {
      if (biari_decode_symbol (dep_dp, &mb_type_contexts[5]))
      {
        act_sym = 12;
        if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
          act_sym += 8;
        if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
          act_sym += 4;
        if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
          act_sym += 2;

        if      (act_sym == 24)  
          act_sym=11;
        else if (act_sym == 26)  
          act_sym = 22;
        else
        {
          if (act_sym == 22)     
            act_sym = 23;
          if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
            act_sym += 1;
        }
      }
      else
      {
        act_sym = 3;
        if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
          act_sym += 4;
        if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
          act_sym += 2;
        if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
          act_sym += 1;
      }
    }
    else
    {
      if (biari_decode_symbol (dep_dp, &mb_type_contexts[6])) 
        act_sym=2;
      else
        act_sym=1;
    }
  }
  else
  {
    act_sym = 0;
  }


  if (act_sym <= 23)
  {
    curr_mb_type = act_sym;
  }
  else  // additional info for 16x16 Intra-mode
  {
    mode_sym = biari_decode_final(dep_dp);
    if( mode_sym == 1 )
    {
      curr_mb_type = 48;
    }
    else
    {
      mb_type_contexts = ctx->mb_type_contexts[1];
      act_ctx = 8;
      mode_sym =  biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx ); // decoding of AC/no AC
      act_sym += mode_sym*12;

      // decoding of cbp: 0,1,2
      act_ctx = 9;
      mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
      if (mode_sym != 0)
      {
        act_sym+=4;
        mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
        if (mode_sym != 0)
          act_sym+=4;
      }

      // decoding of I pred-mode: 0,1,2,3
      act_ctx = 10;
      mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
      act_sym += mode_sym*2;
      mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
      act_sym += mode_sym;
      curr_mb_type = act_sym;
    }
  }

  se->value1 = curr_mb_type;
}


int getSE(Macroblock *currMB, int type)
{
    Slice *currSlice = currMB->p_Slice;  
    const byte *partMap = assignSE2partition[currSlice->dp_mode];
    bool isCabac = currSlice->p_Vid->active_pps->entropy_coding_mode_flag;

    if (type == CTX_MB_TYPE) {
    	DataPartition *dP = &currSlice->partArr[partMap[SE_MBTYPE]];
    	SyntaxElement currSE;
    	currSE.type = SE_MBTYPE;
        if (!isCabac || dP->bitstream->ei_flag)   
            currSE.mapping = linfo_ue;
        else {
        	switch (currSlice->slice_type) {
        	case I_SLICE:
        	case SI_SLICE:
        		currSE.reading = readMB_typeInfo_CABAC_i_slice;
        		break;
        	case P_SLICE:
        	case SP_SLICE:
            	currSE.reading = readMB_typeInfo_CABAC_p_slice;
            	break;
        	case B_SLICE:
            	currSE.reading = readMB_typeInfo_CABAC_b_slice;
            	break;
            }
		}
        dP->readSyntaxElement(currMB, &currSE, dP);
	    if (!dP->bitstream->ei_flag)
	        currMB->ei_flag = 0;
        return currSE.value1;
    }

    return 0;
}
