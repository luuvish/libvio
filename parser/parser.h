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
 *  File      : parser.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 15, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _PARSER_H_
#define _PARSER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bitstream.h"

struct macroblock_dec;
struct syntaxelement_dec;

//! struct to characterize the state of the arithmetic coding engine
typedef struct {
    unsigned int  Drange;
    unsigned int  Dvalue;
    int           DbitsLeft;
    byte         *Dcodestrm;
    int          *Dcodestrm_len;
} DecodingEnvironment;

typedef DecodingEnvironment *DecodingEnvironmentPtr;

//! Syntaxelement
typedef struct syntaxelement_dec {
    int           type;                  //!< type of syntax element for data part.
    int           value1;                //!< numerical value of syntax element
    int           value2;                //!< for blocked symbols, e.g. run/level
    int           len;                   //!< length of code
    int           inf;                   //!< info part of CAVLC code
    unsigned int  bitpattern;            //!< CAVLC bitpattern
    int           context;               //!< CABAC context
    int           k;                     //!< CABAC context for coeff_count,uv

#if TRACE
    #define       TRACESTRING_SIZE 100           //!< size of trace string
    char          tracestring[TRACESTRING_SIZE]; //!< trace string
#endif

    //! for mapping of CAVLC to syntaxElement
    void (*mapping)(int len, int info, int *value1, int *value2);
    //! used for CABAC: refers to actual coding method of each individual syntax element type
    void (*reading)(struct macroblock_dec *currMB, struct syntaxelement_dec *, DecodingEnvironmentPtr);
} SyntaxElement;

//! DataPartition
typedef struct datapartition_dec {
    Bitstream           *bitstream;
    DecodingEnvironment  de_cabac;

    int (*readSyntaxElement)(struct macroblock_dec *currMB, struct syntaxelement_dec *, struct datapartition_dec *);
          /*!< virtual function;
               actual method depends on chosen data partition and
               entropy coding method  */
} DataPartition;

extern DataPartition *AllocPartition(int n);
extern void FreePartition(DataPartition *dp, int n);

#ifdef __cplusplus
}
#endif

#endif /* _PARSER_H_ */
