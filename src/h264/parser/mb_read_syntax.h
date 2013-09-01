#ifndef _MB_READ_SYNTAX_H_
#define _MB_READ_SYNTAX_H_

#include "global.h"

uint32_t parse_mb_skip_run            (mb_t *mb);
bool     parse_mb_skip_flag           (mb_t *mb);
bool     parse_mb_field_decoding_flag (mb_t *mb);
uint32_t parse_mb_type                (mb_t *mb);
uint8_t  parse_sub_mb_type            (mb_t *mb);

bool     parse_transform_size_8x8_flag(mb_t *mb);
int8_t   parse_intra_pred_mode        (mb_t *mb);
uint8_t  parse_intra_chroma_pred_mode (mb_t *mb);
uint8_t  parse_ref_idx                (mb_t *mb, uint8_t list);
int16_t  parse_mvd                    (mb_t *mb, uint8_t xy, uint8_t list);
uint8_t  parse_coded_block_pattern    (mb_t *mb);
int8_t   parse_mb_qp_delta            (mb_t *mb);


#endif /* _MB_READ_SYNTAX_H_ */
