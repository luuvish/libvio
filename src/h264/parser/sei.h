#ifndef _SEI_H_
#define _SEI_H_

// tone mapping information
#define MAX_CODED_BIT_DEPTH 12
#define MAX_SEI_BIT_DEPTH   12
#define MAX_NUM_PIVOTS       (1<<MAX_CODED_BIT_DEPTH)

typedef struct tone_mapping_struct_s {
    bool          seiHasTone_mapping;
    unsigned int  tone_map_repetition_period;
    unsigned char coded_data_bit_depth;
    unsigned char sei_bit_depth;
    unsigned int  model_id;
    unsigned int  count;

    px_t lut[1<<MAX_CODED_BIT_DEPTH];       //<! look up table for mapping the coded data value to output data value

    int        payloadSize;
} ToneMappingSEI;

struct slice_t;

void parse_sei(byte *payload, int size, VideoParameters *p_Vid, struct slice_t *pSlice);

void tone_map               (px_t **imgX, px_t *lut, int size_x, int size_y);
void init_tone_mapping_sei  (ToneMappingSEI *seiToneMapping);
void update_tone_mapping_sei(ToneMappingSEI *seiToneMapping);

#endif /* _SEI_H_ */
