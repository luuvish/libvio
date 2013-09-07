#ifndef _QUANTIZATION_H_
#define _QUANTIZATION_H_


//! Dequantization coefficients
static const int dequant_coef[6][4][4] = {
    {{ 10, 13, 10, 13 },
     { 13, 16, 13, 16 },
     { 10, 13, 10, 13 },
     { 13, 16, 13, 16 }},
    {{ 11, 14, 11, 14 },
     { 14, 18, 14, 18 },
     { 11, 14, 11, 14 },
     { 14, 18, 14, 18 }},
    {{ 13, 16, 13, 16 },
     { 16, 20, 16, 20 },
     { 13, 16, 13, 16 },
     { 16, 20, 16, 20 }},
    {{ 14, 18, 14, 18 },
     { 18, 23, 18, 23 },
     { 14, 18, 14, 18 },
     { 18, 23, 18, 23 }},
    {{ 16, 20, 16, 20 },
     { 20, 25, 20, 25 },
     { 16, 20, 16, 20 },
     { 20, 25, 20, 25 }},
    {{ 18, 23, 18, 23 },
     { 23, 29, 23, 29 },
     { 18, 23, 18, 23 },
     { 23, 29, 23, 29 }}
};

static const int quant_coef[6][4][4] = {
    {{ 13107,  8066, 13107,  8066 },
     {  8066,  5243,  8066,  5243 },
     { 13107,  8066, 13107,  8066 },
     {  8066,  5243,  8066,  5243 }},
    {{ 11916,  7490, 11916,  7490 },
     {  7490,  4660,  7490,  4660 },
     { 11916,  7490, 11916,  7490 },
     {  7490,  4660,  7490,  4660 }},
    {{ 10082,  6554, 10082,  6554 },
     {  6554,  4194,  6554,  4194 },
     { 10082,  6554, 10082,  6554 },
     {  6554,  4194,  6554,  4194 }},
    {{  9362,  5825,  9362,  5825 },
     {  5825,  3647,  5825,  3647 },
     {  9362,  5825,  9362,  5825 },
     {  5825,  3647,  5825,  3647 }},
    {{  8192,  5243,  8192,  5243 },
     {  5243,  3355,  5243,  3355 },
     {  8192,  5243,  8192,  5243 },
     {  5243,  3355,  5243,  3355 }},
    {{  7282,  4559,  7282,  4559 },
     {  4559,  2893,  4559,  2893 },
     {  7282,  4559,  7282,  4559 },
     {  4559,  2893,  4559,  2893 }}
};

struct slice_t;
struct macroblock_t;

// For Q-matrix
void assign_quant_params(slice_t* slice);


struct quantization_t {
    void coeff_luma_dc  (macroblock_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void coeff_luma_ac  (macroblock_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void coeff_chroma_dc(macroblock_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void coeff_chroma_ac(macroblock_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);

    int  inverse_quantize(macroblock_t* mb, bool uv, ColorPlane pl, int i0, int j0, int levarr);
};

extern quantization_t quantization;


#endif /* _QUANTIZATION_H_ */
