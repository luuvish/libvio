//
//  ocl_trquant.cl
//  HM
//
//  Created by Injo Hwang on 12. 9. 29..
//
//

#define MAX_CU_DEPTH 7                   // log2(LCUSize)
#define MAX_CU_SIZE  (1<<(MAX_CU_DEPTH)) // maximum allowable size of CU


#define REG_DCT 65535

#define SHIFT_INV_1ST  7 // Shift after first inverse transform stage
#define SHIFT_INV_2ND 12 // Shift after second inverse transform stage

#define Clip3(min, max, a) ((min) > (a) ? (min) : (max) < (a) ? (max) : (a))

__constant char g_aucConvertToBit[MAX_CU_SIZE+1] = { // from width to log2(width)-2
    -1,
    -1,-1,-1, 0,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 2, // i=  4 -> c=0
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, // i=  8 -> c=1
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // i= 16 -> c=2
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4, // i= 32 -> c=3
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // i= 64 -> c=4
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5  // i=128 -> c=5
};
__constant int g_uiBitIncrement = 0;


__constant short g_aiT4[4][4] = {
    { 64, 64, 64, 64},
    { 83, 36,-36,-83},
    { 64,-64,-64, 64},
    { 36,-83, 83,-36}
};

__constant short g_aiT8[8][8] = {
    { 64, 64, 64, 64, 64, 64, 64, 64},
    { 89, 75, 50, 18,-18,-50,-75,-89},
    { 83, 36,-36,-83,-83,-36, 36, 83},
    { 75,-18,-89,-50, 50, 89, 18,-75},
    { 64,-64,-64, 64, 64,-64,-64, 64},
    { 50,-89, 18, 75,-75,-18, 89,-50},
    { 36,-83, 83,-36,-36, 83,-83, 36},
    { 18,-50, 75,-89, 89,-75, 50,-18}
};

__constant short g_aiT16[16][16] = {
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
    { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
    { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
    { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
    { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
    { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
    { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
    { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
    { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
    { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
    { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
    { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
    { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
    { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
    { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
    {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

__constant short g_aiT32[32][32] = {
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
    { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4,
      -4,-13,-22,-31,-38,-46,-54,-61,-67,-73,-78,-82,-85,-88,-90,-90},
    { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90,
     -90,-87,-80,-70,-57,-43,-25, -9,  9, 25, 43, 57, 70, 80, 87, 90},
    { 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13,
      13, 38, 61, 78, 88, 90, 85, 73, 54, 31,  4,-22,-46,-67,-82,-90},
    { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89,
      89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
    { 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22,
     -22,-61,-85,-90,-73,-38,  4, 46, 78, 90, 82, 54, 13,-31,-67,-88},
    { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87,
     -87,-57, -9, 43, 80, 90, 70, 25,-25,-70,-90,-80,-43,  9, 57, 87},
    { 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31,
      31, 78, 90, 61,  4,-54,-88,-82,-38, 22, 73, 90, 67, 13,-46,-85},
    { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83,
      83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
    { 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38,
     -38,-88,-73, -4, 67, 90, 46,-31,-85,-78,-13, 61, 90, 54,-22,-82},
    { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80,
     -80, -9, 70, 87, 25,-57,-90,-43, 43, 90, 57,-25,-87,-70,  9, 80},
    { 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46,
      46, 90, 38,-54,-90,-31, 61, 88, 22,-67,-85,-13, 73, 82,  4,-78},
    { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75,
      75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
    { 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54,
     -54,-85,  4, 88, 46,-61,-82, 13, 90, 38,-67,-78, 22, 90, 31,-73},
    { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70,
     -70, 43, 87, -9,-90,-25, 80, 57,-57,-80, 25, 90,  9,-87,-43, 70},
    { 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61,
      61, 73,-46,-82, 31, 88,-13,-90, -4, 90, 22,-85,-38, 78, 54,-67},
    { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64,
      64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
    { 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67,
     -67,-54, 78, 38,-85,-22, 90,  4,-90, 13, 88,-31,-82, 46, 73,-61},
    { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57,
     -57, 80, 25,-90,  9, 87,-43,-70, 70, 43,-87, -9, 90,-25,-80, 57},
    { 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73,
      73, 31,-90, 22, 78,-67,-38, 90,-13,-82, 61, 46,-88,  4, 85,-54},
    { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50,
      50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
    { 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78,
     -78, -4, 82,-73,-13, 85,-67,-22, 88,-61,-31, 90,-54,-38, 90,-46},
    { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43,
     -43, 90,-57,-25, 87,-70, -9, 80,-80,  9, 70,-87, 25, 57,-90, 43},
    { 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82,
      82,-22,-54, 90,-61,-13, 78,-85, 31, 46,-90, 67,  4,-73, 88,-38},
    { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36,
      36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
    { 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85,
     -85, 46, 13,-67, 90,-73, 22, 38,-82, 88,-54, -4, 61,-90, 78,-31},
    { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25,
     -25, 70,-90, 80,-43, -9, 57,-87, 87,-57,  9, 43,-80, 90,-70, 25},
    { 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88,
      88,-67, 31, 13,-54, 82,-90, 78,-46,  4, 38,-73, 90,-85, 61,-22},
    { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18,
      18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
    { 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90,
     -90, 82,-67, 46,-22, -4, 31,-54, 73,-85, 90,-88, 78,-61, 38,-13},
    {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9,
      -9, 25,-43, 57,-70, 80,-87, 90,-90, 87,-80, 70,-57, 43,-25,  9},
    {  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90,
      90,-90, 88,-85, 82,-78, 73,-67, 61,-54, 46,-38, 31,-22, 13, -4}
};


// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm
// give identical results
void fastForwardDst_gl(__global short *block, __local short *coeff, int shift)  // input block, output coeff
{
    int i, c[4];
    int rnd_factor = 1<<(shift-1);

    for (i = 0; i < 4; i++) {
        // Intermediate Variables
        c[0] = block[4*i+0] + block[4*i+3];
        c[1] = block[4*i+1] + block[4*i+3];
        c[2] = block[4*i+0] - block[4*i+1];
        c[3] = 74 * block[4*i+2];
        
        coeff[   i] = (29 * c[0] + 55 * c[1]         + c[3]             + rnd_factor) >> shift;
        coeff[ 4+i] = (74 * (block[4*i+0]+ block[4*i+1] - block[4*i+3]) + rnd_factor) >> shift;
        coeff[ 8+i] = (29 * c[2] + 55 * c[0]         - c[3]             + rnd_factor) >> shift;
        coeff[12+i] = (55 * c[2] - 29 * c[1]         + c[3]             + rnd_factor) >> shift;
    }
}
void fastForwardDst_lg(__local short *block, __global short *coeff, int shift)  // input block, output coeff
{
    int i, c[4];
    int rnd_factor = 1<<(shift-1);
    
    for (i = 0; i < 4; i++) {
        // Intermediate Variables
        c[0] = block[4*i+0] + block[4*i+3];
        c[1] = block[4*i+1] + block[4*i+3];
        c[2] = block[4*i+0] - block[4*i+1];
        c[3] = 74 * block[4*i+2];
        
        coeff[   i] = (29 * c[0] + 55 * c[1]         + c[3]             + rnd_factor) >> shift;
        coeff[ 4+i] = (74 * (block[4*i+0]+ block[4*i+1] - block[4*i+3]) + rnd_factor) >> shift;
        coeff[ 8+i] = (29 * c[2] + 55 * c[0]         - c[3]             + rnd_factor) >> shift;
        coeff[12+i] = (55 * c[2] - 29 * c[1]         + c[3]             + rnd_factor) >> shift;
    }
}

void fastInverseDst_gl(__global short *tmp, __local short *block, int shift)  // input tmp, output block
{
    int i, c[4];
    int rnd_factor = 1<<(shift-1);

    for (i = 0; i < 4; i++) {
        // Intermediate Variables
        c[0] = tmp[  i] + tmp[ 8+i];
        c[1] = tmp[8+i] + tmp[12+i];
        c[2] = tmp[  i] - tmp[12+i];
        c[3] = 74 * tmp[4+i];
        
        block[4*i+0] = Clip3(-32768, 32767, (29 * c[0] + 55 * c[1]     + c[3]     + rnd_factor) >> shift);
        block[4*i+1] = Clip3(-32768, 32767, (55 * c[2] - 29 * c[1]     + c[3]     + rnd_factor) >> shift);
        block[4*i+2] = Clip3(-32768, 32767, (74 * (tmp[i] - tmp[8+i] + tmp[12+i]) + rnd_factor) >> shift);
        block[4*i+3] = Clip3(-32768, 32767, (55 * c[0] + 29 * c[2]     - c[3]     + rnd_factor) >> shift);
    }
}
void fastInverseDst_lg(__local short *tmp, __global short *block, int shift)  // input tmp, output block
{
    int i, c[4];
    int rnd_factor = 1<<(shift-1);
    
    for (i = 0; i < 4; i++) {
        // Intermediate Variables
        c[0] = tmp[  i] + tmp[ 8+i];
        c[1] = tmp[8+i] + tmp[12+i];
        c[2] = tmp[  i] - tmp[12+i];
        c[3] = 74 * tmp[4+i];
        
        block[4*i+0] = Clip3(-32768, 32767, (29 * c[0] + 55 * c[1]     + c[3]     + rnd_factor) >> shift);
        block[4*i+1] = Clip3(-32768, 32767, (55 * c[2] - 29 * c[1]     + c[3]     + rnd_factor) >> shift);
        block[4*i+2] = Clip3(-32768, 32767, (74 * (tmp[i] - tmp[8+i] + tmp[12+i]) + rnd_factor) >> shift);
        block[4*i+3] = Clip3(-32768, 32767, (55 * c[0] + 29 * c[2]     - c[3]     + rnd_factor) >> shift);
    }
}

void partialButterfly4_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j;
    int E[2], O[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O */
        E[0] = src[0] + src[3];
        O[0] = src[0] - src[3];
        E[1] = src[1] + src[2];
        O[1] = src[1] - src[2];
        
        dst[0*line] = (g_aiT4[0][0]*E[0] + g_aiT4[0][1]*E[1] + add) >> shift;
        dst[2*line] = (g_aiT4[2][0]*E[0] + g_aiT4[2][1]*E[1] + add) >> shift;
        dst[1*line] = (g_aiT4[1][0]*O[0] + g_aiT4[1][1]*O[1] + add) >> shift;
        dst[3*line] = (g_aiT4[3][0]*O[0] + g_aiT4[3][1]*O[1] + add) >> shift;
        
        src += 4;
        dst++;
    }
}
void partialButterfly4_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j;
    int E[2], O[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O */
        E[0] = src[0] + src[3];
        O[0] = src[0] - src[3];
        E[1] = src[1] + src[2];
        O[1] = src[1] - src[2];
        
        dst[0*line] = (g_aiT4[0][0]*E[0] + g_aiT4[0][1]*E[1] + add) >> shift;
        dst[2*line] = (g_aiT4[2][0]*E[0] + g_aiT4[2][1]*E[1] + add) >> shift;
        dst[1*line] = (g_aiT4[1][0]*O[0] + g_aiT4[1][1]*O[1] + add) >> shift;
        dst[3*line] = (g_aiT4[3][0]*O[0] + g_aiT4[3][1]*O[1] + add) >> shift;
        
        src += 4;
        dst++;
    }
}

void partialButterflyInverse4_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j;
    int E[2], O[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        O[0] = g_aiT4[1][0]*src[1*line] + g_aiT4[3][0]*src[3*line];
        O[1] = g_aiT4[1][1]*src[1*line] + g_aiT4[3][1]*src[3*line];
        E[0] = g_aiT4[0][0]*src[0*line] + g_aiT4[2][0]*src[2*line];
        E[1] = g_aiT4[0][1]*src[0*line] + g_aiT4[2][1]*src[2*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        dst[0] = Clip3(-32768, 32767, (E[0] + O[0] + add) >> shift);
        dst[1] = Clip3(-32768, 32767, (E[1] + O[1] + add) >> shift);
        dst[2] = Clip3(-32768, 32767, (E[1] - O[1] + add) >> shift);
        dst[3] = Clip3(-32768, 32767, (E[0] - O[0] + add) >> shift);
        
        src++;
        dst += 4;
    }
}
void partialButterflyInverse4_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j;
    int E[2], O[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        O[0] = g_aiT4[1][0]*src[1*line] + g_aiT4[3][0]*src[3*line];
        O[1] = g_aiT4[1][1]*src[1*line] + g_aiT4[3][1]*src[3*line];
        E[0] = g_aiT4[0][0]*src[0*line] + g_aiT4[2][0]*src[2*line];
        E[1] = g_aiT4[0][1]*src[0*line] + g_aiT4[2][1]*src[2*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        dst[0] = Clip3(-32768, 32767, (E[0] + O[0] + add) >> shift);
        dst[1] = Clip3(-32768, 32767, (E[1] + O[1] + add) >> shift);
        dst[2] = Clip3(-32768, 32767, (E[1] - O[1] + add) >> shift);
        dst[3] = Clip3(-32768, 32767, (E[0] - O[0] + add) >> shift);
        
        src++;
        dst += 4;
    }
}

void partialButterfly8_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O*/
        for (k = 0; k < 4; k++) {
            E[k] = src[k] + src[7-k];
            O[k] = src[k] - src[7-k];
        }
        /* EE and EO */
        EE[0] = E[0] + E[3];
        EO[0] = E[0] - E[3];
        EE[1] = E[1] + E[2];
        EO[1] = E[1] - E[2];
        
        dst[0*line] = (g_aiT8[0][0]*EE[0] + g_aiT8[0][1]*EE[1] + add) >> shift;
        dst[4*line] = (g_aiT8[4][0]*EE[0] + g_aiT8[4][1]*EE[1] + add) >> shift;
        dst[2*line] = (g_aiT8[2][0]*EO[0] + g_aiT8[2][1]*EO[1] + add) >> shift;
        dst[6*line] = (g_aiT8[6][0]*EO[0] + g_aiT8[6][1]*EO[1] + add) >> shift;
        
        dst[1*line] = (g_aiT8[1][0]*O[0] + g_aiT8[1][1]*O[1]
                     + g_aiT8[1][2]*O[2] + g_aiT8[1][3]*O[3] + add) >> shift;
        dst[3*line] = (g_aiT8[3][0]*O[0] + g_aiT8[3][1]*O[1]
                     + g_aiT8[3][2]*O[2] + g_aiT8[3][3]*O[3] + add) >> shift;
        dst[5*line] = (g_aiT8[5][0]*O[0] + g_aiT8[5][1]*O[1]
                     + g_aiT8[5][2]*O[2] + g_aiT8[5][3]*O[3] + add) >> shift;
        dst[7*line] = (g_aiT8[7][0]*O[0] + g_aiT8[7][1]*O[1]
                     + g_aiT8[7][2]*O[2] + g_aiT8[7][3]*O[3] + add) >> shift;
        
        src += 8;
        dst++;
    }
}
void partialButterfly8_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O*/
        for (k = 0; k < 4; k++) {
            E[k] = src[k] + src[7-k];
            O[k] = src[k] - src[7-k];
        }
        /* EE and EO */
        EE[0] = E[0] + E[3];
        EO[0] = E[0] - E[3];
        EE[1] = E[1] + E[2];
        EO[1] = E[1] - E[2];
        
        dst[0*line] = (g_aiT8[0][0]*EE[0] + g_aiT8[0][1]*EE[1] + add) >> shift;
        dst[4*line] = (g_aiT8[4][0]*EE[0] + g_aiT8[4][1]*EE[1] + add) >> shift;
        dst[2*line] = (g_aiT8[2][0]*EO[0] + g_aiT8[2][1]*EO[1] + add) >> shift;
        dst[6*line] = (g_aiT8[6][0]*EO[0] + g_aiT8[6][1]*EO[1] + add) >> shift;
        
        dst[1*line] = (g_aiT8[1][0]*O[0] + g_aiT8[1][1]*O[1]
                       + g_aiT8[1][2]*O[2] + g_aiT8[1][3]*O[3] + add) >> shift;
        dst[3*line] = (g_aiT8[3][0]*O[0] + g_aiT8[3][1]*O[1]
                       + g_aiT8[3][2]*O[2] + g_aiT8[3][3]*O[3] + add) >> shift;
        dst[5*line] = (g_aiT8[5][0]*O[0] + g_aiT8[5][1]*O[1]
                       + g_aiT8[5][2]*O[2] + g_aiT8[5][3]*O[3] + add) >> shift;
        dst[7*line] = (g_aiT8[7][0]*O[0] + g_aiT8[7][1]*O[1]
                       + g_aiT8[7][2]*O[2] + g_aiT8[7][3]*O[3] + add) >> shift;
        
        src += 8;
        dst++;
    }
}

void partialButterflyInverse8_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 4; k++) {
            O[k] = g_aiT8[ 1][k]*src[1*line] + g_aiT8[ 3][k]*src[3*line]
                 + g_aiT8[ 5][k]*src[5*line] + g_aiT8[ 7][k]*src[7*line];
        }
        
        EO[0] = g_aiT8[2][0]*src[2*line] + g_aiT8[6][0]*src[6*line];
        EO[1] = g_aiT8[2][1]*src[2*line] + g_aiT8[6][1]*src[6*line];
        EE[0] = g_aiT8[0][0]*src[0*line] + g_aiT8[4][0]*src[4*line];
        EE[1] = g_aiT8[0][1]*src[0*line] + g_aiT8[4][1]*src[4*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        E[0] = EE[0] + EO[0];
        E[3] = EE[0] - EO[0];
        E[1] = EE[1] + EO[1];
        E[2] = EE[1] - EO[1];
        for (k = 0; k < 4; k++) {
            dst[k+0] = Clip3(-32768, 32767, (E[k+0] + O[k+0] + add) >> shift);
            dst[k+4] = Clip3(-32768, 32767, (E[3-k] - O[3-k] + add) >> shift);
        }
        src++;
        dst += 8;
    }
}
void partialButterflyInverse8_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 4; k++) {
            O[k] = g_aiT8[ 1][k]*src[1*line] + g_aiT8[ 3][k]*src[3*line]
            + g_aiT8[ 5][k]*src[5*line] + g_aiT8[ 7][k]*src[7*line];
        }
        
        EO[0] = g_aiT8[2][0]*src[2*line] + g_aiT8[6][0]*src[6*line];
        EO[1] = g_aiT8[2][1]*src[2*line] + g_aiT8[6][1]*src[6*line];
        EE[0] = g_aiT8[0][0]*src[0*line] + g_aiT8[4][0]*src[4*line];
        EE[1] = g_aiT8[0][1]*src[0*line] + g_aiT8[4][1]*src[4*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        E[0] = EE[0] + EO[0];
        E[3] = EE[0] - EO[0];
        E[1] = EE[1] + EO[1];
        E[2] = EE[1] - EO[1];
        for (k = 0; k < 4; k++) {
            dst[k+0] = Clip3(-32768, 32767, (E[k+0] + O[k+0] + add) >> shift);
            dst[k+4] = Clip3(-32768, 32767, (E[3-k] - O[3-k] + add) >> shift);
        }
        src++;
        dst += 8;
    }
}

void partialButterfly16_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O*/
        for (k = 0; k < 8; k++) {
            E[k] = src[k] + src[15-k];
            O[k] = src[k] - src[15-k];
        }
        /* EE and EO */
        for (k = 0; k < 4; k++) {
            EE[k] = E[k] + E[7-k];
            EO[k] = E[k] - E[7-k];
        }
        /* EEE and EEO */
        EEE[0] = EE[0] + EE[3];
        EEO[0] = EE[0] - EE[3];
        EEE[1] = EE[1] + EE[2];
        EEO[1] = EE[1] - EE[2];
        
        dst[ 0*line] = (g_aiT16[ 0][0]*EEE[0] + g_aiT16[ 0][1]*EEE[1] + add) >> shift;
        dst[ 8*line] = (g_aiT16[ 8][0]*EEE[0] + g_aiT16[ 8][1]*EEE[1] + add) >> shift;
        dst[ 4*line] = (g_aiT16[ 4][0]*EEO[0] + g_aiT16[ 4][1]*EEO[1] + add) >> shift;
        dst[12*line] = (g_aiT16[12][0]*EEO[0] + g_aiT16[12][1]*EEO[1] + add) >> shift;
        
        for (k = 2; k < 16; k += 4) {
            dst[k*line] = (g_aiT16[k][0]*EO[0] + g_aiT16[k][1]*EO[1]
                         + g_aiT16[k][2]*EO[2] + g_aiT16[k][3]*EO[3] + add) >> shift;
        }
        
        for (k = 1; k < 16; k += 2) {
            dst[k*line] = (g_aiT16[k][0]*O[0] + g_aiT16[k][1]*O[1]
                         + g_aiT16[k][2]*O[2] + g_aiT16[k][3]*O[3]
                         + g_aiT16[k][4]*O[4] + g_aiT16[k][5]*O[5]
                         + g_aiT16[k][6]*O[6] + g_aiT16[k][7]*O[7] + add) >> shift;
        }
        
        src += 16;
        dst ++;
        
    }
}
void partialButterfly16_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O*/
        for (k = 0; k < 8; k++) {
            E[k] = src[k] + src[15-k];
            O[k] = src[k] - src[15-k];
        }
        /* EE and EO */
        for (k = 0; k < 4; k++) {
            EE[k] = E[k] + E[7-k];
            EO[k] = E[k] - E[7-k];
        }
        /* EEE and EEO */
        EEE[0] = EE[0] + EE[3];
        EEO[0] = EE[0] - EE[3];
        EEE[1] = EE[1] + EE[2];
        EEO[1] = EE[1] - EE[2];
        
        dst[ 0*line] = (g_aiT16[ 0][0]*EEE[0] + g_aiT16[ 0][1]*EEE[1] + add) >> shift;
        dst[ 8*line] = (g_aiT16[ 8][0]*EEE[0] + g_aiT16[ 8][1]*EEE[1] + add) >> shift;
        dst[ 4*line] = (g_aiT16[ 4][0]*EEO[0] + g_aiT16[ 4][1]*EEO[1] + add) >> shift;
        dst[12*line] = (g_aiT16[12][0]*EEO[0] + g_aiT16[12][1]*EEO[1] + add) >> shift;
        
        for (k = 2; k < 16; k += 4) {
            dst[k*line] = (g_aiT16[k][0]*EO[0] + g_aiT16[k][1]*EO[1]
                           + g_aiT16[k][2]*EO[2] + g_aiT16[k][3]*EO[3] + add) >> shift;
        }
        
        for (k = 1; k < 16; k += 2) {
            dst[k*line] = (g_aiT16[k][0]*O[0] + g_aiT16[k][1]*O[1]
                           + g_aiT16[k][2]*O[2] + g_aiT16[k][3]*O[3]
                           + g_aiT16[k][4]*O[4] + g_aiT16[k][5]*O[5]
                           + g_aiT16[k][6]*O[6] + g_aiT16[k][7]*O[7] + add) >> shift;
        }
        
        src += 16;
        dst ++;
        
    }
}

void partialButterflyInverse16_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 8; k++) {
            O[k] = g_aiT16[ 1][k]*src[ 1*line] + g_aiT16[ 3][k]*src[ 3*line]
                 + g_aiT16[ 5][k]*src[ 5*line] + g_aiT16[ 7][k]*src[ 7*line]
                 + g_aiT16[ 9][k]*src[ 9*line] + g_aiT16[11][k]*src[11*line]
                 + g_aiT16[13][k]*src[13*line] + g_aiT16[15][k]*src[15*line];
        }
        for (k = 0; k < 4; k++) {
            EO[k] = g_aiT16[ 2][k]*src[ 2*line] + g_aiT16[ 6][k]*src[ 6*line]
                  + g_aiT16[10][k]*src[10*line] + g_aiT16[14][k]*src[14*line];
        }
        EEO[0] = g_aiT16[4][0]*src[4*line] + g_aiT16[12][0]*src[12*line];
        EEE[0] = g_aiT16[0][0]*src[0*line] + g_aiT16[ 8][0]*src[ 8*line];
        EEO[1] = g_aiT16[4][1]*src[4*line] + g_aiT16[12][1]*src[12*line];
        EEE[1] = g_aiT16[0][1]*src[0*line] + g_aiT16[ 8][1]*src[ 8*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        for (k = 0; k < 2; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k+2] = EEE[1-k] - EEO[1-k];
        }
        for (k = 0; k < 4; k++) {
            E[k+0] = EE[k+0] + EO[k+0];
            E[k+4] = EE[3-k] - EO[3-k];
        }
        for (k = 0; k < 8; k++) {
            dst[k+0] = Clip3(-32768, 32767, (E[k+0] + O[k+0] + add) >> shift);
            dst[k+8] = Clip3(-32768, 32767, (E[7-k] - O[7-k] + add) >> shift);
        }
        src++;
        dst += 16;
    }
}
void partialButterflyInverse16_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 8; k++) {
            O[k] = g_aiT16[ 1][k]*src[ 1*line] + g_aiT16[ 3][k]*src[ 3*line]
            + g_aiT16[ 5][k]*src[ 5*line] + g_aiT16[ 7][k]*src[ 7*line]
            + g_aiT16[ 9][k]*src[ 9*line] + g_aiT16[11][k]*src[11*line]
            + g_aiT16[13][k]*src[13*line] + g_aiT16[15][k]*src[15*line];
        }
        for (k = 0; k < 4; k++) {
            EO[k] = g_aiT16[ 2][k]*src[ 2*line] + g_aiT16[ 6][k]*src[ 6*line]
            + g_aiT16[10][k]*src[10*line] + g_aiT16[14][k]*src[14*line];
        }
        EEO[0] = g_aiT16[4][0]*src[4*line] + g_aiT16[12][0]*src[12*line];
        EEE[0] = g_aiT16[0][0]*src[0*line] + g_aiT16[ 8][0]*src[ 8*line];
        EEO[1] = g_aiT16[4][1]*src[4*line] + g_aiT16[12][1]*src[12*line];
        EEE[1] = g_aiT16[0][1]*src[0*line] + g_aiT16[ 8][1]*src[ 8*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        for (k = 0; k < 2; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k+2] = EEE[1-k] - EEO[1-k];
        }
        for (k = 0; k < 4; k++) {
            E[k+0] = EE[k+0] + EO[k+0];
            E[k+4] = EE[3-k] - EO[3-k];
        }
        for (k = 0; k < 8; k++) {
            dst[k+0] = Clip3(-32768, 32767, (E[k+0] + O[k+0] + add) >> shift);
            dst[k+8] = Clip3(-32768, 32767, (E[7-k] - O[7-k] + add) >> shift);
        }
        src++;
        dst += 16;
    }
}

void partialButterfly32_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O*/
        for (k = 0; k < 16; k++) {
            E[k] = src[k] + src[31-k];
            O[k] = src[k] - src[31-k];
        }
        /* EE and EO */
        for (k = 0; k < 8; k++) {
            EE[k] = E[k] + E[15-k];
            EO[k] = E[k] - E[15-k];
        }
        /* EEE and EEO */
        for (k = 0; k < 4; k++) {
            EEE[k] = EE[k] + EE[7-k];
            EEO[k] = EE[k] - EE[7-k];
        }
        /* EEEE and EEEO */
        EEEE[0] = EEE[0] + EEE[3];
        EEEO[0] = EEE[0] - EEE[3];
        EEEE[1] = EEE[1] + EEE[2];
        EEEO[1] = EEE[1] - EEE[2];
        
        dst[ 0*line] = (g_aiT32[ 0][0]*EEEE[0] + g_aiT32[ 0][1]*EEEE[1] + add) >> shift;
        dst[16*line] = (g_aiT32[16][0]*EEEE[0] + g_aiT32[16][1]*EEEE[1] + add) >> shift;
        dst[ 8*line] = (g_aiT32[ 8][0]*EEEO[0] + g_aiT32[ 8][1]*EEEO[1] + add) >> shift;
        dst[24*line] = (g_aiT32[24][0]*EEEO[0] + g_aiT32[24][1]*EEEO[1] + add) >> shift;
        for (k = 4; k < 32; k += 8) {
            dst[k*line] = (g_aiT32[k][0]*EEO[0] + g_aiT32[k][1]*EEO[1]
                         + g_aiT32[k][2]*EEO[2] + g_aiT32[k][3]*EEO[3] + add) >> shift;
        }
        for (k = 2; k < 32; k += 4) {
            dst[k*line] = (g_aiT32[k][0]*EO[0] + g_aiT32[k][1]*EO[1]
                         + g_aiT32[k][2]*EO[2] + g_aiT32[k][3]*EO[3]
                         + g_aiT32[k][4]*EO[4] + g_aiT32[k][5]*EO[5]
                         + g_aiT32[k][6]*EO[6] + g_aiT32[k][7]*EO[7] + add) >> shift;
        }
        for (k = 1; k < 32; k += 2) {
            dst[k*line] = (g_aiT32[k][ 0]*O[ 0] + g_aiT32[k][ 1]*O[ 1]
                         + g_aiT32[k][ 2]*O[ 2] + g_aiT32[k][ 3]*O[ 3]
                         + g_aiT32[k][ 4]*O[ 4] + g_aiT32[k][ 5]*O[ 5]
                         + g_aiT32[k][ 6]*O[ 6] + g_aiT32[k][ 7]*O[ 7]
                         + g_aiT32[k][ 8]*O[ 8] + g_aiT32[k][ 9]*O[ 9]
                         + g_aiT32[k][10]*O[10] + g_aiT32[k][11]*O[11]
                         + g_aiT32[k][12]*O[12] + g_aiT32[k][13]*O[13]
                         + g_aiT32[k][14]*O[14] + g_aiT32[k][15]*O[15] + add)>>shift;
        }
        src += 32;
        dst++;
    }
}
void partialButterfly32_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* E and O*/
        for (k = 0; k < 16; k++) {
            E[k] = src[k] + src[31-k];
            O[k] = src[k] - src[31-k];
        }
        /* EE and EO */
        for (k = 0; k < 8; k++) {
            EE[k] = E[k] + E[15-k];
            EO[k] = E[k] - E[15-k];
        }
        /* EEE and EEO */
        for (k = 0; k < 4; k++) {
            EEE[k] = EE[k] + EE[7-k];
            EEO[k] = EE[k] - EE[7-k];
        }
        /* EEEE and EEEO */
        EEEE[0] = EEE[0] + EEE[3];
        EEEO[0] = EEE[0] - EEE[3];
        EEEE[1] = EEE[1] + EEE[2];
        EEEO[1] = EEE[1] - EEE[2];
        
        dst[ 0*line] = (g_aiT32[ 0][0]*EEEE[0] + g_aiT32[ 0][1]*EEEE[1] + add) >> shift;
        dst[16*line] = (g_aiT32[16][0]*EEEE[0] + g_aiT32[16][1]*EEEE[1] + add) >> shift;
        dst[ 8*line] = (g_aiT32[ 8][0]*EEEO[0] + g_aiT32[ 8][1]*EEEO[1] + add) >> shift;
        dst[24*line] = (g_aiT32[24][0]*EEEO[0] + g_aiT32[24][1]*EEEO[1] + add) >> shift;
        for (k = 4; k < 32; k += 8) {
            dst[k*line] = (g_aiT32[k][0]*EEO[0] + g_aiT32[k][1]*EEO[1]
                           + g_aiT32[k][2]*EEO[2] + g_aiT32[k][3]*EEO[3] + add) >> shift;
        }
        for (k = 2; k < 32; k += 4) {
            dst[k*line] = (g_aiT32[k][0]*EO[0] + g_aiT32[k][1]*EO[1]
                           + g_aiT32[k][2]*EO[2] + g_aiT32[k][3]*EO[3]
                           + g_aiT32[k][4]*EO[4] + g_aiT32[k][5]*EO[5]
                           + g_aiT32[k][6]*EO[6] + g_aiT32[k][7]*EO[7] + add) >> shift;
        }
        for (k = 1; k < 32; k += 2) {
            dst[k*line] = (g_aiT32[k][ 0]*O[ 0] + g_aiT32[k][ 1]*O[ 1]
                           + g_aiT32[k][ 2]*O[ 2] + g_aiT32[k][ 3]*O[ 3]
                           + g_aiT32[k][ 4]*O[ 4] + g_aiT32[k][ 5]*O[ 5]
                           + g_aiT32[k][ 6]*O[ 6] + g_aiT32[k][ 7]*O[ 7]
                           + g_aiT32[k][ 8]*O[ 8] + g_aiT32[k][ 9]*O[ 9]
                           + g_aiT32[k][10]*O[10] + g_aiT32[k][11]*O[11]
                           + g_aiT32[k][12]*O[12] + g_aiT32[k][13]*O[13]
                           + g_aiT32[k][14]*O[14] + g_aiT32[k][15]*O[15] + add)>>shift;
        }
        src += 32;
        dst++;
    }
}

void partialButterflyInverse32_gl(__global short *src, __local short *dst, int shift, int line)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 16; k++) {
            O[k] = g_aiT32[ 1][k]*src[ 1*line] + g_aiT32[ 3][k]*src[ 3*line]
                 + g_aiT32[ 5][k]*src[ 5*line] + g_aiT32[ 7][k]*src[ 7*line]
                 + g_aiT32[ 9][k]*src[ 9*line] + g_aiT32[11][k]*src[11*line]
                 + g_aiT32[13][k]*src[13*line] + g_aiT32[15][k]*src[15*line]
                 + g_aiT32[17][k]*src[17*line] + g_aiT32[19][k]*src[19*line]
                 + g_aiT32[21][k]*src[21*line] + g_aiT32[23][k]*src[23*line]
                 + g_aiT32[25][k]*src[25*line] + g_aiT32[27][k]*src[27*line]
                 + g_aiT32[29][k]*src[29*line] + g_aiT32[31][k]*src[31*line];
        }
        for (k = 0; k < 8; k++) {
            EO[k] = g_aiT32[ 2][k]*src[ 2*line] + g_aiT32[ 6][k]*src[ 6*line]
                  + g_aiT32[10][k]*src[10*line] + g_aiT32[14][k]*src[14*line]
                  + g_aiT32[18][k]*src[18*line] + g_aiT32[22][k]*src[22*line]
                  + g_aiT32[26][k]*src[26*line] + g_aiT32[30][k]*src[30*line];
        }
        for (k = 0; k < 4; k++) {
            EEO[k] = g_aiT32[ 4][k]*src[ 4*line] + g_aiT32[12][k]*src[12*line]
                   + g_aiT32[20][k]*src[20*line] + g_aiT32[28][k]*src[28*line];
        }
        EEEO[0] = g_aiT32[8][0]*src[8*line] + g_aiT32[24][0]*src[24*line];
        EEEO[1] = g_aiT32[8][1]*src[8*line] + g_aiT32[24][1]*src[24*line];
        EEEE[0] = g_aiT32[0][0]*src[0*line] + g_aiT32[16][0]*src[16*line];
        EEEE[1] = g_aiT32[0][1]*src[0*line] + g_aiT32[16][1]*src[16*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        EEE[0] = EEEE[0] + EEEO[0];
        EEE[3] = EEEE[0] - EEEO[0];
        EEE[1] = EEEE[1] + EEEO[1];
        EEE[2] = EEEE[1] - EEEO[1];    
        for (k = 0; k < 4; k++) {
            EE[k+0] = EEE[k+0] + EEO[k+0];
            EE[k+4] = EEE[3-k] - EEO[3-k];
        }    
        for (k = 0; k < 8; k++) {
            E[k+0] = EE[k+0] + EO[k+0];
            E[k+8] = EE[7-k] - EO[7-k];
        }    
        for (k = 0; k < 16; k++) {
            dst[k+ 0] = Clip3(-32768, 32767, (E[k +0] + O[k +0] + add) >> shift);
            dst[k+16] = Clip3(-32768, 32767, (E[15-k] - O[15-k] + add) >> shift);
        }
        src++;
        dst += 32;
    }
}
void partialButterflyInverse32_lg(__local short *src, __global short *dst, int shift, int line)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1<<(shift-1);
    
    for (j = 0; j < line; j++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 16; k++) {
            O[k] = g_aiT32[ 1][k]*src[ 1*line] + g_aiT32[ 3][k]*src[ 3*line]
            + g_aiT32[ 5][k]*src[ 5*line] + g_aiT32[ 7][k]*src[ 7*line]
            + g_aiT32[ 9][k]*src[ 9*line] + g_aiT32[11][k]*src[11*line]
            + g_aiT32[13][k]*src[13*line] + g_aiT32[15][k]*src[15*line]
            + g_aiT32[17][k]*src[17*line] + g_aiT32[19][k]*src[19*line]
            + g_aiT32[21][k]*src[21*line] + g_aiT32[23][k]*src[23*line]
            + g_aiT32[25][k]*src[25*line] + g_aiT32[27][k]*src[27*line]
            + g_aiT32[29][k]*src[29*line] + g_aiT32[31][k]*src[31*line];
        }
        for (k = 0; k < 8; k++) {
            EO[k] = g_aiT32[ 2][k]*src[ 2*line] + g_aiT32[ 6][k]*src[ 6*line]
            + g_aiT32[10][k]*src[10*line] + g_aiT32[14][k]*src[14*line]
            + g_aiT32[18][k]*src[18*line] + g_aiT32[22][k]*src[22*line]
            + g_aiT32[26][k]*src[26*line] + g_aiT32[30][k]*src[30*line];
        }
        for (k = 0; k < 4; k++) {
            EEO[k] = g_aiT32[ 4][k]*src[ 4*line] + g_aiT32[12][k]*src[12*line]
            + g_aiT32[20][k]*src[20*line] + g_aiT32[28][k]*src[28*line];
        }
        EEEO[0] = g_aiT32[8][0]*src[8*line] + g_aiT32[24][0]*src[24*line];
        EEEO[1] = g_aiT32[8][1]*src[8*line] + g_aiT32[24][1]*src[24*line];
        EEEE[0] = g_aiT32[0][0]*src[0*line] + g_aiT32[16][0]*src[16*line];
        EEEE[1] = g_aiT32[0][1]*src[0*line] + g_aiT32[16][1]*src[16*line];
        
        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
        EEE[0] = EEEE[0] + EEEO[0];
        EEE[3] = EEEE[0] - EEEO[0];
        EEE[1] = EEEE[1] + EEEO[1];
        EEE[2] = EEEE[1] - EEEO[1];
        for (k = 0; k < 4; k++) {
            EE[k+0] = EEE[k+0] + EEO[k+0];
            EE[k+4] = EEE[3-k] - EEO[3-k];
        }
        for (k = 0; k < 8; k++) {
            E[k+0] = EE[k+0] + EO[k+0];
            E[k+8] = EE[7-k] - EO[7-k];
        }
        for (k = 0; k < 16; k++) {
            dst[k+ 0] = Clip3(-32768, 32767, (E[k +0] + O[k +0] + add) >> shift);
            dst[k+16] = Clip3(-32768, 32767, (E[15-k] - O[15-k] + add) >> shift);
        }
        src++;
        dst += 32;
    }
}

/** MxN forward transform (2D)
 *  \param block input data (residual)
 *  \param coeff output data (transform coefficients)
 *  \param iWidth input data (width of transform)
 *  \param iHeight input data (height of transform)
 */
__kernel void xTrMxN(__global short *block, __global short *coeff,
                     const int iWidth, const int iHeight, const unsigned int uiMode)
{
    int shift_1st = g_aucConvertToBit[iWidth] + 1 + g_uiBitIncrement; // log2(iWidth) - 1 + g_uiBitIncrement
    int shift_2nd = g_aucConvertToBit[iHeight] + 8;                   // log2(iHeight) + 6
    
    __local short tmp[64 * 64];
    
    if (iWidth == 4 && iHeight == 4) {
        if (uiMode != REG_DCT) {
            fastForwardDst_gl(block, tmp, shift_1st); // Forward DST BY FAST ALGORITHM, block input, tmp output
            fastForwardDst_lg(tmp, coeff, shift_2nd); // Forward DST BY FAST ALGORITHM, tmp input, coeff output
        }
        else {
            partialButterfly4_gl(block, tmp, shift_1st, iHeight);
            partialButterfly4_lg(tmp, coeff, shift_2nd, iWidth);
        }
    }
    else if (iWidth == 8 && iHeight == 8) {
        partialButterfly8_gl(block, tmp, shift_1st, iHeight);
        partialButterfly8_lg(tmp, coeff, shift_2nd, iWidth);
    }
    else if (iWidth == 16 && iHeight == 16) {
        partialButterfly16_gl(block, tmp, shift_1st, iHeight);
        partialButterfly16_lg(tmp, coeff, shift_2nd, iWidth);
    }
    else if (iWidth == 32 && iHeight == 32) {
        partialButterfly32_gl(block, tmp, shift_1st, iHeight);
        partialButterfly32_lg(tmp, coeff, shift_2nd, iWidth);
    }
}

/** MxN inverse transform (2D)
 *  \param coeff input data (transform coefficients)
 *  \param block output data (residual)
 *  \param iWidth input data (width of transform)
 *  \param iHeight input data (height of transform)
 */
__kernel void xITrMxN(__global short *coeff, __global short *block,
                      const int iWidth, const int iHeight, const unsigned int uiMode)
{
    int shift_1st = SHIFT_INV_1ST;
    int shift_2nd = SHIFT_INV_2ND - g_uiBitIncrement;
    
    __local short tmp[64 * 64];

    if (iWidth == 4 && iHeight == 4) {
        if (uiMode != REG_DCT) {
            fastInverseDst_gl(coeff, tmp, shift_1st); // Inverse DST by FAST Algorithm, coeff input, tmp output
            fastInverseDst_lg(tmp, block, shift_2nd); // Inverse DST by FAST Algorithm, tmp input, coeff output
        }
        else {
            partialButterflyInverse4_gl(coeff, tmp, shift_1st, iWidth);
            partialButterflyInverse4_lg(tmp, block, shift_2nd, iHeight);
        }
    }
    else if (iWidth == 8 && iHeight == 8) {
        partialButterflyInverse8_gl(coeff, tmp, shift_1st, iWidth);
        partialButterflyInverse8_lg(tmp, block, shift_2nd, iHeight);
    }
    else if (iWidth == 16 && iHeight == 16) {
        partialButterflyInverse16_gl(coeff,tmp,shift_1st,iWidth);
        partialButterflyInverse16_lg(tmp,block,shift_2nd,iHeight);
    }
    else if (iWidth == 32 && iHeight == 32) {
        partialButterflyInverse32_gl(coeff, tmp, shift_1st, iWidth);
        partialButterflyInverse32_lg(tmp, block, shift_2nd, iHeight);
    }
}

