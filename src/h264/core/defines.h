#ifndef _DEFINES_H_
#define _DEFINES_H_


#define ENABLE_OUTPUT_TONEMAPPING 1    //!< enable tone map the output if tone mapping SEI present
#define DISABLE_ERC               0    //!< Disable any error concealment processes

#define MVC_EXTENSION_ENABLE      1    //!< enable support for the Multiview High Profile

#define MVC_INIT_VIEW_ID          -1
#define MAX_VIEW_NUM              1024   

typedef uint8_t  byte;
typedef uint16_t imgpel;

#define MAX_NUM_SLICES         50
#define MAX_REFERENCE_PICTURES 32               //!< H.264 allows 32 fields
#define MAX_CODED_FRAME_SIZE   8000000         //!< bytes for one frame
#define MAX_NUM_DECSLICES      16
#define MCBUF_LUMA_PAD_X        32
#define MCBUF_LUMA_PAD_Y        12
#define MCBUF_CHROMA_PAD_X      16
#define MCBUF_CHROMA_PAD_Y      8
#define MAX_NUM_DPB_LAYERS      2

//AVC Profile IDC definitions
typedef enum {
  FREXT_CAVLC444 = 44,       //!< YUV 4:4:4/14 "CAVLC 4:4:4"
  BASELINE       = 66,       //!< YUV 4:2:0/8  "Baseline"
  MAIN           = 77,       //!< YUV 4:2:0/8  "Main"
  EXTENDED       = 88,       //!< YUV 4:2:0/8  "Extended"
  FREXT_HP       = 100,      //!< YUV 4:2:0/8  "High"
  FREXT_Hi10P    = 110,      //!< YUV 4:2:0/10 "High 10"
  FREXT_Hi422    = 122,      //!< YUV 4:2:2/10 "High 4:2:2"
  FREXT_Hi444    = 244,      //!< YUV 4:4:4/14 "High 4:4:4"
  MVC_HIGH       = 118,      //!< YUV 4:2:0/8  "Multiview High"
  STEREO_HIGH    = 128       //!< YUV 4:2:0/8  "Stereo High"
} ProfileIDC;

#define BLOCK_SHIFT            2
#define BLOCK_SIZE             4
#define MB_BLOCK_SIZE         16


enum {
  EOS = 1,    //!< End Of Sequence
  SOP = 2,    //!< Start Of Picture
  SOS = 3,     //!< Start Of slice_t
  SOS_CONT = 4
};


#define MAX_PLANE 3


#endif
