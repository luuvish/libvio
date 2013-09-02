#ifndef _DEFINES_H_
#define _DEFINES_H_


#define JM                  "18 (FRExt)"
#define VERSION             "18.5"
#define EXT_VERSION         "(FRExt)"

#define ENABLE_OUTPUT_TONEMAPPING 1    //!< enable tone map the output if tone mapping SEI present
#define DISABLE_ERC               0    //!< Disable any error concealment processes
#define SIMULCAST_ENABLE          0    //!< to test the decoder

#define MVC_EXTENSION_ENABLE      1    //!< enable support for the Multiview High Profile

#define MVC_INIT_VIEW_ID          -1
#define MAX_VIEW_NUM              1024   

typedef unsigned char  byte;     //!< byte type definition
typedef unsigned char  uint8;    //!< type definition for unsigned char (same as byte, 8 bits)
typedef unsigned short uint16;   //!< type definition for unsigned short (16 bits)
typedef unsigned int   uint32;   //!< type definition for unsigned int (32 bits)

typedef          char  int8;
typedef          short int16;
typedef          int   int32;

typedef uint16 imgpel;

enum {
  FALSE,
  TRUE
};
#define Boolean int

#define MAX_NUM_SLICES     50
#define MAX_REFERENCE_PICTURES 32               //!< H.264 allows 32 fields
#define MAX_CODED_FRAME_SIZE 8000000         //!< bytes for one frame
#define MAX_NUM_DECSLICES  16
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

#define FILE_NAME_SIZE  255

#define BLOCK_SHIFT            2
#define BLOCK_SIZE             4
#define BLOCK_SIZE_8x8         8
#define SMB_BLOCK_SIZE         8
#define BLOCK_PIXELS          16
#define MB_BLOCK_SIZE         16
#define BLOCK_MULTIPLE         4 // (MB_BLOCK_SIZE/BLOCK_SIZE)

//  Available MB modes
typedef enum {
  PSKIP        =  0,
  BSKIP_DIRECT =  0,
  P16x16       =  1,
  P16x8        =  2,
  P8x16        =  3,
  SMB8x8       =  4,
  SMB8x4       =  5,
  SMB4x8       =  6,
  SMB4x4       =  7,
  P8x8         =  8,
  I4MB         =  9,
  I16MB        = 10,
  IBLOCK       = 11,
  SI4MB        = 12,
  I8MB         = 13,
  IPCM         = 14,
  MAXMODE      = 15
} MBModeTypes;

enum {
  EOS = 1,    //!< End Of Sequence
  SOP = 2,    //!< Start Of Picture
  SOS = 3,     //!< Start Of slice_t
  SOS_CONT = 4
};


#define MAX_PLANE 3


#endif
