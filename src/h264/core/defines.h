#ifndef _DEFINES_H_
#define _DEFINES_H_


#define DISABLE_ERC               0    //!< Disable any error concealment processes

#define MVC_EXTENSION_ENABLE      1    //!< enable support for the Multiview High Profile
#define MVC_INIT_VIEW_ID          -1
#define MAX_VIEW_NUM              1024   

typedef uint8_t  byte;
typedef uint16_t imgpel;

#define MAX_CODED_FRAME_SIZE   8000000         //!< bytes for one frame
#define MCBUF_LUMA_PAD_X        32
#define MCBUF_LUMA_PAD_Y        12
#define MCBUF_CHROMA_PAD_X      16
#define MCBUF_CHROMA_PAD_Y      8


#endif
