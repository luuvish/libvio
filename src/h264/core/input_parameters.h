#ifndef _CONFIGFILE_H_
#define _CONFIGFILE_H_


#define JM          "18 (FRExt)"
#define VERSION     "18.5"
#define EXT_VERSION "(FRExt)"

#define FILE_NAME_SIZE  255

// input parameters from configuration file
struct InputParameters {
    char        infile [FILE_NAME_SIZE]; //!< H.264 inputfile
    char        outfile[FILE_NAME_SIZE]; //!< Decoded YUV 4:2:0 output
    char        reffile[FILE_NAME_SIZE]; //!< Optional YUV 4:2:0 reference file for SNR measurement

    int         FileFormat;       //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
    int         ref_offset;
    int         poc_scale;
    int         write_uv;
    int         silent;

    // Input/output sequence format related variables

    int         DecodeAllLayers;

    // picture error concealment
    int         conceal_mode;
    int         ref_poc_gap;
    int         poc_gap;
  
    int         iDecFrmNum;

    int         bDisplayDecParams;
    int         dpb_plus[2];

    void        ParseCommand(int ac, char* av[]);
};


#endif
