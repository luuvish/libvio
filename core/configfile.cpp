
/*!
 ***********************************************************************
 * \file
 *    configfile.c
 * \brief
 *    Configuration handling.
 * \author
 *  Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger           <stewe@cs.tu-berlin.de>
 * \note
 *    In the future this module should hide the Parameters and offer only
 *    Functions for their access.  Modules which make frequent use of some parameters
 *    (e.g. picture size in macroblocks) are free to buffer them on local variables.
 *    This will not only avoid global variable and make the code more readable, but also
 *    speed it up.  It will also greatly facilitate future enhancements such as the
 *    handling of different picture sizes in the same sequence.                         \n
 *                                                                                      \n
 *    For now, everything is just copied to the inp_par structure (gulp)
 *
 **************************************************************************************
 * \par Configuration File Format
 **************************************************************************************
 * Format is line oriented, maximum of one parameter per line                           \n
 *                                                                                      \n
 * Lines have the following format:                                                     \n
 * \<ParameterName\> = \<ParameterValue\> # Comments \\n                                    \n
 * Whitespace is space and \\t
 * \par
 * \<ParameterName\> are the predefined names for Parameters and are case sensitive.
 *   See configfile.h for the definition of those names and their mapping to
 *   cfgparams->values.
 * \par
 * \<ParameterValue\> are either integers [0..9]* or strings.
 *   Integers must fit into the wordlengths, signed values are generally assumed.
 *   Strings containing no whitespace characters can be used directly.  Strings containing
 *   whitespace characters are to be inclosed in double quotes ("string with whitespace")
 *   The double quote character is forbidden (may want to implement something smarter here).
 * \par
 * Any Parameters whose ParameterName is undefined lead to the termination of the program
 * with an error message.
 *
 * \par Known bug/Shortcoming:
 *    zero-length strings (i.e. to signal an non-existing file
 *    have to be coded as "".
 *
 * \par Rules for using command files
 *                                                                                      \n
 * All Parameters are initially taken from DEFAULTCONFIGFILENAME, defined in configfile.h.
 * If an -f \<config\> parameter is present in the command line then this file is used to
 * update the defaults of DEFAULTCONFIGFILENAME.  There can be more than one -f parameters
 * present.  If -p <ParameterName = ParameterValue> parameters are present then these
 * override the default and the additional config file's settings, and are themselves
 * overridden by future -p parameters.  There must be whitespace between -f and -p commands
 * and their respective parameters
 ***********************************************************************
 */

#define INCLUDED_BY_CONFIGFILE_C

#include <sys/stat.h>

#include "global.h"
#include "memalloc.h"
#include "bitstream_cabac.h"
#include "configfile.h"


#define MAX_ITEMS_TO_PARSE  10000
#define DEFAULTCONFIGFILENAME "decoder.cfg"

InputParameters cfgparams;

//! Maps parameter name to its address, type etc.
typedef struct {
  const char *TokenName;    //!< name
  void *Place;        //!< address
  int Type;           //!< type:  0-int, 1-char[], 2-double
  double Default;     //!< default value
  int param_limits;   //!< 0: no limits, 1: both min and max, 2: only min (i.e. no negatives), 3: special case for QPs since min needs bitdepth_qp_scale
  double min_limit;
  double max_limit;
  int    char_size;   //!< Dimension of type char[]
} Mapping;

// Mapping_Map Syntax:
// {NAMEinConfigFile,  &cfgparams.VariableName, Type, InitialValue, LimitType, MinLimit, MaxLimit, CharSize}
// Types : {0:int, 1:text, 2: double}
// LimitType: {0:none, 1:both, 2:minimum, 3: QP based}
// We could separate this based on types to make it more flexible and allow also defaults for text types.
Mapping Map[] = {
    {"InputFile",                &cfgparams.infile,                       1,   0.0,                       0,  0.0,              0.0,             FILE_NAME_SIZE, },
    {"OutputFile",               &cfgparams.outfile,                      1,   0.0,                       0,  0.0,              0.0,             FILE_NAME_SIZE, },
    {"RefFile",                  &cfgparams.reffile,                      1,   0.0,                       0,  0.0,              0.0,             FILE_NAME_SIZE, },
    {"WriteUV",                  &cfgparams.write_uv,                     0,   1.0,                       1,  0.0,              1.0,                             },
    {"FileFormat",               &cfgparams.FileFormat,                   0,   0.0,                       1,  0.0,              1.0,                             },
    {"RefOffset",                &cfgparams.ref_offset,                   0,   0.0,                       1,  0.0,              256.0,                             },
    {"POCScale",                 &cfgparams.poc_scale,                    0,   2.0,                       1,  1.0,              10.0,                            },
    {"DisplayDecParams",         &cfgparams.bDisplayDecParams,            0,   1.0,                       1,  0.0,              1.0,                             },
    {"ConcealMode",              &cfgparams.conceal_mode,                 0,   0.0,                       1,  0.0,              2.0,                             },
    {"RefPOCGap",                &cfgparams.ref_poc_gap,                  0,   2.0,                       1,  0.0,              4.0,                             },
    {"POCGap",                   &cfgparams.poc_gap,                      0,   2.0,                       1,  0.0,              4.0,                             },
    {"Silent",                   &cfgparams.silent,                       0,   0.0,                       1,  0.0,              1.0,                             },
    {"DecFrmNum",                &cfgparams.iDecFrmNum,                   0,   0.0,                       2,  0.0,              0.0,                             },
#if (MVC_EXTENSION_ENABLE)
    {"DecodeAllLayers",          &cfgparams.DecodeAllLayers,              0,   0.0,                       1,  0.0,              1.0,                             },
#endif
    {"DPBPLUS0",                 &cfgparams.dpb_plus[0],                  0,   1.0,                       1,  -16.0,            16.0,                             },
    {"DPBPLUS1",                 &cfgparams.dpb_plus[1],                  0,   0.0,                       1,  -16.0,            16.0,                             },
    {NULL,                       NULL,                                   -1,   0.0,                       0,  0.0,              0.0,                             },
};

/*!
 ***********************************************************************
 * \brief
 *    allocates memory buf, opens file Filename in f, reads contents into
 *    buf and returns buf
 * \param Filename
 *    name of config file
 * \return
 *    if successfull, content of config file
 *    NULL in case of error. Error message will be set in errortext
 ***********************************************************************
 */
static char *GetConfigFileContent (const char *Filename)
{
  long FileSize;
  FILE *f;
  char *buf;

  if (NULL == (f = fopen (Filename, "r")))
  {
      snprintf (errortext, ET_SIZE, "Cannot open configuration file %s.", Filename);
      return NULL;
  }

  if (0 != fseek (f, 0, SEEK_END))
  {
    snprintf (errortext, ET_SIZE, "Cannot fseek in configuration file %s.", Filename);
    return NULL;
  }

  FileSize = ftell (f);

  if (FileSize < 0 || FileSize > 150000)
  {
    snprintf (errortext, ET_SIZE, "\nUnreasonable Filesize %ld reported by ftell for configuration file %s.", FileSize, Filename);
    return NULL;
  }
  if (0 != fseek (f, 0, SEEK_SET))
  {
    snprintf (errortext, ET_SIZE, "Cannot fseek in configuration file %s.", Filename);
    return NULL;
  }

  if ((buf = (char *)malloc (FileSize + 1))==NULL) no_mem_exit("GetConfigFileContent: buf");

  // Note that ftell() gives us the file size as the file system sees it.  The actual file size,
  // as reported by fread() below will be often smaller due to CR/LF to CR conversion and/or
  // control characters after the dos EOF marker in the file.

  FileSize = (long) fread (buf, 1, FileSize, f);
  buf[FileSize] = '\0';


  fclose (f);
  return buf;
}


/*!
 ***********************************************************************
 * \brief
 *    Returns the index number from Map[] for a given parameter name.
 * \param Map
 *    Mapping structure
 * \param s
 *    parameter name string
 * \return
 *    the index number if the string is a valid parameter name,         \n
 *    -1 for error
 ***********************************************************************
 */
static int ParameterNameToMapIndex (Mapping *Map, char *s)
{
  int i = 0;

  while (Map[i].TokenName != NULL)
    if (0==strcasecmp (Map[i].TokenName, s))
      return i;
    else
      i++;
  return -1;
}

/*!
 ***********************************************************************
 * \brief
 *    Parses the character array buf and writes global variable input, which is defined in
 *    configfile.h.  This hack will continue to be necessary to facilitate the addition of
 *    new parameters through the Map[] mechanism (Need compiler-generated addresses in map[]).
 * \param p_Inp
 *    InputParameters of configuration
 * \param Map
 *    Mapping structure to specify the name and value mapping relation
 * \param buf
 *    buffer to be parsed
 * \param bufsize
 *    buffer size of buffer
 ***********************************************************************
 */
static void ParseContent (InputParameters *p_Inp, Mapping *Map, char *buf, int bufsize)
{
  char *items[MAX_ITEMS_TO_PARSE] = {NULL};
  int MapIdx;
  int item = 0;
  int InString = 0, InItem = 0;
  char *p = buf;
  char *bufend = &buf[bufsize];
  int IntContent;
  double DoubleContent;
  int i;

  // Stage one: Generate an argc/argv-type list in items[], without comments and whitespace.
  // This is context insensitive and could be done most easily with lex(1).

  while (p < bufend)
  {
    switch (*p)
    {
    case 13:
      ++p;
      break;
    case '#':                 // Found comment
      *p = '\0';              // Replace '#' with '\0' in case of comment immediately following integer or string
      while (*p != '\n' && p < bufend)  // Skip till EOL or EOF, whichever comes first
        ++p;
      InString = 0;
      InItem = 0;
      break;
    case '\n':
      InItem = 0;
      InString = 0;
      *p++='\0';
      break;
    case ' ':
    case '\t':              // Skip whitespace, leave state unchanged
      if (InString)
        p++;
      else
      {                     // Terminate non-strings once whitespace is found
        *p++ = '\0';
        InItem = 0;
      }
      break;

    case '"':               // Begin/End of String
      *p++ = '\0';
      if (!InString)
      {
        items[item++] = p;
        InItem = ~InItem;
      }
      else
        InItem = 0;
      InString = ~InString; // Toggle
      break;

    default:
      if (!InItem)
      {
        items[item++] = p;
        InItem = ~InItem;
      }
      p++;
    }
  }

  item--;

  for (i=0; i<item; i+= 3)
  {
    if (0 > (MapIdx = ParameterNameToMapIndex (Map, items[i])))
    {
      //snprintf (errortext, ET_SIZE, " Parsing error in config file: Parameter Name '%s' not recognized.", items[i]);
      //error (errortext, 300);
      printf ("\n\tParsing error in config file: Parameter Name '%s' not recognized.", items[i]);
      i -= 2 ;
      continue;
    }
    if (strcasecmp ("=", items[i+1]))
    {
      snprintf (errortext, ET_SIZE, " Parsing error in config file: '=' expected as the second token in each line.");
      error (errortext, 300);
    }

    // Now interpret the Value, context sensitive...

    switch (Map[MapIdx].Type)
    {
    case 0:           // Numerical
      if (1 != sscanf (items[i+2], "%d", &IntContent))
      {
        snprintf (errortext, ET_SIZE, " Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
        error (errortext, 300);
      }
      * (int *) (Map[MapIdx].Place) = IntContent;
      printf (".");
      break;
    case 1:
      if (items[i + 2] == NULL)
        memset((char *) Map[MapIdx].Place, 0, Map[MapIdx].char_size);
      else
        strncpy ((char *) Map[MapIdx].Place, items [i+2], Map[MapIdx].char_size);
      printf (".");
      break;
    case 2:           // Numerical double
      if (1 != sscanf (items[i+2], "%lf", &DoubleContent))
      {
        snprintf (errortext, ET_SIZE, " Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
        error (errortext, 300);
      }
      * (double *) (Map[MapIdx].Place) = DoubleContent;
      printf (".");
      break;
    default:
      error ((char *)"Unknown value type in the map definition of configfile.h",-1);
    }
  }
  *p_Inp = cfgparams;
}

/*!
 ***********************************************************************
 * \brief
 *    Sets initial values for encoding parameters.
 * \return
 *    -1 for error
 ***********************************************************************
 */
static int InitParams(Mapping *Map)
{
  int i = 0;

  while (Map[i].TokenName != NULL)
  {
    if (Map[i].Type == 0)
        * (int *) (Map[i].Place) = (int) Map[i].Default;
    else if (Map[i].Type == 2)
    * (double *) (Map[i].Place) = Map[i].Default;
      i++;
  }
  return -1;
}

/*!
 ***********************************************************************
 * \brief
 *    Validates encoding parameters.
 * \return
 *    -1 for error
 ***********************************************************************
 */
static int TestParams(Mapping *Map, int bitdepth_qp_scale[3])
{
  int i = 0;

  while (Map[i].TokenName != NULL)
  {
    if (Map[i].param_limits == 1)
    {
      if (Map[i].Type == 0)
      {
        if ( * (int *) (Map[i].Place) < (int) Map[i].min_limit || * (int *) (Map[i].Place) > (int) Map[i].max_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, (int) Map[i].min_limit,(int)Map[i].max_limit );
          error (errortext, 400);
        }

      }
      else if (Map[i].Type == 2)
      {
        if ( * (double *) (Map[i].Place) < Map[i].min_limit || * (double *) (Map[i].Place) > Map[i].max_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%.2f, %.2f] range.", Map[i].TokenName,Map[i].min_limit ,Map[i].max_limit );
          error (errortext, 400);
        }
      }
    }
    else if (Map[i].param_limits == 2)
    {
      if (Map[i].Type == 0)
      {
        if ( * (int *) (Map[i].Place) < (int) Map[i].min_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should not be smaller than %d.", Map[i].TokenName, (int) Map[i].min_limit);
          error (errortext, 400);
        }
      }
      else if (Map[i].Type == 2)
      {
        if ( * (double *) (Map[i].Place) < Map[i].min_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should not be smaller than %2.f.", Map[i].TokenName,Map[i].min_limit);
          error (errortext, 400);
        }
      }
    }
    else if (Map[i].param_limits == 3) // Only used for QPs
    {
      
      if (Map[i].Type == 0)
      {
        int cur_qp = * (int *) (Map[i].Place);
        int min_qp = (int) (Map[i].min_limit - (bitdepth_qp_scale? bitdepth_qp_scale[0]: 0));
        int max_qp = (int) Map[i].max_limit;
        
        if (( cur_qp < min_qp ) || ( cur_qp > max_qp ))
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, min_qp, max_qp );
          error (errortext, 400);
        }
      }
    }

    i++;
  }
  return -1;
}



/*!
 ***********************************************************************
 * \brief
 *    Outputs encoding parameters.
 * \return
 *    -1 for error
 ***********************************************************************
 */
static int DisplayParams(Mapping *Map, char *message)
{
  int i = 0;

  printf("******************************************************\n");
  printf("*               %s                   *\n", message);
  printf("******************************************************\n");
  while (Map[i].TokenName != NULL)
  {
    if (Map[i].Type == 0)
      printf("Parameter %s = %d\n",Map[i].TokenName,* (int *) (Map[i].Place));
    else if (Map[i].Type == 1)
      printf("Parameter %s = ""%s""\n",Map[i].TokenName,(char *)  (Map[i].Place));
    else if (Map[i].Type == 2)
      printf("Parameter %s = %.2f\n",Map[i].TokenName,* (double *) (Map[i].Place));
      i++;
  }
  printf("******************************************************\n");
  return i;
}

/*!
 ***********************************************************************
 * \brief
 *   print help message and exit
 ***********************************************************************
 */
static void JMDecHelpExit (void)
{
  fprintf( stderr, "\n   ldecod [-h] [-d defdec.cfg] {[-f curenc1.cfg]...[-f curencN.cfg]}"
    " {[-p EncParam1=EncValue1]..[-p EncParamM=EncValueM]}\n\n"
    "## Parameters\n\n"

    "## Options\n"
    "   -h :  prints function usage\n"
    "   -d :  use <defdec.cfg> as default file for parameter initializations.\n"
    "         If not used then file defaults to encoder.cfg in local directory.\n"
    "   -f :  read <curencM.cfg> for reseting selected encoder parameters.\n"
    "         Multiple files could be used that set different parameters\n"
    "   -p :  Set parameter <DecParamM> to <DecValueM>.\n"
    "         See default decoder.cfg file for description of all parameters.\n\n"

    "## Examples of usage:\n"
    "   ldecod\n"
    "   ldecod  -h\n"
    "   ldecod  -d default.cfg\n"
    "   ldecod  -f curenc1.cfg\n"
    "   ldecod  -f curenc1.cfg -p InputFile=\"e:\\data\\container_qcif_30.264\" -p OutputFile=\"dec.yuv\" -p RefFile=\"Rec.yuv\"\n");

  exit(-1);
}

/*!
 ***********************************************************************
 * \brief
 *    Checks the input parameters for consistency.
 ***********************************************************************
 */
static void PatchInp (InputParameters *p_Inp)
{
  //int i;
  //int storedBplus1;
  TestParams(Map, NULL);
}

/*!
 ***********************************************************************
 * \brief
 *    Parse the command line parameters and read the config files.
 * \param p_Vid
 *    VideoParameters structure for encoding
 * \param p_Inp
 *    InputParameters structure as input configuration
 * \param ac
 *    number of command line parameters
 * \param av
 *    command line parameters
 ***********************************************************************
 */
void ParseCommand(InputParameters *p_Inp, int ac, char *av[])
{
  char *content = NULL;
  int CLcount;
  const char *filename=DEFAULTCONFIGFILENAME;

  if (ac==2)
  {
    if (0 == strncmp (av[1], "-v", 2))
    {
      printf("JM " JM ": compiled " __DATE__ " " __TIME__ "\n");
      exit(-1);
    }

    if (0 == strncmp (av[1], "-h", 2))
    {
      JMDecHelpExit();
    }
  }

  memcpy (&cfgparams, p_Inp, sizeof (InputParameters));
  //Set default parameters.
  printf ("Setting Default Parameters...\n");
  InitParams(Map);

  *p_Inp = cfgparams;
  // Process default config file
  CLcount = 1;

  if (ac>=3)
  {
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMDecHelpExit();
    }
  }
  if(filename)
  {
    printf ("Parsing Configfile %s\n", filename);
    content = GetConfigFileContent (filename);
    if (NULL != content)
    {
      ParseContent (p_Inp, Map, content, (int) strlen(content));
      printf ("\n");
      free (content);
    }
  }
  // Parse the command line

  while (CLcount < ac)
  {
    if (0 == strncmp (av[CLcount], "-h", 2))
    {
      JMDecHelpExit();
    }

    if (0 == strncmp (av[CLcount], "-i", 2) || 0 == strncmp (av[CLcount], "-I", 2))  // A file parameter?
    {
      strncpy(p_Inp->infile, av[CLcount+1], FILE_NAME_SIZE);
      CLcount += 2;
    } 
    else if (0 == strncmp (av[CLcount], "-o", 2) || 0 == strncmp (av[CLcount], "-O", 2))  // A file parameter?
    {
      strncpy(p_Inp->outfile, av[CLcount+1], FILE_NAME_SIZE);
      CLcount += 2;
    } 
    else
    {
      snprintf (errortext, ET_SIZE, "Error in command line, ac %d, around string '%s', missing -f or -p parameters?", CLcount, av[CLcount]);
      error (errortext, 300);
    }
  }
  printf ("\n");

  PatchInp(p_Inp);
  cfgparams = *p_Inp;
  if (p_Inp->bDisplayDecParams)
    DisplayParams(Map, (char *)"Decoder Parameters");
}
