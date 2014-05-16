# -*- coding: utf-8 -*-

'''
================================================================================

  This confidential and proprietary software may be used only
 as authorized by a licensing agreement from Thumb o'Cat Inc.
 In the event of publication, the following notice is applicable:

      Copyright (C) 2013 - 2014 Thumb o'Cat
                    All right reserved.

  The entire notice above must be reproduced on all authorized copies.

================================================================================

 File      : jm.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 June 15, 2013    first release
     2.0 May  12, 2014    Executor classify

================================================================================
'''

__all__ = ('JM', )

__version__ = '2.0.0'

from . import ModelExecutor


class JM(ModelExecutor):

    model   = 'jm-18.6'
    codecs  = ('h264', )
    actions = ('decode', 'digest_by_frames', 'compare')

    def __init__(self, codec='h264', **kwargs):
        from os.path import join, normpath, dirname

        super(JM, self).__init__(codec, **kwargs)

        root   = normpath(join(dirname(__file__), '../../..'))
        binary = 'ref/dist/jm-18.6/bin/ldecod.exe'

        self._execute = join(root, binary)

        self._config = '''\
# This is a file containing input parameters to the JVT H.264/AVC decoder.
# The text line following each parameter is discarded by the decoder.
#
# For bug reporting and known issues see:
# https://ipbt.hhi.fraunhofer.de
#
# New Input File Format is as follows
# <ParameterName> = <ParameterValue> # Comment
#
##########################################################################################
# Files
##########################################################################################
InputFile              = "%s"             # H.264/AVC coded bitstream
OutputFile             = "%s"             # Output file, YUV/RGB
RefFile                = "%s"             # Ref sequence (for SNR)
WriteUV                = 1                # Write 4:2:0 chroma components for monochrome streams
FileFormat             = 0                # NAL mode (0=Annex B, 1: RTP packets)
RefOffset              = 0                # SNR computation offset
POCScale               = 2                # Poc Scale (1 or 2)
##########################################################################################
# HRD parameters
##########################################################################################
#R_decoder             = 500000           # Rate_Decoder
#B_decoder             = 104000           # B_decoder
#F_decoder             = 73000            # F_decoder
#LeakyBucketParamFile  = "leakybucketparam.cfg" # LeakyBucket Params
##########################################################################################
# decoder control parameters
##########################################################################################
DisplayDecParams       = 0                # 1: Display parameters;
ConcealMode            = 0                # Err Concealment(0:Off,1:Frame Copy,2:Motion Copy)
RefPOCGap              = 2                # Reference POC gap (2: IPP (Default), 4: IbP / IpP)
POCGap                 = 2                # POC gap (2: IPP /IbP/IpP (Default), 4: IPP with frame skip = 1 etc.)
Silent                 = 0                # Silent decode
IntraProfileDeblocking = 1                # Enable Deblocking filter in intra only profiles (0=disable, 1=filter according to SPS parameters)
DecFrmNum              = 0                # Number of frames to be decoded (-n)
##########################################################################################
# MVC decoding parameters
##########################################################################################
DecodeAllLayers        = 0                 # Decode all views (-mpr)
'''

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-f', source]

    def decode(self, source, target):
        from subprocess import call
        from os import remove
        from os.path import splitext, basename

        srcname, srcext = splitext(basename(source))
        outname, outext = splitext(basename(target))
        optname = srcname + '.opt'
        with open(optname, 'wt') as f:
            f.write(self._config % (source, target, outname+'.rec'))

        execute = self.execute()
        options = self.options(optname, target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
            remove(optname)
            remove('dataDec.txt')
            remove('log.dec')
        except:
            raise Exception('decode error: %s' % basename(source))
