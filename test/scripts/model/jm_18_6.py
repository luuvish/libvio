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
    actions = ('decode', 'encode', 'digest_by_frames', 'compare')

    ext     = '264'

    def __init__(self, codec='h264', **kwargs):
        from os.path import join, normpath, dirname

        super(JM, self).__init__(codec, **kwargs)

        root   = normpath(join(dirname(__file__), '../../..'))
        decode = 'ref/dist/jm-18.6/bin/ldecod.exe'
        encode = 'ref/dist/jm-18.6/bin/lencod.exe'
        cfgdec = 'ref/dist/jm-18.6/bin/decoder.cfg'
        cfgenc = 'ref/dist/jm-18.6/bin/HM-like/encoder_JM_RA_B_HE.cfg'

        self._decode = join(root, decode)
        self._encode = join(root, encode)
        self._cfgdec = join(root, cfgdec)
        self._cfgenc = join(root, cfgenc)

    def execute(self):
        return self._decode

    def options(self, source, target):
        return ['-d', self._cfgdec, '-p', 'InputFile=' + source, '-p' ,'OutputFile=' + target]

    def decode(self, source, target):
        from os import remove

        super(JM, self).decode(source, target)

        remove('dataDec.txt')
        remove('log.dec')

    def encode(self, source, target, option):
        from subprocess import call
        from os.path import exists, basename

        execute = self._encode
        options = [
            '-d', self._cfgenc,
            '-p', 'IntraPeriod=32', '-p', 'NumberBFrames=7',
            '-p', 'FramesToBeEncoded=%d' % option['frames'],
            '-p', 'SourceWidth=%d' % option['width'], '-p', 'SourceHeight=%d' % option['height'],
            '-p', 'InputFile=' + source, '-p', 'OutputFile=' + target
        ]

        self.mkdir(target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
        except Exception as e:
            raise e
        if not exists(target):
            raise Exception('encode error: %s' % basename(source))
