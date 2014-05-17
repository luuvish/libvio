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

 File      : hm.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 Jun 30, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('HM', )

__version__ = '2.0.0'

from . import ModelExecutor


class HM(ModelExecutor):

    model   = 'hm-14.0'
    codecs  = ('hevc', )
    actions = ('decode', 'encode', 'digest_by_frames', 'compare')

    ext     = '265'

    def __init__(self, codec='hevc', **kwargs):
        from os.path import join, normpath, dirname

        super(HM, self).__init__(codec, **kwargs)

        root   = normpath(join(dirname(__file__), '../../..'))
        decode = 'ref/dist/hm-14.0/bin/TAppDecoderStatic'
        encode = 'ref/dist/hm-14.0/bin/TAppEncoderStatic'
        config = 'ref/dist/hm-14.0/cfg/encoder_randomaccess_main.cfg'

        self._decode = join(root, decode)
        self._encode = join(root, encode)
        self._config = join(root, config)

    def execute(self):
        return self._decode

    def options(self, source, target):
        return ['--OutputBitDepth=8', '--BitstreamFile=' + source, '--ReconFile=' + target]

    def encode(self, source, target, option):
        from subprocess import call
        from os.path import exists, basename

        execute = self._encode
        options = [
            '-c', self._config, '--SEIpictureDigest=1',
            '--IntraPeriod=32', '--GOPSize=8', '--FrameRate=30',
            '--FramesToBeEncoded=%d' % option['frames'],
            '--SourceWidth=%d' % option['width'], '--SourceHeight=%d' % option['height'],
            '--InputFile=' + source, '--BitstreamFile=' + target
        ]

        self.mkdir(target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
        except Exception as e:
            raise e
        if not exists(target):
            raise Exception('encode error: %s' % basename(source))
