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

 File      : ffmpeg.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 Aug 28, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('FFmpeg', )

__version__ = '2.0.0'

from . import ModelExecutor


class FFmpeg(ModelExecutor):

    model   = 'ffmpeg'
    codecs  = ('h264', 'hevc', 'vc1', 'vp8', 'vp9')
    actions = ('decode', 'digest_by_frames', 'compare')

    def __init__(self, codec, **kwargs):
        from os.path import join, normpath, dirname

        super(FFmpeg, self).__init__(codec, **kwargs)

        root   = normpath(join(dirname(__file__), '../../..'))
        binary = 'ref/dist/ffmpeg-2.2.1.bin/dist/bin/ffmpeg'

        self._execute = join(root, binary)

        self.defaults['decode'] += ['-codec', codec]

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-i', source, target]
