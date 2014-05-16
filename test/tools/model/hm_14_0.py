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
    actions = ('decode', 'digest_by_frames', 'compare')

    def __init__(self, codec='hevc', **kwargs):
        from os.path import join, normpath, dirname

        super(HM, self).__init__(codec, **kwargs)

        root   = normpath(join(dirname(__file__), '../../..'))
        binary = 'ref/dist/hm-14.0/bin/TAppDecoderStatic'

        self._execute = join(root, binary)

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-b', source, '-o', target]
