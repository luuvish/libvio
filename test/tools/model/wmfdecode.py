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

 File      : wmfdecode.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 May 20, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('WmfDecode', )

__version__ = '2.0.0'

from . import ModelExecutor


class WmfDecode(ModelExecutor):

    model   = 'wmf-decode'
    codecs  = ('vp1', )
    actions = ('decode', 'digest_by_frames', 'compare')

    def __init__(self, codec='vc1', **kwargs):
        from os.path import join, normpath, dirname

        super(WmfDecoder, self).__init__(codec, **kwargs)

        self._execute = normpath(join(dirname(__file__), 'wmfdecode_vc.exe'))

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-v', 'IYUV_WMV', '-cut1stframe', '1', '-ignoreaudio', '1', '-i', source, '-o', target]
