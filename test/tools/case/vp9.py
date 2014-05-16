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

 File      : vp9.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 13, 2014    Executor classify

================================================================================
'''

__all__ = ('models', 'cases')

__version__ = '2.0.0'

from os.path import normpath, join, dirname

from ..model.ffmpeg import FFmpeg
from ..model.libvpx import LibVpx


root = normpath(join(dirname(__file__), '../..'))

models = (FFmpeg, LibVpx)

cases = (
    {
        'case'  : 'libvpx-vp9-digest',
        'model' : 'libvpx',
        'codec' : 'vp9',
        'action': 'digest',
        'stdout': 'vp9-libvpx.log',
        'srcdir': join(root, 'streams/vp9'),
        'outdir': join(root, 'digests/vp9'),
        'includes': ('*.ivf', ),
        'excludes': ()
    },
)
