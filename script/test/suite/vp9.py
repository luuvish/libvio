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

__all__ = ('models', 'suites')

__version__ = '2.0.0'

from os.path import join

from . import rootpath
from ..model.ffmpeg import FFmpeg
from ..model.libvpx import LibVpx


models = (FFmpeg, LibVpx)

suites = (
    {
        'suite' : 'decode-vp9-libvpx',
        'model' : 'libvpx',
        'codec' : 'vp9',
        'action': 'decode',
        'stdout': 'vp9-libvpx.log',
        'srcdir': join(rootpath, 'test/stream/vp9'),
        'outdir': join(rootpath, 'test/image/vp9'),
        'includes': ('*.ivf', '*.webm'),
        'excludes': ('vp91-2-04-yv444.webm', )
    },
    {
        'suite' : 'digest-vp9-libvpx',
        'model' : 'libvpx',
        'codec' : 'vp9',
        'action': 'digest',
        'stdout': 'vp9-libvpx.log',
        'srcdir': join(rootpath, 'test/stream/vp9'),
        'outdir': join(rootpath, 'test/digest/vp9'),
        'includes': ('*.ivf', '*.webm'),
        'excludes': ('vp91-2-04-yv444.webm', )
    },
    {
        'suite' : 'compare-vp9-libvpx',
        'model' : 'libvpx',
        'codec' : 'vp9',
        'action': 'compare',
        'stdout': 'vp9-libvpx.log',
        'srcdir': join(rootpath, 'test/stream/vp9'),
        'outdir': join(rootpath, 'test/digest/vp9'),
        'includes': ('*.ivf', '*.webm'),
        'excludes': ('vp91-2-04-yv444.webm', )
    }
)
