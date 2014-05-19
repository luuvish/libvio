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

 File      : vp8.py
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
from ..model.coda960 import Coda960
from ..model.ffmpeg import FFmpeg
from ..model.libvpx import LibVpx


models = (Coda960, FFmpeg, LibVpx)

suites = (
    {
        'suite' : 'decode-vp8-libvpx',
        'model' : 'libvpx',
        'codec' : 'vp8',
        'action': 'decode',
        'stdout': 'vp8-libvpx.log',
        'srcdir': join(rootpath, 'test/stream/vp8'),
        'outdir': join(rootpath, 'test/image/vp8'),
        'includes': ('*.ivf', ),
        'excludes': ()
    },
    {
        'suite' : 'digest-vp8-libvpx',
        'model' : 'libvpx',
        'codec' : 'vp8',
        'action': 'digest',
        'stdout': 'vp8-libvpx.log',
        'srcdir': join(rootpath, 'test/stream/vp8'),
        'outdir': join(rootpath, 'test/digest/vp8'),
        'includes': ('*.ivf', ),
        'excludes': ()
    },
    {
        'suite' : 'compare-vp8-libvpx',
        'model' : 'libvpx',
        'codec' : 'vp8',
        'action': 'compare',
        'stdout': 'vp8-libvpx.log',
        'srcdir': join(rootpath, 'test/stream/vp8'),
        'outdir': join(rootpath, 'test/digest/vp8'),
        'includes': ('*.ivf', ),
        'excludes': ()
    },
    {
        'suite' : 'compare-vp8-coda960',
        'model' : 'coda960',
        'codec' : 'vp8',
        'action': 'compare',
        'stdout': 'vp8-coda960.log',
        'srcdir': join(rootpath, 'test/stream/vp8'),
        'outdir': join(rootpath, 'test/digest/vp8'),
        'includes': ('*.ivf', ),
        'excludes': ()
    }
)
