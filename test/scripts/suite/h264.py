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

 File      : h264.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 13, 2014    Executor classify

================================================================================
'''

__all__ = ('models', 'suites')

__version__ = '2.0.0'

from os.path import normpath, join, dirname

from ..model.allegro_h264 import AllegroH264
from ..model.coda960 import Coda960
from ..model.ffmpeg import FFmpeg
from ..model.jm_18_6 import JM
from ..model.violib import VioLib


root = normpath(join(dirname(__file__), '../..'))

models = (AllegroH264, Coda960, FFmpeg, JM, VioLib)

suites = (
    {
        'suite' : 'encode-h264-jm-18.6',
        'model' : 'jm-18.6',
        'codec' : 'h264',
        'action': 'encode',
        'stdout': 'h264-jm-18.6.log',
        'srcdir': join(root, 'images/yuv'),
        'outdir': join(root, 'streams/h264'),
        'includes': (('phantom.yuv', {'width':1920, 'height':1080, 'frames':600}), ),
        'excludes': ()
    },
    {
        'suite' : 'decode-h264-allegro',
        'model' : 'allegro-h264',
        'codec' : 'h264',
        'action': 'decode',
        'stdout': 'h264-allegro.log',
        'srcdir': join(root, 'streams/h264'),
        'outdir': join(root, 'images/h264'),
        'includes': (
            'jvt/bp/*',
            'jvt/mp/*',
            'jvt/hp/*',
            'allegro/bp/*',
            'allegro/mp/*',
            'allegro/hp/*'
        ),
        'excludes': (
            'jvt/bp/MR2_TANDBERG_B.264',
        )
    },
    {
        'suite' : 'compare-h264-coda960',
        'model' : 'coda960',
        'codec' : 'h264',
        'action': 'compare',
        'stdout': 'h264-coda960.log',
        'srcdir': join(root, 'streams/h264'),
        'outdir': join(root, 'digests/h264'),
        'includes': (
            'jvt/bp/*',
            'jvt/mp/*',
            'jvt/hp/*',
            'allegro/bp/*',
            'allegro/mp/*',
            'allegro/hp/*'
        ),
        'excludes': (
            'jvt/bp/MR2_TANDBERG_B.264',
            'allegro/mp/Allegro_Inter_CABAC_10_L42_HD2K@60Hz_10.0.26l',
            'allegro/mp/Allegro_Inter_CAVLC_10_L42_HD2K@60Hz_10.0.26l'
        )
    },
    {
        'suite' : 'compare-h264-jm-18.6',
        'model' : 'jm-18.6',
        'codec' : 'h264',
        'action': 'compare',
        'stdout': 'h264-jm-18.6.log',
        'srcdir': join(root, 'streams/h264'),
        'outdir': join(root, 'digests/h264'),
        'includes': (
            'jvt/bp/*',
            'jvt/mp/*',
            'jvt/hp/*'
        ),
        'excludes': (
            'jvt/bp/FM1_FT_B.264',
            'jvt/bp/MR2_TANDBERG_B.264',
            'jvt/bp/sp1_bt_a.h264'
        )
    },
    {
        'suite' : 'compare-h264-violib',
        'model' : 'violib',
        'codec' : 'h264',
        'action': 'compare',
        'stdout': 'h264-violib.log',
        'srcdir': join(root, 'streams/h264'),
        'outdir': join(root, 'digests/h264'),
        'includes': (
            'jvt/bp/*',
            'jvt/mp/*',
            'jvt/hp/*'
        ),
        'excludes': (
            'jvt/bp/FM1_FT_B.264',
            'jvt/bp/MR2_TANDBERG_B.264',
            'jvt/bp/sp1_bt_a.h264'
        )
    }
)
