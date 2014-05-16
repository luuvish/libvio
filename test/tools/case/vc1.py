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

 File      : vc1.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 13, 2014    Executor classify

================================================================================
'''

__all__ = ('models', 'cases')

__version__ = '2.0.0'

from os.path import normpath, join, dirname

from ..model.coda960 import Coda960
from ..model.ffmpeg import FFmpeg
from ..model.smpte_vc1 import SMPTEVc1
from ..model.wmfdecode import WmfDecode


root = normpath(join(dirname(__file__), '../..'))

models = (Coda960, FFmpeg, SMPTEVc1, WmfDecode)

cases = (
    {
        'case'  : 'wmfdecode-vc1-decode',
        'model' : 'wmfdecode',
        'codec' : 'vc1',
        'action': 'decode',
        'stdout': 'vc1-wmfdecode.log',
        'srcdir': join(root, 'streams/vc1'),
        'outdir': join(root, 'images/vc1'),
        'includes': (
            'artificial6/*_Main_Progressive_*',
            'artificial7/*',
            'smpte-vc1/*',
            'conf-vc1/*',
            'conf-wmv9/*',
            'stress/*'
        ),
        'excludes': ()
    },
    {
        'case'  : 'smpte-vc1-decode',
        'model' : 'smpte-vc1',
        'codec' : 'vc1',
        'action': 'decode',
        'stdout': 'vc1-smpte.log',
        'srcdir': join(root, 'streams/vc1'),
        'outdir': join(root, 'images/vc1'),
        'includes': (
            'artificial6/*_Main_Progressive_*',
            'artificial7/*',
            'smpte-vc1/*',
            'conf-vc1/*',
            'conf-wmv9/*',
            'stress/*'
        ),
        'excludes': ()
    },
    {
        'case'  : 'coda960-vc1-compare',
        'model' : 'coda960',
        'codec' : 'vc1',
        'action': 'compare',
        'stdout': 'vc1-coda960.log',
        'srcdir': join(root, 'streams/vc1'),
        'outdir': join(root, 'digests/vc1'),
        'includes': (
            'artificial6/*_Main_Progressive_*',
            'artificial7/*',
            'smpte-vc1/*',
            'conf-vc1/*',
            'conf-wmv9/*',
            'stress/*'
        ),
        'excludes': ()
    }
)
