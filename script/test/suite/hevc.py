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

 File      : hevc.py
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
from ..model.hm_14_0 import HM


bbc = (('bbc/' + profile + '/' + stream[0], stream[1])
    for profile in (
        'i_main',
        'i_main10',
        'ld_main',
        'ld_main10',
        'lp_main',
        'lp_main10',
        'ra_main',
        'ra_main10'
    )
    for stream in (
        ('BasketballDrill_832x480_50_*.bin',                   500),
        ('BasketballDrillText_832x480_50_*.bin',               500),
        ('BasketballDrive_1920x1080_50_*.bin',                 500),
        ('BasketballPass_416x240_50_*.bin',                    500),
        ('BlowingBubbles_416x240_50_*.bin',                    500),
        ('BQMall_832x480_60_*.bin',                            600),
        ('BQSquare_416x240_60_*.bin',                          600),
        ('BQTerrace_1920x1080_60_*.bin',                       600),
        ('Cactus_1920x1080_50_*.bin',                          500),
        ('ChinaSpeed_1024x768_30_*.bin',                       500), # not 300
        ('FourPeople_1280x720_60_*.bin',                       600),
        ('Johnny_1280x720_60_*.bin',                           600),
        ('Kimono1_1920x1080_24_*.bin',                         240),
        ('KristenAndSara_1280x720_60_*.bin',                   600),
        ('NebutaFestival_2560x1600_60_10bit_crop_*.bin',       300), # not 600
        ('ParkScene_1920x1080_24_*.bin',                       240),
        ('PartyScene_832x480_50_*.bin',                        500),
        ('PeopleOnStreet_2560x1600_30_crop_*.bin',             150), # not 300
        ('RaceHorses_416x240_30_*.bin',                        300),
        ('RaceHorses_832x480_30_*.bin',                        300),
        ('SlideEditing_1280x720_30_*.bin',                     300),
        ('SlideShow_1280x720_20_*.bin',                        500), # not 200
        ('SteamLocomotiveTrain_2560x1600_60_10bit_crop_*.bin', 300), # not 600
        ('Traffic_2560x1600_30_crop_*.bin',                    150)  # not 300
    )
)


models = (FFmpeg, HM)

suites = (
    {
        'suite' : 'encode-hevc-hm-14.0',
        'model' : 'hm-14.0',
        'codec' : 'hevc',
        'action': 'encode',
        'stdout': 'hevc-hm-14.0.log',
        'srcdir': join(rootpath, 'test/image/yuv'),
        'outdir': join(rootpath, 'test/stream/hevc'),
        'includes': (('phantom.yuv', {'width':1920, 'height':1080, 'frames':600}), ),
        'excludes': ()
    },
    {
        'suite' : 'decode-hevc-hm-14.0',
        'model' : 'hm-14.0',
        'codec' : 'hevc',
        'action': 'decode',
        'stdout': 'hevc-hm-14.0.log',
        'srcdir': join(rootpath, 'test/stream/hevc'),
        'outdir': join(rootpath, 'test/image/hevc'),
        'includes': bbc,
        'excludes': ()
    },
    {
        'suite' : 'digest-hevc-hm-14.0',
        'model' : 'hm-14.0',
        'codec' : 'hevc',
        'action': 'digest_by_frames',
        'stdout': 'hevc-hm-14.0.log',
        'srcdir': join(rootpath, 'test/stream/hevc'),
        'outdir': join(rootpath, 'test/digest/hevc'),
        'includes': bbc,
        'excludes': ()
    },
    {
        'suite' : 'compare-hevc-hm-14.0',
        'model' : 'hm-14.0',
        'codec' : 'hevc',
        'action': 'compare',
        'stdout': 'hevc-hm-14.0.log',
        'srcdir': join(rootpath, 'test/stream/hevc'),
        'outdir': join(rootpath, 'test/digest/hevc'),
        'includes': bbc,
        'excludes': ()
    },
    {
        'suite' : 'compare-hevc-ffmpeg',
        'model' : 'ffmpeg',
        'codec' : 'hevc',
        'action': 'compare',
        'stdout': 'hevc-ffmpeg.log',
        'srcdir': join(rootpath, 'test/stream/hevc'),
        'outdir': join(rootpath, 'test/digest/hevc'),
        'includes': bbc,
        'excludes': ()
    }
)
