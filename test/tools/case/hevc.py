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

__all__ = ('models', 'cases')

__version__ = '2.0.0'

from os.path import normpath, join, dirname

from ..model.ffmpeg import FFmpeg
from ..model.hm_14_0 import HM


root = normpath(join(dirname(__file__), '../..'))

models = (FFmpeg, HM)
cases  = ()
'''
srcdir = join(root, 'streams/hevc')
outdir = join(root, 'reports/hevc')


bbc = {
    'profiles': ('lp_main', 'lp_main10', 'ra_main', 'ra_main10'),
    'streams': (
        ('BasketballDrill_*.yuv',      500),
        ('BasketballDrillText_*.yuv',  500),
        ('BasketballDrive_*.yuv',      500),
        ('BasketballPass_*.yuv',       500),
        ('BlowingBubbles_*.yuv',       500),
        ('BQMall_*.yuv',               600),
        ('BQSquare_*.yuv',             600),
        ('BQTerrace_*.yuv',            600),
        ('Cactus_*.yuv',               500),
        ('ChinaSpeed_*.yuv',           500), # not 300
        ('FourPeople_*.yuv',           600),
        ('Johnny_*.yuv',               600),
        ('Kimono1_*.yuv',              240),
        ('KristenAndSara_*.yuv',       600),
        ('NebutaFestival_*.yuv',       300), # not 600
        ('ParkScene_*.yuv',            240),
        ('PartyScene_*.yuv',           500),
        ('PeopleOnStreet_*.yuv',       150), # not 300
        ('RaceHorses_*.yuv',           300),
        ('SlideEditing_*.yuv',         300),
        ('SlideShow_*.yuv',            500), # not 200
        ('SteamLocomotiveTrain_*.yuv', 300), # not 600
        ('Traffic_*.yuv',              150)  # not 300
    )
}


def ffmpeg():

    model = FFmpeg('hevc', **{
        'stdout': open('hevc-ffmpeg.stdout', 'wt'),
        'stderr': open('hevc-ffmpeg.stderr', 'wt')
    })

    globber = Globber(
        srcdir=srcdir,
        outdir=outdir,
        includes=('bbc/' + p + '/' + s[0] for p in bbc['profiles'] for s in bbc['streams'])
    )

    model.bind(globber).digest()


def hm_14_0():

    model = HM(**{
        'stdout': open('hevc-hm-14.0.stdout', 'wt'),
        'stderr': open('hevc-hm-14.0.stderr', 'wt')
    })

    globber = Globber(
        srcdir=srcdir,
        outdir=outdir,
        includes=('bbc/' + p + '/' + s[0] for p in bbc['profiles'] for s in bbc['streams'])
    )

    model.bind(globber).digest()


if __name__ == '__main__':
    from argparse import ArgumentParser
    from sys import argv, stderr

    parser = ArgumentParser(prog='hevc reports generator',
            description='reports.hevc.py',
            epilog='Copyright (C) 2014 Thumb \'o Cat luuvish (luuvsih@gmail.com)')

    parser.add_argument('-m', '--mode', help='generator mode',
            dest='mode', action='store', choices=['ffmpeg', 'hm_14_0'])

    options = vars(parser.parse_args(args=argv[1:]))

    actions = {
        'ffmpeg' : ffmpeg,
        'hm_14_0': hm_14_0
    }

    actions[options.get('mode')]()
'''
