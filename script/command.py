#!/usr/bin/env python
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

 File      : command.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 12, 2014    Executor classify

================================================================================
'''

from sys import path, argv, stderr
from os.path import dirname

path.append(dirname(__file__))

from test.globber import Globber
from test.model import ModelExecutor
from test.model.allegro_h264 import AllegroH264
from test.model.coda960 import Coda960
from test.model.ffmpeg import FFmpeg
from test.model.hm_14_0 import HM
from test.model.jm_18_6 import JM
from test.model.libvio import LibVio
from test.model.libvpx import LibVpx
from test.model.smpte_vc1 import SMPTEVc1
from test.model.wmfdecode import WmfDecode


models = (AllegroH264, Coda960, FFmpeg, HM, JM, LibVio, LibVpx, SMPTEVc1, WmfDecode)


class Command(object):

    program     = 'command.py'
    description = 'shell command helper for model executors'
    epilog      = 'Copyright (C) 2014 Thumb \'o Cat luuvish (luuvsih@gmail.com)'

    def __init__(self, *largs):

        self.executors = {}
        self.models    = []
        self.codecs    = []
        self.actions   = []
        self.encode_by_option = False
        self.digest_by_frames = False

        for executor in largs:
            if not (issubclass(executor, ModelExecutor) or isinstance(executor, ModelExecutor)):
                raise Exception('require class inherited from ModelExecutor')
            if executor.model in self.models:
                raise Exception('models must be identical each others')

            self.executors[executor.model] = executor
            self.models.append(executor.model)
            for codec in executor.codecs:
                if codec not in self.codecs:
                    self.codecs.append(codec)
            for action in executor.actions:
                if action == 'encode':
                    self.encode_by_option = True
                if action == 'digest_by_frames':
                    self.digest_by_frames = True
                    action = 'digest'
                if action not in self.actions:
                    self.actions.append(action)

        self.models.sort()
        self.codecs.sort()
        self.actions.sort()

    def parse(self, *largs):
        from argparse import ArgumentParser

        parser = ArgumentParser(prog=self.program,
                description=self.description, epilog=self.epilog)

        if len(self.models) > 1:
            parser.add_argument('-m', '--model', help='codec model',
                    dest='model', action='store', choices=self.models)
        if len(self.codecs) > 1:
            parser.add_argument('-c', '--codec', help='video codec',
                    dest='codec', action='store', choices=self.codecs)
        if len(self.actions) > 1:
            parser.add_argument('-a', '--action', help='action',
                    dest='action', action='store', choices=self.actions)

        if self.encode_by_option:
            parser.add_argument('-x', '--width', help='encoded picture width',
                    dest='width', action='store', type=int, default=0)
            parser.add_argument('-y', '--height', help='encoded picture height',
                    dest='height', action='store', type=int, default=0)
        if self.encode_by_option or self.digest_by_frames:
            parser.add_argument('-f', '--frames', help='frames',
                    dest='frames', action='store', type=int, default=0)

        parser.add_argument('files', help='input filenames',
                action='store', metavar='FILENAMES', nargs='*')

        parser.parse_args(args=largs, namespace=self)

        return vars(self)

    def main(self, *largs):
        from glob import glob

        options = self.parse(*(largs[1:] if len(largs) > 1 else ['-h']))

        model  = options.get('model', self.models[0] if len(self.models) == 1 else None)
        codec  = options.get('codec', self.codecs[0] if len(self.codecs) == 1 else None)
        action = options.get('action', self.actions[0] if len(self.actions) == 1 else None)
        width  = options.get('width', 0)
        height = options.get('height', 0)
        frames = options.get('frames', 0)
        files  = options.get('files', [])

        executor = self.executors[model]
        if issubclass(executor, ModelExecutor):
            executor = executor(codec)

        globber = Globber(stdout=None).bind(executor)

        if action == 'encode' and (width <= 0 or height <= 0 or frames <= 0):
            raise Exception('encode option must be set width, height, frames option')
        if action == 'digest' and frames > 0:
            action = 'digest_by_frames'

        method = getattr(globber, action)
        for file in files:
            if action == 'encode':
                file = (file, {'width':width, 'height':height, 'frames':frames})
            if action == 'digest_by_frames':
                file = (file, frames)
            method(file)


if __name__ == '__main__':
    try:
        Command(*models).main(*argv)
    except Exception as e:
        stderr.write('%s\n' % e.message)
