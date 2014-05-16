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

 File      : model.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('ModelExecutor', )

__version__ = '2.0.0'

from sys import stdout, stderr
from os.path import exists, splitext, basename


class Executor(object):

    codecs  = ('h264', 'hevc', 'vc1', 'vp8', 'vp9')
    actions = ('decode', 'encode', 'digest', 'digest_by_frames', 'compare')

    def __init__(self, **kwargs):

        self.stdout = kwargs.get('stdout', stdout)
        self.stderr = kwargs.get('stderr', stderr)

        self.defaults = {'decode':[], 'encode':[], 'digest':[]}

    def execute(self):
        return ''

    def options(self, source, target):
        return []

    def mkdir(self, target):
        from os import makedirs
        from os.path import normpath, dirname

        outdir = normpath(dirname(target))
        if not exists(outdir):
            makedirs(outdir)

    def decode(self, source, target):
        from subprocess import call

        execute = self.execute()
        options = self.defaults['decode'] + self.options(source, target)

        self.mkdir(target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
        except Exception as e:
            raise e
        if not exists(target):
            raise Exception('decode error: %s' % basename(source))

    def encode(self, source, target):
        from subprocess import call

        execute = self.execute()
        options = self.defaults['encode'] + self.options(source, target)

        self.mkdir(target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
        except Exception as e:
            raise e
        if not exists(target):
            raise Exception('encode error: %s' % basename(source))

    def digest(self, source, target=None):
        from subprocess import call
        from os import remove

        outname, outext = splitext(basename(source))
        output = '.' + outname + '.yuv.md5' if target is None else target

        execute = self.execute()
        options = self.defaults['digest'] + self.options(source, output)

        if target is not None:
            self.mkdir(target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
        except Exception as e:
            raise e
        if not exists(output):
            raise Exception('digest error: %s' % basename(source))

        lines = []
        with open(output, 'rt') as f:
            lines = [line.rstrip().lower() for line in f]
        if target is None:
            remove(output)

        return lines

    def digest_by_frames(self, source, target=None, frames=0):
        from hashlib import md5
        from os import remove
        from os.path import getsize

        if frames <= 0:
            raise Exception('digest error: %s' % basename(source))

        outname, outext = splitext(basename(source))
        output = '.' + outname + '.yuv'

        try:
            self.decode(source, output)
        except Exception as e:
            raise e
        if not exists(output):
            raise Exception('digest error: %s' % basename(source))

        lines = []
        with open(output, 'rb') as f:
            size = getsize(output) / frames
            while True:
                data = f.read(size)
                if len(data) <= 0:
                    break
                md5hash = md5()
                md5hash.update(data)
                lines.append(md5hash.hexdigest().lower())
        remove(output)

        if target is not None:
            with open(target, 'wt') as f:
                for line in lines:
                    f.write('%s\n' % line)

        return lines

    def compare(self, source, target):

        if not exists(target):
            raise Exception('digest no exists: %s' % basename(target))
        hashs = []
        with open(target, 'rt') as f:
            hashs = [line.rstrip().lower() for line in f]
        nhash = len(hashs)

        lines = []
        try:
            if 'digest' in self.actions:
                lines = self.digest(source, None)
            else:
                lines = self.digest_by_frames(source, None, nline)
        except Exception as e:
            raise e
        nline = len(lines)

        if nline != nhash:
            raise Exception('decoded frames is different: %s' % basename(source))
        for i in xrange(nline):
            line = lines[i]
            hash = hashs[i]
            if line != hash:
                raise Exception('mismatch %d %s: %s != %s' % (i, basename(source), line, hash))

        return lines


class ModelExecutor(Executor):

    model = ''

    def __init__(self, codec, **kwargs):
        super(ModelExecutor, self).__init__(**kwargs)

        if codec not in self.codecs:
            raise Exception('codec must be one of %s' % self.codecs)

        for action in Executor.actions:
            if action not in self.actions:
                setattr(self, action, self.noslot)

    def noslot(self, *largs, **kwargs):
        raise Exception('action must be one of %s' % self.actions)

    def bind(self, globber):
        from .. import Globber

        if not isinstance(globber, Globber):
            raise Exception('bind require Globber class')
        return globber.bind(self)
