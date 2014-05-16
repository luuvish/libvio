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

 File      : tools.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 May 19, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('Globber', )

__version__ = '2.0.0'

from functools import partial
from sys import stdout, stderr
from os.path import join, normpath, splitext, basename

from .model import ModelExecutor


class Globber(object):

    def __init__(self, executor=None, **kwargs):

        self.stdout = kwargs.get('stdout', stdout)
        self.stderr = kwargs.get('stderr', stderr)

        self.srcdir = kwargs.get('srcdir', '.')
        self.outdir = kwargs.get('outdir', '.')

        self.includes = kwargs.get('includes', [])
        self.excludes = kwargs.get('excludes', [])

        self.bind(executor)

    def glob(self, source):
        from glob import glob
        from os import getcwd, chdir

        curdir = getcwd()
        chdir(self.srcdir)
        items = glob(source)
        items.sort()
        chdir(curdir)

        for item in items:
            if item not in self.excludes:
                source = normpath(join(self.srcdir, item))
                target = normpath(join(self.outdir, item))
                target = basename(target) if self.outdir == '.' else target
                yield source, target

    def write(self, *largs, **kwargs):
        if type(self.stdout) is str:
            self.stdout = open(self.stdout, 'wt')
        if type(self.stdout) is file:
            self.stdout.write(*largs, **kwargs)

    def error(self, *largs, **kwargs):
        if type(self.stderr) is str:
            self.stderr = open(self.stderr, 'wt')
        if type(self.stderr) is file:
            self.stderr.write(*largs, **kwargs)

    def bind(self, executor):

        def decode(executor, source, target):
            outname, outext = splitext(target)
            target = outname + '.yuv'
            try:
                self.error('decode -i %s -o %s\n' % (basename(source), basename(target)))
                executor.decode(source, target)
            except Exception as e:
                self.error('# %s\n' % e.message)

        def encode(executor, source, target):
            outname, outext = splitext(target)
            target = outname + '.' + executor.ext
            try:
                self.error('encode -i %s -o %s\n' % (basename(source), basename(target)))
                executor.encode(source, target)
            except Exception as e:
                self.error('# %s\n' % e.message)

        def digest(executor, source, target):
            outname, outext = splitext(target)
            target = outname + '.yuv.md5'
            try:
                self.error('digest -i %s -o %s\n' % (basename(source), basename(target)))
                lines = executor.digest(source, target)
                self.write('%8d frames digest: %s\n' % (len(lines), basename(source)))
            except Exception as e:
                self.error('# %s\n' % e.message)

        def digest_by_frames(executor, source, frames):
            if source in self.excludes:
                return

            source = normpath(join(self.srcdir, source))
            target = normpath(join(self.outdir, source))
            target = basename(target) if self.outdir == '.' else target

            outname, outext = splitext(target)
            target = outname + '.yuv.md5'
            try:
                self.error('digest -i %s -o %s\n' % (basename(source), basename(target)))
                lines = executor.digest_by_frames(source, target, frames)
                self.write('%8d frames digest: %s\n' % (len(lines), basename(source)))
            except Exception as e:
                self.error('# %s\n' % e.message)

        def compare(executor, source, target):
            outname, outext = splitext(target)
            target = outname + '.yuv.md5'
            try:
                self.error('digest -i %s -o %s\n' % (basename(source), basename(target)))
                lines = executor.compare(source, target)
                self.write('%8d frames passed: %s\n' % (len(lines), basename(source)))
            except Exception as e:
                self.write('       # failed: %s\n' % e.message)
                self.error('       # failed: %s\n' % e.message)

        def noslot(executor, *largs, **kwargs):
            raise Exception('action must be one of %s' % executor.actions)

        def globber(method, executor, *largs):
            method = partial(method, executor)
            sources = largs if len(largs) > 0 else self.includes
            for file in sources:
                for source, target in self.glob(file):
                    method(source, target)

        actions = {'decode':decode, 'encode':encode, 'digest':digest, 'compare':compare}

        for action, method in actions.iteritems():
            if not isinstance(executor, ModelExecutor):
                method = partial(globber, method)
            elif action in executor.actions:
                method = partial(globber, method, executor)
            else:
                method = partial(noslot, executor)
            setattr(self, action, method)

        if not isinstance(executor, ModelExecutor):
            method = digest_by_frames
        elif 'digest_by_frames' in executor.actions:
            method = partial(digest_by_frames, executor)
        else:
            method = partial(noslot, executor)
        setattr(self, 'digest_by_frames', method)

        return self

    def unbind(self):
        self.bind(None)
        return self
