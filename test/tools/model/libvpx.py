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

 File      : libvpx.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 May 21, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('LibVpx', )

__version__ = '2.0.0'

from . import ModelExecutor


class LibVpx(ModelExecutor):

    model   = 'libvpx'
    codecs  = ('vp8', 'vp9')
    actions = ('decode', 'digest', 'digest_by_frames', 'compare')

    def __init__(self, codec, **kwargs):
        from os.path import join, normpath, dirname

        super(LibVpx, self).__init__(codec, **kwargs)

        root   = normpath(join(dirname(__file__), '../../..'))
        binary = 'ref/dist/libvpx-1.3.0.bin/dist/bin/vpxdec'

        self._execute = join(root, binary)

        self.defaults['decode'] += ['--codec=' + codec, '--i420']
        self.defaults['digest'] += ['--codec=' + codec, '--i420', '--md5']

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-o', target, source]

    def digest(self, source, target=None):
        from subprocess import call
        from os import remove
        from os.path import exists, splitext, basename

        outname, outext = splitext(basename(source))
        output = '.' + outname + '.yuv.md5'

        execute = self.execute()
        options = self.defaults['digest'] + self.options(source, '%s-%%wx%%h-%%4.i420' % outname)

        try:
            f = open(output, 'wt')
            call([execute] + options, stdout=f, stderr=self.stderr)
            f.close()
        except Exception as e:
            raise e
        if not exists(output):
            raise Exception('digest error: %s' % basename(source))

        lines = []
        with open(output, 'rt') as f:
            lines = [line.rstrip().split()[0].lower() for line in f]
        remove(output)

        return lines
