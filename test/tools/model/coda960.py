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

 File      : coda960.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 May 19, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('Coda960', )

__version__ = '2.0.0'

from . import ModelExecutor


class Coda960(ModelExecutor):

    model   = 'coda960'
    codecs  = ('h264', 'vc1', 'vp8')
    actions = ('decode', 'digest', 'digest_by_frames', 'compare')

    def __init__(self, codec, **kwargs):
        from os.path import join, normpath, dirname

        super(Coda960, self).__init__(codec, **kwargs)

        root     = normpath(join(dirname(__file__), '../../..'))
        binary   = 'ref/dist/coda960-v1.0.0/design/ref_c/bin/darwin'
        binaries = {'h264':'avcdecoder', 'vc1':'vc1decoder', 'vp8':'vpxdecoder'}
        executes = {k: join(root, binary, v) for k, v in binaries.iteritems()}

        self._execute = executes[codec]

        self.defaults['digest'] += ['-5']
        if codec == 'vp8':
            self.defaults['decode'] += ['--std', '2']
            self.defaults['digest'] += ['--std', '2']

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-i', source, '-o', target]
