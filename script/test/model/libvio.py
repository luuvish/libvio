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

 File      : libvio.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 Aug 28, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('LibVio', )

__version__ = '2.0.0'

from . import rootpath, ModelExecutor


class LibVio(ModelExecutor):

    model   = 'libvio'
    codecs  = ('h264', )
    actions = ('decode', 'digest_by_frames', 'compare')

    def __init__(self, codec, **kwargs):
        from os.path import join

        super(LibVio, self).__init__(codec, **kwargs)

        binary = 'bin/libvio'

        self._execute = join(rootpath, binary)

        self.defaults['digest'] += ['-5']

    def execute(self):
        return self._execute

    def options(self, source, target):
        return ['-i', source, '-o', target]

    def decode(self, source, target):
        from os import remove

        super(LibVio, self).decode(source, target)

        remove('dataDec.txt')
        remove('log.dec')
