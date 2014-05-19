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

 File      : smpte_vc1.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     1.0 May 20, 2013    first release
     2.0 May 12, 2014    Executor classify

================================================================================
'''

__all__ = ('SMPTEVc1', )

__version__ = '2.0.0'

from . import rootpath, ModelExecutor


class SMPTEVc1(ModelExecutor):

    model   = 'smpte-vc1'
    codecs  = ('vc1', )
    actions = ('decode', 'digest_by_frames', 'compare')

    def __init__(self, codec='vc1', **kwargs):
        from os.path import join

        super(SMPTEVc1, self).__init__(codec, **kwargs)

        binary = 'tool/3rd-party/smpte-vc1-2008/decoder/decoder'

        self._execute = join(rootpath, binary)

        self._config = '''\
# Options for SMPTE VC-1 reference decoder.
# Comments in this file can start with # or ; at any point in the line
; Blank lines are allowed.

#if decoder
DebugMask      : 0x00000000     # See debug.h for bits to set/clear
#endif

Bitstream File : %s
Output YUV     : %s
'''

    def execute(self):
        return self._execute

    def options(self, source, target):
        return [source]

    def decode(self, source, target):
        from subprocess import call
        from os import remove
        from os.path import splitext, basename

        srcname, srcext = splitext(basename(source))
        optname = srcname + '.opt'
        with open(optname, 'wt') as f:
            f.write(self._config % (source, target))

        execute = self.execute()
        options = self.options(optname, target)

        self.mkdir(target)

        try:
            call([execute] + options, stdout=self.stdout, stderr=self.stderr)
            remove(optname)
        except:
            raise Exception('decode error: %s' % basename(source))
