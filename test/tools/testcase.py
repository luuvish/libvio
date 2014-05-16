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

 File      : testcase.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 13, 2014    Executor classify

================================================================================
'''

from sys import path, argv, stdout, stderr
from os import remove
from os.path import join, dirname

path.append(join(dirname(__file__), '..'))

from tools import Globber
from tools.model import ModelExecutor
from tools.case import h264, hevc, vc1, vp8, vp9


codecs = (h264, hevc, vc1, vp8, vp9)
models = (model for codec in codecs for model in codec.models)
cases  = (case for codec in codecs for case in codec.cases)


class TestCase(object):

    program     = 'testcase.py'
    description = 'shell command helper for testcase of video codecs'
    epilog      = 'Copyright (C) 2014 Thumb \'o Cat luuvish (luuvsih@gmail.com)'

    def __init__(self, *largs):

        self.executors = {}
        self.models    = []
        self.testcases = {}
        self.cases     = []

        for arg in largs:
            if type(arg) is dict:
                case = arg['case']
                if case not in self.testcases:
                    self.testcases[case] = arg
                    self.cases.append(case)
            elif issubclass(arg, ModelExecutor):
                model = arg.model
                if model not in self.models:
                    self.executors[model] = arg
                    self.models.append(model)

        self.models.sort()
        self.cases.sort()

    def parse(self, *largs):
        from argparse import ArgumentParser

        parser = ArgumentParser(prog=self.program,
                description=self.description, epilog=self.epilog)

        parser.add_argument('case', help='test case',
                action='store', choices=self.cases)
        parser.print_help = self.help

        parser.parse_args(args=largs, namespace=self)

        return vars(self)

    def help(self):
        message = '''\
usage: %s [-h] testcases...

%s

positional arguments:   testcases

%s

optionsal arguments:
  -h, --help            show this help message and exit

%s
'''
        stderr.write(message % (self.program, self.description, '\n'.join(self.cases), self.epilog))

    def main(self, *largs):

        options = self.parse(*(largs[1:] if len(largs) > 1 else ['-h']))

        case = options.get('case')

        testcase = self.testcases[case]
        model  = testcase.get('model')
        codec  = testcase.get('codec')
        action = testcase.get('action')

        if model is None:
            raise Exception('testcase model is none')
        if codec is None:
            raise Exception('testcase codec is none')
        if action is None:
            raise Exception('testcase action is none')

        logout = open('.testcase.log', 'wt')

        executor = model
        if type(executor) is str:
            executor = self.executors[executor]
        if not issubclass(executor, ModelExecutor):
            raise Exception('Executor is not subclass of ModelExecutor')
        executor = executor(codec, stdout=logout, stderr=logout)

        globber = Globber(**{
            'stdout': testcase.get('stdout', stdout),
            'stderr': testcase.get('stderr', stderr),
            'srcdir': testcase.get('srcdir', '.'),
            'outdir': testcase.get('outdir', '.'),
            'includes': testcase.get('includes', []),
            'excludes': testcase.get('excludes', [])
        }).bind(executor)

        method = getattr(globber, action)
        method()

        remove('.testcase.log')


if __name__ == '__main__':
    try:
        TestCase(*(list(models) + list(cases))).main(*argv)
    except Exception as e:
        stderr.write('%s\n' % e.message)
