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

 File      : testsuite.py
 Author(s) : Luuvish
 Version   : 2.0
 Revision  :
     2.0 May 13, 2014    Executor classify

================================================================================
'''

from sys import path, argv, stdout, stderr
from os import remove
from os.path import dirname

path.append(dirname(__file__))

from test.globber import Globber
from test.model import ModelExecutor
from test.suite import h264, hevc, vc1, vp8, vp9


codecs = (h264, hevc, vc1, vp8, vp9)
models = (model for codec in codecs for model in codec.models)
suites = (suite for codec in codecs for suite in codec.suites)


class TestSuite(object):

    program     = 'testsuite.py'
    description = 'shell command helper for testsuite of video codecs'
    epilog      = 'Copyright (C) 2014 Thumb \'o Cat luuvish (luuvsih@gmail.com)'

    def __init__(self, *largs):

        self.executors  = {}
        self.models     = []
        self.testsuites = {}
        self.suites     = []

        for arg in largs:
            if type(arg) is dict:
                suite = arg['suite']
                if suite not in self.testsuites:
                    self.testsuites[suite] = arg
                    self.suites.append(suite)
            elif issubclass(arg, ModelExecutor):
                model = arg.model
                if model not in self.models:
                    self.executors[model] = arg
                    self.models.append(model)

        self.models.sort()
        self.suites.sort()

    def parse(self, *largs):
        from argparse import ArgumentParser

        parser = ArgumentParser(prog=self.program,
                description=self.description, epilog=self.epilog)

        parser.add_argument('suite', help='test suite',
                action='store', choices=self.suites)
        parser.print_help = self.help

        parser.parse_args(args=largs, namespace=self)

        return vars(self)

    def help(self):
        message = '''\
usage: %s [-h] testsuites...

%s

positional arguments:   testsuites

%s

optionsal arguments:
  -h, --help            show this help message and exit

%s
'''
        stderr.write(message % (self.program, self.description, '\n'.join(self.suites), self.epilog))

    def main(self, *largs):
        import time

        options = self.parse(*(largs[1:] if len(largs) > 1 else ['-h']))

        suite  = options.get('suite')
        testsuite = self.testsuites[suite]
        model  = testsuite.get('model')
        codec  = testsuite.get('codec')
        action = testsuite.get('action')

        if model is None:
            raise Exception('testsuite model is none')
        if codec is None:
            raise Exception('testsuite codec is none')
        if action is None:
            raise Exception('testsuite action is none')

        logout = open('.testsuite.log', 'wt')

        executor = model
        if type(executor) is str:
            executor = self.executors[executor]
        if not issubclass(executor, ModelExecutor):
            raise Exception('Executor is not subclass of ModelExecutor')
        executor = executor(codec, stdout=logout, stderr=logout)

        globber = Globber(**{
            'stdout': testsuite.get('stdout', stdout),
            'stderr': testsuite.get('stderr', stderr),
            'srcdir': testsuite.get('srcdir', '.'),
            'outdir': testsuite.get('outdir', '.'),
            'includes': testsuite.get('includes', []),
            'excludes': testsuite.get('excludes', [])
        }).bind(executor)

        start = time.time()

        method = getattr(globber, action)
        method()

        end = time.time()
        duration = time.gmtime(end - start)
        globber.write('complete in %s.\n' % time.strftime('%H:%M:%Ss', duration))

        remove('.testsuite.log')


if __name__ == '__main__':
    try:
        TestSuite(*(list(models) + list(suites))).main(*argv)
    except Exception as e:
        stderr.write('%s\n' % e.message)
