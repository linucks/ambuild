#!/usr/bin/env python
'''
Created on 14 May 2016

@author: jmht
'''

import os
import sys
import unittest
from paths import BUILDER_DIR

VERBOSITY = 2
suite = unittest.TestLoader().discover(BUILDER_DIR)
if int(suite.countTestCases()) <= 0:
    msg = 'Could not find any tests to run in directory: {0}'.format(BUILDER_DIR) + os.linesep
    sys.stderr.write(msg)
    sys.exit(1)

unittest.TextTestRunner(verbosity=VERBOSITY, buffer=False).run(suite)
