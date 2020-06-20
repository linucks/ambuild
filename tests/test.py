#!/usr/bin/env python
"""
Created on 14 May 2016

@author: jmht
"""
import os
import sys
import unittest

TEST_DIR = "."
VERBOSITY = 2
suite = unittest.TestLoader().discover(TEST_DIR)
if int(suite.countTestCases()) <= 0:
    msg = (
        "Could not find any tests to run in directory: {0}".format(TEST_DIR)
        + os.linesep
    )
    sys.stderr.write(msg)
    sys.exit(1)

unittest.TextTestRunner(verbosity=VERBOSITY, buffer=False).run(suite)
