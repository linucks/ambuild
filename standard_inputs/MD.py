#!/usr/bin/env python

import sys
sys.path.append("/opt/ambuild/builder")

import util
import cell

mycell = util.cellFromPickle("step_628.pkl")

mycell.runMD(doDihedral=True, rCut=5)

