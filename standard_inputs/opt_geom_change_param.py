#!/usr/bin/env python
import sys
sys.path.append("/opt/ambuild/builder")

import util
import cell

mycell = util.cellFromPickle("step_628.pkl")
mycell.optimiseGeometry(doDihedral=True, optCycles = 50000, dt=0.001) 
mycell.optimiseGeometry(doDihedral=True, finc=1.2, fdec=1.0)
