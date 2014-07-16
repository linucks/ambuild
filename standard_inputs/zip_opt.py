#!/usr/bin/env python
import sys
sys.path.append("/opt/ambuild/builder")

import cell

mycell = util.cellFromPickle("step_628.pkl")

added = mycell.zipBlocks(bondMargin=4.0, bondAngleMargin=40)
if added > 0:
      mycell.optimiseGeometry(doDihedral=True, optCycles = 50000, dt=0.001) 
      mycell.optimiseGeometry(doDihedral=True, finc=1.2, fdec=1.0)
      mycell.dump()#4
