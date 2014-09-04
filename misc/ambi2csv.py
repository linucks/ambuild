#!/usr/bin/env python

import sys
import os
import glob

for a in glob.glob("*.ambi"):
    n=os.path.splitext(a)[0]
    with open(n+".csv",'w') as w,open(a) as f:
        for i,line in enumerate(f):
            if i==0:
                w.write("ntype,endgroup,capatom,dihedral,delatom\n")
                continue
            tokens=line.strip().split()
            w.write(",".join(tokens)+"\n")
        w.write("\n")
