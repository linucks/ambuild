#!/usr/bin/env python3
import csv
import os
import glob
import re
import sys
sys.path.append("/opt/ambuild.git")
from ambuild import ab_util

paramsDir = '/Users/jmht/Dropbox/Ambuild/Ambuild_Files/Parameters'

ignore = ['DMF', 'TEA', 'cat']
data = []
for p in glob.glob("*.pkl*"):
    nstep  = re.match('step_([0-9]+)\.pkl.*',p).group(1)
    mycell = ab_util.cellFromPickle(p, paramsDir=paramsDir)
    nblocks = 0
    max_mass = 0
    for i, b in enumerate(mycell.blocks.values()):
        if len(b.fragments) == 1 and b.fragments[0].fragmentType in ignore:
            continue
        nblocks += 1
        max_mass = max(max_mass, b.blockMass())
    data.append([nstep,nblocks,max_mass])


with open('blocks.csv', 'w') as csvfile:
    w = csv.writer(csvfile, delimiter=',')
    w.writerow(['nstep','nblocks','max_mass'])
    for d in data: 
        w.writerow(d)
