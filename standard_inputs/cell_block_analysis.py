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
prev_block_data = {}
for p in sorted(glob.glob("*.pkl*")):
    print("Processing: {}".format(p))
    nstep  = re.match('step_([0-9]+)\.pkl.*',p).group(1)
    mycell = ab_util.cellFromPickle(p, paramsDir=paramsDir)
    nblocks = 0
    max_mass = 0
    zip_bonds = 0
    current_block_data = {}
    for k, b in mycell.blocks.items():
        if len(b.fragments) == 1 and b.fragments[0].fragmentType in ignore:
            continue
        current_block_data[k] = len(b.blockBonds())
        if k in prev_block_data and prev_block_data[k] < current_block_data[k]:
            zip_bonds += current_block_data[k] - prev_block_data[k]
        nblocks += 1
        max_mass = max(max_mass, b.blockMass())
    #print("SAME ",set(current_block_data.keys()).intersection(set(prev_block_data.keys())))
    prev_block_data = current_block_data
    data.append([nstep,nblocks,max_mass, zip_bonds])

with open('blocks.csv', 'w') as csvfile:
    w = csv.writer(csvfile, delimiter=',')
    w.writerow(['nstep','nblocks','max_mass', 'zip_bonds'])
    for d in data: 
        w.writerow(d)
