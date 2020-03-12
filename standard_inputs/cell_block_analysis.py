#!/usr/bin/env python2
import csv
import os
import glob
import re
import sys
sys.path.append("/opt/ambuild.git")
from ambuild import ab_util

def pkl_step_num(pkl_filename):
    return re.match('step_([0-9]+)\.pkl.*',pkl_filename).group(1)

paramsDir = '/home/pierre/Dropbox/Ambuild_Files/Parameters'
ignore = ['DMF', 'TEA', 'cat']

prev_block_data = {}
with open('blocks.csv', 'w') as csvfile:
    w = csv.writer(csvfile, delimiter=',')
    w.writerow(['nstep','nblocks','max_mass', 'zip_bonds'])
    for p in sorted(glob.glob("*.pkl*"), key=pkl_step_num):
        print("Processing: {}".format(p))
        nstep  = pkl_step_num(p)
        try:
            mycell = ab_util.cellFromPickle(p, paramsDir=paramsDir)
        except EOFError:
            print("Picklefile: {} was empty!".format(p))
            continue
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
        prev_block_data = current_block_data
        w.writerow([nstep,nblocks,max_mass, zip_bonds])

