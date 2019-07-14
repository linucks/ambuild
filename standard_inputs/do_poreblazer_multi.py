#!/usr/bin/env python
import glob
import sys
sys.path.append("/opt/ambuild.git")
from ambuild import ab_util

paramsDir = '/home/pierre/Dropbox/Ambuild_Files/Parameters'
poreblazerExe = '/opt/poreblazer/src/poreblazer.exe'

for spkl in sorted(glob.glob('step_*.pkl.gz')):
    print('Processing pkl file: {} '.format(spkl))
    mycell = ab_util.cellFromPickle(spkl, paramsDir=paramsDir)
    mycell.poreblazer(poreblazerExe)
