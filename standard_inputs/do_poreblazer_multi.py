#!/usr/bin/env python
import glob
import multiprocessing
import sys
sys.path.append("/opt/ambuild.git")
from ambuild import ab_util

# This script should be run in a directory containing ambuild pkl files (step_X.pkl.gz)
# Each pkl file needs to be numbered differently (renaming the file won't help).

# Edit these to suit your case
num_processors = 6
paramsDir = '/home/pierre/Dropbox/Ambuild_Files/Parameters'
poreblazerExe = '/opt/poreblazer/src/poreblazer.exe'
#
# Don't change anything below here
#
def run_poreblazer(pkl_file):
    print('Processing pkl file: {} '.format(pkl_file))
    mycell = ab_util.cellFromPickle(pkl_file, paramsDir=paramsDir)
    mycell.poreblazer(poreblazerExe)

pklfiles = glob.glob('step_*.pkl.gz')
pool = multiprocessing.Pool(processes=num_processors)
pool.map(run_poreblazer, pklfiles)
