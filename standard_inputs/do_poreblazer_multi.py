#!/usr/bin/env python
import glob
import multiprocessing
import os
import sys
sys.path.append("/opt/ambuild.git")
from ambuild import ab_util

# This script should be run in a directory containing ambuild pkl files (step_X.pkl.gz)
# Each pkl file needs to be numbered differently (renaming the file won't help).

# Edit these to suit your case
num_processors = 6
params_dir = '/home/patrick/Dropbox/Ambuild_Files/Parameters'
poreblazer_exe = '/opt/poreblazer/src/poreblazer.exe'
#
# Don't change anything below here
#
def run_poreblazer(pkl_file):
    print('Processing pkl file: {} '.format(pkl_file))
    mycell = ab_util.cellFromPickle(pkl_file, paramsDir=params_dir)
    mycell.poreblazer(poreblazer_exe)

pklglob = 'step_*.pkl.gz'
pklfiles = glob.glob(pklglob)
# Abbie-proof
if not ab_util.is_exe(poreblazer_exe):
    sys.stderr.write("Cannot find poreblazer executable: {}\n".format(poreblazer_exe))
    sys.exit(1)
if not pklfiles:
    sys.stderr.write("Cannot find any files named {} in directory: {}\n".format(pklglob, os.path.abspath(os.getcwd())))
    sys.exit(1)

pool = multiprocessing.Pool(processes=num_processors)
pool.map(run_poreblazer, pklfiles)
