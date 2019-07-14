#!/usr/bin/env python
import sys
sys.path.append("/opt/ambuild.git")
from ambuild import ab_util

paramsDir = '/home/pierre/Dropbox/Ambuild_Files/Parameters'
mycell = ab_util.cellFromPickle('step_1.pkl.gz', paramsDir=paramsDir)

poreblazerExe = '/opt/poreblazer/src/poreblazer.exe'
mycell.poreblazer(poreblazerExe)
