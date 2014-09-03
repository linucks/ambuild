#!/usr/bin/env python

import sys
import os
import shutil
import subprocess

if len(sys.argv) != 3:
    print "Usage is: {0} <number_of_times_to_run> <script_to_run>".format(sys.argv[0])
    sys.exit(1)

toRun=int(sys.argv[1])
script=os.path.abspath(sys.argv[2])

if not os.path.isfile(script):
    print "Cannot find script file: {0}".format(script)
    sys.exit(1)

if not os.access(script, os.X_OK):
    print "Script {0} is not executable!".format(script)
    sys.exit(1)

# get basename of the files
name=os.path.splitext(os.path.basename(script))[0]

runDir=os.getcwd()

for i in range(toRun):
    jobNum=i+1
    dir=os.path.join(runDir,"{0}_{1}".format(name,jobNum))
    log=os.path.join(dir,"{0}_{1}.log".format(name,jobNum))
    os.mkdir(dir)
    shutil.copy(script,dir)
    os.chdir(dir)
    print "Running job {0} in directory: {1}. Output will be sent to the file: {2}".format(jobNum,dir,log)
    try:
        f=open(log,'w')
        subprocess.check_call(script,stdout=f,stderr=subprocess.STDOUT)
        f.close()
    except Exception,e:
        print "ERROR RUNNING JOB {0}\n{1}".format(jobNum,e)
    os.chdir(runDir)



