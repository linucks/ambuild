#!/usr/bin/env python3

# Our imports
from ambuild import ab_cell
from ambuild import ab_util

import sys, getopt

inputfile = ''
try:
   opts, args = getopt.getopt(sys.argv[1:],"hi:",["ifile="])
except getopt.GetoptError:
   print (sys.argv[0], ' -i <restartfile>')
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print (sys.argv[0], ' -i <restartfile>')
      sys.exit(0)
   elif opt in ("-i", "--ifile"):
      inputfile = arg
sys.argv = [sys.argv[0]]
#print ('Restart file is', inputfile)

paramsDir='./params'

   #Create Cell and seed it with the blocks
mycell = ab_util.cellFromPickle(inputfile, paramsDir=paramsDir )
