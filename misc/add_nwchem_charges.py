#!/usr/bin/env python
'''
Created on 24 Sep 2015

@author: jmht
'''

import os,sys

def read_nwchem_charges(nwchem_file):
    atoms, charges = [],[]
    with open(nwchem_file) as f:
        line = f.readline()
        # Read to start of charges
        while line:
            line = f.readline()
            if line.startswith(' charge analysis on each atom'):
                break
        if not line.startswith(' charge analysis on each atom'):
            raise RuntimeError, 'Could not find charges in file: {0}'.format(nwchem_file)
        
        for _ in range(4): f.readline() # skip 4 lines
        while True:
            line = f.readline()
            fields = line.split()
            if not fields or fields[0]=='Total': break
            if not len(fields)==5: raise RuntimeError,"Error parsing line: {0}".format(line)
            atoms.append(fields[1])
            charges.append(float(fields[4]))
    if not len(atoms) and len(atoms)==len(charges): raise RuntimeError,"Error reading charges from: {0}".format(nwchem_file)
    return atoms,charges

if not len(sys.argv)==3 or not os.path.isfile(sys.argv[1]) or not os.path.isfile(sys.argv[2]):
    print "Usage: {0} <nwchem_output_file> <car_file>"
    sys.exit(1)

nwchem_file=sys.argv[1]
car_file=sys.argv[2]
               
atoms,charges = read_nwchem_charges(nwchem_file)

newlines=[]
with open(car_file) as f: lines = f.readlines()
istart=-1
for i, line in enumerate(lines):
    newlines.append(line)
    if line.startswith("!DATE"):
        if lines[i+1].startswith('PBC'): istart=i+2
        else: istart=i+1
if istart==-1: print "Cannot read car file: {0}".format(car_file)

print atoms
for i, line in enumerate(lines[istart:]):
    print "line ",i, atoms[i],line
    fields=line.split()
    if len(fields) != 9:
        assert fields[0].startswith('end'),"Error reading car file: {0}".format(car_file)
        newlines.append('end\n')
        newlines.append('end\n')
        break
    
    runin=line[:73]
    atom=fields[7]
    charge=float(fields[8])
    assert atoms[i]==atom,"Mismatching atoms: {0}: {1}".format(atoms[i],line)
