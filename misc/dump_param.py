'''
Created on 28 May 2015

@author: jmht
'''

import csv
import math
import os
import sys

thisd =  os.path.abspath(os.path.dirname(__file__ ))
builder_dir = os.path.abspath(os.path.join(thisd,'..','builder'))

print builder_dir
sys.path.insert(0,builder_dir)

from opt import FfieldParameters

ffield = FfieldParameters()

# bonds
# 'cp-cp'   : { 'k' : 1550.0, 'r0' : 1.384 }

pfile = os.path.join(builder_dir,'bond_params.csv')
header = ['bond','k','r0','comments']
with open(pfile,'w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(header)
    for bond in ffield.bonds.keys():
        k = ffield.bonds[bond]['k']
        r0 = ffield.bonds[bond]['r0']
        writer.writerow([bond,k,r0,""])

# angles
# 'h-c-c'    : { 'k' : 330.0, 't0' : math.radians(109) }
pfile = os.path.join(builder_dir,'angle_params.csv')
header = ['angle','k','t0','comments']
with open(pfile,'w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(header)
    for angle in ffield.angles.keys():
        k = ffield.angles[angle]['k']
        t0 = round(math.degrees(ffield.angles[angle]['t0']))
        writer.writerow([angle,k,t0,""])

# dihedrals
# 'h-c-c-h'     : { 'k' : 70, 'd' : 1, 'n' : 1 }
pfile = os.path.join(builder_dir,'dihedral_params.csv')
header = ['dihedral','k','d','n','comments']
with open(pfile,'w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(header)
    for dihedral in ffield.dihedrals.keys():
        k = ffield.dihedrals[dihedral]['k']
        d = ffield.dihedrals[dihedral]['d']
        n = ffield.dihedrals[dihedral]['n']
        writer.writerow([dihedral,k,d,n,""])

# impropers
# 'cp-cp-cp-np'       : { 'k' : 200, 'chi' : 0.0 }
pfile = os.path.join(builder_dir,'improper_params.csv')
header = ['improper','k','chi','comments']
with open(pfile,'w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(header)
    for improper in ffield.impropers.keys():
        k = ffield.impropers[improper]['k']
        chi = ffield.impropers[improper]['chi']
        writer.writerow([improper,k,chi,""])

# pairs
# ('c',  'h' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  }
pfile = os.path.join(builder_dir,'pair_params.csv')
header = ['atom1','atom2','epsilon','sigma','comments']
with open(pfile,'w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(header)
    for pair in ffield.pairs.keys():
        atom1, atom2 = pair
        epsilon = ffield.pairs[pair]['epsilon']
        sigma = ffield.pairs[pair]['sigma']
        writer.writerow([atom1, atom2, epsilon, sigma, ""])

