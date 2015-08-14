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
header = ['A','B','k','r0','comments']
with open(pfile,'w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(header)
    for bond in ffield.bonds.keys():
        a,b = bond.split('-')
        k = ffield.bonds[bond]['k']
        r0 = ffield.bonds[bond]['r0']
        writer.writerow([a,b,k,r0,""])

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


# Init code for FfieldParameters

"""
        thisd =  os.path.abspath(os.path.dirname(__file__ ))
        
        # REM = bonds are in alphabetical order as in the order of the atoms in a bond
        # e.g. c-h NOT h-c
        # All parameters calculated to fit PCFF adapted forcefield Holden et al j.phys.chem.c 2012
        # from combining rules calculated using dlpoly-prep (Willock)
        # bonds adapted from quartic PCFF
        bond_file = os.path.join(thisd,'bond_params.csv')
        header = ['bond','k','r0','comments']
        self.bonds = {}
        with open(bond_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header: raise RuntimeError,"Header for {0} file, should be: {1}".format(bond_file,",".join(header))
                    else: continue
                bond = row[0]
                k = float(row[1])
                r0 = float(row[2])
                self.bonds[bond] = {'k' : k, 'r0' : r0}
        
        # REM: angles are stored in degrees in the parameter file, but we convert to radians for our uses
        angle_file = os.path.join(thisd,'angle_params.csv')
        header = ['angle','k','t0','comments']
        self.angles = {}
        with open(angle_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header: raise RuntimeError,"Header for {0} file, should be: {1}".format(angle_file,",".join(header))
                    else: continue
                angle = row[0]
                k = float(row[1])
                t0 = math.radians(float(row[2]))
                self.angles[angle] = {'k' : k, 't0' : t0}
        
        dihedral_file = os.path.join(thisd,'dihedral_params.csv')
        header = ['dihedral','k','d','n','comments']
        self.dihedrals = {}
        with open(dihedral_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header: raise RuntimeError,"Header for {0} file, should be: {1}".format(dihedral_file,",".join(header))
                    else: continue
                dihedral = row[0]
                k = int(row[1])
                d = int(row[2])
                n = int(row[3])
                self.dihedrals[dihedral] = {'k' : k, 'd' : d, 'n' : n}

        improper_file = os.path.join(thisd,'improper_params.csv')
        header = ['improper','k','chi','comments']
        self.impropers = {}
        with open(improper_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header: raise RuntimeError,"Header for {0} file, should be: {1}".format(improper_file,",".join(header))
                    else: continue
                improper = row[0]
                k = int(row[1])
                chi = float(row[2])
                self.impropers[improper] = {'k' : k, 'chi' : chi}
                
        pairs_file = os.path.join(thisd,'pair_params.csv')
        header = ['atom1','atom2','epsilon','sigma','comments']
        self.pairs = {}
        with open(pairs_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header: raise RuntimeError,"Header for {0} file, should be: {1}".format(pairs_file,",".join(header))
                    else: continue
                atom1 = row[0]
                atom2 = row[1]
                epsilon = float(row[2])
                sigma = float(row[3])
                self.pairs[(atom1,atom2)] = {'epsilon' : epsilon, 'sigma' : sigma}

"""

