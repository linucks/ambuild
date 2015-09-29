#!/usr/bin/env python
'''
Created on 24 Sep 2015

@author: jmht
'''

import os,sys
import time

def read_cell(fh):
    dim=[None,None,None]
    dima=[None,None,None]
    tag,val=fh.readline().split()
    assert tag=='lat_a'
    dim[0]=float(val)
    tag,val=fh.readline().split()
    assert tag=='lat_b'
    dim[1]=float(val)
    tag,val=fh.readline().split()
    assert tag=='lat_c'
    dim[2]=float(val)
    tag,val=fh.readline().split()
    assert tag=='alpha'
    dima[0]=float(val)
    tag,val=fh.readline().split()
    assert tag=='beta'
    dima[1]=float(val)
    tag,val=fh.readline().split()
    assert tag=='gamma'
    dima[2]=float(val)
    fh.readline()
    return dim,dima

def read_charges(fh):
    atoms,charges=[],[]
    for _ in range(4): fh.readline() # skip 4 lines
    while True:
        line = fh.readline()
        fields = line.split()
        if not fields or fields[0]=='Total': break
        if not len(fields)==5: raise RuntimeError,"Error parsing line: {0}".format(line)
        atoms.append(fields[1])
        charges.append(float(fields[4]))
    if not len(atoms) and len(atoms)==len(charges): raise RuntimeError,"Error reading charges from: {0}".format(fh.name)
    return atoms,charges
                
def read_geometry(fh):
    symbols,coords=[],[]
    line = fh.readline()
    i=0
    while line:
        i+=1
        line=fh.readline()
        if line.strip().startswith('No.       Tag          Charge          X              Y              Z'):
            fh.readline()
            break
        assert i < 10,"Error reading geometry"
    i=0
    while True:
        line=fh.readline()
        fields=line.split()
        if len(fields) != 6: break
        count,tag,charge,x,y,z=fields
        if i==0:
            assert int(count)==1
            i=1
        symbols.append(tag)
        x,y,z=float(x),float(y),float(z)
        coords.append((x,y,z))
    return symbols,coords
        
def read_nwchem_output(nwchem_file):
    symbols,coords,charges,dim,dima = None,None,None,None,None
    with open(nwchem_file) as f:
        line = f.readline()
        # Read to start of charges
        while line:
            line = f.readline()
            if line.strip().lower().startswith('system crystal'):
                dim,dima=read_cell(f)
            if line.strip().startswith('Geometry "geometry" ->'):
                symbols,coords=read_geometry(f)
            if line.startswith(' charge analysis on each atom'):
                _,charges=read_charges(f)
                break
    assert symbols and coords and dim and dima and charges
    return symbols,coords,charges,dim,dima

def write_car_file(car_file,symbols,coords,charges,dim,dima):
    car = "!BIOSYM archive 3\n"
    car += "PBC=ON\n"
    car += "ambuild generated car file\n"
    tstr = time.strftime("%a %b %d %H:%M:%S %Y", time.gmtime())
    car += "!DATE {0}\n".format(tstr)
    car += "PBC  {0: < 9.4F} {1: < 9.4F} {2: < 9.4F} {3: < 9.4F} {4: < 9.4F} {5: < 9.4F}(P1)\n".format(dim[0],
                                                                                                       dim[1],
                                                                                                       dim[2],
                                                                                                       dima[0],
                                                                                                       dima[1],
                                                                                                       dima[2],)
# 1-4    Atom name.
# 6-20    x Cartesian coordinate for the atom (angstrom).
# 21-35    y Cartesian coordinate for the atom (angstrom).
# 36-50    z Cartesian coordinate for the atom (angstrom).
# 52-55    Name of residue containing atom.
# 56-60    Residue sequence number relative to the  beginningof the current molecule.
# 62-65    Potential function atom type (left justified) (ignored; see Molecular Data File).
# 71-72    Element type.
# 74-79    Partial charge on the atom.

    for i, coord in enumerate(coords):
            x = coord[0]
            y = coord[1]
            z = coord[2]
            atype = 'xx'
            charge = charges[i]
            symbol = symbols[i]
            label = '{0}{1}'.format(symbol,i+1)
            #car += "{0: <5} {1: >15.10} {2: >15.10} {3: >15.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format(label, x, y, z, atype, symbol, charge)
            car += "{0: <4} {1: >15.10}{2: >15.10}{3: >15.10} XXXX    1 {4: <4}      {5: <2} {6: > 2.3f}\n".format(label,
                                                                                                                     x,
                                                                                                                     y,
                                                                                                                     z,
                                                                                                                     atype,
                                                                                                                     symbol,
                                                                                                                     charge)
    car += "end\nend\n\n"
    
    with open(car_file, 'w') as f:
        fpath = os.path.abspath(f.name)
        f.writelines(car)
    
    print "Wrote car file: {0}".format(fpath)
    return

if not len(sys.argv)==2 or not os.path.isfile(sys.argv[1]):
    print "Usage: {0} <nwchem_output_file>"
    sys.exit(1)

nwchem_file=sys.argv[1]
name=os.path.splitext(os.path.basename(nwchem_file))[0]
car_file=os.path.abspath(name+'.car')
               
symbols,coords,charges,dim,dima = read_nwchem_output(nwchem_file)
write_car_file(car_file,symbols,coords,charges,dim,dima)

sys.exit()

# 1-4    Atom name.
# 6-20    x Cartesian coordinate for the atom (angstrom).
# 21-35    y Cartesian coordinate for the atom (angstrom).
# 36-50    z Cartesian coordinate for the atom (angstrom).
# 52-55    Name of residue containing atom.
# 56-60    Residue sequence number relative to the  beginningof the current molecule.
# 62-65    Potential function atom type (left justified) (ignored; see Molecular Data File).
# 71-72    Element type.
# 74-79    Partial charge on the atom.
with open(car_file) as f:
    for i,line in enumerate(f):
        if i != 5: continue
        print "Name ",line[0:4]
        print "x ",line[5:20]
        print "y ",line[20:35]
        print "z ",line[35:50]
        print "rname ",line[51:55]
        print "rsn ",line[55:60]
        print "rsn ",line[61:65]
        print "elemen ",line[70:72]
        print "charge ",line[73:79]

