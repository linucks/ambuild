#!/usr/bin/env python
'''
Created on Feb 3, 2013

@author: jmht

Utility functions
'''
import sys
PYTHONFLAVOUR = sys.version_info[0]
import gzip
import logging
import os
import numpy
if PYTHONFLAVOUR < 3:
    import cPickle as pickle
else:
    import pickle

from ab_paths import PARAMS_DIR
import xyz_util

HOOMDVERSION = None
try:
    import hoomd
    HOOMDVERSION = hoomd.__version__
    del hoomd
except Exception:
    pass

PKL_SUFFIX = '.pkl'
GZIP_PKL_SUFFIX = '.pkl.gz'

logger = logging.getLogger()

def cellFromPickle(pickleFile, paramsDir=None):
    """Recreate a cell from a pickled file and apply any hacks so that we can work with older versions"""
    def fixFragment(fragment):
        for e in fragment._endGroups:
            # More horrible hacks for old versions
            if hasattr(e,'_isBonded'):
                e.bonded = e._isBonded
            if not hasattr(e,'blocked'):
                e.blocked = False
        # Need to make sure coords and masses are numpy arrays
        if type(fragment._coords) is list:
            fragment._coords = numpy.array(fragment._coords)
            fragment._masses = numpy.array(fragment._masses)
        if not hasattr(fragment,'_atomTypes'):
            fragment._atomTypes = fragment._types
            fragment._sharedAttrs['_atomTypes'] = None
        if not hasattr(fragment,'solvent'):
            # Solvent is a new attribute so we set to false
            fragment.solvent = False
            fragment._sharedAttrs['solvent'] = None
        if not hasattr(fragment,'onbondFunction'):
            fragment.onbondFunction = None
            fragment._sharedAttrs['onbondFunction'] = None
        if not hasattr(fragment,'markBonded'):
            fragment.markBonded = None
            fragment._sharedAttrs['markBonded'] = None
        if not hasattr(fragment,'unBonded'):
            fragment.unBonded = [ False ] * len(fragment._coords)
            fragment._individualAttrs['unBonded'] = None
        if not hasattr(fragment,'catalyst'):
            fragment.catalyst = False
            fragment._sharedAttrs['catalyst'] = False
        if not hasattr(fragment,'config'):
            fragment.config = [True if eg.free() else False for eg in fragment._endGroups]
            fragment._individualAttrs['config'] = fragment.config
        if not hasattr(fragment,'configStr'):
            fragment.configStr = 'XX' # Not sure if worth recalculating these 
            fragment._individualAttrs['configStr'] = fragment.configStr
        return
    
    if not os.path.isfile(pickleFile):
        raise RuntimeError("Cannot find file: {}".format(pickleFile))
    if pickleFile.endswith(GZIP_PKL_SUFFIX):
        compressed = True
        popen = gzip.open
    elif pickleFile.endswith(PKL_SUFFIX):
        compressed = False
        popen = open
    else:
        raise RuntimeError('Unrecognised pkl file suffix for file {}. Use pkl or pklz.'.format(pickleFile))
    mode = 'r'
    if compressed or PYTHONFLAVOUR == 3:
        mode += 'b'
    # This is another terrible hack for dealing with old pkl files that used the old class names
    import ab_block
    import ab_bond
    import ab_cell
    import ab_fragment
    # Patch the buildingBlock module so it has a reference to Bond
    ab_block.Bond = ab_bond.Bond
    sys.modules['buildingBlock'] = ab_block
    sys.modules['cell'] = ab_cell
    sys.modules['fragment'] = ab_fragment
    # Renamed cell class so need to alias here for old files
    with popen(pickleFile, mode) as f:
        myCell = pickle.load(f)
    del ab_block.Bond
    del sys.modules['buildingBlock']
    del sys.modules['cell']
    del sys.modules['fragment']
    # Need to hack to work with older versions
    if not hasattr(myCell, 'dim'):
        myCell.dim = numpy.array([myCell.A, myCell.B, myCell.C])
        myCell.numBoxes = [myCell.numBoxA, myCell.numBoxB, myCell.numBoxC]
    if not type(myCell.dim) is numpy.ndarray:
        myCell.dim = numpy.array(myCell.dim)
    if not hasattr(myCell, 'pbc'):
        myCell.pbc = [True, True, True]
        myCell.walls = [False, False, False]
    # Set parameter Directory
    if paramsDir is None:
        if hasattr(myCell, 'paramsDir'):
            paramsDir = myCell.paramsDir
        else:
            paramsDir = PARAMS_DIR
    if not os.path.isdir(paramsDir):
        raise RuntimeError("Cannot find cell paramsDir: {0}".format(paramsDir))
    logger.info("Getting parameter files from directory: {0}".format(paramsDir))
    xyz_util.setModuleBondLength(os.path.join(paramsDir,'bond_params.csv'))
    myCell.setMdEngine(HOOMDVERSION, paramsDir)
    # Fix all the fragments
    for fragment in myCell._fragmentLibrary.values():
        fixFragment(fragment)
    for block in myCell.blocks.values():
        for fragment in block.fragments:
            fixFragment(fragment)
    return myCell


def dumpPkl(pickleFile, split=None, nonPeriodic=False):

    fpath = os.path.abspath(pickleFile)
    logger.info("Dumping pkl file: {0}".format(fpath))
    _, fname = os.path.split(fpath)
    prefix = os.path.splitext(fname)[0]

    mycell = cellFromPickle(pickleFile)
    if split == "fragments":
        for t in mycell.fragmentTypes().keys():
            data = mycell.dataDict(fragmentType=t)
            mycell.writeXyz("{0}_{1}_P.xyz".format(prefix, t),
                            data=data,
                            periodic=True)
            mycell.writeCml("{0}_{1}_PV.cml".format(prefix, t),
                            data=data,
                            periodic=True,
                            pruneBonds=True)
    elif split == "blocks":
        periodic = True
        for i, b in enumerate(mycell.blocks.values()):
            # Write out each block to a separate file
            if periodic:
                b.writeXyz("{0}_block{1}.xyz".format(prefix, i),
                           cell=mycell.dim)
                b.writeCml("{0}_block{1}.cml".format(prefix, i),
                           cell=mycell.dim)
            else:
                b.writeXyz("{0}_block{1}.xyz".format(prefix, i))
                b.writeCml("{0}_block{1}.cml".format(prefix, i))
    else:
        if nonPeriodic:
            data = mycell.dataDict(rigidBody=False, periodic=False)
            mycell.writeCml(prefix + "_NP.cml", data,
                            periodic=False, pruneBonds=False)
            mycell.writeXyz(prefix + "_NP.xyz", data=data, periodic=False)
            mycell.writeXyz(prefix + "_NP_types.xyz", data=data, periodic=False, atomTypes=True)
        else:
            data = mycell.dataDict(rigidBody=False)
            mycell.writeXyz(prefix + "_P.xyz", data=data, periodic=True)
            mycell.writeXyz(prefix + "_P_types.xyz", data=data, periodic=True, atomTypes=True)
            # self.writeCar(prefix+"_P.car",data=data,periodic=True)
            mycell.writeCml(prefix + "_PV.cml", data=data, periodic=True, pruneBonds=True)
    return


def dumpDLPOLY(pickleFile, rigidBody=False, skipDihedrals=False):
    fpath = os.path.abspath(pickleFile)
    logger.info("Dumping DLPOLY files from pkl file: {0}".format(fpath))
    mycell = cellFromPickle(pickleFile)

    # Set parameter Directory
    if hasattr(mycell, 'paramsDir'):
        paramsDir = mycell.paramsDir
    else:
        paramsDir = PARAMS_DIR
    if not os.path.isdir(paramsDir):
        raise RuntimeError("Cannot find cell paramsDir: {0}".format(paramsDir))
    logger.info("Getting parameter files from directory: {0}".format(paramsDir))

    # Need to do this here or else hoomdblue gets the command line arguments on import of the module
    import dlpoly
    d = dlpoly.DLPOLY(paramsDir=paramsDir)
    d.writeCONTROL()
    d.writeFIELDandCONFIG(mycell, rigidBody=rigidBody, skipDihedrals=skipDihedrals)
    os.unlink(mycell.logfile)
    os.unlink(mycell.logcsv)
    return


def frange(start, stop, step):
    """
    Range function that works with floating points
    http://stackoverflow.com/questions/4189766/python-range-with-step-of-type-float
    """
    while start < stop:
        yield start
        start += step


def newFilename(filename, separator="_"):
    dname, name = os.path.split(filename)
    name, suffix = os.path.splitext(name)
    try:
        basename, num = name.split(separator)
    except ValueError:
        # No separator so assume is an un-numbered file
        return os.path.join(dname, name + separator + "1" + suffix)
    num = int(num) + 1
    return os.path.join(dname, basename + separator + str(num) + suffix)


def pickleObj(obj, fileName, compress=True):
    """Pickle an object - required as we can't pickle in the cell as otherwise the open filehandle
    is within the cell which is the object we are trying to pickle..."""
    if compress:
        assert fileName.endswith(GZIP_PKL_SUFFIX), 'Wrong suffix for gzip pkl file'
    else:
        assert fileName.endswith(PKL_SUFFIX), 'Wrong suffix for pkl file'
    mode = 'w'
    if compress:
        popen = gzip.open
    else:
        popen = open
        if PYTHONFLAVOUR == 3:
            mode += 'b'
    with popen(fileName, mode) as pfile:
            pickle.dump(obj, pfile)
    return fileName

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('pkl_file', type=str, metavar='pickle_file.pkl', help='Ambuild pickle file')
    p.add_argument('-f', '--split_fragments', action='store_true', default=False, help="Split the cell into fragments")
    p.add_argument('-b', '--split_blocks', action='store_true', default=False, help="Split the cell into blocks")
    p.add_argument('-da', '--dlpoly_allatom', action='store_true', default=False, help="Create all-atom DLPOLY CONFIG and FIELD files")
    p.add_argument('-dr', '--dlpoly_rigid', action='store_true', default=False, help="Create rigid-body DLPOLY CONFIG and FIELD files")
    p.add_argument('-np', '--non_periodic', action='store_true', default=False, help="Dump a non-periodic system")
    p.add_argument('-id', '--ignore_dihedrals', action='store_true', default=False, help="Ignore missing dihedrals when creating DL-POLY files.")
    args=p.parse_args()
    
    split = None
    dlpoly = False
    rigid = True
    nonPeriodic = False
    skipDihedrals = False
    
    if args.split_fragments:
        split='fragments'
    elif args.split_blocks:
        split='blocks'
        
    if args.dlpoly_allatom:
        dlpoly = True
        rigid = False
    elif args.dlpoly_rigid:
        dlpoly = True
        rigid = True
        
    if args.non_periodic:
        nonPeriodic = True
    if args.ignore_dihedrals:
        skipDihedrals = True     

    # Need to reset sys.argv as otherwise hoomdblue eats it and complains
    sys.argv = [sys.argv[0]]

    if dlpoly:
        dumpDLPOLY(args.pkl_file, rigidBody=rigid, skipDihedrals=skipDihedrals)
    else:
        dumpPkl(args.pkl_file, split=split, nonPeriodic=nonPeriodic)
