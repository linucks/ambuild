'''
Created on Jan 15, 2013

@author: abbietrewin
'''
VERSION = "a842475a64c9"
import collections
import copy
import logging
import math
import os
import random as _random
import sys
import time

# External modules
import numpy as np

# Our modules
import ab_analyse
import ab_block
import ab_rigidparticle
import ab_bond
import ab_celldata
import ab_endgroup
import ab_fragment
import ab_subunit
import ab_util
from ab_paths import PARAMS_DIR
import xyz_core
import xyz_util

BONDTYPESEP = "-"  # Character for separating bonds
ENDGROUPSEP = ":"  # Character for separating endGroups in bonds

logger = logging.getLogger(__name__)


class Cell():
    '''
    classdocs
    '''
    def __init__(self,
                 boxDim=None,
                 filePath=None,
                 atomMargin=0.5,
                 bondMargin=0.5,
                 bondAngleMargin=15,
                 debugLog=False,
                 paramsDir=None):
        '''Construct an empty cell:

        Args:
        boxDim - a list with three numbers specifying the size of the cell A,B and C dimensions (angstroms)
                 Dimensions are from 0 - A, B or C
        filePath - a (.car) file with the cell dimensions and the coordinates of molecules that will be kept
                   static throught the simulation.
        atomMargin - the additional distance that will be added on to the VdW radii of two atoms
                     to determine if they are close enough to clash.
        bondMargin - two atoms are considered close enough to bond if they are within the bond length
                     defined for the two atoms +/- the bondMargin
        bondAngleMargin - the tolerance (in degrees) from the ideal of 180 that defines an acceptable bond
        debugLog - True/False - specifies if a log will be created - not recommended as it generates lots of data
                and slows the program.
        paramsDir - path to the directory holding the forcefield parameter csv files (default ../params)
        '''
        # For time being origin always 0,0,0
        self.origin = np.array([0, 0, 0], dtype=np.float64)
        self.dim = None # The cell dimensions
        self.pbc = [True, True, True]
        self.walls = [False, False, False]
        self.wallAtomType = None
        self.wallRadius = None # (determind from the covalent radius of the wallAtomType)
        # additional distance to add on to the characteristic bond length
        # when checking whether 2 atoms are close enough to bond
        self.bondMargin = bondMargin
        # additional distance to add on to the atom covalent radii when checking if two atoms
        # are close enough to clash
        self.atomMargin = atomMargin
        # convert bondAngle and bondMargin to angstroms
        # Could check values are in degrees and not radians?
        self.bondAngleMargin = math.radians(bondAngleMargin)
        self.targetDensity = 10
        self.targetEndGroups = 100  # number of free endgroups left
        # Dict mapping box key to a list of tuples defining the atoms in that box
        self.box1 = {}
        # Dict mapping key to list of boxes surrounding the keyed box
        self.box3 = {}
        # max atom radius - used to calculate box size
        self.boxSize = None
        self.maxAtomRadius = -1
        self.rCut = 5.0
        # number of boxes in A,B,C axes - used for calculating PBC
        self.numBoxes = [None, None, None]
        self._fragmentLibrary = {}  # fragmentType -> parentFragment
        self._endGroup2LibraryFragment = {}  # endGroup type -> parentFragment
        self.bondTypes = []
        self._bondTable = {}  # endGroup type -> set of endGroups it can bond to
        # dictionary mapping id of the block to the block - can't use a list and indices
        # as we add and remove blocks and need to keep track of them
        self.blocks = collections.OrderedDict()
        self.lastAdded = None # Tracks the id of the last block added to the cell
        self.newBonds = [] # Tracks recently added bonds
        self.mdEngineCls = None
        self._possibleBonds = [] # Holds possible bond after checkMove is run
        # Logging functions
        self.logfile = None
        self.logcsv = None
        self.setupLogging(debugLog=debugLog)
        self._setupAnalyse()
        self._fileCount = 0  # for naming output files
        self._deterministicState = 0 # For adding blocks in a non-random manner (for testing)
        self.rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        if paramsDir is not None:
            if not os.path.isdir(paramsDir):
                msg = "Cannot find parameter directory: {0}".format(paramsDir)
                logger.critical(msg)
                raise RuntimeError(msg)
        else:
            paramsDir = PARAMS_DIR
        self.paramsDir = paramsDir
        # Use the parameters to set the bond lengts in the util module
        xyz_util.setModuleBondLength(os.path.join(paramsDir,'bond_params.csv'))
        self.version = VERSION # Save as attribute so we can query pickle files
        logger.info("AMBUILD version: {0}".format(VERSION))
        self.setMdEngineCls(ab_util.HOOMDVERSION)
        if filePath: # Init from a car file
            self.setStaticBlock(filePath)
        if boxDim:
            self.setBoxSize(boxDim)
        assert self.dim[0] > 0 and self.dim[1] > 0 and self.dim[2] > 0
        return

    def addBlock(self, block, idxBlock=None):
        """
        Add the block and put all atoms in their cells
        """
        if idxBlock is None:
            idxBlock = block.id
        self.blocks[idxBlock] = block

        # Each atom has its own cell - X atoms remain as Nones
        block.atomCell = [None] * block.numAtoms()

        # print "nbox ",self.numBoxA,self.numBoxB,self.numBoxC
        for idxCoord, coord in enumerate(block.iterCoord()):
            key = self._getBox(coord)
            block.atomCell[idxCoord] = key
            try:
                # Add the atom to the cell
                self.box1[key].append((idxBlock, idxCoord))
            except KeyError:
                # Add the cell to the list and then add the atom
                self.box1[key] = [(idxBlock, idxCoord)]
                # Map the cells surrounding this one
                self.box3[key] = self.haloCells(key)
        self.lastAdded = idxBlock
        return idxBlock

    def addBlocks(self, blocks):
        """Try and add a list of blocks to the cell"""
        added = 0
        for block in blocks:
            idxBlock = self.addBlock(block)
            if self.checkMove(idxBlock):
                added += 1
                if self.processBonds() > 0:
                    logger.info("Added bond while adding blocks!")
            else:
                self.delBlock(idxBlock)
        return added

    def addBondType(self, bondType):
        """Allow bonds between the two endGroups specified in the bondType.

        endGroups are defined by the fragmentType they belong to (which is set by the fragmentType argument
        to libraryAddFragment), together with the identifier for that endGroup (which is specified by the first column
        in the .csv file). These are separated by a {0}, so an endGroup identifier is of the form:

        FRAGMENT1{0}ENDGROUP1

        A bond is defined by two endGroups, separated by a {1}, so a bond identifier has the form:

        FRAGMENT1{0}ENDGROUP1{1}FRAGMENT2{0}ENDGROUP2

        Args:
        bondType - a string specifying the two endGroups separated by a "{1}"
        """.format(ENDGROUPSEP,BONDTYPESEP)

        try:
            b1EndGroupType, b2EndGroupType = bondType.split(BONDTYPESEP)
            b1FragmentType = b1EndGroupType.split(ENDGROUPSEP)[0]
            b2FragmentType = b2EndGroupType.split(ENDGROUPSEP)[0]
        except ValueError:
            raise RuntimeError("Error adding BondType {0} - string needs to be of form 'A:a-B:b'".format(bondType))

        # Checks
        # HACK FOR ADDING * TO SHOW BONDED GROUPS
        assert b1FragmentType in self._fragmentLibrary, \
        "No fragment type {0} in fragmentLibrary".format(b1FragmentType)
        if not b1EndGroupType.endswith(ab_endgroup.ENDGROUPBONDED):
            assert b1EndGroupType in self._fragmentLibrary[ b1FragmentType ].endGroupTypes(), \
            "Fragment type {0} has no endGroup type {1}".format(b1FragmentType, b1EndGroupType)

            assert b1EndGroupType in self._endGroup2LibraryFragment.keys(), \
            "No endGroup type {0} in fragmentLibrary!".format(b1EndGroupType)
        assert b2FragmentType in self._fragmentLibrary, "No fragment type {0} in fragmentLibrary".format(b2FragmentType)
        if not b2EndGroupType.endswith(ab_endgroup.ENDGROUPBONDED):
            assert b2EndGroupType in self._fragmentLibrary[ b2FragmentType ].endGroupTypes(), \
            "Fragment type {0} has no endGroup type {1}".format(b2FragmentType, b2EndGroupType)

            assert b2EndGroupType in self._endGroup2LibraryFragment.keys(), \
            "No endGroup type {0} in fragmentLibrary!".format(b2EndGroupType)

        bt = (b1EndGroupType, b2EndGroupType)
        if bt in self.bondTypes:
            raise RuntimeError("Adding an existing bond type: {0}".format(bt))
        self.bondTypes.append(bt)
        # Now recalculate the bond table
        self._updateBondTable()
        return

    def _updateBondTable(self):
        """Recalculate the dicionary of what endGroups can bond to which"""
        self._bondTable = {}
        for bondA, bondB in self.bondTypes:
            if bondA not in self._bondTable:
                self._bondTable[bondA] = set()
            if bondB not in self._bondTable:
                self._bondTable[bondB] = set()
            self._bondTable[bondA].add(bondB)
            self._bondTable[bondB].add(bondA)
        return

    def angle(self, c1, c2, c3):
        """Return the angle in radians c1---c2---c3
        where c are the coordinates in a numpy array
        """
        return xyz_core.angle(c1, c2, c3, dim=self.dim, pbc=self.pbc)

    def attachBlock(self, growEndGroup, staticEndGroup, dihedral=None):
        """
        Position growBlock so it can bond to blockS, using the given _endGroups

        Arguments:

        We take responsibility for adding and removing the growBlock from the cell on
        success or failure
        """

        growBlock = growEndGroup.block()
        staticBlock = staticEndGroup.block()
        idxStaticBlock = staticBlock.id

        staticBlock.positionGrowBlock(staticEndGroup, growEndGroup, dihedral=dihedral)

        # Now add growBlock to the cell so we can check for clashes
        blockId = self.addBlock(growBlock)
        #logger.debug("GOT {0} {1}".format(staticEndGroup, growEndGroup))

        # Check it doesn't clash
        if self.checkMove(blockId) and self.processBonds() > 0:
            logger.debug("attachBlock first checkMove returned True")
            return True
        else:
            logger.debug("attachBlock first checkMove returned False")

        # Only attempt rotation if we're not worried about the dihedral
        # NOTE! - should only bother with the rotation if the cell is relatively crowded - otherwise
        # the multiple rotations are the slowest step
        if not dihedral:
            # Didn't work so try rotating the growBlock about the bond to see if that lets it fit

            # Determine the axis and center to rotate about
            blockEndGroup = growBlock.coord(growEndGroup.endGroupIdx())
            blockS = self.blocks[ idxStaticBlock ]
            blockSEndGroup = blockS.coord(staticEndGroup.endGroupIdx())
            axis = blockSEndGroup - blockEndGroup
            center = blockSEndGroup
            # step = math.pi/18 # 10 degree increments
            step = math.pi / 9  # 20 degree increments
            for angle in ab_util.frange(step, math.pi * 2, step):
                # print "attachBlock rotating as clash: {}".format(angle*util.RADIANS2DEGREES)
                # remove the growBlock from the cell
                self.delBlock(blockId)
                # rotate by increment
                logger.debug("attachBlock rotating by {0} to {1}".format(math.degrees(step),
                                                                         math.degrees(angle)))
                growBlock.rotate(axis, angle, center=center)
                self.addBlock(growBlock) # add it and check
                if self.checkMove(blockId) and self.processBonds() > 0:
                    logger.debug("attachBlock rotate worked")
                    return True
        # remove the growBlock from the cell
        self.delBlock(blockId)
        return False

    def bondAllowed(self, endGroup1, endGroup2):
        """Check if the given bond is permitted from the types of the two fragments
        """
        #    logger.debug( "checking bondAllowed {0} {1}".format( ftype1, ftype2 ) )
        if endGroup1.fragmentType() == "cap" or endGroup2.fragmentType() == "cap":
            assert False, "NEED TO FIX CAPS!"
            return True

        # Get the two fragment types for this bond
        eg1 = endGroup1.type()
        eg2 = endGroup2.type()
        for type1, type2 in self.bondTypes:
            # logger.debug("bondTypes {0}".format((type1, type2)) )
            if (eg1 == type1 and eg2 == type2) or (eg1 == type2 and eg2 == type1):
                # logger.debug( "BOND RETURN TRUE" )
                return True
        # logger.debug( "BOND RETURN FALSE" )
        return False

    def bondBlock(self, bond, selfBond=True):
        """ Bond the second block1 to the first and update the data structures
        """
        logger.debug("cell bondBlock: {0}".format(bond))
        # logger.debug("before bond: {0} - {1}".format( bond.idxBlock1, bond.block1._bondObjects) )
        selfBond = True# HACK
        # We need to remove the block even if we are bonding to self as we need to recalculate the atomCell list
        if bond.isInternalBond() and not selfBond:
            logger.info("bondBlock skipped self-bonded Block")
            return False
        self.delBlock(bond.rootId)
        if not bond.isInternalBond():
            self.delBlock(bond.targetId)
        else:
            logger.info("self-bonded block1: {0}".format(bond))
        bond.engage()
        self.addBlock(bond.endGroup1.block())
        return True

    def bondClash(self, bond, clashDist):
        """Check if any atoms are clashDist from bond.
        """
        if  clashDist > self.boxSize:
            raise RuntimeError("clashDist needs to be less than the boxSize: {0}".format(self.boxSize))
        # Create bond coordinates method
        idxAtom1 = bond.endGroup1.endGroupIdx()
        idxCap1 = bond.endGroup1.capIdx()
        idxAtom2 = bond.endGroup2.endGroupIdx()
        idxCap2 = bond.endGroup2.capIdx()
        p1 = bond.endGroup1.block().coord(idxAtom1)
        p2 = bond.endGroup2.block().coord(idxAtom2)
        idxBlock1 = bond.endGroup1.block().id
        idxBlock2 = bond.endGroup2.block().id

        # Get a list of the cells the bond passes through (excluding endpoints)
        cells = self._intersectedCells(p1, p2)
        # Now build up a list of the cells surrounding the cells on the bond vector
        allcells = set(cells)
        for cell in cells:
            try: surround = self.box3[cell]
            except KeyError: continue
            allcells.update(surround)

        # loop through all the atoms in the cells
        # & check if any are too close to the bond vector
        for cell in list(allcells):
            try:
                atomList = self.box1[cell]
            except KeyError:
                continue
            for (idxBlock3, idxAtom3) in atomList:
                # Dont' check the bond atoms or the cap atoms
                if idxBlock3 == idxBlock1 and idxAtom3 == idxAtom1 or \
                idxBlock3 == idxBlock2 and idxAtom3 == idxAtom2 or \
                idxBlock3 == idxBlock1 and idxAtom3 == idxCap1 or \
                idxBlock3 == idxBlock2 and idxAtom3 == idxCap2:
                    continue
                block3 = self.blocks[idxBlock3]
                if block3.solvent():
                    continue
                p3 = block3.coord(idxAtom3)
                # Distance of a point from a line:
                #  http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
                # http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment

                # First need to find whether the point is inside the segment
                p3p1 = self.vecDiff(p3, p1)
                p3p2 = self.vecDiff(p3, p2)
                p2p1 = self.vecDiff(p2, p1)
                t = np.dot((p3p1), (p2p1)) / np.square(np.linalg.norm(p2p1))
                if 0 < t < 1:
                    dist = np.linalg.norm(np.cross(p3p1, p3p2)) / np.linalg.norm(p2p1)
                    if dist < clashDist:
                        return True
        return False

    def _cat1Paf2(self, bond, fragmentTypes, dt=0.00001, optCycles=1000000):
        """A CAT bonded to two PAF goups"""
        # See if either of the blocks connected is the cat
        cfrag = None
        if bond.rootFragment.catalyst and bond.targetFragment.fragmentType in fragmentTypes:
            cfrag = bond.rootFragment
        elif bond.targetFragment.catalyst and bond.rootFragment.fragmentType in fragmentTypes:
            cfrag = bond.targetFragment
        else: return False # Nothing to do

        # The block may have multiple fragments, but we are only interested in this cat fragment
        # and if this has two bonds made to it and both are to PAF
        # For time being assume only two allowed bonds to cat
        endGroups = cfrag.endGroups()
        if not all([eg.bonded for eg in endGroups]):
            return False
        assert len(endGroups) == 2, "Assumption is CAT only has 2 endGroups!"
        catEG1, catEG2 = endGroups
        logger.info("_cat1Paf2 processing bond %s" % bond)
        logger.debug("_cat1Paf2 %s %s" % (catEG1, catEG2))

        # Now get the two paf endGroups that are bonded to the cat EndGroups
        cat = catEG1.block()
        bond1, bond2 = None, None
        for bond in cat._blockBonds:
            if bond.endGroup1 in [catEG1, catEG2] and bond.targetFragment.fragmentType in fragmentTypes:
                if bond1:
                    assert not bond2
                    bond2 = bond
                else:
                    bond1 = bond
            elif bond.endGroup2 in [catEG1, catEG2] and bond.rootFragment.fragmentType in fragmentTypes:
                if bond1:
                    assert not bond2
                    bond2 = bond
                else:
                    bond1 = bond
        #assert paf1EG and paf2EG,"Could not find two PAF endGroups: {0}".format([str(b) for b in cat._blockBonds])
        if  bond1 is None or bond2 is None:
            logger.info("Could not find two PAF-cat bonds : {0}".format([str(b) for b in cat._blockBonds]))
            return
        #return self._joinPaf(catEG1, paf1EG, paf2EG)
        logger.critical("Got two bonds {0} {1}".format(bond1,bond2))
        return self._joinPaf(fragmentTypes, bond1, bond2, dt=dt, optCycles=optCycles)

    def _cat2Paf2(self, cc_bond, fragmentTypes, dt=0.00001, optCycles=1000000):
        """Function to unbond a Ni-catalyst bonded to two PAF groups

        check if have made a cat:a*-cat:a* bond

        select the two cat endGroups at either end of the bond

        reset the cat endGroup types (cat is now free)

        select the two cat fragments

        break the cat-cat bond

        select the PAF groups attached to each cat blocks.

        select (randomly) one of the cat blocks

        select that PAF

        unbond the PAF unit from the cat

        bond the PAF to the bound cat

        do as before
        """

        # See if this bond is a made between two catalysts which both have PAF bonded to them
        cat1EG = cc_bond.endGroup1
        cat2EG = cc_bond.endGroup2

        # First check if this bond is between 2 cat groups
        if not (cat1EG.bondedCatalyst() and cat2EG.bondedCatalyst()):
            return False # Nothing to do
        logger.info("_cat2Paf2 processing bond: {0}".format(cc_bond))

        # Select the block that is to contain all the fragments - we call it cat1 as we're going to break
        # the bond to the other cat
        cat1 = cat1EG.block()

        # Remove the block from the cell
        self.delBlock(cat1.id)

        # Find the two PAF endgroups that are bonded to the two cat fragments
        paf1EG = None
        paf2EG = None
        for bond in cat1._blockBonds:
            if bond.rootFragment.fragmentType in fragmentTypes and bond.endGroup2.fragment in [cat1EG.fragment, cat2EG.fragment]:
                if bond.targetFragment == cat1EG.fragment:
                    paf1EG = bond.endGroup1
                    bond1 = bond
                else:
                    paf2EG = bond.endGroup1
                    bond2 = bond
            elif bond.targetFragment.fragmentType in fragmentTypes and bond.rootFragment in [cat1EG.fragment,cat2EG.fragment]:
                if bond.rootFragment == cat1EG.fragment:
                    paf1EG = bond.endGroup2
                    bond1 = bond
                else:
                    paf2EG = bond.endGroup2
                    bond2 = bond

        assert paf1EG and paf2EG,"Could not find PAF endGroups"

        # Break the bond between the two cat blocks
        cat2 = cat1.deleteBond(cc_bond, root=cat1EG.fragment)

        # The initial assumption was that breaking the cat-cat bond would create two blocks, but if the PAF bonded to the cat
        # blocks are part of a chain, then this won't happen
        # Now pick a PAF-CAT bond and break it to create a lone PAF block that will be joined to the other cat
        # If the bonds were circular we create the separate cat2 block at this stage. Otherwise we are just splitting
        # off the PAF block
        if cat2:
            # We ignore paf2 as it's only the container for the paf fragment, which is joined to cat1 below using the endGroup
            paf2 = cat2.deleteBond(bond2,root=cat2EG.fragment)
        else:
            logger.info("Breaking bond cat-cat did not create two blocks so PAF blocks are joined in a chain")
            cat2 = cat1.deleteBond(bond2,root=paf2EG.fragment)

        # Put the separated off catalyst back in the cell
        self.addBlock(cat2)

        # Now bond the PAF to the first catalyst
        cp_bond1 = ab_bond.Bond(cat1EG,paf2EG)
        cat1.bondBlock(cp_bond1)

        # Put the newly bonded cat block back in the cell
        self.addBlock(cat1)

        # Run optimisation to move CAT away
        logger.info("_cat2Paf2 Optimisation")
        #self.dump()
        self.optimiseGeometry(rigidBody=True, quiet=False, dt=dt, optCycles=optCycles)

        # Now dealing with a CAT bonded to two PAF groups
        # Need to select the other cat-paf bond
        cp_bond2 = None
        for bond in cat1._blockBonds:
            if bond.rootFragment.catalyst and bond.endGroup2 == paf1EG:
                cp_bond2 = bond
                break
            if bond.targetFragment.catalyst and bond.endGroup1 == paf1EG:
                cp_bond2 = bond
                break
        assert cp_bond2,"Could not find second cat/PAF bond"
        self._joinPaf(fragmentTypes, cp_bond1, cp_bond2, dt=dt, optCycles=optCycles)
        return True

    def _joinPaf(self, fragmentTypes, bond1, bond2, dt=0.00001, optCycles=1000000):
        """Given a cat bonded to two PAF groups, break the PAF bonds and form a PAF-PAF bond"""

        logger.info("Entering _joinPaf: {0} {1}".format(bond1, bond2))
        # Get the PAF endgroups and the block
        catBlock = bond1.endGroup1.block() # should all be part of the same block
        catEG, paf1EG, paf2EG = None, None, None
        if bond1.rootFragment.fragmentType in fragmentTypes and bond1.targetFragment.catalyst:
            paf1EG = bond1.endGroup1
            catEG = bond1.endGroup2
        elif bond1.rootFragment.catalyst and bond1.targetFragment.fragmentType in fragmentTypes:
            paf1EG = bond1.endGroup2
            catEG = bond1.endGroup1

        if bond2.rootFragment.fragmentType in fragmentTypes and bond2.targetFragment.catalyst:
            paf2EG = bond2.endGroup1
            #catEG = bond2.endGroup2
        elif bond2.rootFragment.catalyst and bond2.targetFragment.fragmentType in fragmentTypes:
            paf2EG = bond2.endGroup2
            #catEG = bond2.endGroup1

        assert catEG and paf1EG and paf2EG, "Could not find endGroups: {0} {1}".format(bond1, bond2)

        self.delBlock(catBlock.id)
        paf1 = catBlock.deleteBond(bond1,root=catEG.fragment)
        paf2 = catBlock.deleteBond(bond2,root=catEG.fragment)

        # Add the unbonded blocks back to the cell
        self.addBlock(catBlock)
        # Possibly need to think more about what happens when blocks are parts of chains
        if paf1:
            self.addBlock(paf1)
        if paf2:
            self.addBlock(paf2)

        # We now need to bond the two PAF groups
        assert paf1EG.free() and paf2EG.free(),"PAF endgroups aren't free!"
        bond = ab_bond.Bond(paf1EG, paf2EG)
        self.bondBlock(bond)
        logger.info("_joinPaf Optimisation")
        #self.dump()
        self.optimiseGeometry(rigidBody=True, dt=dt, optCycles=optCycles, max_tries=1, retries_on_error=0, dump=True)
        self.clearUnbonded()
        return True

    def clearUnbonded(self):
        "Run after we have optimised to unset the unBonded flag and add capatoms back into the cell"
        for block in self.blocks.values():
            for frag in block.fragments:
                frag.clearUnbonded()
        return

    def cat1Paf2(self, fragmentTypes, dt=0.00001, optCycles=1000000):
        """Function to unbond a Ni-catalyst bonded to two PAF groups"""
        assert type(fragmentTypes) is list and len(fragmentTypes) > 0 and all([type(f) is str for f in fragmentTypes]),\
            "Need a list of fragmentTypes"
        if not len(self.newBonds):
            return False
        #logger.info("ca1tPaf2 got new bonds %s" % [str(b) for b in self.newBonds ])
        return any([self._cat1Paf2(b, fragmentTypes, dt=dt, optCycles=optCycles) for b in self.newBonds])

    def cat2Paf2(self, fragmentTypes, dt=0.00001, optCycles=1000000):
        """Function to unbond a Ni-catalyst bonded to two PAF groups"""
        assert type(fragmentTypes) is list and len(fragmentTypes) > 0 and all([type(f) is str for f in fragmentTypes]),"Need a list of fragmentTypes"
        if not len(self.newBonds): return False
        #logger.info("cat2Paf2 got new bonds %s" % [str(b) for b in self.newBonds ])
        return any([self._cat2Paf2(b, fragmentTypes, dt=dt, optCycles=optCycles) for b in self.newBonds])

    def canBond(self,
                 staticBlock,
                 idxStaticAtom,
                 addBlock,
                 idxAddAtom,
                 distance,
                 bondMargin,
                 bondAngleMargin
                ):
        # The check should have been made before this is called on whether the two atoms are endGroup
        # Check length
        bond_length = xyz_util.bondLength(addBlock.symbol(idxAddAtom), staticBlock.symbol(idxStaticAtom))
        if bond_length < 0:
            raise RuntimeError("Missing bond distance for: {0}-{1}".format(addBlock.symbol(idxAddAtom),
                                                                            staticBlock.symbol(idxStaticAtom)))

        # See if the distance between them is acceptable
        # print "CHECKING BOND ATOMS ",bond_length,self.distance( addCoord, staticCoord )
        if  not (max(0.1, bond_length - bondMargin) < distance < bond_length + bondMargin):
            logger.debug("Cannot bond due to distance: {0}".format(distance))
            return False
        # Now loop over all endGroups seeing if any of the angles are satisfied
        for staticEndGroup in staticBlock.atomEndGroups(idxStaticAtom):
            for addEndGroup in addBlock.atomEndGroups(idxAddAtom):
                # EndGroups in the same fragment can never bond
                if staticEndGroup.fragment == addEndGroup.fragment:
                    break
                assert staticEndGroup.free() and addEndGroup.free()
                # We need to check that we've not already added these endGroups as possible bonds
                if self._endGroupsInPossibleBonds([ staticEndGroup, addEndGroup ]):
                    continue
                if not self.bondAllowed(staticEndGroup, addEndGroup):
                    logger.debug("Bond disallowed by bonding rules: {0} : {1}".format(staticEndGroup,
                                                                                      addEndGroup))
                    continue

                # print "Possible bond for {0} {1} {2} dist {3}".format( idxAddAtom,
                #                                                       idxStaticBlock,
                #                                                       idxStaticAtom,
                #                                                       self.distance( addCoord, staticCoord ) )
                addCoord = addBlock.coord(idxAddAtom)
                staticCoord = staticBlock.coord(idxStaticAtom)
                addCapAtom = addBlock.coord(addEndGroup.capIdx())
                angle1 = self.angle(addCapAtom, addCoord, staticCoord)
                staticCapAtom = staticBlock.coord(staticEndGroup.capIdx())
                angle2 = self.angle(staticCapAtom, staticCoord, addCoord)

                # print "CHECKING ANGLE BETWEEN: {0} | {1} | {2}".format( angleAtom, addCoord, staticCoord )
                # print "{} < {} < {}".format( (self.bondAngle-bondAngleMargin) * util.RADIANS2DEGREES,
                #                             angle * util.RADIANS2DEGREES,
                #                             (self.bondAngle+bondAngleMargin) * util.RADIANS2DEGREES  )

                # Check if atoms are in line (zero degrees) within margin
                if not ((0.0 - bondAngleMargin < angle1 < 0.0 + bondAngleMargin and \
                           0.0 - bondAngleMargin < angle2 < 0.0 + bondAngleMargin)):
                    logger.debug("Cannot bond due to angles: {0} : {1}".format(math.degrees(angle1),
                                                                               math.degrees(angle2)))
                    continue

                logger.debug("Acceptable bond with angles: {0} : {1} | distance {2}".format(math.degrees(angle1),
                                                                                            math.degrees(angle2),
                                                                                            distance))

                # Create bond object and set the parameters
                bond = ab_bond.Bond(staticEndGroup, addEndGroup)
                self._possibleBonds.append(bond)
                logger.debug("canBond returning True with bonds: {0}".format([str(b) for b in self._possibleBonds]))
                return True

        logger.debug("canBond returning False")
        return False

    def capBlocks(self, fragmentType=None, filename=None):

        # Create the cap block
        capBlock = ab_block.Block(filePath=filename, fragmentType='cap')

        # The endgroup is always the first only endGroup
        capEndGroup = capBlock.freeEndGroups()[0]

        # Get a list of blocks - need to do it this way or else we are changing the list of blocks while
        # we iterate through them
        cblocks = self.blocks.copy()
        for blockId, block in cblocks.items():

            if block.numFreeEndGroups() == 0:
                logger.info("capBlocks block {0} already has no free endGroups".format(blockId))
                continue
            else:
                logger.info("capBlocks capping {0} endGroups of block {1}".format(block.numFreeEndGroups(),
                                                                                        blockId))
            # If there are no free this won't loop
            for endGroup in block.freeEndGroups(endGroupTypes=[fragmentType]):
                cblock = capBlock.copy()
                block.positionGrowBlock(endGroup, capEndGroup)
                idxBlock = self.addBlock(cblock)
                # Test for Clashes with other blocks
                if self.checkMove(idxBlock) and self.processBonds() > 0:
                    logger.info("Capped block {0} endGroup {1}".format(blockId, endGroup))
                else:
                    logger.critical("Failed to cap block {0} endGroup {1}".format(blockId, endGroup))
                    # Unsuccessful so remove the block from cell
                    self.delBlock(idxBlock)
        return

    def checkMove(self, idxAddBlock):
        clashing = self._checkMove(idxAddBlock)
        if clashing > 0:
            return False
        return True

    def _checkMove(self, idxAddBlock):
        """
        get pairs of close atoms

        """
        close, wallClashes = self.closeAtoms(idxAddBlock, walls=self.walls) # Get a list of the close atoms
        if wallClashes:
            logger.debug("_checkMove got clash with wall")
            return 1
        if len(close) is 0:
            logger.debug("_checkMove no close contacts")
            return 0
        addBlock = self.blocks[ idxAddBlock ]
        clashAtoms = []
        self._possibleBonds = []
        for (idxAddAtom, staticBlock, idxStaticAtom, distance) in close:
#             logger.debug("CHECKING  ATOMS {}:{} {} -> {}:{} {}= {}".format( idxAddBlock,
#                                                               idxAddAtom,
#                                                               addBlock.symbol(idxAddAtom),
#                                                               staticBlock.id,
#                                                               idxStaticAtom,
#                                                               staticBlock.symbol( idxStaticAtom ),
#                                                               distance ))
            if not (addBlock.isEndGroup(idxAddAtom) and \
                     staticBlock.isEndGroup(idxStaticAtom) and \
                     self.canBond(staticBlock,
                                   idxStaticAtom,
                                   addBlock,
                                   idxAddAtom,
                                   distance,
                                   self.bondMargin,
                                   self.bondAngleMargin)
                    ):
                # No bond so just check if the two atoms are close enough for a clash
                d = addBlock.radius(idxAddAtom) + staticBlock.radius(idxStaticAtom) + self.atomMargin
                if distance <= d:
                    #logger.debug("CLASH {}{}->{}{} = {} > {}".format( addBlock.type(idxAddAtom), addBlock.coord(idxAddAtom),\
                    #                                           staticBlock.type(idxStaticAtom), staticBlock.coord(idxStaticAtom), \
                    #                                       d, distance  ))
                    clashAtoms.append((staticBlock, idxStaticAtom, addBlock, idxAddAtom))

        # Now have list of possible bonds and clashes
        # Nothing so return True
        if not len(self._possibleBonds) and not len(clashAtoms):
            logger.debug("NO BONDS AND NO ATOMS")
            return 0
        # no bonds but a clash - return False
        if not len(self._possibleBonds) and len(clashAtoms):
            logger.debug("No bonds and clashing atoms {0}".format(clashAtoms))
            return len(clashAtoms)
        s = ""
        for b in self._possibleBonds:
            s += str(b) + " | "
        logger.debug("Bonds {0} and clashing atoms {1}".format(s, clashAtoms))

        for bond in self._possibleBonds:
            b1Block = bond.endGroup1.block()
            b2Block = bond.endGroup2.block()
            # We need to remove any clashes with the cap atoms - the added block isn't bonded
            # so the cap atoms aren't excluded
            b1Cap = bond.endGroup1.capIdx()
            b2Cap = bond.endGroup2.capIdx()
            b1EndGroup = bond.endGroup1.endGroupIdx()
            b2EndGroup = bond.endGroup2.endGroupIdx()
            # Also need to remove any clashes of the endGroups with atoms directly bonded to the
            # opposite endGroup
            b1BondAtoms = b1Block.atomBonded1(b1EndGroup)
            b1BondAtoms.remove(b1Cap)
            b2BondAtoms = b2Block.atomBonded1(b2EndGroup)
            b2BondAtoms.remove(b2Cap)
            toGo = []  # Need to remember indices as we can't remove from a list while we cycle through it
            for i, (cellBlock, idxCellAtom, addBlock, idxAddAtom) in enumerate(clashAtoms):
                # logger.debug("CHECKING {0} {1} {2} {3}".format( idxAddBlock, idxAddAtom, idxStaticBlock, idxStaticAtom ) )
                # Remove any clashes with the cap atoms
                if (b1Block == cellBlock and idxCellAtom == b1Cap) or \
                   (b2Block == addBlock  and idxAddAtom == b2Cap):
                    # logger.debug("REMOVING {0} {1} {2} {3}".format( idxAddBlock, idxAddAtom, idxStaticBlock, idxStaticAtom ) )
                    toGo.append(i)
                    continue
                # remove any clashes with directly bonded atoms
                if (b1Block == cellBlock and idxCellAtom == b1EndGroup and idxAddAtom  in b2BondAtoms) or \
                   (b2Block == b2Block   and idxAddAtom == b2EndGroup and idxCellAtom in b1BondAtoms):
                    # logger.info( "REMOVING BOND ATOMS FROM CLASH TEST" )
                    toGo.append(i)
            # End of loop so now remove the items
            for i, idx in enumerate(toGo):
                # Need to subtract i as otherwise the indexing is out
                clashAtoms.pop(idx - i)
        # If there are any clashing atoms remaining this move failed
        if len(clashAtoms):
            logger.debug("Got clash atoms {0}".format(clashAtoms))
            return len(clashAtoms)
        # Got bonds and no clashes
        logger.debug("Checkmove no clashes")
        return 0

    def clear(self):
        """Empty the cell of blocks and reset any data structures"""
        # Remove all blocks from their cells
        self.box1 = {}
        self.box3 = {}
        self.blocks.clear()  # Delete block list
        return

    def closeAtoms(self, idxBlock1, walls=None):
        """
        Find all atoms that are close to the atoms in the given block.

        This version is optimised to use numpy to do the PBC distance calculations, so builds up a vector of points and then
        calculates the distances in one go with the distanceV function

        Args:
        idxBlock1: index of the block in self.blocks

        Returns:
        A tuple with:
        * a list of tuples: (thisAtomIndex, otherBlockIndex, otherAtomIndex, distance) or None if there we no close contacts
          or a clash with the wall
        * True/False if there is a clash with the wall
        """
        # Build up list of the coordinates and an initial contacts array
        c1 = []
        c2 = []
        allContacts = []
        count = 0
        block1 = self.blocks[ idxBlock1 ]
        wallCheck = False if walls is None else any(walls)
        for idxAtom1, coord1 in enumerate(block1.iterCoord()):

            # First check if the atom is close to a wall
            if wallCheck:
                radius = block1.radius(idxAtom1) + self.wallRadius + self.atomMargin
                if walls[0] and (coord1[0] < radius or coord1[0] > self.dim[0] - radius) or \
                   walls[1] and (coord1[1] < radius or coord1[1] > self.dim[1] - radius) or \
                   walls[2] and (coord1[2] < radius or coord1[2] > self.dim[2] - radius):
                    # Got a clash with a wall, so we can stop all other checks
                    return None, True

            # Get the box this atom is in
            key = block1.atomCell[ idxAtom1 ]
            # Skip checking dummy atoms
            if key is None:
                continue
            # For each box loop through all its atoms chekcking for clashes
            surrounding = self.box3[key]
            for sbox in surrounding:
                # For each box, get the list of the atoms as (block,coord1) tuples
                # Check if we have a box with anything in it
                try:
                    alist = self.box1[ sbox ]
                except KeyError:
                    continue
                for (idxBlock2, idxAtom2) in alist:
                    # Check we are not checking ourself - need to check block index too!
                    if idxBlock1 == idxBlock2:
                        continue
                    block2 = self.blocks[idxBlock2]
                    c1.append(coord1)
                    c2.append(block2.coord(idxAtom2))
                    allContacts.append((idxAtom1, block2, idxAtom2))
                    count += 1
        # Have list so now calculate distances
        if count == 0:
            return [], False
        # Calculate array of distances for all coordinates
        distances = self.distance(np.array(c1), np.array(c2))
        # prune contacts array according to distances
        return [ (c[0], c[1], c[2], distances[i]) for i, c in enumerate(allContacts) if distances[i] < self.boxSize ], False

    def cellEndGroupPair(self, cellEndGroups=None):
        """Return two free endGroups from two different blocks in the cell"""
        # Get a list of available free endGroups in the cell
        endGroupTypes2Block = self.endGroupTypes2Block()
        if len(endGroupTypes2Block.keys()) == 0:
            logger.critical("cellEndGroupPair: No available endGroups for: {0}".format(cellEndGroups))
            return None, None
        # If the user supplied a list of cellEndGroups we use this to determine what can bond -
        # other wise we use all available endGroups
        allTypes = endGroupTypes2Block.keys()  # just for informing the user
        if cellEndGroups is not None:
            if isinstance(cellEndGroups, str):
                cellEndGroups = [ cellEndGroups ]
            cellEndGroups = set(cellEndGroups)
        else:
            cellEndGroups = set(endGroupTypes2Block.keys())
        # We create a dictionary mapping which cell endGroups can bond to which
        # This is basically a truncated _bondTable with any endGroups removed that aren't present
        cell2cell = copy.copy(self._bondTable)
        for eg in self._bondTable:
            if eg not in cellEndGroups:
                # First remove all references in the keys
                del cell2cell[ eg ]
            else:
                # The key is valid so remove any missing endGroups from the values
                ceg = cellEndGroups.intersection(cell2cell[ eg ])
                # Need to check there are some other types that this can bond to, but also there is at least
                # one other block if both are of the same type
                # NB the strange next/iter thing is the only way to get an item from the set
                if len(ceg) == 0 or \
                   (len(ceg) == 1 and eg == next(iter(ceg)) and len(endGroupTypes2Block[ eg ]) < 2):
                    # No valid blocks that could bond with this endGroup so delete the entire key
                    del cell2cell[ eg ]
                else:
                    cell2cell[ eg ] = ceg
        if len(cell2cell.keys()) == 0:
            logger.critical("cellEndGroupPair: No endGroups of types {0} are available to bond from {1}".format(cellEndGroups, allTypes))
            return None, None

        # At this point we assume we can definitely bond at least 2 blocks
        logger.debug("cellEndGroupPair got cell/library endGroups: {0}".format(cell2cell))
        # Select a random block/endGroup from the list
        eg1Type = _random.choice(list(cell2cell.keys()))
        block1 = _random.choice(list(endGroupTypes2Block[ eg1Type ]))
        endGroup1 = block1.selectEndGroup(endGroupTypes=[eg1Type])
        # Pick a random endGroup type that can bond to this
        eg2Type = _random.sample(cell2cell[ eg1Type ], 1)[0]
        # Select a random block/endGroup of that type
        # (REM: need to remove the first block from the list of possibles hence the difference thing
        # XXX Also need to convert to list as sets don't support random.choice
        # Using choice and list as sample was giving error "ValueError: sample larger than population"
        # block2 = random.sample( endGroupTypes2Block[ eg2Type ].difference(set([block1])), 1 )[0]
        try:
            block2 = _random.choice(list(endGroupTypes2Block[ eg2Type ].difference(set([block1]))))
            # This will trigger an IndexError if there isn't a free block of the given type
        except IndexError:
            logger.critical("cellEndGroupPair: No 2nd block available for cellEndGroups".format(cellEndGroups))
            return None, None
        endGroup2 = block2.selectEndGroup(endGroupTypes=[eg2Type])

        logger.debug("cellEndGroupPair returning: {0} {1}".format(endGroup1.type(), endGroup2.type()))
        return endGroup1, endGroup2

    def cellData(self, rigidBody=True, periodic=True, center=False, fragmentType=None):
        RIGIDPARTICLES = rigidBody and ab_util.HOOMDVERSION and ab_util.HOOMDVERSION[0] > 1
        # Object to hold the cell data
        d = ab_celldata.CellData()
        d.cell = self.dim
        if fragmentType is not None:
            if RIGIDPARTICLES:
                assert False,"Need to update for HOOMD2"
            # Only returning data for one type of fragment
            assert fragmentType in self.fragmentTypes(), "FragmentType {0} not in cell!".format(fragmentType)
            atomIdx = 0
            for b in self.blocks.values():
                coords, symbols, bonds = b.dataByFragment(fragmentType)
                d.coords += coords
                d.symbols += symbols
                d.bonds += [ (b1 + atomIdx, b2 + atomIdx) for (b1, b2) in bonds ]
                atomIdx += len(coords)
            return d
        atomIdx = 0
        bodyIdx = 0
        for block in self.blocks.values():
            # Do bonds first as counting starts from atomIdx and goes up through the blocks
            if not rigidBody:
                # add all bonds, angles and dihederals throughout the whole block
                d.bonds += [(a1 + atomIdx, a2 + atomIdx) for a1, a2 in block.bonds()]
                d.bondLabels += ["{0}-{1}".format(block.type(a1), block.type(a2)) for a1, a2 in block.bonds()]
                angles, propers, impropers = block.anglesAndDihedrals()
                d.angles += [(a1 + atomIdx, a2 + atomIdx, a3 + atomIdx) for a1, a2, a3 in angles]
                d.angleLabels += ["{0}-{1}-{2}".format(block.type(a1), block.type(a2), block.type(a3)) for a1, a2, a3 in angles]
                d.propers += [(a1 + atomIdx, a2 + atomIdx, a3 + atomIdx, a4 + atomIdx) \
                              for a1, a2, a3, a4 in propers]
                d.properLabels += ["{0}-{1}-{2}-{3}".format(block.type(a1),
                                                            block.type(a2),
                                                            block.type(a3),
                                                            block.type(a4)
                                                            ) for a1, a2, a3, a4 in propers]
                d.impropers += [(a1 + atomIdx, a2 + atomIdx, a3 + atomIdx, a4 + atomIdx) \
                                for a1, a2, a3, a4 in impropers]
                d.improperLabels += ["{0}-{1}-{2}-{3}".format(block.type(a1),
                                                              block.type(a2),
                                                              block.type(a3),
                                                              block.type(a4)
                                                              ) for a1, a2, a3, a4 in impropers]
            else:
                # Just add the bonds between blocks. Also add angles for all atoms connected to the bonds
                # we do this so that we can exclude them from VdW interactions in MD codes
                for b1, b2 in block.blockBonds():
                    # The bonds themselves
                    d.bonds.append((b1 + atomIdx, b2 + atomIdx))
                    d.bondLabels.append("{0}-{1}".format(block.type(b1), block.type(b2)))
                    _angles = set()
                    # Atoms connected to the endGroup that we need to specify as connected so we add as angles
                    for batom in block.atomBonded1(b1):
                        # The opposite endGroup is included in the list bonded to an endGroup so skip
                        if batom == b2:
                            continue
                        _angles.add((batom, b1, b2))
                    for batom in block.atomBonded1(b2):
                        # The opposite endGroup is included in the list bonded to an endGroup so skip
                        if batom == b1:
                            continue
                        _angles.add((b1, b2, batom))
                    # Add to overall lists
                    for a1, a2, a3 in _angles:
                        d.angles.append((a1 + atomIdx, a2 + atomIdx, a3 + atomIdx))
                        l = "{0}-{1}-{2}".format(block.type(a1),
                                                  block.type(a2),
                                                  block.type(a3))
                        d.angleLabels.append(l)
                    # Dihedrals
                    for dindices in block.dihedrals(b1, b2):
                        dihedral = (dindices[0] + atomIdx,
                                     dindices[1] + atomIdx,
                                     dindices[2] + atomIdx,
                                     dindices[3] + atomIdx)
                        d.propers.append(dihedral)
                        dlabel = "{0}-{1}-{2}-{3}".format(block.type(dindices[0]),
                                                           block.type(dindices[1]),
                                                           block.type(dindices[2]),
                                                           block.type(dindices[3]) )
                        d.properLabels.append(dlabel)
            # Now loop through fragments and coordinates
            for frag in block.fragments:
                for body in frag.bodies():
                    d.natoms += body.natoms
                    atomIdx += body.natoms
                    if RIGIDPARTICLES:
                        d.rigidParticles.append(self.rigidParticleMgr.createParticle(body))
                    else:
                        coords = body.coords
                        coords, images = xyz_core.wrapCoord3(coords, dim=self.dim, center=center)
                        d.coords += list(coords)
                        d.images += list(images)
                        d.atomTypes += body.atomTypes
                        d.bodies += [bodyIdx] * body.natoms
                        d.charges += body.charges
                        d.diameters += body.diameters
                        d.masked += body.masked
                        d.masses += list(body.masses)
                        d.static += body.static
                        d.symbols += body.symbols
                bodyIdx += 1
        if RIGIDPARTICLES:
            self.rigidParticleMgr.checkConfigStrClashes(d.atomTypes)
            d.rigidParticleMgr = self.rigidParticleMgr
            
        return d

    def delBlock(self, blockId):
        """
        Remove the block with the given index from the cell
        """
        block = self.blocks[ blockId ]
        # Remove each atom from the list
        keys = []
        for iatom, key in enumerate(block.atomCell):
            # Skip dummy atoms
            if key == None:
                continue
            keys.append(key)
            self.box1[key].remove((blockId, iatom))
        # Now remove any empty keys and corresponding surrounding boxes
        # Might think about keeping the surrounding boxes as they could be use by other atoms?
        # print "now checking keys"
        for key in keys:
            if key in self.box1 and len(self.box1[key]) == 0:
                del self.box1[key]
                del self.box3[key]
        del self.blocks[blockId]
        return

    def deleteBlocks(self, fragmentTypes=None, indices=None, maxFrags=1, numBlocks=0):
        """Remove blocks from the cell.

        Arguments:
        fragmentTypes: the fragmentType, or list of fragmentTypes of the blocks to remove
        indices: list of indices of the blocks in the overall list of blocks in the cell to be deleted
        maxFrags: only blocks containing <= this number of fragments will be removed (default = 1)
        numBlocks: optional - the number of blocks to remove, otherwise all blocks of the specified fragmentTypes will be removed

        Returns:
        A list of the blocks that were removed - suitable for re-adding with restoreBlocks
        """
        if fragmentTypes is not None:
            return self._deleteBlocksType(fragmentTypes, numBlocks=numBlocks, maxFrags=maxFrags)
        elif indices is not None:
            return self._deleteBlocksIndices(indices)

    def _deleteBlocksIndices(self, indices):
        """Delete the index-th(s) block from the cell

        Arguments:
        indices: list of indices of the blocks to be deleted
        """
        if type(indices) is float:
            indices = [indices]
        elif type(indices) is list:
            pass
        else:
            raise RuntimeError("indices needs to be a int or list of int")
        for i, idx in enumerate(indices):
            if not (type(idx) is int and (0 <= idx < len(self.blocks))):
                raise RuntimeError("Bad value for index {0} : {1}".format(i,idx))
        toRemove = [ (blockId, block) for i, (blockId, block) in enumerate(self.blocks.items())  if i in indices ]
        removed = []
        for blockId, block in toRemove:
            self.delBlock(blockId)
            removed.append(block)
        count = len(toRemove)
        logger.info("deleteBlocksIndices deleted {0} blocks. Cell now contains {1} blocks.".format(count,len(self.blocks)))
        return removed

    def _deleteBlocksType(self, fragmentTypes, numBlocks=0, maxFrags=1):
        """Remove numBlocks of fragmentType from the cell.

        Arguments:
        fragmentType: the fragmentType, or list of fragmentTypes of the blocks to remove
        numBlocks: optional - the number of blocks to remove, otherwise all blocks of the specified types will be removed
        multiple: if True remove blocks that contain > 1 fragment, else only single-fragment blocks
        """
        if type(fragmentTypes) is str:
            fragmentTypes = [fragmentTypes]
        fragmentTypes = set(fragmentTypes)
        allBlocks = []
        # First get a list of all the blocks that contain fragments that match the fragment type
        for blockId, block in self.blocks.items():
            # Check if is multiple or not
            if maxFrags > 0 and len(block.fragments) > maxFrags: continue
            # Add this block if it contains any of the fragmentTypes
            if set(block.fragmentTypes()).intersection(fragmentTypes):
                allBlocks.append((blockId,block))
        # We have a list of valid blocks, we now randomly remove them - randomly in case we aren't deleting all of them
        if numBlocks > 0:
            toRemove = min(numBlocks, len(allBlocks))
        else:
            toRemove = len(allBlocks)
        removed = []
        if toRemove > 0:
            for _ in range(toRemove):
                blockId, block = _random.choice(allBlocks)
                self.delBlock(blockId)
                allBlocks.remove((blockId,block))
                removed.append(block)
            logger.info("Delete removed {0} blocks. Cell now contains {1} blocks".format(len(removed), len(self.blocks)))
        else:
            logger.info("Could not remove any blocks of type(s) {0}".format(fragmentTypes))
        return removed

    def deleteBondType(self, bondType):
        """Remove the given bondType from the list of acceptable bonds"""
        b1EndGroupType, b2EndGroupType = bondType.split(BONDTYPESEP)
        bt = (b1EndGroupType, b2EndGroupType)
        assert bt in self.bondTypes,"{0} is not in the list of acceptable bonds {1}".format(bondType,self.bondTypes)
        self.bondTypes.remove(bt)
        self._updateBondTable()
        return

    def deleteFragment(self, frag, block=None):
        logger.info("Deleting fragment: {0}".format(frag.fragmentType))
        if block is None:
            # Find the fekker
            for b in self.blocks.values():
                for f in b.fragments:
                    if f is frag:
                        block = b
                        break

        assert block,"Need to know which block to remove the fragment from!"

        # deleteFragment will return a list of blocks to be added back to the cell
        # The original block may have been deleted during the fragment deletion processs so we always need to
        # remove it and add the blocks returned by deleteFragment
        self.delBlock(block.id)
        blocks = block.deleteFragment(frag)
        if len(blocks):
            logger.info("Deleting fragment resulted in {0} blocks.".format(len(blocks)))
            for b in blocks: self.addBlock(b)
        else:
            logger.info("Deleting fragment deleted a block")
            # Nothing to do here as we've already removed the block from the cell
        return

    def density(self):
        """The density of the cell"""
        d = (sum([ b.blockMass() for b in self.blocks.values() ]) / (self.dim[0] * self.dim[1] * self.dim[2]))
        return d * (10 / 6.022)

    def dihedral(self, p1, p2, p3, p4):
        return xyz_core.dihedral(p1, p2, p3, p4, dim=self.dim, pbc=self.pbc)

    def distance(self, v1, v2):
        """Distance with numpy taking PBC into account
        This works either with 2 points or a vector of any number of points
        Adapted from: http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co
        Changed so that it can cope with distances across more than one cell
        """
        return xyz_core.distance(v1, v2, dim=self.dim, pbc=self.pbc)

    def dump(self, prefix="step", addCount=True):
        """Write out our current state"""

        if addCount:
            self._fileCount += 1
            prefix = prefix + "_{0}".format(self._fileCount)
        # self.writeXyz(prefix+".xyz",data=data, periodic=False)
        # self.writeXyz(prefix+"_P.xyz",data=data, periodic=True)
        # self.writeCar(prefix+"_P.car",data=data,periodic=True)
        # self.writeCml(prefix+"_PV.cml", data=data, allBonds=True, periodic=True, pruneBonds=True)
        # self.writeCml(prefix+".cml", data=data, allBonds=True, periodic=False, pruneBonds=False)
        pklFile = self.writePickle(prefix)
        return pklFile

    def endGroupConfig(self, fragmentType):
        """Return the name of the last pkl file and the endGroupConfig (number of bonded endGroups) for
        blocks containing the fragmentType"""
        s = ""
        for b in self.blocks.values():
            c = b.endGroupConfig(fragmentType)
            if c is not None:
                if len(s) == 0:
                    s = c
                else:
                    s += "|" + c
        return[str(self._fileCount), s]

    def _endGroupsInPossibleBonds(self, endGroups):
        """Check if any of the endGroups are already in the list of possible bonds"""
        eg = set()
        for b in self._possibleBonds:
            eg.update([ b.endGroup1, b.endGroup2 ])
        return bool(eg.intersection(frozenset(endGroups)))

    def endGroupTypes2Block(self):
        """Return a dictionary mapping free endGroup types to a list of the blocks

        We don't check if any are available just return an empty dictionary if not
        """
        endGroupTypes2Block = {}
        for block in self.blocks.values():
            # numFreeEndGroups += block.numFreeEndGroups()
            # Get a list of the free endGroup types in the block
            for endGroupType in block.freeEndGroupTypes():
                if endGroupType not in endGroupTypes2Block:
                    endGroupTypes2Block[ endGroupType ] = set()
                endGroupTypes2Block[ endGroupType ].add(block)
        return endGroupTypes2Block

    def fragMaxEnergy(self,
                  rigidBody=True,
                  doDihedral=False,
                  doImproper=False,
                  xmlFilename="hoomdCalc.xml",
                  **kw):
        if not self.mdEngineCls:
            raise RuntimeError("No mdEngine defined - cannot run.")
        mdEngine = self.mdEngineCls(self.paramsDir)
        data = self.cellData(periodic=True, center=True, rigidBody=rigidBody)
        if 'rCut' in kw:
            self.rCut = kw['rCut']
        else:
            self.rCut = self.MDENGINE.rCut
        logger.info("Running fragMaxEnergy")
        maxe, idxBlock, idxFragment = self.mdEngine.fragMaxEnergy(data,
                                                                  xmlFilename,
                                                                  rigidBody=rigidBody,
                                                                  doDihedral=doDihedral,
                                                                  doImproper=doImproper,
                                                                  **kw)
        return

    def fragmentTypes(self):
        ft = collections.defaultdict(list)
        for b in self.blocks.values():
            d = b.fragmentTypeDict()
            for k, v in d.items():
                try:
                    ft[ k ] += v
                except:
                    ft[ k ] = v
        return ft

    def fragmentTypeFromEndGroupType(self, endGroupType):
        return endGroupType.split(self.ENDGROUPSEP)[0]

    def _getBox(self, coord):
        """Return the box that the coord is in under periodic boundaries"""
        return xyz_util.getCell(coord, self.boxSize, dim=self.dim, pbc=self.pbc)

    def getLibraryBlock(self, fragmentType=None, random=True):
        """Return an initBlock"""
        if fragmentType is None:
            # Need to determine the fragmentType from the endGroupType
            if random:
                fragmentType = _random.choice(list(self._fragmentLibrary.keys()))
            else:
                i = self._deterministicState % len(self._fragmentLibrary.keys())
                fragmentType = sorted(self._fragmentLibrary.keys())[i]
                self._deterministicState += 1
        # sanity check
        if not fragmentType in self._fragmentLibrary:
            raise RuntimeError("Asking for a non-existing initBlock type: {0}".format(fragmentType))
        # Copy the init fragment
        f = self._fragmentLibrary[ fragmentType ].copy()
        return ab_block.Block(initFragment=f)

    def growBlocks(self,
                   toGrow,
                   cellEndGroups=None,
                   libraryEndGroups=None,
                   endGroupType=None,
                   dihedral=None,
                   maxTries=50,
                   random=True):
        """
        Add toGrow new blocks to the cell.

        Args:
        toGrow: number of blocks to add
        cellEndGroups: a list of the endGroup types in the cell that the new blocks will be bonded to.
                      If more than one endGroup type is supplied, the endGroup will be randomly chosen from
                      that list.
        libraryEndGroups: a list of the endGroup types from the library that will be used form the bonds.
                      If more than one endGroup type is supplied, the endGroup will be randomly chosen from
                      that list.
        dihedral: the dihedral angle about the bond (3rd column in csv file)
        maxTries: number of attempts to make before giving up
        """
        logger.info("Growing {0} new blocks".format(toGrow))
        assert len(self.blocks), "Need to seed blocks before growing!"
        if endGroupType:
            logger.warn('endGroupType is deprecated! Use cellEndGroups/libraryEndGroups instead')
            assert not cellEndGroups or libraryEndGroups
            libraryEndGroups = endGroupType

        if dihedral: dihedral = math.radians(dihedral)
        added = 0
        tries = 0
        attemptedPairs = set()
        while added < toGrow:
            if tries >= maxTries:
                logger.critical("growBlocks - exceeded maxtries {0} when joining blocks!".format(maxTries))
                return added
            if self.numFreeEndGroups() == 0:
                logger.critical("growBlocks got no free endGroups!")
                return added
            # Select two random blocks that can be bonded
            try:
                endGroupPair = self.libraryEndGroupPair(cellEndGroups=cellEndGroups,
                                                                         libraryEndGroups=libraryEndGroups,
                                                                         random=random)
            except RuntimeError as e:
                logger.critical("growBlocks cannot grow more blocks: {0}".format(e))
                return added

            # See if we've seen this pair before
            if endGroupPair in attemptedPairs:
                logger.debug("growBlocks got endGroupPair again")
                tries += 1
                continue
            else:
                attemptedPairs.add(endGroupPair)
            cellEndGroup, libraryEndGroup = endGroupPair
            # Apply random rotation in 3 axes to randomise the orientation before we align
            libraryBlock = libraryEndGroup.block()
            if random:
                libraryBlock.randomRotate(origin=self.origin)
            # Try and attach it
            ok = self.attachBlock(libraryEndGroup, cellEndGroup, dihedral=dihedral)
            if ok:
                added += 1
                logger.info("growBlocks added block {0} after {1} tries.".format(added, tries))
                self.analyse.stop('grow', d={'num_tries':tries})
                # print "GOT BLOCK ",[ b  for b in self.blocks.values() ][0]
                tries = 0
            else:
                # del initBlock?
                tries += 1
        logger.info("After growBlocks numBlocks: {0}".format(len(self.blocks)))
        return added

    def growPolymer(self, monomers, ratio, length, random=False, center=False):
        """Create a linear polymer.

        Arguments:
        monomers - list of ambuild fragmentTypes (as specified when adding the fragment to the cell with libraryAddFragment)
                   that will be joined to create a subunit
        ratio - list of integers specifying the number of monomers of each type (needs to be the same length as monomers) in
                the subunit
        length - the number of subumits that will be created
        random - True/False - whether to build up the polymer deterministically (as per monomers/ratio) or stochastically, but
               - where the final ratio of fragments will approach that specified in ratio.
        center - place the first monomer polymer in the center of the cell (if possible)

        The fragment and ratio lists define the construction of the subunit. E.g.:
        monomers = ['A', 'B']
        ratio = [1,2]
        gives: [ABB]

        We could even have:
        monomers = ['A', 'B', 'A', 'B']
        ratio = [1,2,3,1]
        gives: [ABBAAAB]

        """
        logger.info("growPolymer, building polymer with monomers: {0}, ratio: {1}, length: {2}, random={3}".format(monomers,
                                                                                                                    ratio,
                                                                                                                    length,
                                                                                                                    random))
        totalTally = [0] * len(monomers) # Trakcs the overall number of fragments in the polymer
        # Determine the first endGroup type
        if random:
            monomer = _random.choice(monomers)
        else:
            monomer = monomers[0]
        # Seed the first block
        self.seed(nblocks=1, fragmentType=monomer, center=center)
        # Get the polymer block
        idxPolymer = self.lastAdded
        polymer = self.blocks[idxPolymer]
        fragment = polymer.fragments[0]

        # Create subunit going forward
        subunit = ab_subunit.subUnit(monomers=monomers,
                                     ratio=ratio,
                                     polymer=polymer,
                                     totalTally=totalTally,
                                     direction=1,
                                     fragment=fragment,
                                     random=random)
        switched = False
        for _ in range(length-1):
            if not subunit.addMonomer(self):
                if switched:
                    logger.critical("Failed to complete polymer!")
                    break
                logger.debug("Failed to add Monomer to intial endPoint - switching end")
                # Switch to other end of the chain -set switched flag so we know wev've done this
                subunit = ab_subunit.subUnit(monomers=monomers,
                                             ratio=ratio,
                                             polymer=polymer,
                                             totalTally=totalTally,
                                             direction=-1,
                                             random=random)
                switched = True
        logger.info("growPolymer finalTally ({0} % {1}): {2}".format(ratio, length, subunit.totalTally))
        # Bit of a hack - currently just so we can test the results
        self.polymerTally = subunit.totalTally
        return

    def haloCells(self, key):
        return xyz_util.haloCells(key, self.numBoxes, pbc=self.pbc)

    def _intersectedCells(self, p1, p2, endPointCells=True):
        """Return a list of the cells intersected by the vector passing from p1 to p2.

        Filched from:
        http://www.flipcode.com/archives/Raytracing_Topics_Techniques-Part_4_Spatial_Subdivisions.shtml
        http://stackoverflow.com/questions/12367071/how-do-i-initialize-the-t-variables-in-a-fast-voxel-traversal-algorithm-for-ray

        Could potentially speed things up by not returning the endpoints as they will be found by the surroundingCells algorithm
        and we then only search the ends.
        """
        TMAXMAX = max(self.dim[0], self.dim[1], self.dim[2])  # Not 100% sure - just needs to be bigger than any possible value in the cell
        # Wrap into a single cell
        if self.pbc[0]:
            p1[0] % self.dim[0]
            p2[0] % self.dim[0]
        if self.pbc[1]:
            p1[1] % self.dim[1]
            p2[1] % self.dim[1]
        if self.pbc[2]:
            p1[2] % self.dim[2]
            p2[2] % self.dim[2]

        X, Y, Z = self._getBox(p1)  # The cell p1 is in
        outX, outY, outZ = self._getBox(p2)  # The cell p2 is in
        dx, dy, dz = self.vecDiff(p2, p1)  # length components of line
        #
        # X
        #
        # stepX=int(math.copysign(1,dx)) # the direction of travel in the x-direction
        stepX = int(math.copysign(1, dx))  # the direction of travel in the x-direction
        # If there is no direction along a component, we need to make sure tMaxX is > the cell
        if X == outX:
            tDeltaX = 0
            tMaxX = TMAXMAX
        else:
            # tDeltaX - how far we can move along the x-component of the ray for the movement to equal
            # the width of a cell
            tDeltaX = self.boxSize / dx  # How many boxes fit in the x-direction
            # tMaxX - how far we can move along the x-coordinate before we cross a cell boundary
            tMaxX = tDeltaX * (1.0 - math.modf(p1[0] / self.boxSize)[0])
        #
        # Y
        #
        stepY = int(math.copysign(1, dy))
        if Y == outY:
            tDeltaY = 0
            tMaxY = TMAXMAX
        else:
            tDeltaY = self.boxSize / dy
            tMaxY = tDeltaY * (1.0 - math.modf(p1[1] / self.boxSize)[0])
        #
        # Z
        #
        stepZ = int(math.copysign(1, dz))
        if Z == outZ:
            tDeltaZ = 0
            tMaxZ = TMAXMAX
        else:
            tDeltaZ = self.boxSize / dz
            tMaxZ = tDeltaZ * (1.0 - math.modf(p1[2] / self.boxSize)[0])

        # print "dx ",dx,dy,dz
        # print "IN ",X,Y,Z
        # print "OUT ",outX,outY,outZ
        # print "stepX ",stepX
        cells = [(X, Y, Z)]  # Start with this cell
        while True:
            # Stop when we've reached the cell with p2
            if X == outX and Y == outY and Z == outZ: break
            if tMaxX < tMaxY:
                if tMaxX < tMaxZ:
                    # PBC
                    X = X + stepX
                    if self.pbc[0]:
                        if X < 0:
                            X = self.numBoxes[0] - 1
                        elif X > self.numBoxes[0] - 1:
                            X = 0
                    cells.append((X, Y, Z))
                    if X == outX:
                        tMaxX = TMAXMAX
                    else:
                        tMaxX += tDeltaX
                else:
                    Z = Z + stepZ
                    # PBC
                    if self.pbc[2]:
                        if Z < 0:
                            Z = self.numBoxes[2] - 1
                        elif Z > self.numBoxes[2] - 1:
                            Z = 0
                    cells.append((X, Y, Z))
                    if Z == outZ:
                        tMaxZ = TMAXMAX
                    else:
                        tMaxZ += tDeltaZ
            else:
                if tMaxY < tMaxZ:
                    # PBC
                    Y = Y + stepY
                    if self.pbc[1]:
                        if Y < 0:
                            Y = self.numBoxes[1] - 1
                        elif Y > self.numBoxes[1] - 1:
                            Y = 0
                    cells.append((X, Y, Z))
                    if Y == outY:
                        tMaxY = TMAXMAX
                    else:
                        tMaxY += tDeltaY
                else:
                    Z = Z + stepZ
                    # PBC
                    if self.pbc[2]:
                        if Z < 0:
                            Z = self.numBoxes[2] - 1
                        elif Z > self.numBoxes[2] - 1:
                            Z = 0
                    cells.append((X, Y, Z))
                    if Z == outZ:
                        tMaxZ = TMAXMAX
                    else:
                        tMaxZ += tDeltaZ

        # return tDeltaX,tMaxX,X,tDeltaY,tMaxY,Y
        if endPointCells:
            return cells
        else:
            return cells[1:-1]

    def joinBlocks(self, toJoin, cellEndGroups=None, dihedral=None, maxTries=100):
        """
        Bond toJoin blocks together using the endGroup types specified in cellEndGroups

        Args:
        toJoin - number of blocks to join
        cellEndGroups - a list of the different endGroupTypes that should be bonded. If this is None
                        randomly chosen endGroups will be used.
        dihedral: the dihedral angle about the bond (3rd column in csv file)
        maxTries - the maximum number of moves to try when joining
        """
        logger.info("Joining {0} new blocks".format(toJoin))
        if dihedral:
            # Convert dihedral to radians
            dihedral = math.radians(dihedral)
        added = 0
        tries = 0
        while added < toJoin:
            if len (self.blocks) == 1:
                logger.info("joinBlocks Hooray! - no more blocks to join!")
                return added
            if tries > maxTries:
                logger.critical("joinBlocks - exceeded maxtries when joining blocks!")
                return added
            if self.numFreeEndGroups() == 0:
                logger.critical("joinBlocks got no free endGroups!")
                return added
            # Select 2 random blocks that can be joined
            moveEndGroup, staticEndGroup = self.cellEndGroupPair(cellEndGroups=cellEndGroups)
            if moveEndGroup == None or staticEndGroup == None:
                logger.critical("joinBlocks cannot join any more blocks")
                return added
            # Copy the original block so we can replace it if the join fails
            moveBlock = moveEndGroup.block()
            idxMoveBlock = moveBlock.id
            blockCopy = moveBlock.copy()
            # Remove from cell so we don't check against itself and pick a different out
            self.delBlock(idxMoveBlock)
            logger.debug("joinBlocks calling attachBlock: {0} {1}".format(moveEndGroup, staticEndGroup))
            # now attach it
            ok = self.attachBlock(moveEndGroup, staticEndGroup, dihedral=dihedral)
            if ok:
                added += 1
                logger.info("joinBlocks joined block {0} after {1} tries.".format(added, tries))
                tries = 0
            else:
                # Put the original block back in the cell
                self.addBlock(blockCopy)
                tries += 1
        logger.info("After joinBlocks numBlocks: {0}".format(len(self.blocks)))
        return added

    def libraryAddFragment(self, filename, fragmentType='A', solvent=False, markBonded=False, catalyst=False):
        """Add a fragment of type fragmentType defined in the .car file filename

        Args:
        filename - the path to the .car file. There will need to be a corresponding .csv file that defines
                 - the endGroups, capAtoms etc.
        fragmentType - a name that will be used to identify the fragment - cannot contain the ":" character
        solvent - specify that this fragmentType is solvent and so won't be clash-checked in Zip steps
        """
        if fragmentType in self._fragmentLibrary:
            raise RuntimeError("Adding existing ftype {0} again!".format(fragmentType))
        # For now don't allow adding blocks when the cell has blocks in
        # assert not len(self.blocks),"Cannot add library fragments with populated cell!"
        # Make sure the type is valid
        if BONDTYPESEP in fragmentType or ENDGROUPSEP in fragmentType:
            raise RuntimeError("fragmentType cannot containing {0} or {1} characters!".format(BONDTYPESEP, ENDGROUPSEP))
        # Create fragment
        frag = ab_fragment.Fragment(filename,
                                    fragmentType,
                                    solvent=solvent,
                                    markBonded=markBonded,
                                    catalyst=catalyst)
        # Update cell parameters for this fragment
        maxAtomRadius = frag.maxAtomRadius()
        if  maxAtomRadius > self.maxAtomRadius:
            self.updateCellSize(maxAtomRadius=maxAtomRadius)
        # Add to _fragmentLibrary
        self._fragmentLibrary[fragmentType] = frag
        # create dictionary keyed by endGroup types
        for ft in frag.endGroupTypes():
            assert ft not in self._endGroup2LibraryFragment, "Adding existing endGroup type to library: {0}".format(ft)
            self._endGroup2LibraryFragment[ft] = fragmentType
        return

    def  _getCell2Library(self, endGroupTypes2Block, cellEndGroups=None, libraryEndGroups=None):
        if len(endGroupTypes2Block.keys()) == 0: raise RuntimeError("No available endGroups in the cell")
        # We create a dictionary mapping cell endGroups to possible libraryEndGroups
        cell2Library = {}
        for ceg in endGroupTypes2Block.keys():
            # We add those endGroup types that can be bonded to that are also in the library
            if ceg in self._bondTable:
                leg = self._bondTable[ ceg ].intersection(self._endGroup2LibraryFragment.keys())
                if len(leg) > 0: cell2Library[ ceg ] = leg

        # Check that there are some available
        if len(cell2Library.keys()) == 0:
            raise RuntimeError("No library fragments available to bond under the given rules: {0}".format(endGroupTypes2Block.keys()))

        # If the user supplied a list of cellEndGroups we prune the list of cell endGroups to those that are in this list
        if cellEndGroups is not None:
            # Convert to a set - make sure is a list first
            if isinstance(cellEndGroups, str):
                cellEndGroups = [ cellEndGroups ]
            for ceg in list(cell2Library.keys()):
                if ceg not in cellEndGroups:
                    del cell2Library[ ceg ]
            if len(cell2Library.keys()) == 0:
                raise RuntimeError("No free endGroups of types in cellEndGroups: {0} - {1}".format(cellEndGroups,
                                                                                                    endGroupTypes2Block.keys()))

        # If the user supplied a list of libraryEndGroups we remove any library endGroups that aren't in the list
        if libraryEndGroups is not None:
            # Convert to a set - make sure is a list first
            if isinstance(libraryEndGroups, str):
                libraryEndGroups = [ libraryEndGroups ]
            libraryEndGroups = set(libraryEndGroups)  # Save old so we can warn user and also find matching
            tmp = cell2Library
            cell2Library = {}
            for ceg, leg in tmp.items():
                # Only select those that are in the libraryEndGroups
                pleg = libraryEndGroups.intersection(leg)
                if len(pleg) > 0:
                    cell2Library[ ceg ] = pleg
            if not len(cell2Library.keys()):
                raise RuntimeError("No library fragments of type {0} available to bond under the given rules: {0}".format(libraryEndGroups,
                                                                                                                           endGroupTypes2Block.keys()))
        return cell2Library

    def libraryEndGroupPair(self, cellEndGroups=None, libraryEndGroups=None, random=True):
        """Return a fee endGroup from the cell and one from the library that can be bonded to it."""
        endGroupTypes2Block = self.endGroupTypes2Block()
        cell2Library = self._getCell2Library(endGroupTypes2Block,
                                             cellEndGroups=cellEndGroups,
                                             libraryEndGroups=libraryEndGroups)
        logger.debug("libraryEndGroupPair got cell/library endGroups: {0}".format(cell2Library))
        if random:
            # Now we can pick a random endGroup from the cell, get the corresponding library group
            cellEgT = _random.choice(list(cell2Library.keys()))

            # First get a block that contains this type of endGroup
            # Need to use sample as sets don't support random.choice
            cellBlock = _random.sample(endGroupTypes2Block[ cellEgT ], 1)[0]

            # Now select a random endGroup of that type from it
            cellEndGroup = cellBlock.selectEndGroup(endGroupTypes=[cellEgT], random=random)

            # Now get a corresponding library endGroup
            # We need to pick a random one of the types that we can bond to that is also in libraryTypes
            libEgT = _random.sample(cell2Library[cellEgT], 1)[0]

            # Now determine the fragmentType and create the block and fragment
            fragmentType = self._endGroup2LibraryFragment[libEgT]
            libraryBlock = self.getLibraryBlock(fragmentType=fragmentType, random=random)

            # now get the endGroup
            libraryEndGroup = libraryBlock.selectEndGroup(endGroupTypes=[libEgT])
        else:
            # We need to keep selecting a different but deterministic endGroupType each call
            i = self._deterministicState % len(cell2Library.keys())
            cellEgT = sorted(cell2Library.keys())[i]
            i = self._deterministicState % len(endGroupTypes2Block[cellEgT])
            # sort blocks by id - what we use is irrelevant, it just needs to be consistent
            cellBlock = sorted(list(endGroupTypes2Block[cellEgT]),
                               key=lambda block: block.id)[i]
            cellEndGroup = cellBlock.selectEndGroup(endGroupTypes=[cellEgT], random=random)
            i = self._deterministicState % len(cell2Library[cellEgT])
            libEgT = sorted(list(cell2Library[cellEgT]))[i]
            fragmentType = self._endGroup2LibraryFragment[libEgT]
            libraryBlock = self.getLibraryBlock(fragmentType=fragmentType, random=random)
            libraryEndGroup = libraryBlock.selectEndGroup(endGroupTypes=[libEgT], random=random)
            self._deterministicState += 1
        # Return them both - phew!
        logger.debug("libraryEndGroupPair returning: {0} {1}".format(cellEndGroup.type(), libraryEndGroup.type()))
        return (cellEndGroup, libraryEndGroup)

    def numBlocks(self):
        return len(self.blocks)

    def numFragments(self):
        # return sum([len(b.fragments) for b in self.blocks.values()])
        # Hack for old interface
        return sum([len(b.fragments) for b in self.blocks.values()])

    def numFreeEndGroups(self):
        return sum([ b.numFreeEndGroups() for b in self.blocks.values() ])

    def numAtoms(self):
        return sum([ b.numAtoms() for b in self.blocks.values() ])

    def optimiseGeometry(self,
                         rigidBody=True,
                         doDihedral=False,
                         doImproper=False,
                         doCharges=True,
                         xmlFilename="hoomdOpt.xml",
                         **kw):
        """Optimise the geometry with hoomdblue

        Rigid-body optimisation is carried out with the hoomdblue integrator: mode_minimize_rigid_fire

        Args:
        rigidBody - True/False - do rigid body or all-atom optimisation
        doDihderal - True/False - include dihedral terms
        doImproper - True/False - include improper terms
        rCut - the VdW cutoff to use [angstroms]
        optCycles - the number of hoomdblue optimisation cycles to run.
        quiet - True/False - don't print out the normal hoomdblue runtime output to the screen.
        ALL OTHER ARGUMENTS ACCEPTED BY mode_minimize_rigid_fire ARE PASSED TO IT
        """
        logger.info("Running optimisation")
        if not self.mdEngineCls:
            raise RuntimeError("No mdEngine defined - cannot run MD.")
        mdEngine = self.mdEngineCls(self.paramsDir)
        if doDihedral and doImproper:
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.setRcut(rigidBody, mdEngine, kw)
        data = self.cellData(periodic=True, center=True, rigidBody=rigidBody)
        d = {}  # for printing results
        ok = mdEngine.optimiseGeometry(data,
                                        xmlFilename=xmlFilename,
                                        rigidBody=rigidBody,
                                        doDihedral=doDihedral,
                                        doImproper=doImproper,
                                        doCharges=doCharges,
                                        d=d,
                                        walls=self.walls,
                                        wallAtomType=self.wallAtomType,
                                        **kw)
        self.analyse.stop('optimiseGeometry', d)
        if ok:
            logger.info("Optimisation succeeded")
            mdEngine.updateCell(self)
            return True
        else:
            logger.critical("Optimisation Failed")
            return False

    def positionInCell(self, block):
        """Make sure the given block is positioned within the cell"""
        bradius = block.blockRadius()
        oradius = bradius + (2 * self.boxMargin)
        # Check there is enough space
        if oradius >= self.dim[0] or oradius >= self.dim[1] or oradius >= self.dim[2]:
            raise RuntimeError("Cannot fit block with radius {0} into cell {1}!".format(bradius, self.dim))
        # get a range for x, y and z that would fit the block in the cell, pick random values within that
        # and stick the block there
        x = _random.uniform(bradius + self.boxMargin, self.dim[0] - bradius - self.boxMargin)
        y = _random.uniform(bradius + self.boxMargin, self.dim[1] - bradius - self.boxMargin)
        z = _random.uniform(bradius + self.boxMargin, self.dim[2] - bradius - self.boxMargin)
        coord = np.array([x, y, z], dtype=np.float64)
        block.translateCentroid(coord)
        logger.debug("positionInCell block moved to: {0}".format(block.centroid()))
        return

    def processBonds(self, selfBond=True):
        """Make any bonds that were found during checkMove
        return Number of bonds made

        Args:
        addedBlockIdx - the block that was added
        """
        if not len(self._possibleBonds):
            logger.debug("processBonds got no bonds")
            return 0
        # Here no atoms clash and we have a list of possible bonds - so bond'em!
        logger.debug("processBonds got bonds: {0}".format([str(s) for s in self._possibleBonds]))
        logger.debug("processBonds blocks are: {0}".format(sorted(self.blocks)))
        bondsMade = 0
        self.newBonds = []
        for count, bond in enumerate(self._possibleBonds):
            # With the bonding rules some endGroups may become not free when other endGroups in that fragment
            # are involved in bonds so we need to make sure they are free before we do
            if bond.isPossible():
                made = self.bondBlock(bond, selfBond=selfBond)
                if made:
                    logger.debug("Added bond: {0}".format(self._possibleBonds[count]))
                    bondsMade += 1
                    self.newBonds.append(bond)
        # Clear any other possible bonds
        self._possibleBonds = []
        return bondsMade

    def positionBlock(self, block, margin=None, point=None, radius=None, zone=None, random=True):
        """Randomly move the given block
         If margin is given, use this as a buffer from the edges of the cell
         when selecting the coord
        """
        if point:
            assert len(point) == 3,"Point needs to be a list of three floats!"
            if not (0 <= point[0] <= self.dim[0] and 0 <= point[1] <= self.dim[1] and \
                    0 <= point[2] <= self.dim[2] ):
                raise RuntimeError("Point needs to be inside the cell.")
            if radius:
                coord = self.randomSpherePoint(point, radius)
            else:
                coord = np.array(point, dtype=np.float64)
        elif zone:
            if not len(zone) == 6:
                raise RuntimeError("Zone needs to be a list of 6 floats: [x0,x1,y1,y2,z1,z1]")
            if not zone[0] >= 0 and zone[1] <= self.dim[0] and\
                zone[2] >= 0 and zone[3] <= self.dim[1] and \
                zone[4] >= 0 and zone[5] <= self.dim[2]:
                raise RuntimeError("Zone needs to be within the cell")
            #bmargin = block.blockRadius()
            bmargin = 0
            x = _random.uniform(zone[0] + bmargin, zone[1] - bmargin)
            y = _random.uniform(zone[2] + bmargin, zone[3] - bmargin)
            z = _random.uniform(zone[4] + bmargin, zone[5] - bmargin)
            coord = np.array([x, y, z], dtype=np.float64)
        else:
            if margin:
                x = _random.uniform(margin, self.dim[0] - margin)
                y = _random.uniform(margin, self.dim[1] - margin)
                z = _random.uniform(margin, self.dim[2] - margin)
            else:
                x = _random.uniform(0, self.dim[0])
                y = _random.uniform(0, self.dim[1])
                z = _random.uniform(0, self.dim[2])
            coord = np.array([x, y, z], dtype=np.float64)
        # print "Got random coord: {}".format(coord)
        # Move to origin, rotate there and then move to new coord
        # Use the cell axis definitions
        if random:
            block.translateCentroid(self.origin)
            block.randomRotate(origin=self.origin, atOrigin=True)
        block.translateCentroid(coord)
        return

    def randomSpherePoint(self, center, radius):
        """Return a random point on a sphere of radius radius centerd at center.
        From: http://mathworld.wolfram.com/SpherePointPicking.html
        """
        # We use spherical coordinates and calculate theta and phi
        theta = 2 * math.pi * _random.uniform(0,1)
        phi = math.acos(2 * _random.uniform(0,1) - 1)
        rradius = _random.uniform(0,radius)
        # coordinate of point centered at origin is therefore:
        x = rradius * math.cos(theta) * math.sin(phi)
        y = rradius * math.sin(theta) * math.sin(phi)
        z = rradius * math.cos(phi)
        return np.array([x + center[0], y + center[1], z + center[2]], dtype=np.float64)

    def repopulateCells(self, boxShift=None):
        """Add all the blocks to resized cells"""
        blocks = None
        if len(self.blocks):
            # Put all the blocks into the new cells
            # First copy the blocks dictionary
            blocks = copy.copy(self.blocks)
        self.clear()
        if len(blocks):
            logger.debug("repopulateCells, adding blocks into new cells")
            for idxBlock, block in blocks.items():
                if boxShift: block.translate(boxShift)  # Need to move if box has been resizes
                self.addBlock(block, idxBlock=idxBlock)
            del blocks
        return

    def restoreBlocks(self, blocks, fragmentTypes=None):
        assert len(blocks) > 0,"Need blocks to restore!"
        added = 0
        if fragmentTypes is not None:
            if type(fragmentTypes) is str: fragmentTypes = [fragmentTypes]
            toAdd = []
            fragmentTypes = set(fragmentTypes)
            for block in blocks:
                if fragmentTypes.intersection(set(block.fragmentTypes())):
                    toAdd.append(block) #if any of the fragments are within the blocks remove this block
            # Remove the blocks
            if len(toAdd):
                for block in toAdd:
                    blocks.remove(block)
                added = self.addBlocks(toAdd)
        else:
            added = self.addBlocks(blocks)
        logger.info("restoreBlocks re-added {0} blocks to the cell".format(added))
        return added

    def runMD(self,
              xmlFilename="hoomdMD.xml",
              rigidBody=True,
              doDihedral=False,
              doImproper=False,
              doCharges=True,
              **kw):
        """Run a rigidbody molecular dynamics run using the hoomd blue nvt_rigid integrator

        Arguments:
        rigidBody - True/False - do rigid body or all-atom MD
        doDihderal - True/False - include dihedral terms
        doImproper - True/False - include improper terms
        rCut - the VdW cutoff to use [angstroms]
        mdCycles - the number of MD cycles to run.
        quiet - True/False - don't print out the normal hoomdblue runtime output to the screen.
        T - nvt_rigid temperature
        tau - nvt_rigid tau
        dt - nvt_rigid timestep
        """
        logger.info("Running MD")
        if not self.mdEngineCls:
            raise RuntimeError("No MDENGINE defined - cannot run MD.")
        mdEngine = self.mdEngineCls(self.paramsDir)
        if doDihedral and doImproper:
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.setRcut(rigidBody, mdEngine, kw)
        d = {}
        data = self.cellData(periodic=True, center=True, rigidBody=rigidBody)
        ok = mdEngine.runMD(data,
                             xmlFilename=xmlFilename,
                             rigidBody=rigidBody,
                             doDihedral=doDihedral,
                             doImproper=doImproper,
                             doCharges=doCharges,
                             d=d,
                             walls=self.walls,
                             wallAtomType=self.wallAtomType,
                             **kw)
        self.analyse.stop('runMD', d)
        mdEngine.updateCell(self)
        return ok

    def runMDAndOptimise(self,
                         xmlFilename="hoomdMDOpt.xml",
                         rigidBody=True,
                         doDihedral=False,
                         doImproper=False,
                         doCharges=True,
                         **kw):

        """Run an MD simulation followed by a Geometry optimisation.

        Args:
        See runMD and optimiseGeometry for acceptable arguments.
        """
        if not self.mdEngineCls:
            raise RuntimeError("No MDENGINE defined - cannot run MD.")
        mdEngine = self.mdEngineCls(self.paramsDir)
        assert rigidBody, "FIX runMD FOR ALL ATOM!!"
        if doDihedral and doImproper:
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.setRcut(rigidBody, mdEngine, kw)
        d = {}
        data = self.cellData(periodic=True, center=True, rigidBody=rigidBody)
        ok = mdEngine.runMDAndOptimise(data,
                                        xmlFilename=xmlFilename,
                                        rigidBody=rigidBody,
                                        doDihedral=doDihedral,
                                        doImproper=doImproper,
                                        doCharges=doCharges,
                                        d=d,
                                        walls=self.walls,
                                        wallAtomType=self.wallAtomType,
                                        **kw)

        if ok:
            logger.info("runMDAndOptimise succeeded")
        self.analyse.stop('runMDAndOptimise', d)
        mdEngine.updateCell(self)
        return ok

    def seed(self,
             nblocks,
             fragmentType=None,
             maxTries=500,
             center=False,
             point=None,
             radius=None,
             zone=None,
             random=True):
        """ Seed a cell with nblocks of type fragmentType.

        Args:
        nblocks - the number of blocks to add.
        fragmentType - the type of blocks to add. If fragment is None, or omitted, then blocks of a randomly
                       chosen type will be added.
        maxTries - the number of attempts to make when adding a block before the seed step is fails and returns.
        center - True/False - if True, place the first block in the center of the cell.
        point - a list of three floats defiting a point around which the centroids of the blocks will be seeded
                (requires the radius argument).
        radius - a single float specifying the radius of a sphere around point, within which the centroids of
                 the blocks will be seeded (requires point argument).
        zone - a list of 6 floats specifying a box within the cell within which the centroids of the blocks will
               be seeded.
        random - randomly rotate the blocks on seeding
        Returns:
        the number of blocks added
        """
        if self.dim[0] is None or self.dim[1] is None or self.dim[2] is None:
            raise RuntimeError("Need to specify cell before seeding")

        # if not len( self._fragmentLibrary ) or not len( self.bondTypes):
        #    raise RuntimeError,"Must have set an initBlock and bondType before seeding."
        if not len(self._fragmentLibrary):
            raise RuntimeError("Must have set an initBlock before seeding.")
        logger.info("seed adding {0} block of type {1}".format(nblocks, fragmentType))
        numBlocksAdded = 0
        # Loop through the nblocks adding the blocks to the cell
        for seedCount in range(nblocks):
            newblock = self.getLibraryBlock(fragmentType=fragmentType) # Create new block
            tries = 0
            while True:
                # quit on maxTries
                if tries >= maxTries:
                    logger.critical("Exceeded maxtries when seeding after adding {0}".format(numBlocksAdded))
                    self.analyse.stop('seed', d={'num_tries':tries})
                    return numBlocksAdded
                # if center put the first one in the center of the cell
                # only if this is the first attempt as otherwise we always fail if there is already something there
                if center and seedCount == 0 and tries == 0:
                    newblock.translateCentroid([ self.dim[0] / 2, self.dim[1] / 2, self.dim[2] / 2 ])
                else:
                    # Move the block and rotate it
                    self.positionBlock(newblock, point=point, radius=radius, zone=zone, random=random)
                # Add the block so we can check for clashes/bonds
                idxBlock = self.addBlock(newblock)
                # Test for Clashes with other molecules
                if self.checkMove(idxBlock):
                    if self.processBonds() > 0:
                        logger.info("Added bond in seed!")
                    logger.debug("seed added block {0} after {1} tries.".format(seedCount + 1, tries))
                    self.analyse.stop('seed', d={'num_tries':tries})
                    numBlocksAdded += 1
                    break
                # Unsuccessful so remove the block from cell
                self.delBlock(idxBlock)
                # If seed fails with center need to bail on first one.
                if center and seedCount == 0 and tries == 0:
                    logger.warn("Seed with center failed to place first block in center!")
                tries += 1 # increment tries counter
            # End Clash loop
        # End of loop to seed cell
        logger.info("Seed added {0} blocks. Cell now contains {1} blocks".format(numBlocksAdded, len(self.blocks)))
        return numBlocksAdded

    def setBondingFunction(self, fragmentType, onbondFunction):
        fragment = self._fragmentLibrary[ fragmentType ]
        fragment.onbondFunction = onbondFunction
        return

    def setBoxSize(self, boxDim):
        # A, B and C are cell vectors - for the time being we assume they are orthogonal
        # so they are just floats
        if not type(boxDim) is list and len(boxDim) == 3:
            raise RuntimeError("setBoxSize needs list of 3 numbers setting the cell dimensions!")
        old_dim = None
        if self.dim is not None: old_dim = self.dim
        self.dim = np.array([float(boxDim[0]), float(boxDim[1]), float(boxDim[2])])
        assert self.dim[0] > 0 and self.dim[1] > 0 and self.dim[2] > 0, "Invalid box dimensions: {0}".format(self.dim)
        boxShift = None
        if old_dim is not None and not np.allclose(self.dim, old_dim):
            if not (self.dim[0] >= old_dim[0] and \
                self.dim[1] >= old_dim[1] and
                self.dim[2] >= old_dim[2]):
                raise RuntimeError("Can currently only increase box size!")
            # If we've made the cell bigger we need to shift all blocks so that they sit
            # in the centre again
            boxShift = [ (self.dim[0] - old_dim[0]) / 2, \
                         (self.dim[1] - old_dim[1]) / 2, \
                         (self.dim[2] - old_dim[2]) / 2 ]
        self.updateCellSize(boxShift=boxShift)
        return

    def setStaticBlock(self, filePath, replace=False):
        # Create fragment
        name = os.path.splitext(os.path.basename(filePath))[0]
        f = ab_fragment.Fragment(filePath, fragmentType=name, static=True)
        p = f.cellParameters()
        if not p:
            raise RuntimeError("car file needs to have PBC=ON and a PBC line defining the cell!")
        logger.info("Read cell parameters A={0}, B={1}, C={2} from car file: {3}".format(
            p['A'],p['B'],p['C'],filePath))
        self.setBoxSize([p['A'], p['B'], p['C']])
        # Update cell parameters for this fragment
        maxAtomRadius = f.maxAtomRadius()
        if  maxAtomRadius > self.maxAtomRadius:
            self.updateCellSize(maxAtomRadius=maxAtomRadius)
        block = ab_block.Block(initFragment=f)
        # Check the molecule fits in the cell
        for i, coord in enumerate(block.iterCoord()):
            if coord[0] < 0 or coord[0] > self.dim[0] or \
               coord[1] < 0 or coord[1] > self.dim[1] or \
               coord[2] < 0 or coord[2] > self.dim[2]:
                raise RuntimeError("Static block doesn't fit in the cell! First failing coord is #{0}: {1}".format(i, coord))
        if replace:
            self.delBlock(list(self.blocks.keys())[0]) # Assumes a static block is always the first
        idxBlock = self.addBlock(block)
        if len(self.blocks) > 1:
            # If the cell already has blocks in, we need to remove them all and then add this as the first one
            d = collections.OrderedDict()
            d[idxBlock] = block
            for k, v in self.blocks.items():
                d[k] = v
            del self.blocks
            self.blocks = d
            # Also need to test for Clashes with other molecules
            if not self.checkMove(idxBlock):
                raise RuntimeError("Problem adding static block: got clashes!")
            if self.processBonds() > 0:
                raise RuntimeError("Problem adding static block-we made bonds!")
        logger.info("Added static block to cell")
        return

    def setMaxBond(self, bondType, count):
        """Limit the number of bondType bonds to an individual fragment to count bonds.

        Args:
        bondType - the bondType (FRAGMENT1:ENDGROUP1-FRAGMENT2:ENDGROUP2) as was specified with the call
                   to addBondType
        count - the maximum number of permissible bonds for a single fragment.

        """
        fragmentType = bondType.split(ENDGROUPSEP)[0]
        fragment = self._fragmentLibrary[ fragmentType ]
        return fragment.setMaxBond(bondType, count)

    def setMdEngineCls(self, hoomdVersion):
        """Set the class definition for the MdEngine
        
        We set a class rather then an instance as otherwise when running mutiple optimisations certain variables
        get stored in the instance and aren't deleted on restart, generating errors in hoomd1 along the lines of:
        
        *Warning*: Not all saved variables were cleared before calling reset()
        *Warning*: 17 references to the particle data still exist somewhere
        *Warning*: Going to try and reset anyways, further errors (such as out of memory) may result
        
        Whenever we ditch hoomdblue 1, we can (probably) revert to using an instance.
        """
        if hoomdVersion is None:
            logger.critical("HOOMD-BLUE could not be found! MD functionality will be unavaiable.")
            return
        if hoomdVersion[0] < 2:
            from hoomd1 import Hoomd1
            self.mdEngineCls = Hoomd1
            #self.MDENGINE_CLS = Hoomd1
        else:
            from hoomd2 import Hoomd2
            self.mdEngineCls = Hoomd2
        logger.info("Using HOOMD-BLUE version: {0}.{1}.{2}".format(*ab_util.HOOMDVERSION))
        return
    
    def setRcut(self, rigidBody, mdEngine, kw):
        """if rCut not in kw, for hoomd2 rigidBodies calculate from max block size or use mdEngine default"""
        RIGIDPARTICLES = rigidBody and ab_util.HOOMDVERSION and ab_util.HOOMDVERSION[0] > 1
        if 'rCut' in kw:
            rCut = kw['rCut']
        elif RIGIDPARTICLES:
            rmax =  max([b.blockRadius() for b in self.blocks.values()])
            cell_r = np.max(self.dim) / 2
            rCut = min(rmax, cell_r)
            logger.info("Rigid Particle code setting rCut to: %f", rCut)
        else:
            rCut = mdEngine.rCut
        self.rCut =rCut
        logger.debug("Set rCut to %f", rCut)
        return

    def setWall(self, XOY=False, XOZ=False, YOZ=False, wallAtomType='c'):
        """Create walls along the specified sides.
        """
        self.walls = [XOY, XOZ, YOZ]
        self.pbc = [ not b for b in self.walls ]

        self.wallAtomType = wallAtomType
        symbol = wallAtomType.upper()
        if symbol not in xyz_core.SYMBOL_TO_NUMBER:
            raise RuntimeError("Data for atomType {0} not available. Supported atomTypes are: {1}".format(wallAtomType,sorted(xyz_core.SYMBOL_TO_NUMBER.keys())))
        z = xyz_core.SYMBOL_TO_NUMBER[ symbol ]
        r = xyz_core.COVALENT_RADII[z] * xyz_core.BOHR2ANGSTROM
        self.wallRadius = r
        return
    
    def _setupAnalyse(self, logfile='ambuild.csv'):
        self.analyse = ab_analyse.Analyse(self, logfile=logfile)
        return

    def setupLogging(self, logfile="ambuild.log", mode='w', debugLog=False):
        """
        Set up the various log files/console logging and return the logger
        """
        logger = logging.getLogger()
        # if not debugLog:
        #    logging.disable(logging.DEBUG)

        # First check if there are any handlers active - if so we remove them
        if len(logger.handlers):
            # Need to copy as otherwise the list changes while we're cycling through it
            for hdlr in copy.copy(logger.handlers):
                logger.removeHandler(hdlr)
        # Not entirely sure why this needed - set overall level of the logger to debug
        logger.setLevel(logging.DEBUG)
        # create file handler and set level to debug
        self.logfile = logfile
        fl = logging.FileHandler(self.logfile, mode=mode)
        if debugLog:
            fl.setLevel(logging.DEBUG)
        else:
            fl.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fl.setFormatter(formatter)
        logger.addHandler(fl)
        # Now create console logger for outputting stuff
        # create file handler and set level to debug
        # Seems they changed the api in python 2.6->2.7
        try:
            cl = logging.StreamHandler(stream=sys.stdout)
        except TypeError:
            cl = logging.StreamHandler(strm=sys.stdout)
        cl.setLevel(logging.INFO)
        #cl.setLevel( logging.DEBUG )
        # Always add a blank line after every print
        formatter = logging.Formatter('%(asctime)s: %(message)s\n')
        cl.setFormatter(formatter)
        logger.addHandler(cl)
        logger = logger
        return

    def updateCellSize(self, boxMargin=None, maxAtomRadius=None, MARGIN=0.01, boxShift=None):
        """The cell size is the vdw radius of the largest atom plus the largest of atom or bondMargin
        plus a MARGIN to account for overflows
        """
        if maxAtomRadius != None:
            self.maxAtomRadius = maxAtomRadius
        if not boxMargin is None:
            self.boxMargin = boxMargin
        else:
            assert self.atomMargin and self.bondMargin
            # box margin is largest we'll be looking for plus a bit to account for arithmetic overflow
            self.boxMargin = max(self.atomMargin, self.bondMargin) + MARGIN
        assert self.boxMargin != 0 and self.maxAtomRadius != 0
        self.boxSize = (self.maxAtomRadius * 2) + self.boxMargin
        # jmht - ceil or floor
        self.numBoxes[0] = int(math.ceil(self.dim[0] / self.boxSize))
        self.numBoxes[1] = int(math.ceil(self.dim[1] / self.boxSize))
        self.numBoxes[2] = int(math.ceil(self.dim[2] / self.boxSize))
        logger.debug("updateCellSize: boxSize {0} nboxes: {1} maxR {2} margin {3}".format(self.boxSize,
                                                                                      self.numBoxes,
                                                                                     self.maxAtomRadius,
                                                                                     self.boxMargin))
        # If there are already blocks in the cell, we need to add them to the new boxes
        if len(self.blocks):
            self.repopulateCells(boxShift=boxShift)
        return

    def vecDiff(self, p1, p2):
        return xyz_core.vecDiff(p1, p2, dim=self.dim, pbc=self.pbc)

    def writePickle(self, fileStem, compress=True):
        """Pickle ourselves"""
        # No idea why I can't get the log to close and then reopen with append mode
        # see _get_state_ and _set_state_ in this class
        if False:
            self._clLogHandler.close()
            self._flLogHandler.close()
            logger.removeHandler(self._clLogHandler)
            logger.removeHandler(self._flLogHandler)
            del self._clLogHandler
            del self._flLogHandler
        if compress:
            suffix = ab_util.GZIP_PKL_SUFFIX
        else:
            suffix = ab_util.PKL_SUFFIX
        fileName = os.path.abspath(fileStem + suffix)
        fileName = ab_util.pickleObj(self, fileName, compress=compress)
        # Restart logging with append mode
        # self.setupLogging( mode='a' )
        logger.info("Wrote pickle file: {0}".format(fileName))
        return fileName

    def writeCar(self, ofile="ambuild.car", data=None, periodic=True, skipDummy=False):
        """Car File
        """
        if not data:
            data = self.cellData()
        car = "!BIOSYM archive 3\n"
        if periodic:
            car += "PBC=ON\n"
        else:
            car += "PBC=OFF\n"
        car += "ambuild generated car file\n"
        tstr = time.strftime("%a %b %d %H:%M:%S %Y", time.gmtime())
        car += "!DATE {0}\n".format(tstr)
        if periodic:
            car += "PBC  {0: < 9.4F} {1: < 9.4F} {2: < 9.4F}  90.0000   90.0000   90.0000 (P1)\n".format(data['A'],
                                                                                                          data['B'],
                                                                                                          data['C'])
#         for i, ( idxBlock, block ) in enumerate( self.blocks.items() ):
#             for j, coord in enumerate( block.iterCoord() ):
        for i, coord in enumerate(data['coord']):
                if periodic:
                    coord, _ = xyz_core.wrapCoord3(coord, np.array([data['A'], data['B'], data['C']]), center=False)
                atype = data['type'][ i ]
                charge = data['charge'][ i ]
                label = data['label'][ i ][:5]
                symbol = data['symbol'][ i ]
                car += "{0: <5} {1: >15.10} {2: >15.10} {3: >15.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format(label, coord[0], coord[1], coord[2], atype, symbol, charge)

        car += "end\nend\n\n"

        with open(ofile, 'w') as f:
            fpath = os.path.abspath(f.name)
            f.writelines(car)

        logger.info("Wrote car file: {0}".format(fpath))
        return

    def writeCml(self, cmlFilename, data=None, rigidBody=True, periodic=True, pruneBonds=False, prettyPrint=False):
        if data is None:
            d = self.cellData(periodic=periodic, rigidBody=rigidBody, fragmentType=None)
        else:
            d = data

        cell = None
        if periodic:
            cell = self.dim

        cmlFilename = xyz_util.writeCml(cmlFilename,
                                        d.coords,
                                        d.symbols,
                                        bonds=d.bonds,
                                        atomTypes=d.atomTypes,
                                        cell=cell,
                                        pruneBonds=pruneBonds,
                                        prettyPrint=prettyPrint)

        logger.info("Wrote cml file: {0}".format(cmlFilename))
        return

    def writeXyz(self, ofile, data=None, periodic=False, atomTypes=False):
        """Write out the cell atoms to an xyz file
        If label is true we write out the atom label and block, otherwise the symbol
        """
        if data is None:
            d = self.cellData(periodic=periodic, fragmentType=None)
        else:
            d = data
        if periodic:
            if atomTypes:
                fpath = xyz_util.writeXyz(ofile, d.coords, d.atomTypes, cell=self.dim)
            else:
                fpath = xyz_util.writeXyz(ofile, d.coords, d.symbols, cell=self.dim)
        else:
            if atomTypes:
                fpath = xyz_util.writeXyz(ofile, d.coords,d.atomTypes)
            else:
                fpath = xyz_util.writeXyz(ofile, d.coords, d.symbols)

        logger.info("Wrote cell file: {0}".format(fpath))
        return

    def zipBlocks(self,
                  bondMargin=0.5,
                  bondAngleMargin=15,
                  clashCheck=False,
                  clashDist=1.6,
                  selfBond=True):
        """Join existing blocks in the cell by changing the bondMargin and bondAngleMargin parameters that were
        specified when the cell was created, and then looping over all the free endGroups to see if any can bond
        with the new parameters. The blocks are not moved in this step.

        Args:
        bondMargin - the new bondMargin [degrees] (see __init__ for definition)
        bondAngleMargin - the new bondAngleMargin [degrees] (see __init__ for definition)
        clashCheck - True/False check for clashes between the bond and any atoms that fall within a cylinder
                     of radius clashDist (default=1.6A) centered on the bond axis.
        clashDist  - a float specifying the perpendicular distance from the bond axis that determines if an atom
                     is clashing with the bond. clashDist needs to be < the the cell box size as otherwise we won't
                     see the atom.
        selfBond  - boolean to specify if zip will allow a block to bond to itself (True) or not (False) [default: True]
        """
        if bondMargin > max(self.dim):
            raise RuntimeError("bondMargin is greater then the cell")
        logger.info("Zipping blocks with bondMargin: {0} bondAngleMargin {1}".format(bondMargin, bondAngleMargin))
        # Convert to radians
        bondAngleMargin = math.radians(bondAngleMargin)
        # Calculate the number of boxes - should calculate max possible bond length
        maxBondLength = 2.5
        boxSize = maxBondLength + bondMargin
        numBoxes = [None, None, None]
        numBoxes[0] = int(math.ceil(self.dim[0] / boxSize))
        numBoxes[1] = int(math.ceil(self.dim[1] / boxSize))
        numBoxes[2] = int(math.ceil(self.dim[2] / boxSize))
        # Create empty cells to hold data
        cell1 = {}
        cell3 = {}
        # Get a list of all block, endGroups
        endGroups = []
        for block in self.blocks.values():
            egs = block.freeEndGroups()
            if len(egs):
                block.zipCell = {}
                for endGroup in egs:
                    endGroups.append((block, endGroup.endGroupIdx()))
        if not len(endGroups) > 0:
            logger.warn("zipBlocks found no free endGroups!")
            return 0
        # Add all (block, idxEndGroup) tuples to the cells
        for (block, idxEndGroup) in endGroups:
            coord = block.coord(idxEndGroup)
            # Periodic Boundaries
            x = coord[0] % self.dim[0] if self.pbc[0] else coord[0]
            y = coord[1] % self.dim[1] if self.pbc[1] else coord[1]
            z = coord[1] % self.dim[2] if self.pbc[2] else coord[2]
            # Calculate which cell the atom is in
            a = int(math.floor(x / boxSize))
            b = int(math.floor(y / boxSize))
            c = int(math.floor(z / boxSize))
            key = (a, b, c)
            block.zipCell[idxEndGroup] = key
            try:
                # Add the atom to the cell
                cell1[key].append((block, idxEndGroup))
            except KeyError:
                # Add the cell to the list and then add the atom
                cell1[key] = [(block, idxEndGroup)]
                # Map the cells surrounding this one
                cell3[key] = xyz_util.haloCells(key, numBoxes, pbc=self.pbc)
        # Loop through all the end groups and run canBond
        egPairs = collections.OrderedDict() # Use dict keys as effectively implements an ordered set
        c1 = []
        c2 = []
        for block1, idxEndGroup1 in endGroups:
            key = block1.zipCell[idxEndGroup1] # Get the box this atom is in
            surroundingCells = cell3[key]
            for sbox in surroundingCells:
                # Check if we have a box with anything in it
                # use exception so we don't have to search through the whole list
                try:
                    alist = cell1[sbox]
                except KeyError:
                    continue
                for (block2, idxEndGroup2) in alist:
                    if block1 == block2: 
                        # Self-bonded blocks need special care
                        if not selfBond:
                            continue
                        if idxEndGroup1 == idxEndGroup2: # Don't check endGroups against themselves
                            continue
                        # Make sure the two atoms are separated by at least 3 bonds.
                        # Could put this check in canBond but it would slow the normal bonding down
                        if idxEndGroup2 in block1.atomBonded3(idxEndGroup1):
                            continue
                    pair1 = (block1, idxEndGroup1, block2, idxEndGroup2)
                    pair2 = (block2, idxEndGroup2, block1, idxEndGroup1)
                    if pair1 not in egPairs and pair2 not in egPairs:
                        # Check if it is already present as we loop over all endGroups so we will have both sides twice
                        egPairs[pair1] = None
                        c1.append(block1.coord(idxEndGroup1))
                        c2.append(block2.coord(idxEndGroup2))
        if len(egPairs.keys()) < 1:
            logger.info("zipBlocks: no endGroups close enough to bond")
            return 0
        distances = self.distance(np.array(c1), np.array(c2))
        self._possibleBonds = []
        for i, (block1, idxEndGroup1, block2, idxEndGroup2) in enumerate(egPairs.keys()):
            got = self.canBond(block1,
                                idxEndGroup1,
                                block2,
                                idxEndGroup2,
                                distances[i],
                                bondMargin,
                                bondAngleMargin)
        # Process any bonds
        if len(self._possibleBonds) == 0:
            logger.info("zipBlocks: no acceptable bonds found")
            return 0
        # Check the bonds don't clash with anything
        if clashCheck:
            logger.info("zipBlocks: checking for clashes with bonds...")
            toRemove = [bond for bond in self._possibleBonds if self.bondClash(bond=bond, clashDist=clashDist)]
            if len(toRemove):
                logger.info("zipBlocks: {0} bonds not accepted due to clashes".format(len(toRemove)))
                for b in toRemove: self._possibleBonds.remove(b)
                if not len(self._possibleBonds):
                    logger.info("zipBlocks: No bonds remaining after clash checks")
                    return 0
        logger.info("zipBlocks found %d additional bonds", len(self._possibleBonds))
#         for b in self._possibleBonds:
#             print "Attempting to bond: {0} {1} {2} -> {3} {4} {5}".format( b.block1.id(),
#                                                                    b.endGroup1.blockEndGroupIdx,
#                                                                    b.block1.atomCoord( b.endGroup1.blockEndGroupIdx),
#                                                                    b.block2.id(),
#                                                                    b.endGroup2.blockEndGroupIdx,
#                                                                    b.block2.atomCoord( b.endGroup2.blockEndGroupIdx),
#                                                                 )
#
        todo = len(self._possibleBonds)
        bondsMade = self.processBonds(selfBond=selfBond)
        if bondsMade != todo:
            logger.debug("Made fewer bonds than expected in zip: %d -> %d", todo, bondsMade)
        self.analyse.stop('zip')
        return bondsMade

    def __getstate__(self):
        """Return a dict of objects we want to pickle.

        This is required as we can't pickle objects containing a logger as they have
        file handles (can remove by calling logging.shutdown() ) and lock objects (not removed
        by calling shutdown)."""
        # Return everything bar our logger
        d = dict(self.__dict__)
        d['logcsv'] = d['analyse'].logfile
        del d['analyse']
        if 'MDENGINE' in d:
            del d['MDENGINE'] # Contains a reference to the hood-blue logger
        return d

    def __setstate__(self, d):
        """Called when we are unpickled """
        self.__dict__.update(d)
        if 'logfile' in d:  # Hack for older versions with no logfile attribute
            logfile = ab_util.newFilename(d['logfile'])
        else:
            logfile = 'ambuild_1.log'
        self.setupLogging(logfile=logfile)
        if 'logcsv' in d:
            self.logcsv = ab_util.newFilename(d['logcsv'])
        else:
            self.logcsv = 'ambuild_1.csv'
        self._setupAnalyse(logfile=self.logcsv)
        return
