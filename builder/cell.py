'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import collections
import copy
import cPickle
import csv
import logging
import math
import os
import random
import sys
import time
import unittest
import xml.etree.ElementTree as ET

# External modules
import numpy

# Our modules
import buildingBlock
import fragment
import opt
import util

class Analyse():
    def __init__(self, cell, logfile="ambuild.csv"):

        self.fieldnames = ['step',
                           'type',
                           'tot_time',
                           'time',
                           'num_frags',
                           'num_particles',
                           'num_blocks',
                           'density',
                           'num_free_endGroups',
                           'potential_energy',
                           'num_tries',
                           'fragment_types',
                           'file_count',
                            ]
        self.cell = cell

        self.step = 0
        self._startTime = time.time()
        self._stepTime  = time.time()

        # Need to create an initial entry as we query the previous one for any data we don't have
        d = {}
        for f in self.fieldnames:
            if f == 'time':
                d[ f ] = self._startTime
            #elif f == 'file_count':
            elif f == 'type':
                d[ f ] = 'init'
            else:
                d[ f ] = 0

        self.last = d

        self.logfile = logfile
        self._logWriter = csv.DictWriter( open(self.logfile, 'w'), self.fieldnames )

        self._logWriter.writeheader()

        return

    def start(self):
        """Called whenever we start a step"""
        assert self._stepTime == None
        assert self.last
        self.step += 1
        self._stepTime = time.time()
        return

    def stop(self, stype, d={}):
        """Called at the end of a step with the data to write to the csv file"""

        new = {}

        for f in self.fieldnames:
            if f == 'type':
                new[ f ] = stype
            elif f == 'step':
                new [ f ] = self.step
            elif f == 'time':
                new[ f ] = time.time() - self._stepTime
            elif f == 'tot_time':
                new[ f ] = time.time() - self._startTime
                self._stepTime
            elif f == 'num_blocks':
                new[ f ] = self.cell.numBlocks()
            elif f == 'num_frags':
                new[ f ] = self.cell.numFragments()
            elif f == 'num_particles':
                new[ f ] = self.cell.numAtoms()
            elif f == 'num_free_endGroups':
                new[ f ] = self.cell.numFreeEndGroups()
            elif f == 'density':
                new[ f ] = self.cell.density()
            elif f == 'fragment_types':
                new[ f ] = str( self.cell.fragmentTypes() )
            elif f == 'file_count':
                new[ f ] = self.cell._fileCount
            elif f in d:
                new[ f ] = d[ f ]
            else:
                new[ f ] = self.last[ f ]

        self._logWriter.writerow( new )

        self.last = new
        self._stepTime = None
        self.start()
        return

class CellData(object):
    def __init__(self):

        self.cell                    = []

        self.atomTypes               = []
        self.bodies                  = []
        self.coords                  = []
        self.charges                 = []
        self.diameters               = []
        self.images                  = []
        self.masses                  = []
        self.symbols                 = []

        self.bonds                   = []
        self.bondLabels              = []
        self.angles                  = []
        self.angleLabels             = []
        self.propers                 = []
        self.properLabels            = []
        self.impropers               = []
        self.improperLabels          = []

        # for computing block/fragment enegies
        self.tagIndices              = []
#         self.blockBonds              = []
#         self.blockBondLabels         = []
#         self.blockBondAngles         = []
#         self.blockBondAngleLabels    = []
#         self.blockBondDihedrals      = []
#         self.blockBondDihedralLabels = []
        return

class Cell():
    '''
    classdocs
    '''

    BONDTYPESEP     = "-" # Character for separating bonds
    ENDGROUPSEP     = ":" # Character for separating endGroups in bonds

    def __init__( self, boxDim, atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15, doLog=False ):
        '''Construct an empty cell:

        Args:
        boxDim - a list with three numbers specifying the size of the cell A,B and C dimensions (angstroms)
        atomMargin - the additional distance that will be added on to the VdW radii of two atoms
                     to determine if they are close enough to clash.
        bondMargin - two atoms are considered close enough to bond if they are within the bond length
                     defined for the two atoms +/- the bondMargin
        bondAngleMargin - the tolerance (in degrees) from the ideal of 180 that defines an acceptable bond
        doLog - True/False - specifies if a log will be created - not recommended as it generates lots of data
                and slows the program.
        '''

        # A, B and C are cell vectors - for the time being we assume they are orthogonal
        # so they are just floats
        if not type(boxDim) is list and len(boxDim) == 3:
            raise RuntimeError,"Cell constructor needs the first argument to be list of 3 numbers setting the cell dimensions!"

        self.A = float(boxDim[0])
        self.B = float(boxDim[1])
        self.C = float(boxDim[2])

        assert self.A >0 and self.B > 0 and self.C > 0

        # For time being origin always 0,0,0
        self.origin = numpy.array( [0,0,0], dtype=numpy.float64 )

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
        self.targetEndGroups = 100 # number of free endgroups left

        # Dict mapping box key to a list of tuples defining the atoms in that box
        self.box1={}

        # Dict mapping key to list of boxes surrounding the keyed box
        self.box3={}

        # max atom radius - used to calculate box size
        self.boxSize = None
        self.maxAtomRadius=-1

        # see optimiseGeometry
        self.rCut = 5.0
        self.minCell = False
        self.minCellData = None

        # number of boxes in A,B,C axes - used for calculating PBC
        self.numBoxA = None
        self.numBoxB = None
        self.numBoxC = None

        self._fragmentLibrary = {} # fragmentType -> parentFragment
        self._endGroup2LibraryFragment = {} # endGroup type -> parentFragment
        self.bondTypes = []
        self._bondTable = {} # endGroup type -> set of endGroups it can bond to

        # dictionary mapping id of the block to the block - can't use a list and indices
        # as we add and remove blocks and need to keep track of them
        self.blocks = collections.OrderedDict()

        # Holds possible bond after checkMove is run
        self._possibleBonds = []

        # Logging functions
        self.logger = None
        self.setupLogging( doLog=doLog )

        # For analysis csv
        self._setupAnalyse()

        self._fileCount=0 # for naming output files

        return

    def addBlock( self, block, idxBlock=None ):
        """
        Add the block and put all atoms in their cells
        """

        # The id of the new block
        if idxBlock is None:
            idxBlock = block.id

        # Add to the dict
        self.blocks[ idxBlock ] = block

        # Each atom has its own cell - X atoms remain as Nones
        block.atomCell = [ None ] * block.numAtoms()

        #print "nbox ",self.numBoxA,self.numBoxB,self.numBoxC
        for idxCoord,coord in enumerate( block.iterCoord() ):

            # Periodic Boundaries
            x = coord[0] % self.A
            y = coord[1] % self.B
            z = coord[2] % self.C

            # Calculate which cell the atom is in
            a=int( math.floor( x / self.boxSize ) )
            b=int( math.floor( y / self.boxSize ) )
            c=int( math.floor( z / self.boxSize ) )
            key = (a,b,c)
            block.atomCell[ idxCoord ] = key

            try:
                # Add the atom to the cell
                self.box1[ key ].append( ( idxBlock, idxCoord ) )
            except KeyError:
                # Add the cell to the list and then add the atom
                self.box1[ key ] = [ ( idxBlock, idxCoord ) ]
                # Map the cells surrounding this one
                self.box3[ key ] = self.surroundCells( key )

        return idxBlock

    def addBondType( self, bondType ):
        """Allow bonds between the two endGroups specified in the bondType.

        endGroups are defined by the fragmentType they belong to (which is set by the fragmentType argument
        to libraryAddFragment), together with the identifier for that endGroup (which is specified by the first column
        in the .csv file). These are separated by a colon, so an endGroup identifier is of the form:

        FRAGMENT1:ENDGROUP1

        A bond is defined by two endGroups, separated by a hyphen, so a bond identifier has the form:

        FRAGMENT1:ENDGROUP1-FRAGMENT2:ENDGROUP2

        Args:
        bondType - a string specifying the two endGroups separated by a "-"
        """

        try:
            b1EndGroupType, b2EndGroupType = bondType.split(self.BONDTYPESEP)
            b1FragmentType = b1EndGroupType.split(self.ENDGROUPSEP)[0]
            b2FragmentType = b2EndGroupType.split(self.ENDGROUPSEP)[0]
        except ValueError:
            raise RuntimeError,"Error adding BondType {0} - string needs to be of form 'A:a-B:b'".format( bondType )

        # Checks
        assert b1FragmentType in self._fragmentLibrary,\
        "No fragment type {0} in fragmentLibrary".format( b1FragmentType )
        assert b1EndGroupType in self._fragmentLibrary[ b1FragmentType ].endGroupTypes(),\
        "Fragment type {0} has no endGroup type {1}".format( b1FragmentType, b1EndGroupType )
        assert b2FragmentType in self._fragmentLibrary,"No fragment type {0} in fragmentLibrary".format( b2FragmentType )
        assert b2EndGroupType in self._fragmentLibrary[ b2FragmentType ].endGroupTypes(),\
        "Fragment type {0} has no endGroup type {1}".format( b2FragmentType, b2EndGroupType )

        assert b1EndGroupType in self._endGroup2LibraryFragment.keys(),\
        "No endGroup type {0} in fragmentLibrary!".format(b1EndGroupType)
        assert b2EndGroupType in self._endGroup2LibraryFragment.keys(),\
        "No endGroup type {0} in fragmentLibrary!".format(b2EndGroupType)

        bt = (b1EndGroupType, b2EndGroupType)
        if bt in self.bondTypes:
            raise RuntimeError,"Adding an existing bond type: {0}".format( bt )

        self.bondTypes.append( bt )

        # Now recalculate the bond table
        self._updateBondTable()

        return

    def _updateBondTable(self):
        """Recalculate the dicionary of what endGroups can bond to which"""
        self._bondTable = {}
        for bondA, bondB in self.bondTypes:
            if bondA not in self._bondTable:
                self._bondTable[ bondA ] = set()
            if bondB not in self._bondTable:
                self._bondTable[ bondB ] = set()

            self._bondTable[ bondA ].add( bondB )
            self._bondTable[ bondB ].add( bondA )
        return

    def libraryAddFragment( self, filename, fragmentType='A' ):
        """Add a fragment of type fragmentType defined in the .car file filename

        Args:
        filename - the path to the .car file. There will need to be a corresponding .csv file that defines
                 - the endGroups, capAtoms etc.
        fragmentType - a name that will be used to identify the fragment - cannot contain the ":" character
        """

        if self._fragmentLibrary.has_key( fragmentType ):
            raise RuntimeError,"Adding existing ftype {0} again!".format( fragmentType )

        # For now don't allow adding blocks when the cell has blocks in
        assert not len(self.blocks),"Cannot add library fragments with populated cell!"

        # Make sure the type is valid
        if ":" in fragmentType or "-" in fragmentType:
            raise RuntimeError,"fragmentType cannot contain - or : characters!"

        # Create fragment
        f = fragment.Fragment( filename, fragmentType )

        # Update cell parameters for this fragment
        maxAtomRadius = f.maxAtomRadius()
        if  maxAtomRadius > self.maxAtomRadius:
            self.updateCellSize( maxAtomRadius=maxAtomRadius )

        # Add to _fragmentLibrary
        self._fragmentLibrary[ fragmentType ] = f

        # create dictionary keyed by endGroup types
        for ft in f.endGroupTypes():
            assert ft not in self._endGroup2LibraryFragment,"Adding existing endGroup type to library: {0}".format(ft)
            self._endGroup2LibraryFragment[ ft ] = fragmentType

        return

    def angle(self, c1, c2, c3):
        """Return the angle in radians c1---c2---c3
        where c are the coordinates in a numpy array"""
        return util.angle(c1, c2, c3, cell=numpy.array([self.A,self.B,self.C]))

    def attachBlock(self, growEndGroup, staticEndGroup, dihedral=None ):
        """
        Position growBlock so it can bond to blockS, using the given _endGroups

        Arguments:

        We take responsibility for adding and removing the growBlock from the cell on
        success or failure
        """

        growBlock      = growEndGroup.block()
        staticBlock    = staticEndGroup.block()
        idxStaticBlock = staticBlock.id

        staticBlock.positionGrowBlock( staticEndGroup, growEndGroup, dihedral=dihedral )

        # Now add growBlock to the cell so we can check for clashes
        blockId = self.addBlock( growBlock )

        #print "attachBlock got ",blockId, self.blocks
        #self.dump()
        #sys.exit()

        self.logger.debug("GOT {0} {1}".format(staticEndGroup, growEndGroup))

        # Check it doesn't clash
        if self.checkMove( blockId ) and self.processBonds() > 0:
            self.logger.debug("attachBlock first checkMove returned True")
            return True
        else:
            self.logger.debug("attachBlock first checkMove returned False")
            #self.writeCml("foo.cml", periodic=True, pruneBonds=False)
            #sys.exit()


        # Only attempt rotation if we're not worried about the dihedral
        # NOTE! - should only bother with the rotation if the cell is relatively crowded - otherwise
        # the multiple rotations are the slowest step
        if not dihedral:
            # Didn't work so try rotating the growBlock about the bond to see if that lets it fit

            # Determine the axis and center to rotate about
            blockEndGroup = growBlock.coord( growEndGroup.endGroupIdx() )
            blockS = self.blocks[ idxStaticBlock ]
            blockSEndGroup = blockS.coord( staticEndGroup.endGroupIdx() )
            axis = blockSEndGroup - blockEndGroup
            center = blockSEndGroup

            #step = math.pi/18 # 10 degree increments
            step = math.pi/9 # 20 degree increments
            for angle in util.frange( step, math.pi*2, step):

                #print "attachBlock rotating as clash: {}".format(angle*util.RADIANS2DEGREES)

                # remove the growBlock from the cell
                self.delBlock(blockId)

                # rotate by increment
                self.logger.debug("attachBlock rotating by {0} to {1}".format( math.degrees(step),
                                                                               math.degrees(angle) ) )
                growBlock.rotate( axis, angle, center=center)

                # add it and check
                self.addBlock(growBlock)

                if self.checkMove( blockId ) and self.processBonds() > 0:
                    self.logger.debug("attachBlock rotate worked")
                    return True

        # remove the growBlock from the cell
        self.delBlock(blockId)

        return False

    def bondAllowed( self, endGroup1, endGroup2 ):
        """Check if the given bond is permitted from the types of the two fragments
        """

        #    self.logger.debug( "checking bondAllowed {0} {1}".format( ftype1, ftype2 ) )
        if endGroup1.fragmentType() == "cap" or endGroup2.fragmentType() == "cap":
            assert False,"NEED TO FIX CAPS!"
            return True

        # Get the two fragment types for this bond
        eg1 = endGroup1.type()
        eg2 = endGroup2.type()
        for type1, type2 in self.bondTypes:
            #self.logger.debug("bondTypes {0}".format((type1, type2)) )
            if ( eg1 == type1 and eg2 == type2 ) or ( eg1 == type2 and eg2 == type1 ):
                #self.logger.debug( "BOND RETURN TRUE" )
                return True

        #self.logger.debug( "BOND RETURN FALSE" )
        return False

    def bondBlock(self, bond ):
        """ Bond the second block1 to the first and update the data structures
        """
        self.logger.debug( "cell bondBlock: {0}".format( bond ) )
        #self.logger.debug("before bond: {0} - {1}".format( bond.idxBlock1, bond.block1._bondObjects) )

        # We need to remove the block even if we are bonding to self as we need to recalculate the atomCell list
        self.delBlock( bond.endGroup1.block().id )
        if bond.endGroup1.block() != bond.endGroup2.block():
            self.delBlock( bond.endGroup2.block().id )
        else:
            self.logger.info("self-bonded block1: {0}".format( bond ) )

        #bond.block1.bondBlock( bond )
        bond.endGroup1.block().bondBlock( bond )
        #self.logger.debug("after bond: {0} - {1}".format( idxBlock1, block1._bondObjects) )
        self.addBlock( bond.endGroup1.block() )
        return

    def canBond( self,
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
        bond_length = util.bondLength( addBlock.symbol( idxAddAtom ), staticBlock.symbol( idxStaticAtom ) )
        if bond_length < 0:
            raise RuntimeError,"Missing bond distance for: {0}-{1}".format( addBlock.symbol( idxAddAtom ),
                                                                            staticBlock.symbol( idxStaticAtom ) )

        # See if the distance between them is acceptable
        #print "CHECKING BOND ATOMS ",bond_length,self.distance( addCoord, staticCoord )
        if  not ( max( 0.1, bond_length - bondMargin ) < distance < bond_length + bondMargin ):
            self.logger.debug( "Cannot bond due to distance: {0}".format( distance ) )
            return False

        # Now loop over all endGroups seeing if any of the angles are satisfied
        for staticEndGroup in staticBlock.atomEndGroups( idxStaticAtom ):
            for addEndGroup in addBlock.atomEndGroups( idxAddAtom ):

                # EndGroups in the same fragment can never bond
                if staticEndGroup.fragment == addEndGroup.fragment:
                    break

                assert staticEndGroup.free() and addEndGroup.free()

                # We need to check that we've not already added these endGroups as possible bonds
                if self._endGroupsInPossibleBonds( [ staticEndGroup, addEndGroup ] ):
                    continue

                # First check if endGroups of this type can bond
                if not self.bondAllowed( staticEndGroup, addEndGroup ):
                    self.logger.debug( "Bond disallowed by bonding rules: {0} : {1}".format( staticEndGroup,
                                                                                              addEndGroup
                                                                                             ) )
                    continue

                #print "Possible bond for {0} {1} {2} dist {3}".format( idxAddAtom,
                #                                                       idxStaticBlock,
                #                                                       idxStaticAtom,
                #                                                       self.distance( addCoord, staticCoord ) )
                addCoord      = addBlock.coord( idxAddAtom )
                staticCoord   = staticBlock.coord( idxStaticAtom )
                addCapAtom    = addBlock.coord( addEndGroup.capIdx() )
                angle1        = self.angle( addCapAtom, addCoord, staticCoord )
                staticCapAtom = staticBlock.coord( staticEndGroup.capIdx() )
                angle2        = self.angle( staticCapAtom, staticCoord, addCoord )

                #print "CHECKING ANGLE BETWEEN: {0} | {1} | {2}".format( angleAtom, addCoord, staticCoord )
                #print "{} < {} < {}".format( (self.bondAngle-bondAngleMargin) * util.RADIANS2DEGREES,
                #                             angle * util.RADIANS2DEGREES,
                #                             (self.bondAngle+bondAngleMargin) * util.RADIANS2DEGREES  )

                # Check if atoms are in line (zero degrees) within margin
                if not ( ( 0.0-bondAngleMargin < angle1 < 0.0+bondAngleMargin and \
                           0.0-bondAngleMargin < angle2 < 0.0+bondAngleMargin ) ):
                    self.logger.debug( "Cannot bond due to angles: {0} : {1}".format( math.degrees(angle1),
                                                                                     math.degrees(angle2) ) )
                    continue

                self.logger.debug( "Acceptable bond with angles: {0} : {1} | distance {2}".format( math.degrees(angle1),
                                                                                    math.degrees(angle2),
                                                                                    distance ) )

                # Create bond object and set the parameters
                bond = buildingBlock.Bond(staticEndGroup,addEndGroup)
                self._possibleBonds.append( bond )
                self.logger.debug( "canBond returning True with bonds: {0}".format( [str(b) for b in self._possibleBonds] ) )
                return True

        self.logger.debug( "canBond returning False" )
        return False

    def capBlocks(self, fragmentType=None, filename=None ):

        # Create the cap block
        capBlock = buildingBlock.Block( filePath=filename, fragmentType='cap' )

        # The endgroup is always the first only endGroup
        capEndGroup = capBlock.freeEndGroups()[0]

        # Get a list of blocks - need to do it this way or else we are changing the list of blocks while
        # we iterate through them
        cblocks = self.blocks.copy()
        for blockId, block in cblocks.iteritems():

            if block.numFreeEndGroups() == 0:
                self.logger.info("capBlocks block {0} already has no free endGroups".format( blockId ) )
                continue
            else:
                self.logger.info("capBlocks capping {0} endGroups of block {1}".format( block.numFreeEndGroups(),
                                                                                        blockId  )  )

            # If there are no free this won't loop
            for endGroup in block.freeEndGroupsFromTypes( [ fragmentType ] ):

                cblock = capBlock.copy()
                #print "ADDING CAPBLOCK ",id(cblock)
                block.positionGrowBlock( endGroup, capEndGroup )

                idxBlock = self.addBlock( cblock )

                # Test for Clashes with other blocks
                if self.checkMove( idxBlock ) and self.processBonds() > 0:
                    self.logger.info("Capped block {0} endGroup {1}".format( blockId, endGroup ) )
                else:
                    self.logger.critical("Failed to cap block {0} endGroup {1}".format( blockId, endGroup ) )
                    # Unsuccessful so remove the block from cell
                    self.delBlock(idxBlock)

        return

    def checkMove( self, idxAddBlock ):
        clashing = self._checkMove( idxAddBlock )
        if clashing > 0:
            return False

        return True

    def _checkMove( self, idxAddBlock ):
        """
        get pairs of close atoms

        bond=False
        for each pair:
            if are bonding atoms:
                if distance acceptable:
                    if angle acceptable and bond allowed:
                        bond=True
                        add to possible bonds

            if not bond:
                # Add atoms to list clash atoms

        if possible bonds:
            remove atoms directly bonded to bond atoms from clash atoms

        if there are still clash atoms:
            return False

        if possible bonds:
            process the bonds in order

        """

        # Get a list of the close atoms
        close = self.closeAtoms( idxAddBlock )
        #if close:
        #    print "GOT {0} CLOSE ATOMS ".format( len(close) )

        if not close:
            self.logger.debug("_checkMove no close contacts")
            # Nothing to see so move along
            return 0

        # Get the block1
        addBlock = self.blocks[ idxAddBlock ]

        clashAtoms = []
        self._possibleBonds = []
        for ( idxAddAtom, staticBlock, idxStaticAtom, distance ) in close:

#             print "CHECKING  ATOMS {}:{} {} -> {}:{} {}= {}".format( idxAddBlock,
#                                                               idxAddAtom,
#                                                               addSymbol,
#                                                               idxStaticBlock,
#                                                               idxStaticAtom,
#                                                               staticBlock.atomSymbol( idxStaticAtom ),
#                                                               distance )

            if not ( addBlock.isEndGroup( idxAddAtom ) and \
                     staticBlock.isEndGroup( idxStaticAtom ) and \
                     self.canBond( staticBlock,
                                   idxStaticAtom,
                                   addBlock,
                                   idxAddAtom,
                                   distance,
                                   self.bondMargin,
                                   self.bondAngleMargin,
                                                        )
                    ):

                # No bond so just check if the two atoms are close enough for a clash
                if distance <= addBlock.radius( idxAddAtom ) + staticBlock.radius( idxStaticAtom ) + self.atomMargin:
                    #print "CLASH {}->{} = {} < {}".format( addCoord,staticCoord, d, l  )
                    clashAtoms.append( ( staticBlock, idxStaticAtom, addBlock, idxAddAtom ) )

        # Now have list of possible bonds and clashes

        # Nothing so return True
        if not len( self._possibleBonds ) and not len( clashAtoms ):
            self.logger.debug( "NO BONDS AND NO ATOMS" )
            return 0

        # no bonds but a clash - return False
        if not len( self._possibleBonds ) and len( clashAtoms ):
            self.logger.debug( "No bonds and clashing atoms {0}".format( clashAtoms ) )
            return len( clashAtoms )

        s = ""
        for b in self._possibleBonds:
            s += str(b)+ " | "
        self.logger.debug( "Bonds {0} and clashing atoms {1}".format( s, clashAtoms ) )

        addBlock    = None
        staticBlock = None

        for bond in self._possibleBonds:

            b1Block = bond.endGroup1.block()
            b2Block = bond.endGroup2.block()

            # We need to remove any clashes with the cap atoms - the added block isn't bonded
            # so the cap atoms aren't excluded
            b1Cap      = bond.endGroup1.capIdx()
            b2Cap      = bond.endGroup2.capIdx()
            b1EndGroup = bond.endGroup1.endGroupIdx()
            b2EndGroup = bond.endGroup2.endGroupIdx()

            # Also need to remove any clashes of the endGroups with atoms directly bonded to the
            # opposite endGroup
            b1BondAtoms = b1Block.atomBonded1( b1EndGroup )
            b1BondAtoms.remove( b1Cap )
            b2BondAtoms    = b2Block.atomBonded1( b2EndGroup )
            b2BondAtoms.remove( b2Cap )

            toGo = [] # Need to remember indices as we can't remove from a list while we cycle through it
            for i, ( cellBlock, idxCellAtom, addBlock, idxAddAtom) in enumerate(clashAtoms):
                #self.logger.debug("CHECKING {0} {1} {2} {3}".format( idxAddBlock, idxAddAtom, idxStaticBlock, idxStaticAtom ) )
                # Remove any clashes with the cap atoms
                if ( b1Block == cellBlock and idxCellAtom == b1Cap ) or \
                   ( b2Block == addBlock  and idxAddAtom  == b2Cap ):
                    #self.logger.debug("REMOVING {0} {1} {2} {3}".format( idxAddBlock, idxAddAtom, idxStaticBlock, idxStaticAtom ) )
                    toGo.append( i )
                    continue

                # remove any clashes with directly bonded atoms
                if (b1Block == cellBlock and idxCellAtom == b1EndGroup and idxAddAtom  in b2BondAtoms ) or \
                   (b2Block == b2Block   and idxAddAtom  == b2EndGroup and idxCellAtom in b1BondAtoms ):
                    #self.logger.info( "REMOVING BOND ATOMS FROM CLASH TEST" )
                    toGo.append( i )

            # End of loop so now remove the items
            for i, idx in enumerate(toGo):
                # Need to subtract i as otherwise the indexing is out
                clashAtoms.pop( idx - i )

        # If there are any clashing atoms remaining this move failed
        if len( clashAtoms ):
            self.logger.debug( "Got clash atoms {0}".format( clashAtoms ) )
            return len( clashAtoms )

        # Got bonds and no clashes
        self.logger.debug( "Checkmove no clashes" )
        return 0

    def closeAtoms(self, idxBlock1):
        """
        Find all atoms that are close to the atoms in the given block.

        This version is optimised to use numpy to do the PBC distance calculations, so builds up a vector of points and then
        calculates the distances in one go with the distanceV function

        Args:
        idxBlock1: index of the block in self.blocks

        Returns:
        a list of tuples: (thisAtomIndex, otherBlockIndex, otherAtomIndex, distance) or None if there we no close contacts
        """

        # Build up list of the coordinates and an initial contacts array
        c1 = []
        c2 = []
        allContacts=[]

        count=0
        block1=self.blocks[ idxBlock1 ]
        for idxAtom1,coord1 in enumerate( block1.iterCoord() ):

            # Get the box this atom is in
            key = block1.atomCell[ idxAtom1 ]

            # Skip checking dummy atoms
            if key is None:
                continue

            #print "Close checking [{}] {}: {} : {}".format(key,idxBlock1,idxAtom1,coord1)

            # Get a list of the boxes surrounding this one
            surrounding = self.box3[key]

            #For each box loop through all its atoms chekcking for clashes
            for sbox in surrounding:

                #print "KEY ",i,sbox
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

                    c1.append( coord1 )
                    c2.append( block2.coord( idxAtom2 ) )
                    allContacts.append( ( idxAtom1, block2, idxAtom2 ) )
                    count += 1

        # Have list so now calculate distances
        if count == 0:
            return []

        # Calculate array of distances for all coordinates
        distances = self.distance(c1, c2)

        # prune contacts array according to distances
        return [ ( c[0], c[1], c[2], distances[i] ) for i, c in enumerate( allContacts ) if distances[i] < self.boxSize ]

    def cellEndGroupPair(self, cellEndGroups=None):
        """Return two free endGroups from two different blocks in the cell"""

        #
        # Get a list of available free endGroups in the cell
        #
        endGroupTypes2Block = self.endGroupTypes2Block()
        if len(endGroupTypes2Block.keys()) == 0:
            raise RuntimeError,"No available endGroups in the cell"
        #
        # If the user supplied a list of cellEndGroups we use this to determine what can bond -
        # other wise we use all available endGroups
        allTypes = endGroupTypes2Block.keys() # just for informing the user
        if cellEndGroups is not None:
            if isinstance(cellEndGroups,str):
                cellEndGroups = [ cellEndGroups ]
            cellEndGroups = set( cellEndGroups )
        else:
            cellEndGroups = set(endGroupTypes2Block.keys())
        #
        # We create a dictionary mapping which cell endGroups can bond to which
        # This is basically a truncated _bondTable with any endGroups removed that aren't present
        #
        cell2cell = copy.copy(self._bondTable)
        for eg in self._bondTable:
            if eg not in cellEndGroups:
                # First remove all references in the keys
                del cell2cell[ eg ]
            else:
                # The key is valid so remove any missing endGroups from the values
                ceg = cellEndGroups.intersection( cell2cell[ eg ] )
                # Need to check there are some other types that this can bond to, but also there is at least
                # one other block if both are of the same type
                # NB the strange next/iter thing is the only way to get an item from the set
                if len( ceg ) == 0 or \
                   (len( ceg ) == 1 and eg == next(iter(ceg)) and len( endGroupTypes2Block[ eg ] ) < 2 ):
                    # No valid blocks that could bond with this endGroup so delete the entire key
                    del cell2cell[ eg ]
                else:
                    cell2cell[ eg ] = ceg

        if len(cell2cell.keys()) == 0:
            raise RuntimeError,"No endGroups of types {0} are available to bond from {1}".format( cellEndGroups, allTypes )

        # At this point we assume we can definitely bond at least 2 blocks
        self.logger.debug("cellEndGroupPair got cell/library endGroups: {0}".format(cell2cell) )

        # Select a random block/endGroup from the list
        eg1Type = random.choice( cell2cell.keys() )
        block1 = random.sample( endGroupTypes2Block[ eg1Type ], 1 )[0]
        endGroup1 = block1.randomEndGroup( endGroupTypes=[eg1Type] )

        # Pick a random endGroup type that can bond to this
        eg2Type = random.sample( cell2cell[ eg1Type ], 1 )[0]

        # Select a random block/endGroup of that type
        # (REM: need to remove the first block from the list of possibles hence the difference thing
        block2 = random.sample( endGroupTypes2Block[ eg2Type ].difference(set([block1])), 1 )[0]
        endGroup2 = block2.randomEndGroup( endGroupTypes=[eg2Type] )

        self.logger.debug("cellEndGroupPair returning: {0} {1}".format(endGroup1.type(),endGroup2.type()) )
        return endGroup1, endGroup2

    def dataDict( self, rigidBody=True, periodic=True, center=False, fragmentType=None ):

        # Object to hold the cell data
        d = CellData()

        d.cell = [self.A,self.B,self.C]

        if fragmentType is not None:
            # Only returning data for one type of fragment
            assert fragmentType in self.fragmentTypes(),"FragmentType {0} not in cell!".format(fragmentType)
            atomCount=0
            for b in self.blocks.values():
                coords, symbols, bonds = b.dataByFragment(fragmentType)
                d.coords += coords
                d.symbols += symbols
                d.bonds += [ (b1+atomCount, b2+atomCount) for (b1, b2) in bonds ]
                atomCount += len(coords)

            return d

        atomCount = 0 # Global count in cell
        bodyCount = -1

        #for idxBlock, block in enumerate(self.blocks.itervalues()):
        for idxBlock, block in self.blocks.iteritems():

            # Do bonds first as counting starts from atomCount and goes up through the blocks
            if not rigidBody:
                # add all bonds, angles and dihederals throughout the whole block
                # Add all bonds
                d.bonds += [ (a1+atomCount, a2+atomCount) for a1, a2 in block.bonds() ]
                d.bondLabels += [ "{0}-{1}".format(block.type(a1),block.type(a2)) for a1, a2 in block.bonds() ]
                angles, propers, impropers = block.anglesAndDihedrals()
                # Add all angles
                d.angles += [ (a1+atomCount, a2+atomCount, a3+atomCount) for a1, a2, a3 in angles ]
                d.angleLabels += [ "{0}-{1}-{2}".format(block.type(a1),block.type(a2),block.type(a3)) for a1, a2, a3 in angles ]
                # Add all propers
                d.propers += [ (a1+atomCount, a2+atomCount, a3+atomCount, a4+atomCount) \
                              for a1, a2, a3, a4 in propers ]
                d.properLabels += [ "{0}-{1}-{2}-{3}".format(block.type(a1),
                                                             block.type(a2),
                                                             block.type(a3),
                                                             block.type(a4)
                                                             ) for a1, a2, a3, a4 in propers ]
                # Add all impropers
                d.impropers += [ (a1+atomCount, a2+atomCount, a3+atomCount, a4+atomCount) \
                                for a1, a2, a3, a4 in impropers ]
                d.improperLabels += [ "{0}-{1}-{2}-{3}".format(block.type(a1),
                                                             block.type(a2),
                                                             block.type(a3),
                                                             block.type(a4)
                                                             ) for a1, a2, a3, a4 in impropers ]

            else:
                # Just add the bonds between blocks. Also add angles for all atoms connected to the bonds
                # we do this so that we can exclude them from VdW interactions in MD codes
                for b1, b2 in block.blockBonds():

                    # The bonds themselves
                    d.bonds.append( (b1+atomCount,b2+atomCount) )
                    d.bondLabels.append("{0}-{1}".format( block.type(b1),block.type(b2) ) )

                    _angles=set()
                    # Atoms connected to the endGroup that we need to specify as connected so we add as angles
                    for batom in block.atomBonded1( b1 ):
                        # The opposite endGroup is included in the list bonded to an endGroup so skip
                        if batom == b2:
                            continue
                        _angles.add( ( batom, b1, b2 ) )

                    for batom in block.atomBonded1( b2 ):
                        # The opposite endGroup is included in the list bonded to an endGroup so skip
                        if batom == b1:
                            continue
                        _angles.add( ( b1, b2, batom ) )

                    # Add to overall lists
                    for a1, a2, a3 in _angles:
                        d.angles.append( (a1+atomCount, a2+atomCount, a3+atomCount) )
                        l = "{0}-{1}-{2}".format( block.type( a1 ),
                                                  block.type( a2 ),
                                                  block.type( a3 ) )
                        d.angleLabels.append( l )

                    #
                    # Dihedrals
                    #
                    for dindices in block.dihedrals( b1, b2 ):
                        dihedral = ( dindices[0] + atomCount,
                                     dindices[1] + atomCount,
                                     dindices[2] + atomCount,
                                     dindices[3] + atomCount )
                        d.propers.append( dihedral )
                        dlabel = "{0}-{1}-{2}-{3}".format( block.type( dindices[0] ),
                                                           block.type( dindices[1] ),
                                                           block.type( dindices[2] ),
                                                           block.type( dindices[3] )
                                                        )
                        d.properLabels.append( dlabel )

            # Now loop through fragments and coordinates
            #for i, coord in enumerate(block.iterCoord() ):
            for idxFrag,frag in enumerate(block._fragments): # need index of fragment in block

                # Body count always increments with fragment although it may go up within a fragment too
                bodyCount += 1

                lastBody=frag.body(0)
                for i,coord in enumerate(frag.iterCoord()):
                    if periodic:
                        x, ix = util.wrapCoord( coord[0], self.A, center=center )
                        y, iy = util.wrapCoord( coord[1], self.B, center=center )
                        z, iz = util.wrapCoord( coord[2], self.C, center=center )
                        d.coords.append( numpy.array([x,y,z]))
                        d.images.append([ix,iy,iz])
                    else:
                        d.coords.append(coord)
                        d.images.append([0,0,0])

                    d.atomTypes.append(frag.type(i))
                    d.charges.append(frag.charge(i))
                    d.diameters.append(0.1)
                    d.masses.append(frag.mass(i))
                    d.symbols.append(frag.symbol(i))

                    # Work out which body this is in
                    b = frag.body(i)
                    if b != lastBody:
                        bodyCount +=1
                        lastBody = b
                    d.bodies.append( bodyCount )

                    # Increment global atom count
                    atomCount += 1

                # Work out which fragment this is in
                # REM blockCount is NOT idxBlock in the dict - need to rationalise this.
                d.tagIndices.append((idxBlock,idxFrag,atomCount-i,atomCount))

            #END loop through fragments
            # Only set the count here
            #atomCount += blockAtomCount

            # End block loop

        return d

    def delBlock(self,blockId):
        """
        Remove the block with the given index from the cell
        """

        #print "Deleting block: {} from {}".format(blockId,self.blocks.keys())
        block =  self.blocks[ blockId ]

        # Remove each atom from the list
        keys = []
        for iatom, key in enumerate( block.atomCell ):

            # Skip dummy atoms
            if key == None:
                continue

            #print "removing ",blockId,iatom,key
            keys.append(key)
            #print "B4 ",self.box1[key]
            self.box1[key].remove( ( blockId, iatom ) )
            #print "AFTER ",self.box1[key]

        # Now remove any empty keys and corresponding surrounding boxes
        # Might think about keeping the surrounding boxes as they could be use by other atoms?
        #print "now checking keys"
        for key in keys:
            if self.box1.has_key(key) and len( self.box1[key] ) == 0:
                del self.box1[key]
                del self.box3[key]
        del self.blocks[blockId]
        del block

        return

    def density(self):
        """The density of the cell"""
        d = ( sum( [ b.blockMass() for b in self.blocks.itervalues() ] ) / ( self.A * self.B * self.C ) )
        return d * (10/6.022)

    def dihedral(self, p1, p2, p3, p4):
        return util.dihedral(p1, p2, p3, p4,cell=[self.A,self.B,self.C])

    def distance(self, v1, v2, cell=None ):
        """Distance with numpy taking PBC into account
        This works either with 2 points or a vector of any number of points
        Adapted from: http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co
        Changed so that it can cope with distances across more than one cell
        """
        return util.distance(v1, v2, cell=numpy.array([self.A,self.B,self.C]))

    def dump(self, prefix="step", addCount=True ):
        """Write out our current state"""

        if addCount:
            self._fileCount+=1
            prefix=prefix+"_{0}".format(self._fileCount)

        #self.writeXyz(prefix+".xyz",data=data, periodic=False)
        #self.writeXyz(prefix+"_P.xyz",data=data, periodic=True)
        #self.writeCar(prefix+"_P.car",data=data,periodic=True)
        #self.writeCml(prefix+"_PV.cml", data=data, allBonds=True, periodic=True, pruneBonds=True)
        #self.writeCml(prefix+".cml", data=data, allBonds=True, periodic=False, pruneBonds=False)

        pklFile=os.path.abspath(prefix+".pkl")
        self.writePickle(pklFile)
        return pklFile

    def _endGroupsInPossibleBonds(self, endGroups ):
        """Check if any of the endGroups are already in the list of possible bonds"""
        eg = set()
        for b in self._possibleBonds:
            eg.update( [ b.endGroup1, b.endGroup2 ] )
        return bool( eg.intersection( frozenset( endGroups ) ) )

    def endGroupTypes2Block(self ):
        """Return a dictionary mapping free endGroup types to a list of the blocks

        We don't check if any are available just return an empty dictionary if not
        """
        endGroupTypes2Block = {}
        for block in self.blocks.itervalues():
            #numFreeEndGroups += block.numFreeEndGroups()
            # Get a list of the free endGroup types in the block
            for endGroupType in block.freeEndGroupTypes():
                if endGroupType not in endGroupTypes2Block:
                    endGroupTypes2Block[ endGroupType ] = set()
                endGroupTypes2Block[ endGroupType ].add( block )
        return endGroupTypes2Block

    def fragMaxEnergy(self,
                  rigidBody=True,
                  doDihedral=False,
                  doImproper=False,
                  xmlFilename="hoomdCalc.xml",
                  **kw
                  ):

        # Loop through all blocks, fragments and atoms creating the labels
        # do this in dataDict - create list of labels and start:stop indices
        data = self.dataDict(periodic=True, center=True, rigidBody=rigidBody)

        print "GOT TAGS ",data.tagIndices

        # in hoomdblue. loopt through list of labels and create groups and computes for each one
        # opt must hold the list of groups and computes
        o = opt.HoomdOptimiser()
        if 'rCut' in kw:
            self.rCut = kw['rCut']
        else:
            self.rCut = o.rCut

        self.logger.info( "Running fragMaxEnergy" )

        maxe,idxBlock,idxFragment = o.fragMaxEnergy(data,
                                                    xmlFilename,
                                                    rigidBody=rigidBody,
                                                    doDihedral=doDihedral,
                                                    doImproper=doImproper,
                                                    **kw )

        print "GOT ",maxe,idxBlock,idxFragment

        return

    def fragmentTypes(self):
        ft = {}
        for b in self.blocks.itervalues():
            d = b.fragmentTypeDict()
            for k,v in d.iteritems():
                try:
                    ft[ k ] += v
                except:
                    ft[ k ] = v
        return ft

    def fragmentTypeFromEndGroupType(self, endGroupType):
        return endGroupType.split( self.ENDGROUPSEP)[0]

    def fromHoomdblueSystem(self, system ):
        """Reset the particle positions from hoomdblue system"""        # Should really check HOOMD version but...
        if hasattr(system.box,"Lx"):
            Lx = system.box.Lx
            Ly = system.box.Ly
            Lz = system.box.Lz
        else:
            Lx = system.box[0]
            Ly = system.box[1]
            Lz = system.box[2]

        if self.minCell:
            # Need to take unwrapped coords and put back into
            # the original cell
            assert self.minCellData

        # Read back in the particle positions
        atomCount=0
        for block in self.blocks.itervalues():
            for k in range( block.numAtoms() ):

                p = system.particles[ atomCount ]
                xt, yt, zt  = p.position
                ix, iy, iz = p.image

                if self.minCell:

                    assert False,"FIX MINCELL"
                    assert self.minCellData['A'] == Lx
                    assert self.minCellData['B'] == Ly
                    assert self.minCellData['C'] == Lz

                    # Unwrap the coordinates in the centered cell
                    x = util.unWrapCoord( xt, ix, Lx, centered=False )
                    y = util.unWrapCoord( yt, iy, Ly, centered=False )
                    z = util.unWrapCoord( zt, iz, Lz, centered=False )

                    # Need to take unwrapped coords and put back into
                    # the original cell
                    x = x + Lx/2 + self.minCellData['minA']
                    y = y + Ly/2 + self.minCellData['minB']
                    z = z + Lz/2 + self.minCellData['minC']

                    # Not the case as could be in a different image
                    #assert x >= 0 and x <= self.A
                    #assert y >= 0 and y <= self.B
                    #assert z >= 0 and z <= self.C
                else:
                    x = util.unWrapCoord( xt, ix, Lx, centered=True )
                    y = util.unWrapCoord( yt, iy, Ly, centered=True )
                    z = util.unWrapCoord( zt, iz, Lz, centered=True )

                #block.atomCoord( k )[0] = x
                #block.atomCoord( k )[1] = y
                #block.atomCoord( k )[2] = z
                block.coord(k, [x,y,z])

                atomCount += 1

        if atomCount != len( system.particles ):
            raise RuntimeError,"Read {0} positions but there were {1} particles!".format( atomCount, len( system.particles ) )

        # Now have the new coordinates, so we need to put the atoms in their new cells
        self.repopulateCells()

        return

    def getInitBlock( self, fragmentType=None ):
        """Return an initBlock"""

        if fragmentType is None:
            # Need to determine the fragmentType from the endGroupType
            fragmentType = random.choice( list( self._fragmentLibrary.keys() ) )

        # sanity check
        if not self._fragmentLibrary.has_key( fragmentType ):
            raise RuntimeError, "Asking for a non-existing initBlock type: {0}".format( fragmentType )

        # Copy the init fragment
        f = self._fragmentLibrary[ fragmentType ].copy()

        return buildingBlock.Block( initFragment=f )

    def growBlocks(self,
                   toGrow,
                   cellEndGroups=None,
                   libraryEndGroups=None,
                   endGroupType=None,
                   dihedral=None,
                   maxTries=50 ):
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

        self.logger.info( "Growing {0} new blocks".format( toGrow ) )

        assert len(self.blocks),"Need to seed blocks before growing!"

        if endGroupType:
            self.logger.warn('endGroupType is deprecated! Use cellEndGroups/libraryEndGroups instead')
            assert not cellEndGroups or libraryEndGroups
            libraryEndGroups = endGroupType

        if dihedral:
            # Convert dihedral to radians
            dihedral = math.radians( dihedral )

        added=0
        tries=0
        while added < toGrow:

            if tries >= maxTries:
                self.logger.critical("growBlocks - exceeded maxtries {0} when joining blocks!".format( maxTries ) )
                return added

            if self.numFreeEndGroups() == 0:
                self.logger.critical("growBlocks got no free endGroups!")
                return added

            # Select two random blocks that can be bonded
            #initBlock, newEG, idxStaticBlock, staticEG = self.randomInitAttachments( endGroupType=endGroupType )
            cellEndGroup, libraryEndGroup = self.libraryEndGroupPair(cellEndGroups=cellEndGroups,
                                                                     libraryEndGroups=libraryEndGroups)

            # Apply random rotation in 3 axes to randomise the orientation before we align
            libraryBlock = libraryEndGroup.block()
            libraryBlock.randomRotate( origin=self.origin )

            # Try and attach it
            ok =  self.attachBlock( libraryEndGroup, cellEndGroup, dihedral=dihedral )

            if ok:
                added +=1
                self.logger.info("growBlocks added block {0} after {1} tries.".format( added, tries ) )
                self.analyse.stop('grow', d={'num_tries':tries})
                #print "GOT BLOCK ",[ b  for b in self.blocks.itervalues() ][0]
                tries = 0
            else:
                # del initBlock?
                tries += 1

        self.logger.info("After growBlocks numBlocks: {0}".format( len(self.blocks) ) )

        return added

    def joinBlocks(self, toJoin, cellEndGroups=None, dihedral=None, maxTries=100 ):
        """
        Bond toJoin blocks together using the endGroup types specified in cellEndGroups

        Args:
        toJoin - number of blocks to join
        cellEndGroups - a list of the different endGroupTypes that should be bonded. If this is None
                        randomly chosen endGroups will be used.
        dihedral: the dihedral angle about the bond (3rd column in csv file)
        maxTries - the maximum number of moves to try when joining
        """

        self.logger.info( "Joining {0} new blocks".format( toJoin ) )

        if dihedral:
            # Convert dihedral to radians
            dihedral = math.radians( dihedral )

        added=0
        tries=0
        while added < toJoin:

            if len ( self.blocks ) == 1:
                self.logger.info("joinBlocks Hooray! - no more blocks to join!")
                return added

            if tries > maxTries:
                self.logger.critical("joinBlocks - exceeded maxtries when joining blocks!")
                return added

            if self.numFreeEndGroups() == 0:
                self.logger.critical("joinBlocks got no free endGroups!")
                return added

            # Select 2 random blocks that can be joined
            moveEndGroup, staticEndGroup = self.cellEndGroupPair(cellEndGroups=cellEndGroups)

            # Copy the original block so we can replace it if the join fails
            moveBlock = moveEndGroup.block()
            idxMoveBlock = moveBlock.id
            blockCopy = moveBlock.copy()

            # Remove from cell so we don't check against itself and pick a different out
            self.delBlock( idxMoveBlock )

            self.logger.debug( "joinBlocks calling attachBlock: {0} {1}".format( moveEndGroup, staticEndGroup ) )

            # now attach it
            ok = self.attachBlock( moveEndGroup, staticEndGroup, dihedral=dihedral )
            if ok:
                added+=1
                self.logger.info("joinBlocks joined block {0} after {1} tries.".format( added, tries ) )
                tries=0
            else:
                # Put the original block back in the cell
                self.addBlock( blockCopy )
                tries+=1

        self.logger.info("After joinBlocks numBlocks: {0}".format( len(self.blocks) ) )

        return added

    def libraryEndGroupPair(self, cellEndGroups=None, libraryEndGroups=None ):
        """Return a fee endGroup from the cell and one from the library that can be bonded to it."""

        #
        # Get a list of available free endGroups in the cell
        #
        endGroupTypes2Block = self.endGroupTypes2Block()
        if len(endGroupTypes2Block.keys()) == 0:
            raise RuntimeError,"No available endGroups in the cell"
        #
        # We create a dictionary mapping cell endGroups to possible libraryEndGroups
        #
        cell2Library = {}
        for ceg in endGroupTypes2Block.keys():
            # We add those endGroup types that can be bonded to that are also in the library
            leg = self._bondTable[ ceg ].intersection( self._endGroup2LibraryFragment.keys() )
            if len(leg) > 0:
                cell2Library[ ceg ] = leg

        # Check that there are some available
        if len(cell2Library.keys()) == 0:
            raise RuntimeError,"No library fragments available to bond under the given rules: {0}".format(endGroupTypes2Block.keys())

        #
        # If the user supplied a list of cellEndGroups we prune the list of cell endGroups to those that are in this list
        #
        if cellEndGroups is not None:
            # Convert to a set - make sure is a list first
            if isinstance(cellEndGroups, str):
                cellEndGroups = [ cellEndGroups ]
            for ceg in cell2Library.keys():
                if ceg not in cellEndGroups:
                    del cell2Library[ ceg ]

            if len(cell2Library.keys()) == 0:
                raise RuntimeError,"No free endGroups of types in cellEndGroups: {0} - {1}".format(cellEndGroups,
                                                                                                   endGroupTypes2Block.keys())

        #
        # If the user supplied a list of libraryEndGroups we remove any library endGroups that aren't in the list
        #
        if libraryEndGroups is not None:
            # Convert to a set - make sure is a list first
            if isinstance(libraryEndGroups, str):
                libraryEndGroups = [ libraryEndGroups ]
            libraryEndGroups = set( libraryEndGroups ) # Save old so we can warn user and also find matching
            tmp = cell2Library
            cell2Library = {}
            for ceg, leg in tmp.iteritems():
                # Only select those that are in the libraryEndGroups
                pleg = libraryEndGroups.intersection( leg )
                if len(pleg) > 0:
                    cell2Library[ ceg ] = pleg

            if not len(cell2Library.keys()):
                raise RuntimeError,"No library fragments of type {0} available to bond under the given rules: {0}".format(libraryEndGroups,endGroupTypes2Block.keys())

        self.logger.debug("libraryEndGroupPair got cell/library endGroups: {0}".format(cell2Library) )

        # Now we can pick a random endGroup from the cell, get the corresponding library group
        cellEgT = random.choice( cell2Library.keys() )

        # First get a block that contains this type of endGroup
        # Need to use sample as sets don't support random.choice
        cellBlock = random.sample( endGroupTypes2Block[ cellEgT ], 1 )[0]

        # Now select a random endGroup of that type from it
        cellEndGroup = cellBlock.randomEndGroup( endGroupTypes=[cellEgT] )

        # Now get a corresponding library endGroup
        # We need to pick a random one of the types that we can bond to that is also in libraryTypes
        libEgT = random.sample( cell2Library[ cellEgT ], 1)[0]

        # Now determine the fragmentType and create the block and fragment
        fragmentType = self._endGroup2LibraryFragment[ libEgT ]
        libraryBlock = self.getInitBlock( fragmentType=fragmentType )

        # now get the endGroup
        libraryEndGroup = libraryBlock.randomEndGroup( endGroupTypes=[libEgT] )

        # Return them both - phew!
        self.logger.debug("libraryEndGroupPair returning: {0} {1}".format(cellEndGroup.type(),libraryEndGroup.type()) )
        return cellEndGroup, libraryEndGroup

    def numBlocks(self):
        return len( self.blocks)

    def numFragments(self):
        return sum( [ len(b._fragments) for b in self.blocks.itervalues() ] )

    def numFreeEndGroups(self):
        return sum( [ b.numFreeEndGroups() for b in self.blocks.itervalues() ] )

    def numAtoms(self):
        return sum( [ b.numAtoms() for b in self.blocks.itervalues() ] )

    def optimiseGeometry(self,
                         rigidBody=True,
                         doDihedral=False,
                         doImproper=False,
                         doCharges=True,
                         xmlFilename="hoomdOpt.xml",
                         **kw ):
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

        self.logger.info( "Running optimisation" )

        # HACK
        minCell=False
        self.minCell = minCell

        if doDihedral and doImproper:
            raise RuntimeError,"Cannot have impropers and dihedrals at the same time"

        optimiser = opt.HoomdOptimiser()
        if 'rCut' in kw:
            self.rCut = kw['rCut']
        else:
            self.rCut = optimiser.rCut

        data = self.dataDict(periodic=True, center=True, rigidBody=rigidBody)
        d = {} # for printing results
        ok = optimiser.optimiseGeometry( data,
                                         xmlFilename=xmlFilename,
                                         rigidBody=rigidBody,
                                         doDihedral=doDihedral,
                                         doImproper=doImproper,
                                         doCharges=doCharges,
                                         d=d,
                                         **kw )
        self.analyse.stop('optimiseGeometry', d )

        if ok:
            self.logger.info( "Optimisation succeeded" )
            self.fromHoomdblueSystem( optimiser.system )
            return True
        else:
            self.logger.critical( "Optimisation Failed" )
            return False

    def positionInCell(self, block):
        """Make sure the given block is positioned within the cell"""

        bradius = block.blockRadius()

        oradius = bradius + ( 2 * self.boxMargin )

        # Check there is enough space
        if oradius >= self.A or oradius >= self.B or oradius >= self.C:
            raise RuntimeError, "Cannot fit block with radius {0} into cell [{1}, {2}, {3}]!".format( bradius,
                                                                                                      self.A,
                                                                                                      self.B,
                                                                                                      self.C
                                                                                                       )

        # get a range for x, y and z that would fit the block in the cell, pick random values within that
        # and stick the block there
        x = random.uniform( bradius+self.boxMargin, self.A-bradius-self.boxMargin )
        y = random.uniform( bradius+self.boxMargin, self.B-bradius-self.boxMargin )
        z = random.uniform( bradius+self.boxMargin, self.C-bradius-self.boxMargin )

        coord = numpy.array([x,y,z], dtype=numpy.float64 )

        block.translateCentroid( coord )

        self.logger.debug( "positionInCell block moved to: {0}".format( block.centroid() ) )

        return

    def processBonds(self):
        """Make any bonds that were found during checkMove
        return Number of bonds made

        Args:
        addedBlockIdx - the block that was added
        """

        if not len( self._possibleBonds ):
            self.logger.debug("processBonds got no bonds" )
            return 0

        #self.logger.debug = lambda x: sys.stdout.write(x + "\n")

        # Here no atoms clash and we have a list of possible bonds - so bond'em!
        self.logger.debug("processBonds got bonds: {0}".format( self._possibleBonds ) )
        self.logger.debug("processBonds blocks are: {0}".format( sorted( self.blocks ) ) )

        bondsMade = 0
        for count, bond in enumerate( self._possibleBonds ):
            # With the bonding rules some endGroups may become not free when other endGroups in that fragment
            # are involved in bonds so we need to make sure they are free before we do
            if bond.endGroup1.free() and bond.endGroup2.free():
                self.bondBlock( bond )
                self.logger.debug("Added bond: {0}".format( self._possibleBonds[count] ) )
                bondsMade += 1
                #self.logger.debug( "process bonds after bond self.blocks is  ",sorted(self.blocks) )

        # Clear any other possible bonds
        self._possibleBonds = []

        return bondsMade

    def randomBlockId( self, checkFree=True ):
        """Return numBlocks random block ids"""

        MAXCOUNT=100

        blockId = random.choice( list( self.blocks.keys() ) )

        if checkFree:
            if self.numFreeEndGroups() == 0:
                return False

            count = 0
            while True:
                count += 1
                if count > MAXCOUNT:
                    return False

                block = self.blocks[ blockId ]
                if block.numFreeEndGroups > 0:
                    break

                blockId = random.choice( list( self.blocks.keys() ) )

        return blockId

    def randomMoveBlock(self, block, margin=None ):
        """Randomly move the given block
         If buffer is given, use this as a buffer from the edges of the cell
         when selecting the coord
        """
        # Get _coords of random point in the cell that we will move the block to after
        # we've rotated at the origin
        if margin:
            x = random.uniform(margin,self.A-margin)
            y = random.uniform(margin,self.B-margin)
            z = random.uniform(margin,self.C-margin)
        else:
            x = random.uniform(0,self.A)
            y = random.uniform(0,self.B)
            z = random.uniform(0,self.C)

        coord = numpy.array([x,y,z], dtype=numpy.float64 )

        #print "Got random coord: {}".format(coord)

        # Move to origin, rotate there and then move to new coord
        # Use the cell axis definitions
        block.translateCentroid( self.origin )
        block.randomRotate( origin=self.origin, atOrigin=True )

        # Now move to new coord
        block.translateCentroid( coord )

        return

    def randomMoveAroundCenter(self, move_block, center, radius ):
        """
        Move the move_block to a random point so that the its centroid is withiin
        radius of the center
        """

        # Calculate new coord
        x = random.uniform(-radius,radius)
        y = random.uniform(-radius,radius)
        z = random.uniform(-radius,radius)
        xyz = numpy.array( [x,y,z], dtype=numpy.float64 )
        coord = numpy.add( center, xyz )

        move_block.translateCentroid( coord )

        move_block.randomRotateBlock( origin=self.origin )
        return

    def resizeCells( self, boxMargin ):
        # Change the cell size
        self.updateCellSize( boxMargin=boxMargin )
        self.repopulateCells()
        return

    def repopulateCells(self):
        """Add all the blocks to resized cells"""

        # Remove all blocks from their cells
        self.box1 = {}
        self.box3 = {}

        if len( self.blocks ):
            # Put all the blocks into the new cells
            # First copy the blocks dictionary
            blocks = copy.copy( self.blocks )

            self.blocks.clear() # Delete the old one and reset

            self.logger.debug("repopulateCells, adding blocks into new cells")
            for idxBlock, block in blocks.iteritems():
                self.addBlock( block, idxBlock=idxBlock )

            del blocks
        return

    def runMD(self,
              xmlFilename="hoomdMD.xml",
              rigidBody=True,
              doDihedral=False,
              doImproper=False,
              doCharges=True,
              **kw ):

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

        assert rigidBody,"FIX runMD FOR ALL ATOM!!"

        if doDihedral and doImproper:
            raise RuntimeError,"Cannot have impropers and dihedrals at the same time"

        optimiser = opt.HoomdOptimiser()
        if 'rCut' in kw:
            self.rCut = kw['rCut']
        else:
            self.rCut = optimiser.rCut

        d = {}
        data = self.dataDict(periodic=True, center=True, rigidBody=rigidBody)
        ok = optimiser.runMD(data,
                             xmlFilename=xmlFilename,
                             rigidBody=rigidBody,
                             doDihedral=doDihedral,
                             doImproper=doImproper,
                             doCharges=doCharges,
                             d=d,
                             **kw )

        self.analyse.stop('runMD', d )

        self.fromHoomdblueSystem( optimiser.system )

        return ok

    def runMDAndOptimise(self,
                         xmlFilename="hoomdMDOpt.xml",
                         rigidBody=True,
                         doDihedral=False,
                         doImproper=False,
                         doCharges=True,
                         **kw ):

        """Run an MD simulation followed by a Geometry optimisation.

        Args:
        See runMD and optimiseGeometry for acceptable arguments.
        """

        assert rigidBody,"FIX runMD FOR ALL ATOM!!"

        if doDihedral and doImproper:
            raise RuntimeError,"Cannot have impropers and dihedrals at the same time"

        optimiser = opt.HoomdOptimiser()
        if 'rCut' in kw:
            self.rCut = kw['rCut']
        else:
            self.rCut = optimiser.rCut

        d = {}
        data = self.dataDict(periodic=True, center=True, rigidBody=rigidBody)
        ok = optimiser.runMDAndOptimise(data,
                                        xmlFilename=xmlFilename,
                                        rigidBody=rigidBody,
                                        doDihedral=doDihedral,
                                        doImproper=doImproper,
                                        doCharges=doCharges,
                                        d=d,
                                        **kw )

        if ok:
            self.logger.info( "runMDAndOptimise succeeded" )


        self.analyse.stop('runMDAndOptimise', d )

        self.fromHoomdblueSystem( optimiser.system )

        return ok

    def seed( self, nblocks, fragmentType=None, maxTries=500, center=False ):
        """ Seed a cell with nblocks of type fragmentType.

        Args:
        nblocks - the number of blocks to add.
        fragmentType - the type of blocks to add. If fragment is None, or omitted, then blocks of a randomly
                       chosen type will be added.
        maxTries - the number of attempts to make when adding a block before the seed step is fails and returns.
        center - True/False - if True, place the first block in the center of the cell.

        Returns:
        the number of blocks added
        """

        if self.A is None or self.B is None or self.C is None:
            raise RuntimeError,"Need to specify cell before seeding"

        #if not len( self._fragmentLibrary ) or not len( self.bondTypes):
        #    raise RuntimeError,"Must have set an initBlock and bondType before seeding."
        if not len( self._fragmentLibrary ):
            raise RuntimeError,"Must have set an initBlock before seeding."

        self.logger.info("seed adding {0} block of type {1}".format( nblocks, fragmentType ) )

        numBlocksAdded = 0
        # Loop through the nblocks adding the blocks to the cell
        for seedCount in range( nblocks ):

            # Create new block
            newblock = self.getInitBlock( fragmentType=fragmentType )

            tries = 0
            while True:
                # quit on maxTries
                if tries >= maxTries:
                    self.logger.critical("Exceeded maxtries when seeding after adding {0}".format(numBlocksAdded))
                    self.analyse.stop( 'seed',d={'num_tries':tries} )
                    return numBlocksAdded

                # if center put the first one in the center of the cell
                # only if this is the first attempt as otherwise we always fail if there is already something there
                if center and seedCount==0 and tries==0:
                    newblock.translateCentroid( [ self.A/2, self.B/2, self.C/2 ] )
                else:
                    # Move the block and rotate it
                    self.randomMoveBlock( newblock )

                #Add the block so we can check for clashes/bonds
                idxBlock = self.addBlock( newblock )

                # Test for Clashes with other molecules
                if self.checkMove( idxBlock ):
                    if self.processBonds() > 0:
                        self.logger.info("Added bond in seed!")
                    self.logger.debug("seed added block {0} after {1} tries.".format( seedCount+1, tries ) )
                    self.analyse.stop('seed',d={'num_tries':tries} )
                    numBlocksAdded += 1
                    break

                # Unsuccessful so remove the block from cell
                self.delBlock(idxBlock)

                # If seed fails with center need to bail on first one.
                if center and seedCount==0 and tries==0:
                    self.logger.warn("Seed with center failed to place first block in center!")

                # increment tries counter
                tries += 1

            # End Clash loop
        # End of loop to seed cell

        self.logger.info("Seed added {0} blocks. Cell now contains {1} blocks".format( numBlocksAdded,
                                                                                       len(self.blocks) ) )

        return numBlocksAdded

    def setMaxBond(self, bondType, count ):
        """Limit the number of bondType bonds to an individual fragment to count bonds.

        Args:
        bondType - the bondType (FRAGMENT1:ENDGROUP1-FRAGMENT2:ENDGROUP2) as was specified with the call
                   to addBondType
        count - the maximum number of permissible bonds for a single fragment.

        """
        # Get fragmentType and endGroupType
        fragmentType, endGroupType = bondType.split(self.ENDGROUPSEP)

        # Now get the library fragment to set it's maxBond parameter
        fragment = self._fragmentLibrary[ fragmentType ]
        return fragment.setMaxBond( bondType, count )

    def _setupAnalyse(self, logfile='ambuild.csv'):
        self.analyse = Analyse( self, logfile=logfile )
        return

    def setupLogging( self, logfile="ambuild.log", mode='w', doLog=False ):
        """
        Set up the various log files/console logging and return the logger
        """

        logger = logging.getLogger()
        #if not doLog:
        #    logging.disable(logging.DEBUG)

        # First check if there are any handlers active - if so we remove them
        if len(logger.handlers):
            # Need to copy as otherwise the list changes while we're cycling through it
            for hdlr in copy.copy(logger.handlers):
                logger.removeHandler(hdlr)

        # Not entirely sure why this needed - set overall level of the logger to debug
        logger.setLevel( logging.DEBUG )

        # create file handler and set level to debug
        self.logfile=logfile
        fl = logging.FileHandler( self.logfile, mode=mode )
        if doLog:
            fl.setLevel( logging.DEBUG )
        else:
            fl.setLevel( logging.INFO )

        # create formatter for fl
        formatter = logging.Formatter( '%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

        # add formatter to fl
        fl.setFormatter( formatter )

        # add fl to logger
        logger.addHandler( fl )

        # Now create console logger for outputting stuff
        # create file handler and set level to debug
        # Seems they changed the api in python 2.6->2.7
        try:
            cl = logging.StreamHandler( stream=sys.stdout )
        except TypeError:
            cl = logging.StreamHandler( strm=sys.stdout )
        cl.setLevel( logging.INFO )
        #cl.setLevel( logging.DEBUG )

        # create formatter for fl
        # Always add a blank line after every print
        formatter = logging.Formatter( '%(message)s\n' )

        # add formatter to fl
        cl.setFormatter( formatter  )

        # add fl to logger
        logger.addHandler( cl )

        self.logger = logger
        return

    def surroundCells(self, key ):
        return self._surroundCells( key, self.numBoxA, self.numBoxB, self.numBoxC )

    def _surroundCells(self, key, numBoxA, numBoxB, numBoxC ):
        """
        return a list of the cells surrounding the one with the given key
        """

        #print "box size {} : {},{},{}".format( self.boxSize, self.numBoxA, self.numBoxB, self.numBoxC )

        # REM - max box num is numboxes - 1
        a,b,c = key
        l = []
        for  i in [ 0, -1, +1 ]:
            for j in [ 0, -1, +1 ]:
                for k in [ 0, -1, +1 ]:
                    # Impose periodic boundaries
                    ai = a+i
                    if ai < 0:
                        ai = numBoxA-1
                    elif ai > numBoxA-1:
                        ai = 0
                    bj = b+j
                    if bj < 0:
                        bj = numBoxB-1
                    elif bj > numBoxB-1:
                        bj = 0
                    ck = c+k
                    if ck < 0:
                        ck = numBoxC-1
                    elif ck > numBoxC-1:
                        ck = 0
                    skey = (ai, bj, ck)
                    #print "sKey ({},{},{})->({})".format(a,b,c,skey)
                    if skey not in l:
                        l.append(skey)

        return l

    def updateFromBlock(self, block ):
        """Update cell parameters from the block"""

        if block.maxAtomRadius <= 0:
            raise RuntimeError,"Error updating cell from block"

        if block.maxAtomRadius() > self.maxAtomRadius:

            # What happens to already added blocks if we call this after we've added them?
            # Not sure so being safe...
            if len( self.blocks ) > 0:
                raise RuntimeError,"Adding initblock after blocks have been added - not sure what to do!"

            self.updateCellSize( maxAtomRadius=block.maxAtomRadius() )

        return

    def updateCellSize(self, boxMargin=None, maxAtomRadius=None, MARGIN=0.01 ):
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
            self.boxMargin = max( self.atomMargin, self.bondMargin ) + MARGIN

        assert self.boxMargin != 0 and self.maxAtomRadius != 0

        self.boxSize = ( self.maxAtomRadius * 2 ) + self.boxMargin

        #jmht - ceil or floor
        self.numBoxA = int(math.ceil( self.A / self.boxSize ) )
        self.numBoxB = int(math.ceil( self.B / self.boxSize ) )
        self.numBoxC = int(math.ceil( self.C / self.boxSize ) )

        self.logger.debug( "updateCellSize: boxSize {0} nboxes: {1} maxR {2} margin {3}".format( self.boxSize,
                                                                                      ( self.numBoxA, self.numBoxB, self.numBoxC ),
                                                                                     self.maxAtomRadius,
                                                                                     self.boxMargin)
                          )

        return

    def XupdateIndices(self):
        """Update the index where the data for each block starts in the overall cell list"""

        atomCount = 0 # Global count in cell
        bodyCount = -1
        for idxBlock, block in cell.blocks.iteritems():
            for idxFrag,frag in enumerate(block._fragments): # need index of fragment in block
                bodyCount += frag.numBodies()
                atomCount += frag.lenData()
        return

    def writePickle( self, fileName="cell.pkl" ):
        """Pickle ourselves"""

        # Need to close all open filehandles and the logger handlers
        #for l in self.logger.handlers:
        #    print "GOT HANDLER1 ",l

        # No idea why I can't get the log to close and then reopen with append mode
        # see _get_state_ and _set_state_ in this class
        if False:
            self._clLogHandler.close()
            self._flLogHandler.close()
            self.logger.removeHandler( self._clLogHandler )
            self.logger.removeHandler( self._flLogHandler )
            del self._clLogHandler
            del self._flLogHandler


        fileName = os.path.abspath( fileName )

        # Create the pickle file
        util.pickleObj( self, fileName )

        # Restart logging with append mode
        #self.setupLogging( mode='a' )

        self.logger.info( "Wrote pickle file: {0}".format(fileName) )

        return

    def writeCar( self, ofile="ambuild.car", data=None, periodic=True, skipDummy=False ):
        """Car File
        """

        if not data:
            data = self.dataDict()

        car = "!BIOSYM archive 3\n"
        if periodic:
            car += "PBC=ON\n"
        else:
            car += "PBC=OFF\n"

        car += "ambuild generated car file\n"
        tstr = time.strftime( "%a %b %d %H:%M:%S %Y", time.gmtime() )
        car += "!DATE {0}\n".format( tstr )

        if periodic:
            car += "PBC  {0: < 9.4F} {1: < 9.4F} {2: < 9.4F}  90.0000   90.0000   90.0000 (P1)\n".format( data['A'],
                                                                                                          data['B'],
                                                                                                          data['C'] )
#         for i, ( idxBlock, block ) in enumerate( self.blocks.iteritems() ):
#             for j, coord in enumerate( block.iterCoord() ):
        for i, coord in enumerate( data['coord'] ):
                if periodic:
                    x, ix = util.wrapCoord( coord[0], data['A'], center=False )
                    y, iy = util.wrapCoord( coord[1], data['B'], center=False )
                    z, iz = util.wrapCoord( coord[2], data['C'], center=False )
                else:
                    x = coord[0]
                    y = coord[1]
                    z = coord[2]

                atype = data['type'][ i ]
                charge = data['charge'][ i ]
                label = data['label'][ i ][:5]
                symbol = data['symbol'][ i ]
                car += "{0: <5} {1: >15.10} {2: >15.10} {3: >15.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format( label, x, y, z, atype, symbol, charge )

        car += "end\nend\n\n"

        with open( ofile, 'w' ) as f:
            fpath = os.path.abspath(f.name)
            f.writelines( car )

        self.logger.info( "Wrote car file: {0}".format(fpath) )
        return

    def writeCml(self, cmlFilename, rigidBody=True, data=None, periodic=True, pruneBonds=False):

#         atomTypes = []
#         coords    = []
#         symbols   = []
#         bonds     = []
#         count     = 0
#         for block in self.blocks.itervalues():
#             bonds += [ (b1+count, b2+count) for (b1, b2 ) in block.bonds() ]
#             for i, c in enumerate(block.iterCoord() ):
#                 coords.append(c)
#                 symbols.append(block.symbol(i))
#                 atomTypes.append(block.type(i))
#                 count += 1

        if data is None:
            d = self.dataDict(periodic=periodic, rigidBody=rigidBody, fragmentType=None)
        else:
            d = data

        cell = None
        if periodic:
            cell = [self.A,self.B,self.C]

        cmlFilename = util.writeCml(cmlFilename,
                                    d.coords,
                                    d.symbols,
                                    bonds=d.bonds,
                                    atomTypes=d.atomTypes,
                                    cell=cell,
                                    pruneBonds=pruneBonds)

        self.logger.info( "Wrote cml file: {0}".format(cmlFilename) )
        return

    def writeXyz(self, ofile, data=None, periodic=False ):
        """Write out the cell atoms to an xyz file
        If label is true we write out the atom label and block, otherwise the symbol
        """

        if data is None:
            d = self.dataDict(periodic=periodic, fragmentType=None)
        else:
            d = data

        if periodic:
            fpath = util.writeXyz(ofile,d.coords,d.symbols)
        else:
            fpath = util.writeXyz(ofile,d.coords,d.symbols,cell=[self.A,self.B,self.C])

        self.logger.info( "Wrote cell file: {0}".format(fpath) )
        return

    def zipBlocks(self, bondMargin=None, bondAngleMargin=None):
        """Join existing blocks in the cell by changing the bondMargin and bondAngleMargin parameters that were
        specified when the cell was created, and then looping over all the free endGroups to see if any can bond
        with the new parameters. The blocks are not moved in this step and no check is made as to whether there
        are any other atoms between the two endGroups when a bond is made.

        Args:
        bondMargin - the new bondMargin [degrees] (see __init__ for definition)
        bondAngleMargin - the new bondAngleMargin [degrees] (see __init__ for definition)"""

        self.logger.info("Zipping blocks with bondMargin: {0} bondAngleMargin {1}".format(bondMargin, bondAngleMargin )  )

        # Convert to radians
        bondAngleMargin = math.radians(bondAngleMargin)

        # Calculate the number of boxes

        # Should calculate max possible bond length
        maxBondLength=2.5
        boxSize = maxBondLength + bondMargin
        numBoxA = int(math.ceil( self.A / boxSize ) )
        numBoxB = int(math.ceil( self.B / boxSize ) )
        numBoxC = int(math.ceil( self.C / boxSize ) )

        # Create empty cells to hold data
        cell1 = {}
        cell3 = {}

        # Get a list of all block, endGroups
        endGroups = []
        for block in self.blocks.itervalues():
            egs = block.freeEndGroups()
            if len(egs):
                block.zipCell = {}
                for endGroup in egs:
                    endGroups.append( ( block, endGroup.endGroupIdx() ) )

        if not len(endGroups) > 0:
            self.logger.warn("zipBlocks found no free endGroups!")
            return 0

        # Add all (block, idxEndGroup) tuples to the cells
        for (block, idxEndGroup) in endGroups:

            coord = block.coord( idxEndGroup )

            # Periodic Boundaries
            x = coord[0] % self.A
            y = coord[1] % self.B
            z = coord[2] % self.C

            # Calculate which cell the atom is in
            a=int( math.floor( x / boxSize ) )
            b=int( math.floor( y / boxSize ) )
            c=int( math.floor( z / boxSize ) )
            key = (a,b,c)
            block.zipCell[ idxEndGroup ] = key

            try:
                # Add the atom to the cell
                cell1[ key ].append( ( block, idxEndGroup ) )
            except KeyError:
                # Add the cell to the list and then add the atom
                cell1[ key ] = [ ( block, idxEndGroup ) ]
                # Map the cells surrounding this one
                cell3[ key ] = self._surroundCells( key, numBoxA, numBoxB, numBoxC )

        # Loop through all the end groups and run canBond
        self._possibleBonds = []
        egPairs             = []
        c1                  = []
        c2                  = []
        for block1, idxEndGroup1 in endGroups:

            # Get the box this atom is in
            key = block1.zipCell[ idxEndGroup1 ]

            # Get a list of the boxes surrounding this one
            surrounding = cell3[ key ]

            # For each box loop through all its atoms
            for i, sbox in enumerate( surrounding ):

                # Check if we have a box with anything in it
                # use exception so we don't have to search through the whole list
                try:
                    alist = cell1[ sbox ]
                except KeyError:
                    continue

                for (block2, idxEndGroup2) in alist:

                    # Self-bonded blocks need special care
                    if block1 == block2:
                        # Don't check endGroups against themselves
                        if idxEndGroup1 == idxEndGroup2:
                            continue
                        # Make sure the two atoms are separated by at least 3 bonds - could
                        # probably put this check in canBond but it would slow the normal
                        # bonding down - need to think about best way to do this
                        if idxEndGroup2 in block1.atomBonded3( idxEndGroup1 ):
                            continue

                    # PROBABLY A BETTER WAY OF DOING THIS
                    p1 = ( block1,idxEndGroup1, block2, idxEndGroup2 )
                    p2 = ( block2, idxEndGroup2, block1, idxEndGroup1 )
                    if p1 not in egPairs and p2 not in egPairs:
                        # Need to check if it is already in there as we loop over all endGroups
                        # so we will have both sides twice
                        egPairs.append( p1 )
                        c1.append( block1.coord( idxEndGroup1 ) )
                        c2.append( block2.coord( idxEndGroup2 ) )

        if not len(egPairs) > 0:
            self.logger.info("zipBlocks: no endGroups close enough to bond" )
            return 0

        # Calculate distances between all pairs
        #distances = util.distance(c1, c2)
        distances = self.distance(c1, c2)

        # Now check for bonds
        for i, ( block1, idxEndGroup1, block2, idxEndGroup2 ) in enumerate( egPairs ):
            got = self.canBond( block1,
                                idxEndGroup1,
                                block2,
                                idxEndGroup2,
                                distances[i],
                                bondMargin,
                                bondAngleMargin,
                                )

        # Assumption is that zipBlocks called on a valid structure, so we don't do any checks for clashes
        # just process the bonds

        # Process any bonds
        todo = len( self._possibleBonds )

        self.logger.info("zipBlocks: found {0} additional bonds".format( todo ) )
#         for b in self._possibleBonds:
#             print "Attempting to bond: {0} {1} {2} -> {3} {4} {5}".format( b.block1.id(),
#                                                                    b.endGroup1.blockEndGroupIdx,
#                                                                    b.block1.atomCoord( b.endGroup1.blockEndGroupIdx),
#                                                                    b.block2.id(),
#                                                                    b.endGroup2.blockEndGroupIdx,
#                                                                    b.block2.atomCoord( b.endGroup2.blockEndGroupIdx),
#                                                                 )
#
        bondsMade = 0
        if todo > 0:
            bondsMade = self.processBonds()
            if bondsMade != todo:
                self.logger.debug("Made fewer bonds than expected in zip: {0} -> {1}".format(
                                                                                                todo,
                                                                                                bondsMade ) )

        self.analyse.stop('zip')

        return bondsMade

    def __str__(self):
        """
        """
        s = ""
        s += "InitBlocks: {0}".format( self._fragmentLibrary )
        s += "Blocks: "
        for block in self.blocks.values():
            s+= str(block)

        return s

    def __getstate__(self):
        """
        Return a dict of objects we want to pickle.

        This is required as we can't pickle objects containing a logger as they have
        file handles (can remove by calling logging.shutdown() ) and lock objects (not removed
        by calling shutdown)."""

        # Return everything bar our logger
        d = dict(self.__dict__)
        del d['logger']
        d['analyseLogfile']=d['analyse'].logfile
        del d['analyse']
        return d

    def __setstate__(self, d):
        """Called when we are unpickled """
        self.__dict__.update(d)
        if 'logfile' in d: # Hack for older versions with no logfile attribute
            logfile = util.newFilename( d['logfile'] )
        else:
            logfile = 'ambuild_1.log'
        self.setupLogging( logfile=logfile )

        if 'analyseLogfile' in d:
            logfile = util.newFilename( d['analyseLogfile'] )
        else:
            logfile = 'ambuild_1.csv'
        self._setupAnalyse( logfile=logfile )

        return

class TestCell(unittest.TestCase):

#     def setUp(self):
#
#         thisd =  os.path.abspath( os.path.dirname( __file__ ) )
#         paths = thisd.split( os.sep )
#         self.ambuildDir = os.sep.join( paths[ : -1 ] )
#
#         self.cx4Car = os.path.join( self.ambuildDir, "blocks", "cx4.car" )
#         self.ch4Car = os.path.join( self.ambuildDir, "blocks", "ch4.car" )
#         self.capLinker = os.path.join( self.ambuildDir, "blocks", "cap_linker.car" )
#         self.benzeneCar = os.path.join( self.ambuildDir, "blocks", "benzene.car" )
#         self.benzene2Car = os.path.join( self.ambuildDir, "blocks", "benzene2.car" )
#         self.pafCar = os.path.join( self.ambuildDir, "blocks", "PAF_bb_typed.car" )
#
#         return

    @classmethod
    def setUpClass(cls):

        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ambuildDir = os.sep.join( paths[ : -1 ] )

        cls.cx4Car = os.path.join( cls.ambuildDir, "blocks", "cx4.car" )
        cls.ch4Car = os.path.join( cls.ambuildDir, "blocks", "ch4.car" )
        cls.ch4_1Car = os.path.join( cls.ambuildDir, "blocks", "ch4_1.car" )
        cls.capLinker = os.path.join( cls.ambuildDir, "blocks", "cap_linker.car" )
        cls.benzeneCar = os.path.join( cls.ambuildDir, "blocks", "benzene.car" )
        cls.benzene2Car = os.path.join( cls.ambuildDir, "blocks", "benzene2.car" )
        cls.pafCar = os.path.join( cls.ambuildDir, "blocks", "PAF_bb_typed.car" )
        cls.dcxCar = os.path.join( cls.ambuildDir, "blocks", "DCX.car" )
        cls.ch4Ca2Car = os.path.join( cls.ambuildDir, "blocks", "ch4Ca2.car" )
        cls.amineCar = os.path.join( cls.ambuildDir, "blocks", "amine_typed.car" )
        cls.triquinCar = os.path.join( cls.ambuildDir, "blocks", "triquin_typed.car" )

        print "START TEST CELL"
        if True:
            # Cell dimensions need to be: L > 2*(r_cut+r_buff) and L < 3*(r_cut+r_buff)
            # From code looks like default r_buff is 0.4 and our default r_cut is 5.0
            boxDim=[20,20,20]
            mycell = Cell(boxDim)
            mycell.libraryAddFragment(filename=cls.benzene2Car, fragmentType='A')
            mycell.addBondType( 'A:a-A:a')
            mycell.seed( 5 )
            mycell.growBlocks( 8 )
            print "FINISHED TEST CELL"
            cls.testCell = mycell

        return

    def testCX4(self):
        """First pass"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment( fragmentType='A', filename=self.cx4Car )
        mycell.addBondType( 'A:a-A:a' )
        mycell.seed( 1 )
        mycell.growBlocks( 1 )

        return

#     def XXXtimeCheck(self):
#         """NOT A TEST JUST CODE TO TIME CHECKMOVE"""
#
#         def cellFromPickle(pickleFile):
#             with open(pickleFile) as f:
#                 myCell=cPickle.load(f)
#             return myCell
#
#         mycell = cellFromPickle("step_1.pkl")
#         # Get the last block
#         idxb = mycell.blocks.keys()[-1]
#         b = mycell.blocks[ idxb ]
#         be = b.freeEndGroups()[0]
#
#         i =  mycell.getInitBlock()
#         ie = i.freeEndGroups()[0]
#
#         def run():
#             global b, be, i, ie
#             for _ in xrange( 1000):
#                 b.positionGrowBlock( be, ie, dihedral=None )
#                 blockId = mycell.addBlock( i )
#                 mycell.checkMove( blockId )
#                 mycell.delBlock(blockId)
#
#         cProfile.run('run()','restats')
#         return

    def testAmbody(self):

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.ch4Ca2Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')
        added = mycell.seed( 3 )
        self.assertEqual(added, 3)

        return

    def testBond(self):
        """First pass"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        fragmentType='A'
        mycell.libraryAddFragment( fragmentType=fragmentType, filename=self.benzeneCar )
        mycell.addBondType( 'A:a-A:a' )

        centralBlock = mycell.getInitBlock( fragmentType=fragmentType )

        block1 = mycell.getInitBlock( fragmentType=fragmentType )
        block2 = mycell.getInitBlock( fragmentType=fragmentType )

        centralBlock.positionGrowBlock( centralBlock.freeEndGroups()[ 0 ], block1.freeEndGroups()[ 0 ] )
        centralBlock.positionGrowBlock( centralBlock.freeEndGroups()[ 1 ], block2.freeEndGroups()[ 0 ] )

        # Now add growBlock to the cell so we can check for clashes
        mycell.addBlock( block1 )
        mycell.addBlock( block2 )

        # Now add the central block - it will have 2 bonds to make
        centralId = mycell.addBlock( centralBlock )

        self.assertTrue( mycell.checkMove( centralId ), "checkMove failed!" )
        self.assertEqual( mycell.processBonds(), 2 )

        return

    def testBlockTypes(self):
        """Test we can add a block correctly"""


        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment( fragmentType='A', filename=self.ch4Car )
        mycell.libraryAddFragment( fragmentType='B', filename=self.ch4Car)
        mycell.libraryAddFragment( fragmentType='C', filename=self.ch4Car )
        mycell.libraryAddFragment( fragmentType='D', filename=self.ch4Car )

        mycell.addBondType( 'A:a-B:a' )
        mycell.addBondType( 'A:a-C:a' )
        mycell.addBondType( 'A:a-D:a' )

        mycell.seed( 1, fragmentType='A' )
        toGrow = 3
        grew = mycell.growBlocks( toGrow, endGroupType=None, maxTries=5)
        self.assertEqual( toGrow, grew, "testBlockTypes not ok")
        return

    def XtestCapBlocks(self):
        self.testCell.capBlocks(fragmentType='A', filename=self.capLinker )
        #mycell.dump()
        block = self.testCell.blocks[ self.testCell.blocks.keys()[0] ]
        return

    def testCellIO(self):
        """Check we can write out and then read in a cell
        """

        # Remember a coordinate for checking
        test_coord = self.testCell.blocks[ self.testCell.blocks.keys()[0] ].coord(4)

        outfile = "./testCellIO.pkl"
        self.testCell.writePickle( outfile )
        with open( outfile ) as f:
            newCell = cPickle.load( f )

        self.assertTrue( numpy.allclose( test_coord, self.testCell.blocks[ self.testCell.blocks.keys()[0] ].coord(4),
                                         rtol=1e-9, atol=1e-9 ),
                         msg="Incorrect testCoordinate of cell.")

        self.testCell.growBlocks( 5 )
        os.unlink( outfile )

        return

    def testCloseAtoms1(self):

        mycell = Cell( [2.1,2.1,2.1],atomMargin=0.1, bondMargin=0.1, bondAngleMargin=15 )
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')
        block1 = mycell.getInitBlock('A')

        natoms = block1.numAtoms()

        # block radius is 1.8
        block1.translateCentroid( [ mycell.A/2, mycell.B/2, mycell.C/2 ] )
        block1Idx = mycell.addBlock(block1)

        # Alone in cell but in center
        self.assertFalse( mycell.closeAtoms(block1Idx) )

        # Add second block overlapping first but in other image
        block2 = mycell.getInitBlock('A')
        block2.translateCentroid( [ mycell.A/2 + mycell.A,
                                   mycell.B/2 + mycell.B,
                                   mycell.C/2 + mycell.C ] )
        block2Idx = mycell.addBlock(block2)

        # Make sure every atom overlaps with ever other
        self.assertEquals( natoms*natoms, len(mycell.closeAtoms(block1Idx) ) )

        # See we have enough clashing atoms - NOT CHECKED THIS NUMBER
        self.assertEquals( 15, mycell._checkMove( block2Idx ) )

        return

    def testCloseAtoms(self):

        mycell = Cell([30,30,30], atomMargin=0.1, bondMargin=0.1, bondAngleMargin=15 )
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')
        block1 = mycell.getInitBlock('A')

        block1.translateCentroid( [ 1, 1, 1 ] )
        block1_id = mycell.addBlock(block1)
        block2 = mycell.getInitBlock('A')
        block2.translateCentroid( [ 2.2, 2.2, 2.2 ] )
        block2_id = mycell.addBlock(block2)

        closeList =  mycell.closeAtoms(block1_id)
        # Updated refPairs as now have
        refPairs = [ (0, 3), (1, 3), (2, 3) ] # Also used to check PBC
        closePairs = []
        for iatom, block, ioatom, distance in closeList:
            closePairs.append( (iatom,ioatom) )

        #mycell.writeXyz("close1.xyz", label=False)

        self.assertEqual(closePairs,
                        refPairs,
                         "Many contacts: {0}".format(closePairs))

        # Too far for any contacts
        mycell.delBlock(block2_id)
        block2 = mycell.getInitBlock('A')
        block2.translateCentroid( [10,10,10] )
        block2_id = mycell.addBlock(block2)

        close = mycell.closeAtoms(block1_id)
        self.assertFalse(close, "No contacts: ".format(close))

        # Now check across periodic boundary
        mycell.delBlock(block2_id)
        block2 = mycell.getInitBlock('A')
        x = 2.2 + 2 * mycell.A
        y = 2.2 + 2 * mycell.B
        z = 2.2 + 2 * mycell.C

        block2.translateCentroid( [x,y,z] )
        block2_id = mycell.addBlock(block2)

        #mycell.writeXyz("close2.xyz", label=False)

        closeList =  mycell.closeAtoms(block1_id)
        closePairs = []
        for iatom, block, ioatom, distance in closeList:
            closePairs.append( (iatom,ioatom) )

        self.assertEqual(closePairs, refPairs, "Periodic boundary: {}".format(closePairs))

        return

    def testCloseDistance(self):
        """
        Test distance and close together
        """
        boxDim=[30,30,30]
        mycell = Cell(boxDim)
        #mycell.initCell("../ch4_typed.car",incell=False)
        #block1 = mycell.initBlock.copy()
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')
        block1 = mycell.getInitBlock('A')

        b1 = numpy.array([2,2,2], dtype=numpy.float64 )
        block1.translateCentroid( b1 )
        block1_id = mycell.addBlock(block1)

        #block2=mycell.initBlock.copy()
        block2 = mycell.getInitBlock('A')
        b2 = numpy.array([3,3,3], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        block2_id = mycell.addBlock(block2)

        #mycell.writeXyz("close1.xyz", label=False)

#        closeList =  mycell.closeAtoms(block1_id)
#        for iatom, ioblock, ioatom in closeList:
#            b1c = block1._coords[iatom]
#            b2c = mycell.blocks[ioblock]._coords[ioatom]
#            distance = mycell.distance(b1c,b2c)
#            print "{}: {} -> {}:{} = {}".format(iatom,b1c,ioatom,b2c,distance)

        # Distance measured with Avogadro so it MUST be right...
        refd = 0.673354948616
        distance = mycell.distance( block1.coord(1), mycell.blocks[ block2_id ].coord(3) )
        self.assertAlmostEqual( refd, distance, 12, "Closest atoms: {}".format(distance) )

        return

    def testDLPOLY(self):
        """


        """

        boxDim=[100.0,100.0,100.0]
        mycell = Cell(boxDim)

        ch4Car=self.ch4Car
        #mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b2 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b3 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )


        # b1 in center
        b1.translateCentroid( [ 25, 25, 25 ] )
        endGroup1=b1.freeEndGroups()[ 0 ]
        endGroup2=b2.freeEndGroups()[ 0 ]
        b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )

        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock( bond )
        mycell.addBlock(b1)

        b3.translateCentroid( [ 25, 25, 20 ] )
        mycell.addBlock(b3)

        d = opt.DLPOLY()

        #data = mycell.dataDict(periodic=True, center=True, rigidBody=True)
        d.writeCONFIG(mycell)
        d.writeFIELD(mycell)

        return

    def testDistance(self):
        """Test the distance under periodic boundary conditions"""

        CELLA=CELLB=CELLC=10
        boxDim=[CELLA,CELLB,CELLC]
        mycell = Cell(boxDim)

        v1 = [ 2.46803012, 1.67131881, 1.96745421]
        v2 = [ 1.07988345, 0.10567109, 1.64897769]

        nv1 = numpy.array( v1 )
        nv2 = numpy.array( v2 )

        dc1 = mycell.distance(nv1,nv2)
        dn = numpy.linalg.norm(nv2-nv1)
        self.assertEqual( dc1, dn, "Distance within cell:{} | {}".format(dc1,dn) )

        x = v2[0] + 2 * CELLA
        y = v2[1] + 2 * CELLB
        z = v2[2] + 2 * CELLC
        nv2 = numpy.array([ x, y, z ])
        dc2 = mycell.distance(nv1,nv2)
        self.assertAlmostEqual( dc1, dc2, 11,"Distance across multiple cells +ve: {} | {}".format(dc1,dc2) )

        x = v2[0] - 2 * CELLA
        y = v2[1] - 2 * CELLB
        z = v2[2] - 2 * CELLC
        nv2 = numpy.array([ x, y, z ])
        dc3 = mycell.distance(nv1,nv2)
        self.assertAlmostEqual( dc1, dc3, 11, "Distance across multiple cells -ve: {} | {}".format(dc1,dc3) )

        v1 = numpy.array([ 0.0, 0.0, 0.0 ])
        v2 = numpy.array([ 0.0, 0.0, 8.0 ])
        dc = mycell.distance(v1,v2)
        self.assertEqual( dc, 2.0, "Distance across boundary cell:{}".format(dc) )

        return

    def testDihedral(self):

        CELLDIM=30
        boxDim=[CELLDIM,CELLDIM,CELLDIM]
        mycell = Cell(boxDim)

        p1 = numpy.array([ 0.0, 0.0, 0.0 ])
        p2 = numpy.array([ 10.0, 0.0, 0.0 ])
        p3 = numpy.array([ 10.0, 10.0, 0.0 ])
        p4 = numpy.array([ 20.0, 10.0, 10.0 ])

        ref = util.dihedral( p1, p2, p3, p4)

        self.assertEqual( ref, mycell.dihedral( p1, p2, p3, p4) )

        # Move by a full cell along x-axis - result should be the same
        p3 = numpy.array([ 10.0+CELLDIM, 10.0, 0.0 ])
        p4 = numpy.array([ 20.0+CELLDIM, 10.0, 10.0 ])

        self.assertEqual( ref, mycell.dihedral( p1, p2, p3, p4) )

        return

    def testEndGroupTypes(self):
        """Test we can add a block correctly"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment( fragmentType='A', filename=self.ch4Car )
        mycell.libraryAddFragment( fragmentType='B', filename=self.ch4Car)
        mycell.libraryAddFragment( fragmentType='C', filename=self.ch4Car )
        mycell.libraryAddFragment( fragmentType='D', filename=self.ch4Car )

        # Everything can bond to A (apart from A itself), but nothing can bond to anything else
        mycell.addBondType( 'A:a-B:a' )
        mycell.addBondType( 'A:a-C:a' )
        mycell.addBondType( 'A:a-D:a' )

        mycell.seed( 1, fragmentType='A' )
        mycell.seed( 1, fragmentType='B' )
        mycell.seed( 1, fragmentType='C' )

        banned = ['B:a','C:a','D:a']
        for _ in xrange( 5 ):
            cellEndGroup, libraryEndGroup = mycell.libraryEndGroupPair()
            if cellEndGroup.type() != 'A:a':
                self.assertNotIn( libraryEndGroup.type(),
                                  banned,
                                  "cell: {0} library: {1}".format( cellEndGroup.type(), libraryEndGroup.type() )
                                  )

        # Check it works if we give a cellEndGroup to bond
        cellEndGroup, libraryEndGroup = mycell.libraryEndGroupPair(cellEndGroups='B:a')

        # Check it fails if we ask for a block not in the cell
        self.assertRaises(RuntimeError, mycell.libraryEndGroupPair,{'cellEndGroups','D:a'} )

        # Check it works if we give it a list of libraryBlocks
        cellEndGroup, libraryEndGroup = mycell.libraryEndGroupPair(libraryEndGroups='B:a')
        self.assertEqual(  libraryEndGroup.type(),'B:a' )

        return

    def testFragMaxEnergy(self):
        """

        run calc
        find highest energy fragment => need to link from a compute/tag to a block/fragment

        - need list of block computes
        - need list of fragment computes

        each groupname is of the form X:Y where X is block index in cell and Y is fragment in index in the block

        The list of groups isn't needed, just the labels

        have a list of labels as long as there are groups - same as number of computes

        find highest energy compute - the index links back to the list of labels, and from the label we get the block/framgent

        block compute could be the first in the list

        each block has tags
        each fragment has tags

        tags = [
                 [
                   [s,e],[s,e],[s,e]
                   ]

        createGroups



        """

        boxDim=[100.0,100.0,100.0]
        mycell = Cell(boxDim)

        ch4Car=self.ch4Car
        #mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b2 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b3 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )


        # b1 in center
        b1.translateCentroid( [ 25, 25, 25 ] )
        endGroup1=b1.freeEndGroups()[ 0 ]
        endGroup2=b2.freeEndGroups()[ 0 ]
        b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )

        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock( bond )
        mycell.addBlock(b1)

        b3.translateCentroid( [ 25, 25, 20 ] )
        mycell.addBlock(b3)

        #mycell.writeXyz("foo.xyz")

        #Calculate the energy
        mycell.fragMaxEnergy()

        return

    def testGrowBlocks(self):
        """Test we can add blocks correctly"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        #mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        added = mycell.seed( 1 )
        self.assertEqual( added, 1, 'seed')
        natoms = mycell.numAtoms()
        nblocks=10
        added = mycell.growBlocks( nblocks, endGroupType=None, maxTries=1 )

        #mycell.writeCml("foo.cml", periodic=True, pruneBonds=False)
        self.assertEqual( added, nblocks, "growBlocks did not return ok")
        self.assertEqual(1,len(mycell.blocks), "Growing blocks found {0} blocks".format( len(mycell.blocks) ) )

        natoms2 = mycell.numAtoms()
        nblocks += 1
        # Need to subtract cap atoms
        self.assertEqual( natoms2, (natoms*nblocks)- (nblocks-1)*2 )

        return

    def testGrowLimited(self):
        """Test we can add blocks correctly"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.ch4_1Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:b' )
        mycell.setMaxBond( 'A:a', 1 )

        added = mycell.seed( 1 )
        mycell.growBlocks(10, cellEndGroups=None, maxTries=10)

        #mycell.dump()

        return

    def testGrowLimited2(self):
        """Test we can add blocks correctly"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.dcxCar, fragmentType='A')
        mycell.addBondType( 'A:CH-A:CCl' )
        mycell.setMaxBond( 'A:CH', 1 )

        added = mycell.seed( 1 )
        mycell.growBlocks(10, cellEndGroups=['A:CH'])

        #mycell.dump()

        return

    def testGrowBlocks2(self):
        """Test we can add blocks correctly"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        block1 = mycell.getInitBlock('A')

        # Position block so that it's aligned along x-axis
        # - use two opposing C-atoms 0 & 3
        block1.alignAtoms( 0, 3, [ 1, 0, 0 ] )

        # Now put it at the center of the cell
        block1.translateCentroid( [mycell.A/2, mycell.B/2, mycell.C/2] )

        # Add to the cell
        mycell.addBlock(block1)

        # Try adding 6 blocks - only 5 will fit
        nblocks=6
        added = mycell.growBlocks( nblocks, endGroupType=None, maxTries=5 )

        self.assertEqual( added, 5, "growBlocks did not add 5 blocks")
        self.assertEqual(1,len(mycell.blocks), "Growing blocks found {0} blocks".format( len(mycell.blocks) ) )

        return

    def testGrowBlocksDihedral(self):
        """Test we can add blocks correctly"""

        boxDim=[10,10,10]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        added = mycell.seed( 1, center=True )

        nblocks=2
        dihedral=90
        added = mycell.growBlocks( nblocks, cellEndGroups=None, libraryEndGroups=None, dihedral=dihedral, maxTries=5 )

        # Check angle between two specified dihedrals is correct
        block1 = mycell.blocks[ mycell.blocks.keys()[0] ]
        bond1 = block1._blockBonds[0]

        p1 = block1.coord( bond1.endGroup1.dihedralIdx() )
        p2 = block1.coord( bond1.endGroup1.endGroupIdx() )
        p3 = block1.coord( bond1.endGroup2.endGroupIdx() )
        p4 = block1.coord( bond1.endGroup2.dihedralIdx() )

        self.assertAlmostEqual( math.degrees( util.dihedral( p1, p2, p3, p4) ), dihedral )

        return

    def testGrowBlocksUw(self):
        """Test we can add blocks correctly"""

        boxDim=[10,10,10]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.amineCar, fragmentType='amine')
        mycell.libraryAddFragment(filename=self.triquinCar, fragmentType='triquin')
        mycell.addBondType('amine:a-triquin:b')

        mycell.seed( 1, fragmentType='triquin', center=True )
        mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups=['amine:a'], maxTries=1)

        return


    def testJoinBlocks(self):
        """Test we can add a block correctly"""

        boxDim=[30,30,30]
        mycell = Cell(boxDim)

        #mycell.libraryAddFragment(filename=self.pafCar, fragmentType='A')
        mycell.libraryAddFragment(filename=self.benzene2Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        toAdd = 5
        added = mycell.seed( toAdd, center=True )
        self.assertEqual( added, toAdd, 'seed')
        nc = 0
        for block in mycell.blocks.itervalues():
            nc += block.numAtoms()

        toJoin=4
        added = mycell.joinBlocks( toJoin, cellEndGroups=None, maxTries=1 )
        self.assertEqual( added, toJoin, "joinBlocks did join enough")
        self.assertEqual( 1, len(mycell.blocks), "joinBlocks found {0} blocks".format( len( mycell.blocks ) ) )

        nc2 = 0
        for block in mycell.blocks.itervalues():
            nc2 += block.numAtoms()

        # Need to subtract cap atoms
        self.assertEqual(nc-(toJoin*2), nc2, "Growing blocks found {0} coords".format( nc2 ) )

        return

    def testOptimiseGeometryRigid(self):
        """
        """
        #self.testCell.writeCml("foo.cml")
        self.testCell.optimiseGeometry( rigidBody=True,
                                        doDihedral=True,
                                        quiet=True,
                                        rCut=5.0,
                                        optCycles = 1000000,
                                        dt=0.005,
                                        Etol=1e-5,
                                         )
        os.unlink("hoomdOpt.xml")
        return

    def testOptimiseGeometryAll(self):
        """
        """
        #self.testCell.writeCml("foo.cml")
        self.testCell.optimiseGeometry( rigidBody=False, doDihedral=True, quiet=True )
        os.unlink("hoomdOpt.xml")
        return

    def XtestOptimiseGeometryMinCell(self):
        """
        """

        # Create a large cell and populate it with a block
        CELLA = CELLB = CELLC = 100
        mycell = Cell( )
        mycell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        mycell.libraryAddFragment( filename=self.benzeneCar, fragmentType='A' )
        mycell.addBondType( 'A:a-A:a' )

        mycell.seed( 1 )
        mycell.growBlocks(3,'A')

        # Get the block and put it in the centre of the cell
        block = mycell.blocks[ mycell.blocks.keys()[0] ]
        block.translateCentroid(  [ CELLA/2, CELLB/2, CELLC/2 ] )

        #mycell.dump()
        # Center of mass shouldn't have changed
        com = 0.0
        for blockIdx, block in mycell.blocks.iteritems():
            com += block.centerOfMass()

        mycell.optimiseGeometry(minCell=True, optAttempts=1 )

        comA = 0.0
        for blockIdx, block in mycell.blocks.iteritems():
            comA += block.centerOfMass()

        self.assertTrue( numpy.allclose( com, comA) )

        return

    def testOptimiseGeometryDihedral(self):
        """
        """
        self.testCell.optimiseGeometry( doDihedral=True, quiet=True )
        os.unlink("hoomdOpt.xml")
        return

    def testRunMD(self):
        """
        """
        self.testCell.runMD( doDihedral=True,
                             quiet=True,
                             rCut=5.0,
                             mdCycles=100,
                             T=1.0,
                             tau=0.5,
                             dt=0.0005,
                         )
        os.unlink("hoomdMD.xml")
        return

    def testRunMDAndOptimise(self):
        """
        """
        self.testCell.runMDAndOptimise( doDihedral=True, quiet=True )
        os.unlink("hoomdMDOpt.xml")
        return

    def testPeriodic(self):

        import hoomdblue

        mycell = self.testCell

        # Grab coords
        coords = []
        for block in mycell.blocks.itervalues():
            for i, coord in enumerate( block.iterCoord() ):
                coords.append( coord )

        # Wrap them
        wcoords = []
        images = []
        for c in coords:
            x, xi = util.wrapCoord( c[0], mycell.A, center=False )
            y, yi = util.wrapCoord( c[1], mycell.B, center=False )
            z, zi = util.wrapCoord( c[2], mycell.C, center=False )
            wcoords.append( [ x, y, z ] )
            images.append( [ xi, yi, zi ] )

        # Now umwrap them
        for i, c in enumerate( wcoords ):
            x = util.unWrapCoord( c[0], images[i][0], mycell.A, centered=False )
            y = util.unWrapCoord( c[1], images[i][1], mycell.B, centered=False )
            z = util.unWrapCoord( c[2], images[i][2], mycell.C, centered=False )

            self.assertAlmostEqual( x,
                                    coords[ i ][ 0 ],
                                    places=9,
                                    msg = "{0} {1}".format( x, coords[ i ][ 0 ] ) )
            self.assertAlmostEqual( y,
                                    coords[ i ][ 1 ],
                                    places=9,
                                    msg = "{0} {1}".format( y, coords[ i ][ 1 ] ) )
            self.assertAlmostEqual( z,
                                    coords[ i ][ 2 ],
                                    places=9,
                                    msg = "{0} {1}".format( z, coords[ i ][ 2 ] ) )

        # Now wrap them with centering
        wcoords = []
        images = []
        for c in coords:
            x, xi = util.wrapCoord( c[0], mycell.A, center=True )
            y, yi = util.wrapCoord( c[1], mycell.B, center=True )
            z, zi = util.wrapCoord( c[2], mycell.C, center=True )
            wcoords.append( [ x, y, z ] )
            images.append( [ xi, yi, zi ] )

        # Now umwrap them
        for i, c in enumerate( wcoords ):
            x = util.unWrapCoord( c[0], images[i][0], mycell.A, centered=True )
            y = util.unWrapCoord( c[1], images[i][1], mycell.B, centered=True )
            z = util.unWrapCoord( c[2], images[i][2], mycell.C, centered=True )

            self.assertAlmostEqual( x,
                                    coords[ i ][ 0 ],
                                    places=9,
                                    msg = "{0} {1}".format( x, coords[ i ][ 0 ] ) )

            self.assertAlmostEqual( y,
                                    coords[ i ][ 1 ],
                                    places=9,
                                    msg = "{0} {1}".format( y, coords[ i ][ 1 ] ) )
            self.assertAlmostEqual( z,
                                    coords[ i ][ 2 ],
                                    places=9,
                                    msg = "{0} {1}".format( z, coords[ i ][ 2 ] ) )

        # Now test with HOOMD-Blue
        filename = "periodicTest.xml"
        data = mycell.dataDict(periodic=True, center=True, rigidBody=True)
        o = opt.HoomdOptimiser()
        o.writeXml(data,
                   xmlFilename=filename,
                   rigidBody=True,
                   doDihedral=True,
                   doImproper=False,
                   doCharges=True,
                   )

        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        system = hoomdblue.init.read_xml( filename=filename )

        wcoords = []
        for i, p in enumerate( system.particles ):
            x, y, z  = p.position
            ix, iy, iz = p.image

            x = util.unWrapCoord( x, ix, mycell.A, centered=True )
            y = util.unWrapCoord( y, iy, mycell.B, centered=True )
            z = util.unWrapCoord( z, iz, mycell.C, centered=True )

            self.assertAlmostEqual( x,
                                    coords[ i ][ 0 ],
                                    places=6,
                                    msg = "{0} {1}".format( x, coords[ i ][ 0 ] ) )

            self.assertAlmostEqual( y,
                                    coords[ i ][ 1 ],
                                    places=6,
                                    msg = "{0} {1}".format( y, coords[ i ][ 1 ] ) )
            self.assertAlmostEqual( z,
                                    coords[ i ][ 2 ],
                                    places=6,
                                    msg = "{0} {1}".format( z, coords[ i ][ 2 ] ) )


        os.unlink( filename )

        return

    def testSeed(self):
        """Test we can seed correctly"""


        boxDim=[50,50,50]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment( filename=self.benzeneCar, fragmentType='A' )
        #mycell.libraryAddFragment( filename=self.ch4Car, fragmentType='A' )
        mycell.addBondType( 'A:a-A:a' )

        nblocks = 10
        added = mycell.seed( nblocks )

        self.assertEqual( nblocks, added, "Incorrect number of cell blocks" )
        self.assertEqual( nblocks, len(mycell.blocks), "Incorrect number of cell blocks" )

        # Check no block is outside the unit cell
        # We seed by the centroid so the edges could stick out by radius
        #radius = ch4.radius()

        bad = []
        for i,b in mycell.blocks.iteritems():
            radius = b.blockRadius()
            for c in b.iterCoord():
                if ( 0-radius > c[0] > CELLA[0]+ radius ) or \
                   ( 0-radius > c[1] > CELLB[1]+ radius ) or \
                   ( 0-radius > c[2] > CELLC[2]+ radius ):

                    bad.append( b )

        self.assertEqual( 0, len(bad), "Got {} blocks outside cell: {}".format( len(bad), bad ) )

        #mycell.writeXyz("seedTest.xyz")

        return

    def testSurroundBoxes(self):
        """
        """
        boxDim=[5,5,5]
        mycell = Cell(boxDim)
        # box size=1 - need to set manually as not reading in a block
        mycell.maxAtomRadius = 0.5
        mycell.atomMargin = 0.0

        mycell.boxSize = ( mycell.maxAtomRadius * 2 ) + mycell.atomMargin
        mycell.numBoxA = int(math.floor( mycell.A / mycell.boxSize ) )
        mycell.numBoxB = int(math.floor( mycell.B / mycell.boxSize ) )
        mycell.numBoxC = int(math.floor( mycell.C / mycell.boxSize ) )

        s = [(2, 2, 2), (2, 2, 1), (2, 2, 3), (2, 1, 2), (2, 1, 1), (2, 1, 3), (2, 3, 2), (2, 3, 1),
         (2, 3, 3), (1, 2, 2), (1, 2, 1), (1, 2, 3), (1, 1, 2), (1, 1, 1), (1, 1, 3), (1, 3, 2),
          (1, 3, 1), (1, 3, 3), (3, 2, 2), (3, 2, 1), (3, 2, 3), (3, 1, 2), (3, 1, 1), (3, 1, 3),
          (3, 3, 2), (3, 3, 1), (3, 3, 3)]
        self.assertEqual(s, mycell.surroundCells( (2,2,2) ), "in center")

        sb = mycell.surroundCells( (0,0,0) )
        s = [ (0, 0, 0), (0, 0, 4), (0, 0, 1), (0, 4, 0), (0, 4, 4), (0, 4, 1), (0, 1, 0), (0, 1, 4),
        (0, 1, 1), (4, 0, 0), (4, 0, 4), (4, 0, 1), (4, 4, 0), (4, 4, 4), (4, 4, 1), (4, 1, 0), (4, 1, 4),
         (4, 1, 1), (1, 0, 0), (1, 0, 4), (1, 0, 1), (1, 4, 0), (1, 4, 4), (1, 4, 1), (1, 1, 0), (1, 1, 4), (1, 1, 1)]
        self.assertEqual(s, sb , "periodic: {0}".format( sb ) )
        return

    def testZipBlocks(self):

        boxDim=[12.0,12.0,12.0]
        mycell = Cell(boxDim)

        ch4Car=self.ch4Car
        #mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b2 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b3 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )
        b4 = buildingBlock.Block( filePath=ch4Car, fragmentType='A'  )

        # b1 in center of cell
        b1.translateCentroid( [ mycell.A/2, mycell.B/2, mycell.C/2 ] )
        mycell.addBlock(b1)
        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Position b2 - all blocks will be positioned around this one
        b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )
        mycell.addBlock(b2)

        # Position b3
        endGroup1 = b2.freeEndGroups()[ 1 ]
        endGroup2 = b3.freeEndGroups()[ 0 ]
        b2.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )
        mycell.addBlock(b3)

        # Position b4
        endGroup1 = b2.freeEndGroups()[ 2 ]
        endGroup2 = b4.freeEndGroups()[ 0 ]
        b2.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )
        mycell.addBlock(b4)

        made = mycell.zipBlocks( bondMargin=0.5, bondAngleMargin=0.5 )

        self.assertEqual(made,3)

        return

    def testZipBlocks2(self):

        boxDim=[12.0,12.0,12.0]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A'  )
        b2 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A'  )

        # Position block so that it's aligned along x-axis
        # - use two opposing C-atoms 0 & 3
        b1.alignAtoms( 0, 3, [ 1, 0, 0 ] )

        b1.translateCentroid( [ mycell.A/2, mycell.B/2, mycell.C/2 ] )

        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Position the block
        b1.positionGrowBlock( endGroup1, endGroup2 )

        # and bond
        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock( bond )

        b1_id = mycell.addBlock(b1)

        made = mycell.zipBlocks( bondMargin=5, bondAngleMargin=5 )

        self.assertEqual( made, 1)
        return

    def testZipBlocks3(self):
        """Bonding with margin"""

        boxDim=[12.0,12.0,12.0]
        mycell = Cell(boxDim)

        ch4Car=self.ch4Car
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block( filePath=ch4Car, fragmentType='A' )
        b2 = buildingBlock.Block( filePath=ch4Car, fragmentType='A' )

        # Align bond along x-axis
        b1.alignAtoms( 0, 1, [ 1, 0, 0 ] )

        # b1 in center of cell
        b1.translateCentroid( [ mycell.A/2, mycell.B/2, mycell.C/2 ] )
        mycell.addBlock(b1)

        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Position b2
        b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )
        mycell.addBlock(b2)

        # Rotate on y-axis so that it's slightly off
        b2.rotateT([0,1,0],math.radians(15))

        made = mycell.zipBlocks( bondMargin=0.5, bondAngleMargin=5 )
        self.assertEqual( made, 0 )

        made = mycell.zipBlocks( bondMargin=0.5, bondAngleMargin=16 )
        self.assertEqual( made, 1 )

        return

    def testWriteCml(self):
        """
        write out cml
        """

        boxDim=[100,100,100]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='B')
        mycell.addBondType( 'B:a-B:a')

        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='C')
        mycell.addBondType( 'C:a-C:a')

        # Create block manually - this is so we have reproducible results
        b1 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )
        b2 = b1.copy()
        b3 = b1.copy()
        b4 = b1.copy()

        # b1 in center of cell
        b1.translateCentroid( [ mycell.A/2, mycell.B/2, mycell.C/2 ] )
        mycell.addBlock(b1)

        for b in [ b2, b3, b4]:
            endGroup1 = b1.freeEndGroups()[ 0 ]
            endGroup2 = b.freeEndGroups()[ 0 ]

            # Add b2
            b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )

            # bond them
            bond = buildingBlock.Bond(endGroup1, endGroup2)
            b1.bondBlock( bond )

        # Add another block that's not overlapping with the first
        b1 = buildingBlock.Block( filePath=self.ch4Car, fragmentType='B' )
        b2 = b1.copy()
        b3 = b1.copy()
        b4 = b1.copy()
        b1.translateCentroid( [ 10, 10, 10] )
        mycell.addBlock(b1)

        for b in [ b2, b3, b4]:
            endGroup1 = b1.freeEndGroups()[ 0 ]
            endGroup2 = b.freeEndGroups()[ 0 ]

            # Add b2
            b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )

            # bond them
            bond = buildingBlock.Bond(endGroup1, endGroup2)
            b1.bondBlock( bond )

        # Add another block that's not overlapping with the others
        b1 = buildingBlock.Block( filePath=self.ch4Car, fragmentType='C' )
        b2 = b1.copy()
        b3 = b1.copy()
        b4 = b1.copy()
        b1.translateCentroid( [ 90, 90, 90] )
        mycell.addBlock(b1)

        for b in [ b2, b3, b4]:
            endGroup1 = b1.freeEndGroups()[ 0 ]
            endGroup2 = b.freeEndGroups()[ 0 ]

            # Add b2
            b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )

            # bond them
            bond = buildingBlock.Bond(endGroup1, endGroup2)
            b1.bondBlock( bond )


        fname = "test.cml"
        mycell.writeCml(fname, periodic=False, rigidBody=True)
        # Test is same as reference
        with open(fname) as f:
            test = f.readlines()
        with open(os.path.join( self.ambuildDir, "tests", "testCellRigid.cml" )) as f:
            ref = f.readlines()

        self.assertEqual(test,ref,"cml compare rigid")

        fname = "test.cml"
        mycell.writeCml(fname, periodic=False, rigidBody=False)
        # Test is same as reference
        with open(fname) as f:
            test = f.readlines()
        with open(os.path.join( self.ambuildDir, "tests", "testCellAll.cml" )) as f:
            ref = f.readlines()

        self.assertEqual(test,ref,"cml compare all")
        os.unlink(fname)
        return


    def testWriteHoomdblue(self):
        """
        write out hoomdblue xml
        """
        import hoomdblue

        boxDim=[20,20,20]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType( 'A:a-A:a')

        # Create block manually - this is so we have reproducible results
        b1 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )
        b2 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )

        # b1 in center of cell
        b1.translateCentroid( [ mycell.A/2, mycell.B/2, mycell.C/2 ] )
        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Bond b2 to it
        b1.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi )
        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock( bond )
        mycell.addBlock(b1)

        initcoords = []
        for block in mycell.blocks.itervalues():
            for c in block.iterCoord():
                initcoords.append( c )

        xmlFilename = "testWriteHoomdblue.xml"
        data = mycell.dataDict(periodic=True, center=True, rigidBody=True)
        o = opt.HoomdOptimiser()
        o.writeXml(data,
                   xmlFilename=xmlFilename,
                   rigidBody=True,
                   doDihedral=True,
                   doImproper=False,
                   doCharges=True,
                   )

        # Test what we've written out matches the reference file
        with open(xmlFilename) as f:
            test = f.readlines()
        with open(os.path.join( self.ambuildDir, "tests", xmlFilename )) as f:
            ref = f.readlines()
        self.assertEqual(test,ref,"xml compare")

        # Init the sytem from the file
        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        system = hoomdblue.init.read_xml( filename=xmlFilename )

        # Read it back in to make sure we get the same values
        mycell.fromHoomdblueSystem( system )
        finalcoords = []
        for block in mycell.blocks.itervalues():
            for c in block.iterCoord():
                finalcoords.append( c )

        self.assertTrue( all( map( lambda x : numpy.allclose( x[0], x[1] ),  zip( initcoords, finalcoords ) ) ),
                         "coords don't match")

        #ok = opt.HoomdOptimiser().optimiseGeometry( xmlFilename,doDihedral=doDihedral)

        os.unlink( xmlFilename )

        return

    def XtestWriteHoomdblueMinCell(self):
        """
        write out hoomdblue xml
        """

        import hoomdblue

        # Create a large cell and populate it with a block
        CELLA = CELLB = CELLC = 100
        mycell = Cell( )
        mycell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        mycell.libraryAddFragment( filename=self.ch4Car, fragmentType='A' )
        mycell.addBondType( 'A:a-A:a' )

        mycell.seed( 1 )
        mycell.growBlocks(3,'A')

        # Get the block and put it in the centre of the cell
        block = mycell.blocks[ mycell.blocks.keys()[0] ]
        block.translateCentroid(  [ CELLA/2, CELLB/2, CELLC/2 ] )

        initcoords = []
        for block in mycell.blocks.itervalues():
            initcoords += block._coords

        filename = "testWriteHoomdblue.xml"
        natoms = mycell.writeHoomdXml( xmlFilename=filename )

        # Init the sytem from the file
        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        system = hoomdblue.init.read_xml( filename=filename )

        # Read it back in to make sure we get the same values
        mycell.fromHoomdblueSystem( system, minCell=True )

        finalcoords = []
        for block in mycell.blocks.itervalues():
            finalcoords += block._coords

        self.assertEqual( initcoords, finalcoords, "coords don't match")

        os.unlink( filename )

        return

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
