'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import copy
import math
import os
import random
import sys
import unittest

import numpy

# Our modules
import util
import buildingBlock

class Cell():
    '''
    classdocs
    '''


    def __init__( self, bondMargin=0.2, blockMargin=2.0, atomMargin=0.1, bondAngle=180, bondAngleMargin=15 ):
        '''
        Constructor
        '''
        
        
       
        # A, B and C are cell vectors - for the time being we assume they are orthogonal
        self.A = None
        self.B = None
        self.C = None
        
        # additional distance to add on when checking whether 2 atoms are close enough to bond
        self.bondMargin = bondMargin
        
        # additional distance to add on when checking if two blocks are close enough to require
        # checking their interactions
        self.blockMargin = blockMargin
        
        # additional distance to add on when checking if two atoms are close enough to clash
        self.atomMargin = atomMargin
        
        # The acceptable bond angle
        self.bondAngle = bondAngle
        
        self.bondAngleMargin = bondAngleMargin
        
        # Dict mapping box key to a list of tuples defining the atoms in that box
        self.box1={}
        
        # Dict mapping key to list of boxes surrounding the keyed box
        self.box3={}
        
        #jmht - hack
        self.maxAtomR=None
        
        # first block is kept separate as everything is built up from it
        self.initBlock = None
        
        # dictionary mapping id of the block to the block - can't use a list and indices
        # as we add and remove blocks and need to keep track of them
        self.blocks = {}
    
    def addBlock( self, block ):
        """
        Add the block and the calculate the boxes for all atoms.
        
        Args:
        block -- block to add
        
        Returns:
        index of the new block in the list
        """
        
        boxSize = self.boxSize()
        
        self.blocks[ id(block) ] = block
        
        # The index of the new block
        iblock = id(block)
        for icoord,coord in enumerate(block.coords):
            a=int( math.floor( coord[0] / boxSize ) )
            b=int( math.floor( coord[1] / boxSize ) ) 
            c=int( math.floor( coord[2] / boxSize ) )
            key = (a,b,c)
            block.atomCell[icoord] = key
            try:
                self.box1[key].append( (iblock,icoord) )
                #print "APPENDING KEY ",key,iblock,icoord
            except KeyError:
                #print "ADDING KEY ",key,iblock,icoord
                # Add to main list
                self.box1[key] = [(iblock,icoord)]
                # Map surrounding boxes
                self.box3[key] = self.surroundBoxes(key)
                
        return iblock
                
    def delBlock(self,blockId):
        """
        Remove the block with the given index from the cell
        """
        
        block =  self.blocks[blockId]
        
        # Remove each atom from the list
        keys = []
        for iatom, key in enumerate( block.atomCell ):
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
        
    def alignBlocks(self, block1, block1EndGroupIndex, block2, block2EndGroupIndex):
        """
        Align block2, so that the bond defined by the endGroup on block 2 is aligned with
        the bond defined by the endGroup on block 1
        
        This assumes that the blocks have already been positioned with the contact atoms at the origin
        """
        
        block1EndGroup = block1.coords[ block1EndGroupIndex ]
        
        # Aign with each axis in turn. The axis to align to is just the coord of bock1EndGroup
        self._alignAxis(block2, block2EndGroupIndex, block1EndGroup, axisLabel="x" )
        self._alignAxis(block2, block2EndGroupIndex, block1EndGroup, axisLabel="y" )
        self._alignAxis(block2, block2EndGroupIndex, block1EndGroup, axisLabel="z" )
        
        # Don't know why, but for the time being do it twice
        self._alignAxis(block2, block2EndGroupIndex, block1EndGroup, axisLabel="x" )
        self._alignAxis(block2, block2EndGroupIndex, block1EndGroup, axisLabel="y" )
        self._alignAxis(block2, block2EndGroupIndex, block1EndGroup, axisLabel="z" )
        
        return
    
    def _alignAxis(self, block, endGroupIndex, targetVector, axisLabel=None ):
        """
        Align the given block so that the bond deined by the endGroupIndex
        is aligned with the targetVector along the given axis
        """
        
        if axisLabel.lower() == "x":
            axis = numpy.array( [1,0,0] )
        elif axisLabel.lower() == "y":
            axis = numpy.array( [0,1,0] )
        elif axisLabel.lower() == "z":
            axis = numpy.array( [0,0,1] )
        else:
            raise RuntimeError,"Unrecognised Axis!"
        
        # For comparing angles
        small = 1.0e-8
        
        # Find the angle targetVector makes with the axis
        tvAngle = util.vectorAngle(targetVector, axis)
        #print "tvAngle: {}".format(tvAngle)
        
        # Find angle the block makes with axis
        blockEndGroup = block.coords[ endGroupIndex ]
        bAngle = util.vectorAngle(blockEndGroup, axis)
        #print "bAngle b4: {}".format(bAngle)
        
        # Need to rotate block on X-axis by tvAngle-bAngle
        angle = tvAngle-bAngle
        
        #print "Got angle: {}".format(angle*util.RADIANS2DEGREES)
        
        if not ( ( -small < angle < small ) or  ( numpy.pi-small < angle < numpy.pi+small ) ):
            rotAxis = numpy.cross( blockEndGroup, axis)
            #print "rotAxis {}".format(rotAxis)
            block.rotate( rotAxis, angle )
        
        #blockEndGroup = block.coords[ endGroupIndex ]
        #bAngle = util.vectorAngle(blockEndGroup, axis)
        #print "bAngle after: {}".format(bAngle)
        #print block
        
    def cellAxis(self,A=None,B=None,C=None):
        """
        Get or set the cell axes
        """
        if A and B and C:
            self.A = numpy.array( [A,0.0,0.0], dtype=numpy.float64 )
            self.B = numpy.array( [0.0,B,0.0], dtype=numpy.float64 )
            self.C = numpy.array( [0.0,0.0,C], dtype=numpy.float64 )
        else:
            return (A,B,C)
        
    def checkClash(self, newblock, ignore=[] ):
        """
        See if the given block clashes with any of the others
        return True if there is a clash
        ignore is a list of block indices to ignore
        """
        
        for i,block in self.blocks.iteritems():
            if i in ignore:
                continue
            if block.clash( newblock, blockMargin=self.blockMargin, atomMargin=self.atomMargin ):
                return True
            
        return False
    
    
    def boxSize(self):
        """
        Return the size of the box
        """
        return (self.maxAtomR*2)+self.atomMargin
    
    def surroundBoxes(self, key ):
        """
        return a list of the boxes surrounding the box with the given key
        """
        
        # jmht - why minus 1?
        maxKeyA = int(math.floor( self.A[0] / self.boxSize() ) )
        maxKeyB = int(math.floor( self.B[1] / self.boxSize() ) )
        maxKeyC = int(math.floor( self.C[2] / self.boxSize() ) )
        #print "box size {} : {},{},{}".format(self.boxSize(),maxKeyA,maxKeyB,maxKeyC)
        
        a,b,c = key
        l = []
        for  i in [ 0, -1, +1 ]:
            for j in [ 0, -1, +1 ]:
                for k in [ 0, -1, +1 ]:
                    # Impose periodic boundaries
                    #jmht - why is test minus 1?
                    ai = a+i
                    if ai < 0:
                        ai = maxKeyA
                    elif ai >= maxKeyA:
                        ai = 0
                    bj = b+j
                    if bj < 0:
                        bj = maxKeyB
                    elif bj >= maxKeyB:
                        bj = 0
                        
                    ck = c+k
                    if ck < 0:
                        ck = maxKeyC
                    elif ck >= maxKeyC:
                        ck = 0
                    #print "sKey ({},{},{})->({},{},{})-".format(a,b,c,ai,bj,ck)
                    l.append((ai, bj, ck))
        return l
        
#    def initCellLists(self):
#        """
#        initialise the cell lists
#        
#        
#        Data structure
#        a cell (3 coordinates)
#        each cell has a list of the 26+ itself cells surrounding it - Box3
#        each cell contains a list of the atoms that are within the cell - Box1
#        """
#        
#        boxSize = self.boxSize()
#        
#        # Dict mapping box key to a list of tuples defining the atoms in that box
#        self.box1 = {}
#        
#        # Build up a list of which atoms are in each box
#        for iblock,block in enumerate(self.blocks):
#            for icoord,coord in enumerate(block.coords):
#                a=int( math.floor( coord[0] / boxSize ) )
#                b=int( math.floor( coord[1] / boxSize ) ) 
#                c=int( math.floor( coord[2] / boxSize ) )
#                key = (a,b,c)
#                block.atomCell[icoord] = key
#                try:
#                    self.box1[key].append( (iblock,icoord) )
#                    #print "Adding box1: {},{}".format(iblock,icoord)
#                except KeyError:
#                    self.box1[key] = [ (iblock,icoord) ]
#                    #print "Adding box1: {},{}".format(i,j)
#        
#        # Now for each Box1, map the surrounding boxes
#        self.box3 = {}
#        for key in self.box1.keys():
#            self.box3[key] = self.surroundBoxes(key)
        
    def closeAtoms(self, iblock):
        """
        Find all atoms that are close to the atom in the given block.
        
        Args:
        iblock: index of the block in self.blocks
        
        Returns:
        a list of tuples: (thisAtomIndex, otherBlockIndex, otherAtomIndex) or None if there we no close contacts
        """
        
        contacts=[]
        
        block=self.blocks[iblock]
        for icoord,coord in enumerate(block.coords):
            
            # Get the box this atom is in
            key = block.atomCell[icoord]
            
            # Get a list of the boxes surrounding this one
            surrounding = self.box3[key]
            
            #For each box loop through all its atoms chekcking for clashes
            for sbox in surrounding:
                
                # For each box, get the list of the atoms as (block,coord) tuples
                # Check if we have a box with anything in it
                if not self.box1.has_key(sbox):
                    continue
                    
                for (ioblock, iocoord) in self.box1[ sbox ]:
                    
                    # Check we are not checking ourself - need to check block index too!
                    if iblock == ioblock:
                        continue
                    
                    oblock = self.blocks[ioblock]
                    ocoord = oblock.coords[iocoord]
                    
                    if ( self.distance( ocoord,coord ) < self.boxSize() ):
                        #print "ATOMS {}-{}:({}) and {}-{}:({}) close".format( iblock,icoord,coord,ioblock,iocoord,ocoord )
                        contacts.append( (icoord, ioblock,iocoord) )
                        
        if len(contacts):
            return contacts
        else:
            return None

        
    def checkMove(self, iblock, closeMargin=None, bondAngle=None, blockList=None ):
        """
        See how this move went.
        Arguments
        iblock - the index of the block in the list of blocks for the cell
        closeMargin - the extra bit to add on to the radii of the blocks to see if they are close
        bondAngle - the angle to define a suitable bond
        blockList - a list of tuples of blocks to sample - (index,block)
        
        we return:
        False - nothing was close
        clash - we clashed with something
        bonds - [ [iblock, (i,j] ]- a list of the bonds as the index of the other block and a pair
        of the atom in this and the other block forming the bond
        """
        
        assert closeMargin, bondAngle
        
        # Loop over all blocks
        if not blockList:
            blockList = [ i for i in range( len(self.blocks) ) ]
        
        
        block = self.blocks[iblock]
        bonds=[]
        #for i, oblock in enumerate( self.blocks ):
        for i in blockList:
            
            print "Checking block: {}".format(i)
            # skip the block we are using
            if i == iblock:
                continue
            
            # Get the next block
            oblock = self.blocks[i]
            
            # First see if we are close enough to consider bonding
            # Shouldn't need to use this as we are using findClose, but it seems to be required
            # - need to look into this
            if not block.close( oblock, blockMargin = self.blockMargin ):
                continue
            
            # See if we can bond
            bond = block.canBond( oblock, bondMargin=self.bondMargin, bondAngle=bondAngle )
            
            # No bond so move  to next block
            if not bond:
                continue
            
            # Either a bond or a clash
            if bond == "clash":
                # A clash so break out of the whole loop so we can reject this step
                return "clash"
            else:
                # Got a bond, so add it to the list of possible bonds
                print "Found bond between block {} and {}".format( iblock, i)
                bonds.append( [i, bond] )

            #End Loop Over Blocks
        
        if len(bonds):
            return bonds
        else:
            return False
        
        # Never get here
        assert False

    def directedShimmy(self, nsteps=100, nmoves=50, bondAngle=None ):
        """
        Shimmy by selecting a block and then moving around it.
        Lots of possiblities for speeding up - preseeclting which blocks
        to sample and a box around the target to avoid looping over undeeded atoms
        """
        
        #For writing out our progress
        filename= "SHIMMY_0.xyz"
        
        # For time being ensure we are given a bond angle
        assert bondAngle
        
        CLOSE_MARGIN=4.0 # how close 2 blocks are before we consider checking if they can bond
        BLOCK_MARGIN=2.0 # margin between the ranges that are used to sample the smaller moves
        
        nbonds=0
        for step in range(nsteps):
            
            if len(self.blocks) == 1:
                print "NO MORE BLOCKS TO BOND _ HOORAY!"
                return
            
            if not step % 100:
                print "Step: {}".format(step)
                filename = util.newFilename(filename)
                self.writeXyz( filename )
            
            # Get two blocks to sample round
            iblock, jblock = self.randomBlockId(count=2)
            block = self.blocks[iblock]
            oblock = self.blocks[jblock]
            print "Sampling block {} about {}".format(iblock,jblock)
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
            
            # Get a list of which blocks are close
            # Don't use - it's slower for some unfathomable reason
            closeBlocks = self.findClose( oblock, block, closeMargin=CLOSE_MARGIN, blockMargin=BLOCK_MARGIN )
            #print "Got close blocks: {}".format(closeBlocks)
            #closeBlocks=None
            
            clash=0
            noclash=0
            gotBond=False
            for _ in range( nmoves ):
                
                self.randomMoveAroundBlock( oblock, block, margin=BLOCK_MARGIN )
                
                print "Move: {} about {}".format( oblock.centroid(), block.centroid() )
                
                # Loop over all blocks to see if we can bond or if we clash
                bonds = self.checkMove( iblock, closeMargin=CLOSE_MARGIN, bondAngle=bondAngle, blockList=closeBlocks )
                
                # See what happend
                if not bonds:
                    noclash+=1
                    continue
                    
                if bonds == "clash":
                    clash+=1
                    continue
                
                # Got some bonds so deal with them
                # Need to decrement index if we are removing blocks
                bcount=0
                print "FOUND {} BOND(S)!!!!".format(len(bonds))
                nbonds+=len(bonds)
                for b,bond in bonds:
                    block.bond( self.blocks[b-bcount], bond )
                    # Remove the block from the list as it is now part of the other one
                    self.blocks.pop(b-bcount)
                    bcount+=1
                
                # No need to loop anymore
                gotBond=True
                break
                            
            print "End of moves  clash/noclash: {}/{}".format(clash,noclash)
            
            if not gotBond:
                # If no bonds place the bond back at it's original coord
                self.blocks[iblock] = orig_block
            
            #End move loop
        #End step loop
        
        print "END OF DIRECTED SHIMMY\nMade {} bonds and got {} clusters".format(nbonds,len(self.blocks))
    #End directedShimmy
    
    
    def distance(self, v1, v2 ):
        """
        Calculate the distance between two vectors in the cell
        under periodic boundary conditions
        """
        
        dx = v2[0] - v1[0]
        if math.fabs(dx) > self.A[0] * 0.5:
            dx = dx - math.copysign( self.A[0], dx)
        dy = v2[1] - v1[1]
        if math.fabs(dy) > self.B[1] * 0.5:
            dy = dy - math.copysign( self.B[1], dy)
        dz = v2[2] - v1[2]
        if math.fabs(dz) > self.C[2] * 0.5:
            dz = dz - math.copysign( self.C[2], dz)
        
        return math.sqrt( dx*dx + dy*dy + dz*dz )
    
    def Xdistance(self, x,y):
        """
        Calculate the distance between two vectors - currently just use for testing
        Actual one lives in cell
        """
        return numpy.linalg.norm(y-x)

    def findClose(self, oblock, block):
        """
        Find the blocks around oblock that might be close to block if we were moving
        it by the parameter used for the random move around block.
        return a list of the indexes
        """
        
        # Need to add the central block radius to the diameter of the sample block
        # plus the margins
        rb = block.radius()
        ro = oblock.radius()
        margin = ro + 2*rb + self.blockMargin + self.closeMargin
        
        centroid = oblock.centroid()
        closeBlocks = []
        # This should implicitly include oblock
        for i,b in self.block.iteritems():
            # See if we are close enough
            dist = numpy.linalg.norm( centroid - b.centroid() )
            if dist < margin + b.radius():
                closeBlocks.append( i )
        
        return closeBlocks
        
    def fromXyz(self, xyzFile ):
        """ Read in an xyz file containing a cell and recreate the cell object"""
        
        blocks = self._readXyz(xyzFile)
        # Need to initialise cell from first block
        initBlock = False
        for labels,coords in blocks:
            block = buildingBlock.BuildingBlock()
            block.createFromLabelAndCoords( labels, coords )
            if not initBlock:
                self.setInitBlock(block.copy())
                initBlock=True
            self.addBlock(block)
        
    def _readXyz(self, xyzFile ):
        """"Read in an xyz file containing a cell - cell axes is title line and 
        atoms are written out with their labels.
        This routine sets the axes and returns a list of the labels and coordinates
        """
        
        # Each block is a list of [label, coords
        blocks = []
        labels=[]
        coords=[]
        
        with open( xyzFile ) as f:
            
            # First line is number of atoms
            line = f.readline()
            natoms = int(line.strip())
            
            # Title line contains cell info
            line = f.readline().strip()
            fields = line.split(":")
            
            if fields[0] != "Axes":
                raise RuntimeError, "TITLE LINE NEEDS TO CONTAIN THE AXES"
            
            A,B,C = None, None, None
            # bit hacky but I'm lazy...
            exec( "A = {}".format(fields[1]) )
            exec( "B = {}".format(fields[2]) )
            exec( "C = {}".format(fields[3]) )
            
            self.cellAxis(A=A, B=B, C=C)
            
            lastBlock=-1
            for _ in range(natoms):
                
                line = f.readline().strip()
                fields = line.split()
                label = fields[0]
                
                # Determine block from last bit of label
                labelf = label.split("_")
                iblock = int( labelf[-1].split("#")[1] )
                
                if iblock != lastBlock:
                    # Add new block
                    if lastBlock != -1:
                        blocks.append( (labels,coords) )
                    lastBlock = iblock
                    labels=[]
                    coords=[]
                
                if len(labelf) == 2:
                    label = labelf[0]
                elif len(labelf) == 3:
                    label = labelf[0]+"_"+labelf[1]
                else:
                    raise RuntimeError,"HELLLP"
                
                labels.append(label) 
                coords.append( numpy.array(fields[1:4], dtype=numpy.float64) )
        
        # Add the last block
        blocks.append( (labels,coords) )
        
        return blocks


    def growBlock(self, block2):
        """
        Bond block2 to a randomly chosen block at randomly chosen endGroups in each
        """

        # Select 2 random blocks
        block1id = self.randomBlockId()
        block1 = self.blocks[ block1id ]
        
        # pick random endgroup and contact in target
        block1EndGroupIndex = block1.randomEndGroupIndex()
        # Need to copy this so we can remember it - otherwise it changes with the moves
        block1ContactIndex = block1.endGroupContactIndex( block1EndGroupIndex )
        block1Contact = block1.coords[ block1ContactIndex ].copy()
        
        # pick random endGroup on the new block
        block2EndGroupIndex = block2.randomEndGroupIndex()
        block2ContactIndex = block2.endGroupContactIndex( block2EndGroupIndex )
        block2Contact = block2.coords[ block2ContactIndex ].copy()
        
        # get the coord where the next block should bond
        bondPos = block1.newBondPosition( block1EndGroupIndex, block2, block2EndGroupIndex )
        #print "got bondPos: {}".format( bondPos )
        
        # Move block1Contact to center so that block1Contact -> block1EndGroup defines
        # the alignment we are after
        block1.translate( -block1Contact )
        # Same for block2
        block2.translate( -block2Contact )
        
        # Align along all 3 axes
        self.alignBlocks( block1, block1EndGroupIndex, block2, block2EndGroupIndex)
        
        # Move block1 back to its original coord
        block1.translate( block1Contact )
        
        # now turn the second block around
        
        # Need to get the new vectors as these will have changed due to the moves
        block2EndGroup = block2.coords[ block2EndGroupIndex ]
        
        # Find vector perpendicular to the bond axis
        # Dot product needs to be 0
        # xu + yv + zw = 0 - set u and v to 1, so w = (x + y)/z
        # vector is 1, 1, w
        w =  -1.0 * ( block2EndGroup[0] + block2EndGroup[1] ) / block2EndGroup[2]
        orth = numpy.array( [1.0, 1.0, w] )
        
        # Find axis that we can rotate about
        rotAxis = numpy.cross(block2EndGroup,orth)
        
        # Rotate by 180
        block2.rotate( rotAxis, numpy.pi )
        
        # Now need to coord the endGroup at the bond coord - at the moment
        # The contact atom is on the origin so we need to subtract the vector to the endGroup
        # from the translation vector
        
        #jmht FIX!
        block2.translate( bondPos + block2EndGroup )
        
        
        # Check it doesn't clash - could try rotating about the dihedral
        if not self.checkClash( block2, ignore=[block1id] ):
            # Bond the two blocks
            block1.bond(  block2, (block1EndGroupIndex,block2EndGroupIndex) )
            return True
        
        return False
        
    
    def initCell(self,inputFile, incell=True):
        """
        Read in the first block from the input file and save it
        """
        
        ib = buildingBlock.BuildingBlock()
        ib.fromCarFile( inputFile )
        
        if incell:
            # Make sure it is positioned in the cell - if not keep moving it about randomly till
            # it fits
            while True:
                #print "moving initblock"
                #print "{} | {} | {} | {}".format( ib.centroid()[0], ib.centroid()[1], ib.centroid()[2], ib.radius()  )
                if ib.centroid()[0] - ib.radius() < 0 or \
                   ib.centroid()[0] + ib.radius() > self.A[0] or \
                   ib.centroid()[1] - ib.radius() < 0 or \
                   ib.centroid()[1] + ib.radius() > self.B[1] or \
                   ib.centroid()[2] - ib.radius() < 0 or \
                   ib.centroid()[2] + ib.radius() > self.C[2]:
                    self.randomMove(ib, buffer=ib.radius())
                else:
                    # Its in!
                    break

        self.setInitBlock(ib)
        print "initCell - block radius: ",ib.radius()
        return
    
    def setInitBlock(self,block):
        """
        Set the init block
        """
        self.maxAtomR=block.maxAtomRadius()
        self.initBlock = block
        return
        
        
    def XXgrowBlock(self):
        
        # For comparing angles
        small = 1.0e-8
        
        block1 = self.blocks[ 0 ]
        block1EndGroup = block1.coords[1]
        block1Contact = block1.coords[0]
        
        block2 = self.blocks[ 1 ]
        block2EndGroup = block2.coords[1]
        block2Contact = block2.coords[0]
        
        
        # Move block1Contact to center so that block1Contact -> block1EndGroup defines
        # the alignment we are after
        block1.translate( -block1Contact )
        
        # Same for block2
        block2.translate( -block2Contact )
        
        
        # Find the angle this makes with the target orientation in the XY direction
        newXY = block2EndGroup.copy()
        newXY[2] = 0
        targetXY = block1EndGroup.copy()
        targetXY[2] = 0
        angle = util.vectorAngle(newXY, targetXY)
        if not ( ( -small < angle < small ) or  ( numpy.pi-small < angle < numpy.pi+small ) ):
            rotAxis = numpy.cross(newXY,targetXY,)
            rotAxis = rotAxis/numpy.linalg.norm(rotAxis)
            print "rot XY {}".format(rotAxis)
            print "GOT ANGLE XY: {}".format( angle*util.RADIANS2DEGREES )
            block2.rotate( rotAxis, -angle )
        
        # Check it worked
        # Get new coordinate - can't use old as the object as the rotation changes the objects
        block2EndGroup = block2.coords[ 1 ]
        newXY = block2EndGroup.copy()
        newXY[2] = 0
        angle = util.vectorAngle(newXY, targetXY)
        print "GOT ANGLEXY AFTER: {}".format( angle*util.RADIANS2DEGREES)
        
        
        # Find the angle this makes with the target orientation in the YZ direction
        #block2EndGroup = block2.coords[ 1 ]
        newYZ = block2EndGroup.copy()
        newYZ[1] = 0
        targetYZ = block1EndGroup.copy()
        targetYZ[1] = 0
        angle = util.vectorAngle(newYZ, targetYZ)
        if not ( ( -small < angle < small ) or  ( numpy.pi-small < angle < numpy.pi+small ) ):
            print "GOT ANGLE YZ B4: {}".format( angle*util.RADIANS2DEGREES )
            rotAxis = numpy.cross(newYZ,targetYZ,)
            rotAxis = rotAxis/numpy.linalg.norm(rotAxis)
            print "rot YZ {}".format(rotAxis)
            print "GOT ANGLE YZ: {}".format( angle*util.RADIANS2DEGREES )
            block2.rotate( rotAxis, -angle )
        
        # Check it worked
        # Get new coordinate - can't use old as the object as the rotation changes the objects
        block2EndGroup = block2.coords[ 1 ]
        newYZ = block2EndGroup.copy()
        newYZ[2] = 0
        angle = util.vectorAngle(newYZ, targetYZ)
        print "GOT ANGLEYZ AFTER: {}".format( angle*util.RADIANS2DEGREES)
 
 
        # Find the angle this makes with the target orientation in the XZ direction
        #block2EndGroup = block2.coords[ 1 ]
        newXZ = block2EndGroup.copy()
        newXZ[1] = 0
        targetXZ = block1EndGroup.copy()
        targetXZ[1] = 0
        angle = util.vectorAngle(newXZ, targetXZ)
        if not ( ( -small < angle < small ) or  ( numpy.pi-small < angle < numpy.pi+small ) ):
            print "GOT ANGLE XZ: B4 {}".format( angle*util.RADIANS2DEGREES )
            rotAxis = numpy.cross(newXZ,targetXZ,)
            rotAxis = rotAxis/numpy.linalg.norm(rotAxis)
            print "rot XZ {}".format(rotAxis)
            print "GOT ANGLE XZ: {}".format( angle*util.RADIANS2DEGREES )
            block2.rotate( rotAxis, -angle )
        
        # Check it worked
        # Get new coordinate - can't use old as the object as the rotation changes the objects
        block2EndGroup = block2.coords[ 1 ]
        newXZ = block2EndGroup.copy()
        newXZ[2] = 0
        angle = util.vectorAngle(newXZ, targetXZ)
        print "GOT ANGLE XZ AFTER: {}".format( angle*util.RADIANS2DEGREES)
        
        ##End growBlock
        

    def _recreateBlocks(self, labels, coords):
        """Given a list of labels and coordinates, recreate the blocks
        Uses a stoopid but simple connectivity algorithm
        """
        
        toler=0.5
        
        bonds = []
        
        for i,i_coord in enumerate( coords ):
            i_symbol = util.label2symbol( labels[i] )
            #zi = util.SYMBOL_TO_NUMBER[ i_symbol ]
            #ri = util.COVALENT_RADII[ zi ]
            
            # Start from j so we don't loop over the same atoms twice
            for j in range( i ):
                j_coord = coords[ j ]
                j_symbol = util.label2symbol( labels[j] )
                #zj = util.SYMBOL_TO_NUMBER[ j_symbol ]
                #rj = util.COVALENT_RADII[ zj ]
                bondLength = util.bondLength(i_symbol, j_symbol) + toler
                dist = numpy.linalg.norm( i_coord - j_coord )
                if dist <= bondLength:
                    #print "bondLength({}|{}): {}".format(i_symbol,j_symbol,bondLength)
                    #print "dist: {}".format(dist)
                    print "Added bond {}:{}".format(i,j)
                    # Got bond so add it to the list
                    bonds.append( (i,j) )
                    
        print "GOT {} BONDS".format(bonds)
                    
        #End loop through coordinates
        
        # Add first atom to first block
        blocks=[]
        blocks.append( [bonds[0][0]] )
        
        # Now loop through the bonds and create the blocks
        for (i,j) in bonds:
            found = False
            for k,block in enumerate(blocks):
                #print "Checking block[{}] {}".format(k,block)
                if i in block:
                    #print "{} in block but {} not so adding {}".format(i,j,j)
                    block.append( j )
                    found=True
                    break
                elif j in block:
                    block.append( i )
                    #print "{} in block but {} not so adding {}".format(j,i,i)
                    found=True
                    break
                elif j in block and i in block:
                    raise RuntimeError, "GOT TWO BONDS FOR INDEXES {}:{}".format(i,j)
            #End loop over blocks
        
            if not found:
                # New block
                blocks.append( [i,j] )
                #print "Added atoms {} and {} to NEW BLOCK {}".format(i,j,len(blocks)-1)
 
        print "Got blocks: {}".format(len(blocks))
        
        # Now create the new blocks
        # for block in blocks:
        
        
    def randomBlockId(self,count=1):
        """Return count random block ids"""
        if count == 1:
            return random.choice( list( self.blocks.keys() ) )
        else:
            ids = []
            while len(ids) < count:
                id = random.choice( list( self.blocks.keys() ) )
                if id not in ids:
                    ids.append(id)
            return ids

        
    
    def randomMove(self, block, buffer=None ):
        """Randomly move the given block
         Defintely needs more work on the rotation
         If buffer is given, use this as a buffer from the edges of the cell
         when selecting the coord
        """
        
        # Get coords of random point in the cell
        if buffer:
            x = random.uniform(buffer,self.A[0]-buffer)
            y = random.uniform(buffer,self.B[1]-buffer)
            z = random.uniform(buffer,self.C[2]-buffer)
        else:
            x = random.uniform(0,self.A[0])
            y = random.uniform(0,self.B[1])
            z = random.uniform(0,self.C[2])
            
        coord = numpy.array([x,y,z], dtype=numpy.float64 )
        
        #print "Got random coord: {}".format(coord)
        
        # Move to origin, rotate there and then move to new coord
        # Use the cell axis definitions
        origin = numpy.array([0,0,0], dtype=numpy.float64 )
        block.translateCentroid( origin )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate( self.A, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate( self.B, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate(self.C, angle )
        
        # Now move to new coord
        block.translateCentroid( coord )
        
        #print "After rotate centroid at: {}".format( block.centroid() )
        
        
    def randomMoveAroundBlock(self, oblock, block, margin=1.0 ):
        
        """
        Move the subject block (block) around the centroid of the target oblock
        """
        
        centroid = oblock.centroid()
        #print "centroid at: {}".format(centroid)
        rb = block.radius()
        ro = oblock.radius()
        
        prange = rb + ro + margin
        prange = prange/2
        
        # Calculate new coord
        x = random.uniform(-prange,prange)
        y = random.uniform(-prange,prange)
        z = random.uniform(-prange,prange)
        xyz = numpy.array( [x,y,z], dtype=numpy.float64 )
        coord = numpy.add( centroid, xyz )
        
        # possibly ovverklll - move to origin, rotate there and then move to new coord
        # Use the cell axis definitions
        origin = numpy.array([0,0,0], dtype=numpy.float64 )
        block.translateCentroid( origin )
        
        # Make a random rotation
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate( self.A, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate( self.B, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate(self.C, angle )
        
        # Now move to new coord
        block.translateCentroid( coord )
        #print "block now at: {}".format( block.centroid() )
        


    def seed( self, nblocks, inputFile ):
        """ Seed a cell with nblocks based on the block that will be created
        from the input block
        """
        
        if self.A == None or self.B == None or self.C == None:
            raise RuntimeError,"Need to specify cell before seeding"
        
        self.initCell( inputFile )
        
        # If it's just one, add the given block
        self.addBlock( self.initBlock.copy() )
        if nblocks == 1:
            return
        
        MAXTRIES = 1000
        # Loop through the nblocks adding the blocks to
        # the cell - nblocks-1 as we've already added the first
        for seedCount in range( nblocks-1 ):
            # Create new block
            newblock = self.initBlock.copy()
            tries = 0
            ok = False
            print "Adding block: {}".format(seedCount+1)
            while not ok:
                # quit on maxTries
                if tries >= MAXTRIES:
                    print "EXCEEDED MAXTRIES WHEN SEEDING"
                    sys.exit(1)
                
                # Move the block and rotate it
                self.randomMove(newblock)
                #print "RANDOM TO MOVE TO: {}".format( newblock.centroid() )
                
                #Add the block so we can check for clashes/bonds
                iblock = self.addBlock(newblock)
                
                # Test for Clashes with other molecules
                ok = self.newCheckMove( iblock )
                
                # Break out if no clashes
                if ok:
                    break
                
                # increment tries counter
                tries += 1
            
            # End Clash loop
        # End of loop to seed cell
    # End seed
    
    def newCheckMove(self,iblock):
        """
        See what happened with this move
        """
        # Get a list of the close atoms
        close = self.closeAtoms(iblock)
        
        if not close:
            # Nothing to see so move along
            return True
        
        # Get the block
        block = self.blocks[iblock]
        
        bonds = [] # list of possible bond atoms
        for ( iatom, ioblock, ioatom ) in close:
            symbol = block.symbols[iatom]
            radius = block.atom_radii[iatom]
            oblock = self.blocks[ioblock]
            coord = block.coords[iatom]
            ocoord = oblock.coords[ioatom]
            # First see if both atoms are endGroups
            if iatom in block.endGroups and ioatom in oblock.endGroups:
                # Checking for a bond
                # NB ASSUMPION FOR BOND LENGTH CHECK IS BOTH BLOCKS HAVE SAME ATOM TYPES
                osymbol = oblock.symbols[ioatom]
                
                # This uses -the 'optimised' bond lenght check - probably uneeded
                bond_length = block.bondLength( symbol, osymbol )
                
                # THINK ABOUT BETTER SHORT BOND LENGTH CHECK
                if  bond_length - self.bondMargin < self.distance( coord, ocoord ) < bond_length + self.bondMargin:
                    
                    print "possible bond for ",iatom,ioblock,ioatom
                    # Possible bond so check the angle
                    icontact = block.endGroupContactIndex( iatom )
                    contact = block.coords[icontact]
                    angle = block.angle( contact, coord, ocoord )
                    
                    print "Got bond angle D: {}".format(angle)
                    
                    if ( self.bondAngle-self.bondAngleMargin < angle < self.bondAngle+self.bondAngleMargin ):
                        bonds.append( iatom, ioblock, ioatom )
                    else:
                        print "Cannot bond due to angle"
                        
            #Finished here so check the next ones
            continue
           
            # No bond so just check if the two atoms are close enough for a clash
            oradius = oblock.atom_radii[ioatom]
            if self.distance( coord, ocoord ) < radius+oradius+self.atomMargin:
                #print "ATOMS CLASH ",iatom,ioblock,ioatom
                return False
    
        # Here no atoms clash and we have a list of possible bonds - so bond'em!
        if len(bonds):
            for iatom, ioblock, ioatom in bonds:
                oblock = self.blocks[ioblock]
                block.bond( oblock, (iatom, ioatom) )
                # Remove the block from the list as it is now part of the other one
                self.delBlock(oblock)
            print "ADDED {} bonds".format(len(bonds))
        
        # Either got bonds or no clashes
        return True

        
    def shimmy(self, nsteps=100, nmoves=50):
        """ Shuffle the molecules about making bonds where necessary for nsteps
        minimoves is number of sub-moves to attempt when the blocks are close
        """
        
        for step in range( nsteps ):
            
            if not step % 20:
                print "step {}".format(step)
                
            iblock = self.randomBlockId()
            block = self.blocks[iblock]
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
            
            # Remove the block from the cell so we don't check against itself
            self.delBlock(iblock)
             
            # Make a random move
            self.randomMove( block )
            
            #Add the block so we can check for clashes/bonds
            new_iblock = self.addBlock(block)
            
            # Test for Clashes with other molecules
            ok = self.newCheckMove( new_iblock )
            
            # Break out if no clashes
            if not ok:
                # Put it back where we got it from
                self.delBlock(new_iblock)
                self.addBlock(orig_block)
            
        #End shimmy
            
    

            
    def OLDshimmy(self, nsteps=100, nmoves=50, bondAngle=None ):
        """ Shuffle the molecules about making bonds where necessary for nsteps
        minimoves is number of sub-moves to attempt when the blocks are close
        """
        
        # For time being ensure we are given a bond angle
        assert bondAngle
        
        PRANGE=4.0 # the range within which to make the smaller minimoves
        
        for step in range( nsteps ):
            
            if not step % 20:
                print "step {}".format(step)
                
            iblock = self.randomBlockId()
            
            block = self.blocks[iblock]
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
             
            # Make a random move
            self.randomMove( block )
            
            # Remember the original cog of the block 
            cog = block.centroid()
            
            # Loop over all blocks to see if we can bond or if we clash
            removed = []
            for i in range( len( self.blocks ) ):
                
                # skip the block we are using
                if i == iblock:
                    continue
                
                # Get the next block
                oblock = self.blocks[i]
                
                # First see if we are close enough to consider bonding
                if not block.close( oblock, blockMargin = self.blockMargin ):
                    continue
                
                # See if we can bond and shimmy minimoves times
                # Get the original cog so we sample about it
                clashmove=0
                noclashmove=0
                for j in range( nmoves ):
                    # Try to make a bond
                    bond = block.canBond( oblock, bondAngle=bondAngle )
                    if not bond or bond == "clash":
                        # Try a smaller move
                        self.randomMoveAroundBlock( block, cog, prange=PRANGE )
                        # Just for logging
                        if bond == "clash":
                            clashmove+=1
                        else:
                            noclashmove+=1
                        
                    else:
                        # Bonding worked so break out of loop
                        print "GOT BOND"
                        break
                
                print "Next Block after clash/noclash: {}/{}".format(clashmove,noclashmove)
                if bond:
                    if bond == "clash":
                        # Reject this move and overwrite the changed block with the original one
                        self.blocks[i] = orig_block
                    else:
                        block.bond( oblock, bond )
                        print "Bonded block {} with block {}".format( iblock, i)
                        # Now delete block
                        removed.append(i)
                i+=1
                #End Loop Over Blocks
            
            #
            for r in removed:
                self.blocks.pop(r)
            #End step loop

             
    def writeXyz(self, ofile, label=False ):
        """Write out the cell atoms to an xyz file
        If label is true we write out the atom label and block, otherwise the symbol
        """
        
        natoms=0
        xyz = ""
        for i,block in self.blocks.iteritems():
            for j, c in enumerate( block.coords ):
                if label:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( "{}_block#{}".format(block.labels[j], i), c[0], c[1], c[2] )
                else:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.symbols[j], c[0], c[1], c[2] )
                natoms += 1
        
        # Write out natoms and axes as title
        xyz = "{}\nAxes:{}:{}:{}\n".format(natoms, self.A[0], self.B[1], self.C[2] ) + xyz
        
        with open( ofile, 'w' ) as f:
            fpath = os.path.abspath(f.name)
            f.writelines( xyz )
            
        print "Wrote cell file: {0}".format(fpath)
    
    def __str__(self):
        """
        """
        s = ""
        for block in self.blocks.values():
            s+= str(block)
            
        return s
    
    
class TestCell(unittest.TestCase):
    
    # Import only here as only needed for testing
    import buildingBlock

    def makeCh4(self, cell):
        """Create a CH4 molecule for testing"""
        
        coords = [ numpy.array([  0.000000,  0.000000,  0.000000 ] ),
        numpy.array([  0.000000,  0.000000,  1.089000 ]),
        numpy.array([  1.026719,  0.000000, -0.363000 ]),
        numpy.array([ -0.513360, -0.889165, -0.363000 ]),
        numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]

        #numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]
        labels = [ 'C', 'H', 'H', 'H', 'H' ]
        
        endGroups = [ 1,2,3,4 ]
        endGroupContacts = { 1:0, 2:0, 3:0, 4:0 }
        
        ch4 = self.buildingBlock.BuildingBlock()
        ch4.createFromArgs( coords, labels, endGroups, endGroupContacts )
        return ch4
    
    def makePaf(self, cell):
        """Return the PAF molecule for testing"""
        
        paf = self.buildingBlock.BuildingBlock()
        paf.fromCarFile("../PAF_bb_typed.car")
        return paf
        
    def XtestAlignBlocks(self):
        """Test we can align two blocks correctly"""
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        #block = self.makeCh4( cell )
        cell.seed( 2, "../PAF_bb_typed.car" )
        
        # Get the two blocks
        block1 = cell.blocks[0]
        block2 = cell.blocks[1]
        
        # Get two end Groups
        block1EndGroupIndex=7
        block2EndGroupIndex=17
        
        # Get the atoms that define things
        block1ContactIndex = block1.endGroupContactIndex( block1EndGroupIndex )
        block1Contact = block1.coords[ block1ContactIndex ]
        block2ContactIndex = block2.endGroupContactIndex( block2EndGroupIndex )
        block2Contact = block2.coords[ block2ContactIndex ]
        
        # Move block1Contact to center so that block1Contact -> block1EndGroup defines
        # the alignment we are after
        block1.translate( -block1Contact )
        
        # Same for block2
        block2.translate( -block2Contact )
        
        cell.alignBlocks( block1, block1EndGroupIndex, block2, block2EndGroupIndex )
        
        # Check the relevant atoms are in the right place
        block1 = cell.blocks[0]
        block2 = cell.blocks[1]
        block1EndGroup = block1.coords[ block1EndGroupIndex ]
        block2EndGroup = block2.coords[ block2EndGroupIndex ]
        block1Contact = block1.coords[  block1ContactIndex ]
        block2Contact = block2.coords[  block2ContactIndex ]
        
        print "{}\n{}".format( block1EndGroup, block2EndGroup)
        print "{}\n{}".format( block1Contact, block2Contact)
        
        # Shocking tolerances - need to work out why...
        self.assertTrue( numpy.allclose( block1EndGroup, block2EndGroup, rtol=1e-1, atol=1e-1 ),
                         msg="End Group incorrectly positioned")
        self.assertTrue( numpy.allclose( block1Contact, block2Contact, rtol=1e-7, atol=1e-7 ),
                         msg="Contact atom incorrectly positioned")

    def XtestCellIO(self):
        """Check we can write out and then read in a cell
        """
        
        nblocks = 4
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        cell.seed( nblocks, "../PAF_bb_typed.car" )
        
        # Remember a coordinate for checking
        test_coord = cell.blocks[3].coords[4]
        
        self.assertEqual( nblocks, len(cell.blocks), "Incorrect number of cell blocks at start: {}".format( len(cell.blocks) ))
        
        outfile = "./testCell.xyz"
        cell.writeXyz( outfile, label=True )
        
        newCell = Cell()
        
        newCell.fromXyz( outfile )
        self.assertEqual( nblocks, len(newCell.blocks), "Incorrect number of cell blocks after read: {}".format(len(newCell.blocks) ))
        
        self.assertTrue( numpy.allclose( test_coord, cell.blocks[3].coords[4], rtol=1e-9, atol=1e-9 ),
                         msg="Incorrect testCoordinate of cell.")
        
        
    def XtestDistance(self):
        """Test the distance under periodic boundary conditions"""
        
        CELLA = 10
        CELLB = 10
        CELLC = 10
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        v1 = numpy.array([ 1.0, 1.0, 1.0 ])
        v2 = numpy.array([ 5.0, 5.0, 5.0 ])
        
        dc = cell.distance(v1,v2)
        dn = numpy.linalg.norm(v2-v1)
        self.assertEqual( dc, dn, "Distance within cell:{} | {}".format(dc,dn) )
        
        v1 = numpy.array([ 0.0, 0.0, 0.0 ])
        v2 = numpy.array([ 0.0, 0.0, 8.0 ])
        dc = cell.distance(v1,v2)
        self.assertEqual( dc, 2.0, "Distance within cell:{} | {}".format(dc,dn) )
        
    
    def XtestGrowBlock(self):
        """Test we can add a block correctly"""
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        #block = self.makePaf( cell )
        cell.seed( 1, "../PAF_bb_typed.car" )
        for _ in range( 10 ):
            ok = cell.growBlock( cell.initBlock.copy() )
            if not ok:
                print "Failed to add block"
            
        cell.writeXyz("JENS.xyz", label=False)
        
    def testSeed(self):
        """Test we can seed correctly"""
        
        
        nblocks = 10
        CELLA = 50
        CELLB = 50
        CELLC = 50
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        #ch4 = self.makeCh4( cell )
        cell.seed( nblocks, "../PAF_bb_typed.car" )
        
        self.assertEqual( nblocks, len(cell.blocks), "Incorrect number of cell blocks" )
        
        # Check no block is outside the unit cell
        # We seed by the centroid so the edges could stick out by radius
        #radius = ch4.radius()
        radius = cell.initBlock.radius()
        
        bad = []
        for i,b in cell.blocks.iteritems():
            for c in b.coords:
                if ( 0-radius > c[0] > CELLA[0]+ radius ) or \
                   ( 0-radius > c[1] > CELLB[1]+ radius ) or \
                   ( 0-radius > c[2] > CELLC[2]+ radius ):
                    
                    bad.append( b )
        
        self.assertEqual( 0, len(bad), "Got {} blocks outside cell: {}".format( len(bad), bad ) )
        
        
        
    def XtestCloseAtoms(self):
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        #block1 = self.makePaf( cell )
        cell.initCell("../ch4_typed.car",incell=False)
        block1 = cell.initBlock.copy()
        
        b1 = numpy.array([2,2,2], dtype=numpy.float64 )
        block1.translateCentroid( b1 )
        cell.addBlock(block1)
        
        block2=cell.initBlock.copy()
        b2 = numpy.array([3,3,3], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        cell.addBlock(block2)
        
        cell.initCellLists()
        self.assertEqual(cell.closeAtoms(0),
                         [(0, 1, 3), (0, 1, 4), (0, 1, 2), (0, 1, 0), (0, 1, 1), (1, 1, 3), (1, 1, 4),
                          (1, 1, 0), (1, 1, 1), (1, 1, 2), (2, 1, 2), (2, 1, 0), (2, 1, 1), (2, 1, 3),
                          (2, 1, 4), (3, 1, 3), (4, 1, 4), (4, 1, 3), (4, 1, 2), (4, 1, 0)]
                         ,"Many contacts")
        
        
        b2 = numpy.array([10,10,10], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        cell.initCellLists()
        self.assertEqual(cell.closeAtoms(0), None, "No contacts")
        
        b2 = numpy.array([29,2,2], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        cell.initCellLists()
        self.assertEqual(cell.closeAtoms(0), [(0, 1, 2), (1, 1, 2), (3, 1, 0), (3, 1, 2), (4, 1, 0), (4, 1, 2)], "Periodic boundary")
        
        
        #cell.writeXyz("jens.xyz", label=False)
        
    def XtestCellLists(self):
        
        CELLA = 10
        CELLB = 10
        CELLC = 10
        
        cell = Cell( atomMargin=0.1)
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        b1 = self.buildingBlock.BuildingBlock()
        coords = [ numpy.array([  0.050000,  0.000000,  0.000000 ] ) ]
        labels = [ 'C' ]
        endGroups = [ 1 ]
        endGroupContacts = { 1:1 }
        b1.createFromArgs( coords, labels, endGroups, endGroupContacts )
        
        b2 = self.buildingBlock.BuildingBlock()
        coords = [ numpy.array([  1.000000,  0.000000,  0.000000 ] ) ]
        b2.createFromArgs( coords, labels, endGroups, endGroupContacts )
        
        cell.blocks.append(b1)
        cell.blocks.append(b2)
        
        cell.initCellLists()
        print "first"
        cell.closeAtoms(1)
#        
#        del cell.blocks[1]
#        b2 = self.buildingBlock.BuildingBlock( cell.distance )
#        coords = [ numpy.array([  7.400000,  2.000000,  2.000000 ] ) ]
#        b2.createFromArgs( coords, labels, endGroups, endGroupContacts )
#        cell.blocks.append(b2)
#        cell.initCellLists()
#        print "second"
#        cell.newClash(1)
#        
#        del cell.blocks[1]
#        b2 = self.buildingBlock.BuildingBlock( cell.distance )
#        coords = [ numpy.array([  9.950000,  1.000000,  1.000000 ] ) ]
#        b2.createFromArgs( coords, labels, endGroups, endGroupContacts )
#        cell.blocks.append(b2)
#        cell.initCellLists()
#        print "third"
#        cell.newClash(1)


if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()