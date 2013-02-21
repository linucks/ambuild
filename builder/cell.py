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


    def __init__( self ):
        '''
        Constructor
        '''
        
        # A, B and C are cell vectors - for the time being we assume they are orthogonal
        self.A = None
        self.B = None
        self.C = None
        
        # first block is kept separate as everything is built up from it
        self.initBlock = None
        self.blocks = []
        
            
        
    
    def addBlock( self, block ):
        self.blocks.append( block )
        
        
    def alignBlocks(self, block1, block1EndGroupIndex, block2, block2EndGroupIndex):
        """
        Align block2, so that the bond defined by the endGroup on block 2 is aligned with
        the bond defined by the endGroup on block 1
        
        This assumes that the blocks have already been positioned with the contact atoms at the origin
        """
        
        block1EndGroup = block1.position( block1EndGroupIndex )
        
        # Aign with each axis in turn. The axis to align to is just the position of bock1EndGroup
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
        blockEndGroup = block.position( endGroupIndex )
        bAngle = util.vectorAngle(blockEndGroup, axis)
        #print "bAngle b4: {}".format(bAngle)
        
        # Need to rotate block on X-axis by tvAngle-bAngle
        angle = tvAngle-bAngle
        
        #print "Got angle: {}".format(angle*util.RADIANS2DEGREES)
        
        if not ( ( -small < angle < small ) or  ( numpy.pi-small < angle < numpy.pi+small ) ):
            rotAxis = numpy.cross( blockEndGroup, axis)
            #print "rotAxis {}".format(rotAxis)
            block.rotate( rotAxis, angle )
        
        #blockEndGroup = block.position( endGroupIndex )
        #bAngle = util.vectorAngle(blockEndGroup, axis)
        #print "bAngle after: {}".format(bAngle)
        #print block
        
    def cellAxis(self,A=None,B=None,C=None):
        """
        Get or set the cell axes
        """
        if A and B and C:
            # A, B and C are cell vectors - for the time being we assume they are orthogonal
            assert A[1] == 0 and A[2] == 0, "Cell vectors need to be orthogonal and in correct order"
            assert B[0] == 0 and B[2] == 0, "Cell vectors need to be orthogonal and in correct order"
            assert C[0] == 0 and C[1] == 0, "Cell vectors need to be orthogonal and in correct order"
            self.A = numpy.array( A, dtype=numpy.float64 )
            self.B = numpy.array( B, dtype=numpy.float64 )
            self.C = numpy.array( C, dtype=numpy.float64 )
        else:
            return (A,B,C)
        
    def checkClash(self, newblock, ignore=[] ):
        """
        See if the given block clashes with any of the others
        return True if there is a clash
        ignore is a list of block indices to ignore
        """
        
        for i,block in enumerate(self.blocks):
            if i in ignore:
                continue
            if block.clash( newblock ):
                return True
            
        return False
        
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
            if not block.close( oblock, margin = closeMargin ):
                continue
            
            # See if we can bond
            bond = block.canBond( oblock, bondAngle=bondAngle )
            
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
            iblock, jblock = self.randomBlocks()
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
                # If no bonds place the bond back at it's original position
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

    def findClose(self, oblock, block, closeMargin=None, blockMargin=None ):
        """
        Find the blocks around oblock that might be close to block if we were moving
        it by the parameter used for the random move around block.
        return a list of the indexes
        """
        
        # Need to add the central block radius to the diameter of the sample block
        # plus the margins
        rb = block.radius()
        ro = oblock.radius()
        margin = ro + 2*rb + blockMargin + closeMargin
        
        centroid = oblock.centroid()
        closeBlocks = []
        # This should implicitly include oblock
        for i,b in enumerate( self.blocks ):
            # See if we are close enough
            dist = numpy.linalg.norm( centroid - b.centroid() )
            if dist < margin + b.radius():
                closeBlocks.append( i )
        
        return closeBlocks
        
    def fromXyz(self, xyzFile ):
        """ Read in an xyz file containing a cell and recreate the cell object"""
        
        blocks = self._readXyz(xyzFile)
        for labels,coords in blocks:
            block = buildingBlock.BuildingBlock( self.distance )
            block.createFromLabelAndCoords( labels, coords )
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
        block1index = self.randomBlockIndex()
        block1 = self.blocks[ block1index ]
        
        # pick random endgroup and contact in target
        block1EndGroupIndex = block1.randomEndGroupIndex()
        # Need to copy this so we can remember it - otherwise it changes with the moves
        block1ContactIndex = block1.endGroupContactIndex( block1EndGroupIndex )
        block1Contact = block1.position( block1ContactIndex ).copy()
        
        # pick random endGroup on the new block
        block2EndGroupIndex = block2.randomEndGroupIndex()
        block2ContactIndex = block2.endGroupContactIndex( block2EndGroupIndex )
        block2Contact = block2.position( block2ContactIndex ).copy()
        
        # get the position where the next block should bond
        bondPos = block1.newBondPosition( block1EndGroupIndex, block2, block2EndGroupIndex )
        #print "got bondPos: {}".format( bondPos )
        
        # Move block1Contact to center so that block1Contact -> block1EndGroup defines
        # the alignment we are after
        block1.translate( -block1Contact )
        # Same for block2
        block2.translate( -block2Contact )
        
        # Align along all 3 axes
        self.alignBlocks( block1, block1EndGroupIndex, block2, block2EndGroupIndex)
        
        # Move block1 back to its original position
        block1.translate( block1Contact )
        
        # now turn the second block around
        
        # Need to get the new vectors as these will have changed due to the moves
        block2EndGroup = block2.position( block2EndGroupIndex )
        
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
        
        # Now need to position the endGroup at the bond position - at the moment
        # The contact atom is on the origin so we need to subtract the vector to the endGroup
        # from the translation vector
        
        #jmht FIX!
        block2.translate( bondPos + block2EndGroup )
        
        
        # Check it doesn't clash - could try rotating about the dihedral
        if not self.checkClash( block2, ignore=[block1index] ):
            # Bond the two blocks
            block1.bond(  block2, (block1EndGroupIndex,block2EndGroupIndex) )
            return True
        
        return False
        
    
    def initCell(self,inputFile):
        """
        Read in the first block from the input file, position it so it is in the cell
        and then save it
        """
        
        ib = buildingBlock.BuildingBlock( self.distance )
        ib.fromCarFile( inputFile )
        
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

        self.initBlock = ib
        
        return
    ##End initCell
        
        
    def XXgrowBlock(self):
        
        # For comparing angles
        small = 1.0e-8
        
        block1 = self.blocks[ 0 ]
        block1EndGroup = block1.position(1)
        block1Contact = block1.position(0)
        
        block2 = self.blocks[ 1 ]
        block2EndGroup = block2.position(1)
        block2Contact = block2.position(0)
        
        
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
        block2EndGroup = block2.position( 1 )
        newXY = block2EndGroup.copy()
        newXY[2] = 0
        angle = util.vectorAngle(newXY, targetXY)
        print "GOT ANGLEXY AFTER: {}".format( angle*util.RADIANS2DEGREES)
        
        
        # Find the angle this makes with the target orientation in the YZ direction
        #block2EndGroup = block2.position( 1 )
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
        block2EndGroup = block2.position( 1 )
        newYZ = block2EndGroup.copy()
        newYZ[2] = 0
        angle = util.vectorAngle(newYZ, targetYZ)
        print "GOT ANGLEYZ AFTER: {}".format( angle*util.RADIANS2DEGREES)
 
 
        # Find the angle this makes with the target orientation in the XZ direction
        #block2EndGroup = block2.position( 1 )
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
        block2EndGroup = block2.position( 1 )
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
        
        

            

    def randomBlockIndex(self):
        """Return the index of one of the blocks"""
        return random.randint( 0, len(self.blocks)-1 )
    
    def randomBlocks(self):
        """Return two random block indices"""
        
        # Pick a block to move
        iblock = self.randomBlockIndex()
        
        # Find a different one to sample around
        jblock = iblock
        while jblock == iblock:
            jblock = self.randomBlockIndex()
        
        return (iblock, jblock)

        
    
    def randomMove(self, block, buffer=None ):
        """Randomly move the given block
         Defintely needs more work on the rotation
         If buffer is given, use this as a buffer from the edges of the cell
         when selecting the position
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
            
        position = numpy.array([x,y,z], dtype=numpy.float64 )
        
        #print "Got random position: {}".format(position)
        
        # Move to origin, rotate there and then move to new position
        # Use the cell axis definitions
        origin = numpy.array([0,0,0], dtype=numpy.float64 )
        block.translateCentroid( origin )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate( self.A, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate( self.B, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        block.rotate(self.C, angle )
        
        # Now move to new position
        block.translateCentroid( position )
        
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
        
        # Calculate new position
        x = random.uniform(-prange,prange)
        y = random.uniform(-prange,prange)
        z = random.uniform(-prange,prange)
        xyz = numpy.array( [x,y,z], dtype=numpy.float64 )
        position = numpy.add( centroid, xyz )
        
        # possibly ovverklll - move to origin, rotate there and then move to new position
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
        
        # Now move to new position
        block.translateCentroid( position )
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
            clash = True
            print "Adding block: {}".format(seedCount)
            while clash:
                # quit on maxTries
                if tries >= MAXTRIES:
                    print "EXCEEDED MAXTRIES WHEN SEEDING"
                    sys.exit(1)
                
                # Move the block and rotate it
                self.randomMove(newblock)
                #print "RANDOM TO MOVE TO: {}".format( newblock.centroid() )
                
                # Test for Clashes with other molecules
                clash = self.checkClash( newblock )
                
                # Break out if no clashes
                if not clash:
                    break
                
                # increment tries counter
                tries += 1
            
            # End Seeding loop
            if not clash:
                #print "ADDED BLOCK AT POS: {}".format( newblock.centroid() )
                self.addBlock(newblock)
                #newblock.index(index=seedCount)
            else:
                print "ERROR ADDING BLOCK"
                sys.exit(1)
    
        # End of loop to seed cell
    
    # End seed

        
    def shimmy(self, nsteps=100, nmoves=50, bondAngle=None ):
        """ Shuffle the molecules about making bonds where necessary for nsteps
        minimoves is number of sub-moves to attempt when the blocks are close
        """
        
        # For time being ensure we are given a bond angle
        assert bondAngle

        
        CLOSE_MARGIN=4.0 # how close 2 blocks are before we consider checking if they can bond
        PRANGE=4.0 # the range within which to make the smaller minimoves
        
        
        for step in range( nsteps ):
            
            if not step % 20:
                print "step {}".format(step)
                
            iblock = self.randomBlockIndex()
            block = self.blocks[iblock]
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
             
            # Make a random move
            self.randomMove( block )
            
            # Loop over all blocks to see if we can bond or if we clash
            bonds = self.checkMove( iblock, closeMargin=CLOSE_MARGIN, bondAngle=bondAngle )
            
            # See what happend
            if bonds:
                if bonds == "clash":
                    # Reject this move and overwrite the changed block with the original one
                    self.blocks[iblock] = orig_block
                elif len( bonds ):
                    # Got some bonds so deal with them
                    # Need to decrement index if we are removing blocks
                    bcount=0
                    for b,bond in bonds:
                        block.bond( self.blocks[b-bcount], bond )
                        # Remove the block from the list as it is now part of the other one
                        self.blocks.pop(b-bcount)
                        bcount+=1
            #End step loop
            
        #End shimmy
            
    

            
    def OLDshimmy(self, nsteps=100, nmoves=50, bondAngle=None ):
        """ Shuffle the molecules about making bonds where necessary for nsteps
        minimoves is number of sub-moves to attempt when the blocks are close
        """
        
        # For time being ensure we are given a bond angle
        assert bondAngle

        
        CLOSE_MARGIN=1.5 # how close 2 blocks are before we consider checking if they can bond
        PRANGE=4.0 # the range within which to make the smaller minimoves
        
        for step in range( nsteps ):
            
            if not step % 20:
                print "step {}".format(step)
                
            iblock = self.randomBlockIndex()
            
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
                if not block.close( oblock, margin = CLOSE_MARGIN ):
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
        for i,block in enumerate(self.blocks):
            for j, c in enumerate( block.coords ):
                if label:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( "{}_block#{}".format(block.labels[j], i), c[0], c[1], c[2] )
                else:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.symbols[j], c[0], c[1], c[2] )
                natoms += 1
        
        # Write out natoms and axes as title
        xyz = "{}\nAxes:{}:{}:{}\n".format(natoms, self.A.tolist(), self.B.tolist(), self.C.tolist() ) + xyz
        
        with open( ofile, 'w' ) as f:
            fpath = os.path.abspath(f.name)
            f.writelines( xyz )
            
        print "Wrote cell file: {0}".format(fpath)
    
    def __str__(self):
        """
        """
        s = ""
        for block in self.blocks:
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
        
        ch4 = self.buildingBlock.BuildingBlock( cell.distance )
        ch4.createFromArgs( coords, labels, endGroups, endGroupContacts )
        return ch4
    
    def makePaf(self, cell):
        """Return the PAF molecule for testing"""
        
        paf = self.buildingBlock.BuildingBlock( cell.distance )
        paf.fromCarFile("../PAF_bb_typed.car")
        return paf
        
    def testAlignBlocks(self):
        """Test we can align two blocks correctly"""
        
        CELLA = [ 30,  0,  0 ]
        CELLB = [ 0, 30,  0 ]
        CELLC = [ 0,  0, 30 ]
        
        #b1.symbols[1] = 'F'

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
        block1Contact = block1.position( block1ContactIndex )
        block2ContactIndex = block2.endGroupContactIndex( block2EndGroupIndex )
        block2Contact = block2.position( block2ContactIndex )
        
        # Move block1Contact to center so that block1Contact -> block1EndGroup defines
        # the alignment we are after
        block1.translate( -block1Contact )
        
        # Same for block2
        block2.translate( -block2Contact )
        
        cell.alignBlocks( block1, block1EndGroupIndex, block2, block2EndGroupIndex )
        
        # Check the relevant atoms are in the right place
        block1 = cell.blocks[0]
        block2 = cell.blocks[1]
        block1EndGroup = block1.position( block1EndGroupIndex )
        block2EndGroup = block2.position( block2EndGroupIndex )
        block1Contact = block1.position(  block1ContactIndex )
        block2Contact = block2.position(  block2ContactIndex )
        
        print "{}\n{}".format( block1EndGroup, block2EndGroup)
        print "{}\n{}".format( block1Contact, block2Contact)
        
        # Shocking tolerances - need to work out why...
        self.assertTrue( numpy.allclose( block1EndGroup, block2EndGroup, rtol=1e-1, atol=1e-1 ),
                         msg="End Group incorrectly positioned")
        self.assertTrue( numpy.allclose( block1Contact, block2Contact, rtol=1e-7, atol=1e-7 ),
                         msg="Contact atom incorrectly positioned")

    def testCellIO(self):
        """Check we can write out and then read in a cell
        """
        
        nblocks = 4
        CELLA = [ 30,  0,  0 ]
        CELLB = [ 0, 30,  0 ]
        CELLC = [ 0,  0, 30 ]
        
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        #paf = self.makePaf( cell )
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
        
        
    def testDistance(self):
        """Test the distance under periodic boundary conditions"""
        
        CELLA = [ 10,  0,  0 ]
        CELLB = [ 0, 10,  0 ]
        CELLC = [ 0,  0, 10 ]
        
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
        
    
    def testGrowBlock(self):
        """Test we can add a block correctly"""
        
        CELLA = [ 30,  0,  0 ]
        CELLB = [ 0, 30,  0 ]
        CELLC = [ 0,  0, 30 ]
        
        #block = self.makeCh4()
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
        CELLA = [ 50,  0,  0 ]
        CELLB = [ 0, 50,  0 ]
        CELLC = [ 0,  0, 50 ]
        
        
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
        for b in cell.blocks:
            for c in b.coords:
                if ( 0-radius > c[0] > CELLA[0]+ radius ) or \
                   ( 0-radius > c[1] > CELLB[1]+ radius ) or \
                   ( 0-radius > c[2] > CELLC[2]+ radius ):
                    
                    bad.append( b )
        
        self.assertEqual( 0, len(bad), "Got {} blocks outside cell: {}".format( len(bad), bad ) )
        

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()