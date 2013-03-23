'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import collections
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


    def __init__( self, atomMargin=0.5, boxMargin=1.0, bondMargin=0.5, bondAngle=180, bondAngleMargin=15 ):
        '''
        Constructor
        '''

        # A, B and C are cell vectors - for the time being we assume they are orthogonal
        self.A = None
        self.B = None
        self.C = None
        
        # additional distance to add on to the characteristic bond length
        # when checking whether 2 atoms are close enough to bond
        self.bondMargin = bondMargin
        
        # additional distance to add on when checking if two blocks are close enough to require
        # checking their interactions
        # CURRENTLY UNUSED
        self.blockMargin = 2.0
        
        # additional distance to add on to the atom covalent radii when checking if two atoms 
        # are close enough to clash
        self.atomMargin = atomMargin
        
        # additional distance to add on to the atom covalent radii when checking if two atoms 
        # are close enough to add them to the interatction boxes
        self.boxMargin = boxMargin
        
        # The acceptable bond angle
        self.bondAngle = bondAngle
        
        self.bondAngleMargin = bondAngleMargin
        
        # Dict mapping box key to a list of tuples defining the atoms in that box
        self.box1={}
        
        # Dict mapping key to list of boxes surrounding the keyed box
        self.box3={}
        
        # max atom radius - used to calculate box size
        self.boxSize = None
        self.maxAtomR=None
        # number of boxes in A,B,C axes - used for calculating PBC
        self.numBoxA = None
        self.numBoxB = None
        self.numBoxC = None
        
        # first block is kept separate as everything is built up from it
        self.initBlock = None
        
        # dictionary mapping id of the block to the block - can't use a list and indices
        # as we add and remove blocks and need to keep track of them
        self.blocks = {}
    
    def addBlock( self, block):
        """
        Add the block and the calculate the boxes for all atoms.
        
        Args:
        block -- block to add
        blockid -- the id to use for the block - otherwise we use the python id
        
        Returns:
        index of the new block in the list
        
        jmht  - rather then the check here, we could create a copy of the coords of the 
        block in the cell and then use these when checking distances - this way we wouldn't
        need to use PBC in the distance check, which would possibly be quicker at the expense of
        storing 2 sets of coordinates, which is negligible - test!
        """
        
        # The index of the new block
        blockid = id(block)
        
        # Add to the dict
        self.blocks[ blockid ] = block
        
        #print "nbox ",self.numBoxA,self.numBoxB,self.numBoxC
        for icoord,coord in enumerate(block.coords):
            
            #print "ADDING COORD ",coord
            
            # Periodic Boundaries
            x = coord[0] % self.A[0]
            y = coord[1] % self.B[1]
            z = coord[2] % self.C[2]
            
            a=int( math.floor( x / self.boxSize ) )
            b=int( math.floor( y / self.boxSize ) ) 
            c=int( math.floor( z / self.boxSize ) )
            
            key = (a,b,c)
            #print "NEW KEY ",key
            block.atomCell[icoord] = key
            try:
                self.box1[key].append( (blockid,icoord) )
                #print "APPENDING KEY ",key,blockid,icoord
            except KeyError:
                #print "ADDING KEY ",key,blockid,icoord
                # Add to main list
                self.box1[key] = [(blockid,icoord)]
                # Map surrounding boxes
                self.box3[key] = self.surroundBoxes(key)
                #surround = self.surroundBoxes(key)
                #print "SURROUND ",surround
                #self.box3[key] = surround
        return blockid
    
    def alignBlocks(self, block, iblockEndGroup, refVector ):
        """
        Align block2, so that the bond defined by the endGroup on block 2 is aligned with
        the refVector
        
        This assumes that the block has already been positioned with the contact atoms at the origin
        """
        
        blockEndGroup = block.coords[ iblockEndGroup ]
        
        print "alignBlocks BEFORE"
        print blockEndGroup
        print refVector
        
        # Makes no sense if they are the same already...
        if numpy.array_equal(blockEndGroup,refVector):
            print "alignBlocks - block already aligned along vector. May not be a problem, but you should know..."
            return
        
        # Calculate normlised cross product to find an axis orthogonal 
        # to both that we can rotate about
        cross = numpy.cross( refVector, blockEndGroup )
        ncross = cross / numpy.linalg.norm( cross )
        
        if numpy.array_equal( cross, [0,0,0]):
            raise RuntimeError,"alignBlocks - vectors are parallel or one is zero"
        
        # Find angle
        angle = util.vectorAngle(refVector, blockEndGroup)
        
        # Rotate
        block.rotate( ncross, angle )
        
        blockEndGroup = block.coords[ iblockEndGroup ]
        
        print "alignBlocks AFTER"
        print blockEndGroup
        print refVector
        return
    
    def X_alignAxis(self, block, endGroupIndex, targetVector, axisLabel=None ):
        """
        Align the given block so that the bond defined by the endGroupIndex
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
        return

    def bondBlock(self, bond):
        """ Bond the second block to the first and update the data structures"""
        
        iblock, iatom, ioblock, ioatom = bond
        block = self.blocks[iblock]
        lenbcoords = len(block.coords)
        oblock = self.blocks[ioblock]
        block.bond( oblock, (iatom, ioatom) )
        
        # Update the cells of the new block - WRITE TEST!
        #print self.box1
        #print "iblock, ioblock ",iblock,ioblock
        for icoord, key in enumerate(oblock.atomCell):
            #print "CHECKING ",icoord,key
            i = self.box1[key].index( ( ioblock, icoord )  )
            #print "UPDATING INDEX ",i
            # NEED TO UPDATE INDEX OF NEW COORD
            self.box1[key][i] =  ( iblock, icoord+lenbcoords )
        
        del self.blocks[ioblock]
        
        return 
    
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
    
    def checkMove(self,iblock):
        """
        See what happened with this move
        
        Return:
        """
        # Get a list of the close atoms
        close = self.closeAtoms(iblock)
        
        #print "GOT {} CLOSE ATOMS ".format(len(close))
        
        if not close:
            # Nothing to see so move along
            return True
        
        # convert bondAngle and bondMargin to angstroms
        bondAngleMargin = self.bondAngleMargin / util.RADIANS2DEGREES
        bondAngle = self.bondAngle / util.RADIANS2DEGREES
        
        # Get the block
        block = self.blocks[iblock]
        
        bonds = [] # list of possible bond atoms
        for ( iatom, ioblock, ioatom ) in close:
            
            symbol = block.symbols[iatom]
            radius = block.atom_radii[iatom]
            oblock = self.blocks[ioblock]
            coord = block.coords[iatom]
            ocoord = oblock.coords[ioatom]
            #print "CHECKING  ATOMS {}:{}->{}:{} = {}".format(iatom,iblock, ioatom,ioblock,self.distance( coord, ocoord ) )
            
            # First see if both atoms are endGroups
            if iatom in block.endGroups and ioatom in oblock.endGroups:
                # Checking for a bond
                # NB ASSUMPION FOR BOND LENGTH CHECK IS BOTH BLOCKS HAVE SAME ATOM TYPES
                osymbol = oblock.symbols[ioatom]
                
                # This uses -the 'optimised' bond lenght check - probably uneeded
                bond_length = block.bondLength( symbol, osymbol )
                
                #print "CHECKING BOND ATOMS ",bond_length,self.distance( coord, ocoord )
                
                # THINK ABOUT BETTER SHORT BOND LENGTH CHECK
                if  bond_length - self.bondMargin < self.distance( coord, ocoord ) < bond_length + self.bondMargin:
                    
                    print "Possible bond for ",iatom,ioblock,ioatom
                    # Possible bond so check the angle
                    icontact = block.endGroupContactIndex( iatom )
                    contact = block.coords[icontact]
                    print "CHECKING ANGLE BETWEEN: {0} | {1} | {2}".format( contact, coord, ocoord )
                    angle = util.angle( contact, coord, ocoord )
                    #print "{} < {} < {}".format( bondAngle-bondAngleMargin, angle, bondAngle+bondAngleMargin  )
                    
                    if ( bondAngle-bondAngleMargin < angle < bondAngle+bondAngleMargin ):
                        bonds.append( (iblock, iatom, ioblock, ioatom) )
                    else:
                        print "Cannot bond due to angle: {}".format(angle  * util.RADIANS2DEGREES)
                        return False
                        
                # Finished checking for bonds so move onto the next atoms
                continue
           
            # No bond so just check if the two atoms are close enough for a clash
            oradius = oblock.atom_radii[ioatom]
            #d = self.distance( coord, ocoord )
            #l = radius+oradius+self.atomMargin
            if self.distance( coord, ocoord ) < radius+oradius+self.atomMargin:
                #print "CLASH {}->{} = {} < {}".format( coord,ocoord, d, l  )
                return False
    
        # Here no atoms clash and we have a list of possible bonds - so bond'em!
        if len(bonds):
            if len(bonds) > 1:
                raise RuntimeError,"JENS FIX FOR MULTIPLE BONDS!!!"
            # FIX _ JUST ONE BOND FOR NOW AS WE NEED TO THINK ABOUT DATA STRUCTURES
            #for bond in bonds:
            #    self.bondBlock(bond)
            self.bondBlock(bonds[0])
            print "Added {} Bonds".format(len(bonds))
        
        # Either got bonds or no clashes
        return True

        
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
        
        #print "box1 ",self.box1
        
        block=self.blocks[iblock]
        for icoord,coord in enumerate(block.coords):
            
            # Get the box this atom is in
            key = block.atomCell[icoord]
            #print "Close checking [{}] {}: {} : {}".format(key,iblock,icoord,coord)
            
            # Get a list of the boxes surrounding this one
            surrounding = self.box3[key]
            
            #For each box loop through all its atoms chekcking for clashes
            for i,sbox in enumerate(surrounding):
                
                #print "KEY ",i,sbox
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
                    #print "AGAINST        [{}] {}: {} : {}".format(sbox,ioblock, iocoord,ocoord)
                    #x = ocoord[0] % self.A[0]
                    #y = ocoord[1] % self.B[1]
                    #z = ocoord[2] % self.C[2]
                    #print "PBC: {}                         {}".format(self.distance( ocoord,coord ),[x,y,z] )
                    
                    if ( self.distance( ocoord,coord ) < self.boxSize ):
                        #print "CLOSE {}-{}:({}) and {}-{}:({}): {}".format( iblock,icoord,coord,ioblock,iocoord,ocoord, self.distance( ocoord,coord ))
                        contacts.append( (icoord, ioblock,iocoord) )
                        
        if len(contacts):
            return contacts
        else:
            return None

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
        
        return 

#    def directedShimmy(self, nsteps=100, nmoves=50):
#        """ Shuffle the molecules about making bonds where necessary for nsteps
#        minimoves is number of sub-moves to attempt when the blocks are close
#        """
#        
#        #For writing out our progress
#        filename= "SHIMMY_0.xyz"
#        self.writeXyz( filename )
#        
#        for step in range( nsteps ):
#            
#            if len(self.blocks) == 1:
#                print "NO MORE BLOCKS TO BOND _ HOORAY!"
#                return
#            
#            #if not step % 100:
#            print "Step: {}".format(step)
#            print "BLOCKS",self.blocks.keys()
#            print "KEYS ",self.box1
#            filename = util.newFilename(filename)
#            self.writeXyz( filename )
#            
#            imove_block, istatic_block = self.randomBlockId(count=2)
#            move_block = self.blocks[imove_block]
#            static_block = self.blocks[istatic_block]
#            
#            # Copy the original coordinates so we can reject the move
#            # we copy the whole block so we don't need to recalculate
#            # anything - not sure if this quicker then saving the coords & updating tho
#            orig_block = copy.deepcopy( move_block )
#            
#            # Calculate how far to move
#            circ = move_block.radius() + static_block.radius()
#            radius = (circ/2) + self.atomMargin
#            
#            for move in range( nmoves ):
#                
#                # Remove the block from the cell so we don't check against itself
#                self.delBlock(imove_block)
#                
#                self.randomMoveAroundCenter( move_block, static_block.centroid(), radius )
#                
#                #Add the block so we can check for clashes/bonds
#                imove_block = self.addBlock(move_block)
#                
#                # Test for Clashes with other molecules
#                ok = self.checkMove( imove_block )
#                
#                # Break out if no clashes
#                if ok:
#                    print "Successful move ",move
#                    break
#                
#                # Put it back where we got it from
#                self.delBlock(imove_block)
#                imove_block = self.addBlock(orig_block)
#                
#                #End move loop
#            #End step loop
#        #End shimmy
#        return

    def distance(self, v1, v2):
        """
        my attempt to do PBC
        """
        
        dx = v2[0] % self.A[0] - v1[0] % self.A[0]
        if math.fabs(dx) > self.A[0] * 0.5:
            dx = dx - math.copysign( self.A[0], dx)
            
        dy = v2[1] % self.B[1] - v1[1] % self.B[1]
        if math.fabs(dy) > self.B[1] * 0.5:
            dy = dy - math.copysign( self.B[1], dy)
            
        dz = v2[2] % self.C[2] - v1[2] % self.C[2]
        if math.fabs(dz) > self.C[2] * 0.5:
            dz = dz - math.copysign( self.C[2], dz)
            
        return math.sqrt( dx*dx + dy*dy + dz*dz )
#        
#    def findClose(self, oblock, block):
#        """
#        Find the blocks around oblock that might be close to block if we were moving
#        it by the parameter used for the random move around block.
#        return a list of the indexes
#        """
#        
#        # Need to add the central block radius to the diameter of the sample block
#        # plus the margins
#        rb = block.radius()
#        ro = oblock.radius()
#        margin = ro + 2*rb + self.blockMargin + self.closeMargin
#        
#        centroid = oblock.centroid()
#        closeBlocks = []
#        # This should implicitly include oblock
#        for i,b in self.block.iteritems():
#            # See if we are close enough
#            dist = numpy.linalg.norm( centroid - b.centroid() )
#            if dist < margin + b.radius():
#                closeBlocks.append( i )
#        
#        return closeBlocks
        
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


    def growBlock(self, block, iblockEndGroup, iblockS, iblockSEndGroup):
        """
        Bond block to blockS, using the given endGroups
        
        Arguments:
        block: the block we are attaching
        iblockEndGroup: the index of the endGroup to use for the bond
        iblockS: the id of the block we are bonding to in the dict of cell blocks
        iblockSEndGroup: the index of the endGroup to use for the bond
        """

        print "Growblock Adding to block: {0}".format(iblockS)
        blockS = self.blocks[ iblockS ]
        
        # The vector we want to align along is the vector from the contact
        # to the endGroup
        blockSEndGroup =  blockS.coords[ iblockSEndGroup ]
        iblockSContact =  blockS.endGroupContactIndex( iblockSEndGroup )
        blockSContact  =  blockS.coords[ iblockSContact ]
        refVector      =  blockSEndGroup - blockSContact
        
        # Get the contact - need to copy this so we can remember it -
        # otherwise it changes with the moves
        iblockContact = block.endGroupContactIndex( iblockEndGroup )
        blockContact = block.coords[ iblockContact ].copy()
        
        # get the coord where the next block should bond
        # symbol of endGroup tells us the sort of bond we are making which determines
        # the bond length
        symbol = block.symbols[ iblockEndGroup ]
        bondPos = blockS.newBondPosition( iblockSEndGroup, symbol )
        #print "got bondPos: {}".format( bondPos )
        
        # Shift block so contact at center, so the vector of the endGroup can be
        # aligned with the refVector
        block.translate( -blockContact )
        
        # Align along the blockS bond
        self.alignBlocks( block, iblockEndGroup, refVector )
        
        # Now turn the second block around
        # - need to get the new vectors as these will have changed due to the moves
        blockEndGroup = block.coords[ iblockEndGroup ]
        
        # Find vector perpendicular to the bond axis
        # Dot product needs to be 0
        # xu + yv + zw = 0 - set u and v to 1, so w = (x + y)/z
        # vector is 1, 1, w
        w =  -1.0 * ( blockEndGroup[0] + blockEndGroup[1] ) / blockEndGroup[2]
        orth = numpy.array( [1.0, 1.0, w] )
        
        # Find axis that we can rotate about
        rotAxis = numpy.cross(blockEndGroup,orth)
        
        # Rotate by 180
        block.rotate( rotAxis, numpy.pi )
        
        # Now need to place the endGroup at the bond coord - at the moment
        # the contact atom is on the origin so we need to subtract the vector to the endGroup
        # from the translation vector
        
        #jmht FIX!
        block.translate( bondPos + blockEndGroup )
        
        # Now add block to the cell so we can check for clashes
        block_id = self.addBlock(block)
        
        # Check it doesn't clash - could try rotating about the dihedral
        if not self.checkMove(block_id):
#            # jmht - just for checking
#            for i in range(len(block1.coords)):
#                block1.labels[i] = 'F'
#            for i in range(len(block.coords)):
#                block.labels[i] = 'F'
#            self.writeXyz("failedBond.xyz")
#            self.writeXyz("failedBondP.xyz",periodic=True)
#            sys.exit(1)
            
            self.delBlock(block_id)
            return False
        else:
            return True
        
    
    def initCell(self,inputFile, incell=True):
        """
        Read in the first block from the input file and save it
        as the initBlock
        Also center the block in the cell
        """
        
        ib = buildingBlock.BuildingBlock()
        ib.fromCarFile( inputFile )
        
        if incell:
            # Make sure it is positioned in the cell - if not keep moving it about randomly till
            # it fits
            while True:
                print "moving initblock"
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
        
    def randomBlockId(self,count=1):
        """Return count random block ids"""
        if count == 1:
            return random.choice( list( self.blocks.keys() ) )
        else:
            ids = []
            while len(ids) < count:
                myid = random.choice( list( self.blocks.keys() ) )
                if myid not in ids:
                    ids.append(myid)
            return ids
    
    def randomGrowBlock(self, block):
        """
        Attach the given block to a randomly selected block, at randomly
        selected endGroups in each
        """
        
        # Select a random block to bond to
        iblockS = self.randomBlockId()
        blockS = self.blocks[ iblockS ]
        
        print "Growblock Adding to block: {0}".format( iblockS )
        
        # pick random endgroup and contact in target
        iblockSEndGroup = blockS.randomEndGroupIndex()
        
        # pick random endGroup on the block we are attaching
        iblockEndGroup = block.randomEndGroupIndex()
        
        # now attach it
        self.growBlock(block, iblockEndGroup, iblockS, iblockSEndGroup)
        
    
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
        
        # move to origin, rotate there and then move to new coord
        # Use the cell axis definitions
        origin = numpy.array([0,0,0], dtype=numpy.float64 )
        move_block.translateCentroid( origin )
        
        # Make a random rotation
        angle = random.uniform( 0, 2*numpy.pi)
        move_block.rotate( self.A, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        move_block.rotate( self.B, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        move_block.rotate(self.C, angle )
        
        # Now move to new coord
        move_block.translateCentroid( coord )
        
        return

    def seed( self, nblocks, inputFile ):
        """ Seed a cell with nblocks based on the block that will be created
        from the input block
        """
        
        if self.A == None or self.B == None or self.C == None:
            raise RuntimeError,"Need to specify cell before seeding"
        
        self.initCell( inputFile )
        
        # If it's just one, add the init block - after a random move
        block = self.initBlock.copy()
        self.randomMove( block )
        self.addBlock( block )
        print "Added block 1"
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
                ok = self.checkMove( iblock )
                
                # Break out of try loop if no clashes
                if ok:
                    print "Added block {0} after {1} tries.".format( seedCount+2, tries )
                    break
                
                # Unsuccessful so remove the block from cell
                self.delBlock(iblock)
                
                # increment tries counter
                tries += 1
            
            # End Clash loop
        # End of loop to seed cell
    # End seed
    
    def setInitBlock(self,block):
        """
        Set the init block
        """
        self.maxAtomR=block.maxAtomRadius()
        if self.maxAtomR <= 0:
            raise RuntimeError,"Error setting initBlock"
        
        self.boxSize = ( self.maxAtomR * 2 ) + self.atomMargin
        
        #jmht - ceil or floor
        self.numBoxA = int(math.ceil( self.A[0] / self.boxSize ) )
        self.numBoxB = int(math.ceil( self.B[1] / self.boxSize ) )
        self.numBoxC = int(math.ceil( self.C[2] / self.boxSize ) )
        
        self.initBlock = block
        return
    
    def surroundBoxes(self, key ):
        """
        return a list of the boxes surrounding the box with the given key
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
                        ai = self.numBoxA-1
                    elif ai > self.numBoxA-1:
                        ai = 0
                    bj = b+j
                    if bj < 0:
                        bj = self.numBoxB-1
                    elif bj > self.numBoxB-1:
                        bj = 0
                    ck = c+k
                    if ck < 0:
                        ck = self.numBoxC-1
                    elif ck > self.numBoxC-1:
                        ck = 0
                    skey = (ai, bj, ck)
                    #print "sKey ({},{},{})->({})".format(a,b,c,skey)
                    if skey not in l:
                        l.append(skey)
                        
        return l
    
    def shimmy(self, nsteps=100, nmoves=50, stype=None):
        """ Shuffle the molecules about making bonds where necessary for nsteps
        minimoves is number of sub-moves to attempt when the blocks are close
        """
        
        startBlocks=len(self.blocks)

        filename = "SHIMMY_0.xyz"
        self.writeXyz( filename )
        self.writeXyz( filename+".lab", label=True )
        
        for step in range( nsteps ):
            
            if len(self.blocks) == 1:
                filename = util.newFilename(filename)
                print "NO MORE BLOCKS!\nResult in file: {}".format(filename)
                self.writeXyz( filename )
                self.writeXyz( os.path.splitext(filename)[0]+"incell.xyz", periodic=True )
                self.writeXyz( filename+".lab", label=True )
                break
            
            #if True:
            if not step % 100:
                print "step {}".format(step)
                filename = util.newFilename(filename)
                self.writeXyz( filename )
                self.writeXyz( filename+".lab", label=True )
                self.writeXyz( os.path.splitext(filename)[0]+"incell.xyz", periodic=True )
                
            istatic_block = None
            if stype == "block" or stype == "bond":
                imove_block, istatic_block = self.randomBlockId(count=2)
                static_block = self.blocks[istatic_block]    
            else:
                imove_block = self.randomBlockId()
                
            move_block = self.blocks[imove_block]
            
            # Copy the original coordinates so we can reject the move
            orig_coords = copy.deepcopy( move_block.coords )
            
            center = None
            if stype == "block":
                # Calculate how far to move
                circ = move_block.radius() + static_block.radius()
                radius = (circ/2) + self.atomMargin
                center = static_block.centroid()
            elif stype == "bond":
                staticEndGroupIndex = static_block.randomEndGroupIndex()
                moveEndGroupIndex = move_block.randomEndGroupIndex()
                center = static_block.newBondPosition( staticEndGroupIndex, move_block, moveEndGroupIndex)
                # jmht - this is wrong!
                radius = self.atomMargin
            else:
                nmoves=1
            
            for move in range(nmoves):
                
                #print "move ",move
                
                # Remove the move_block from the cell so we don't check against itself
                self.delBlock(imove_block)
                
                if stype == "block" or stype == "bond":
                    self.randomMoveAroundCenter( move_block, center, radius )
                else:
                    self.randomMove( move_block )
                
                #Add the move_block back so we can check for clashes/bonds
                icheck = self.addBlock(move_block)
                if icheck != imove_block:
                    raise RuntimeError,"BAD ADD IN SHIMMY1"
                
                # Test for Clashes with other molecules
                ok = self.checkMove( imove_block )
                
                # If the move failed, put the move_block back
                if ok:
                    # End the moves and go onto the next step
                    break
                else:
                    # Put it back where we got it from
                    self.delBlock( imove_block )
                    move_block.coords = copy.deepcopy(orig_coords)
                    move_block.update()
                    icheck = self.addBlock(move_block)
                    if icheck != imove_block:
                        raise RuntimeError,"BAD ADD IN SHIMMY1"
        
        # End of shimmy loop
        
        endBlocks = len(self.blocks)
        made = startBlocks-endBlocks
        if made > 0:
            print "Shimmy bonded {} blocks".format(made)
        else:
            print "Shimmy made no bonds"
        
        return
        #End shimmy
             
    def writeXyz(self, ofile, label=False, periodic=False ):
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
                    if periodic:
                        # PBC
                        x = c[0] % self.A[0]
                        y = c[1] % self.B[1]
                        z = c[2] % self.C[2]
                        #xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.symbols[j], c[0], c[1], c[2] )
                        xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format(  util.label2symbol(block.labels[j]), x, y, z )
                    else:
                        xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol(block.labels[j]), c[0], c[1], c[2] )
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
        
    def testAlignBlocks(self):
        """Test we can align two blocks correctly"""
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        #block = self.makeCh4( cell )
        cell.seed( 2, "../PAF_bb_typed.car" )
        
        # Get the two blocks
        ids = cell.blocks.keys()
        blockS = cell.blocks[ ids[0] ]
        block = cell.blocks[ ids[1] ]
        
        # Get the atoms that define things
        iblockSEndGroup = blockS.randomEndGroupIndex()
        blockSEndGroup = blockS.coords[ iblockSEndGroup ]
        iblockSContact = blockS.endGroupContactIndex( iblockSEndGroup )
        blockSContact = blockS.coords[ iblockSContact ]
        
        iblockEndGroup=block.randomEndGroupIndex()
        iblockContact = block.endGroupContactIndex( iblockEndGroup )
        blockContact = block.coords[ iblockContact ]
        
        # we want to align along block1Contact -> block1EndGroup
        refVector = blockSEndGroup - blockSContact 
        
        # Position block so contact is at origin
        block.translate( -blockContact )
        
        cell.alignBlocks( block, iblockEndGroup, refVector )
        
        # Check the relevant atoms are in the right place
        blockEndGroup = block.coords[ iblockEndGroup ]
        blockContact = block.coords[  iblockContact ]
        
        newVector = blockEndGroup - blockContact
        
        # Normalise two vectors so we can compare them
        newNorm = newVector / numpy.linalg.norm(newVector)
        refNorm = refVector / numpy.linalg.norm(refVector)
        
        # Slack tolerances - need to work out why...
        self.assertTrue( numpy.allclose( newNorm, refNorm),
                         msg="End Group incorrectly positioned: {0} | {1}".format(newNorm, refNorm )  )

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
        
        
    def testDistance(self):
        """Test the distance under periodic boundary conditions"""
        
        CELLA = 10
        CELLB = 10
        CELLC = 10
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        
        v1 = [ 2.46803012, 1.67131881, 1.96745421]
        v2 = [ 1.07988345, 0.10567109, 1.64897769]
        
        nv1 = numpy.array( v1 )
        nv2 = numpy.array( v2 )
        
        dc1 = cell.distance(nv1,nv2)
        dn = numpy.linalg.norm(nv2-nv1)
        self.assertEqual( dc1, dn, "Distance within cell:{} | {}".format(dc1,dn) )
        
        x = v2[0] + 2 * CELLA
        y = v2[1] + 2 * CELLB
        z = v2[2] + 2 * CELLC
        nv2 = numpy.array([ x, y, z ])
        dc2 = cell.distance(nv1,nv2)
        self.assertAlmostEqual( dc1, dc2, 11,"Distance across multiple cells +ve: {} | {}".format(dc1,dc2) )
        
        x = v2[0] - 2 * CELLA
        y = v2[1] - 2 * CELLB
        z = v2[2] - 2 * CELLC
        nv2 = numpy.array([ x, y, z ])
        dc3 = cell.distance(nv1,nv2)
        self.assertAlmostEqual( dc1, dc3, 11, "Distance across multiple cells -ve: {} | {}".format(dc1,dc3) )
        
        v1 = numpy.array([ 0.0, 0.0, 0.0 ])
        v2 = numpy.array([ 0.0, 0.0, 8.0 ])
        dc = cell.distance(v1,v2)
        self.assertEqual( dc, 2.0, "Distance across boundary cell:{}".format(dc) )
        
    def testGrowBlock(self):
        """Test we can add a block correctly"""
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        cell.seed( 1, "../PAF_bb_typed.car" )
        #cell.writeXyz("JENS.xyz")

        nblocks=10
        for i in range( nblocks ):
            ok = cell.randomGrowBlock( cell.initBlock.copy() )
            #cell.writeXyz("JENS"+str(i)+".xyz")
            if not ok:
                print "Failed to add block"

        self.assertEqual(1,len(cell.blocks), "Growing blocks found {0} blocks".format( len(cell.blocks) ) )
            
        #cell.writeXyz("JENS.xyz", label=False)
        
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
        
        
        
    def testCloseAtoms(self):
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( atomMargin=0.1, boxMargin=0.1, bondMargin=0.5, bondAngle=180, bondAngleMargin=15 )
        
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        #block1 = self.makePaf( cell )
        cell.initCell("../ch4_typed.car",incell=False)
        block1 = cell.initBlock.copy()
        
        b1 = numpy.array([1,1,1], dtype=numpy.float64 )
        block1.translateCentroid( b1 )
        block1_id = cell.addBlock(block1)
        
        block2=cell.initBlock.copy()
        b2 = numpy.array([2.2,2.2,2.2], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        block2_id = cell.addBlock(block2)
        
        closeList =  cell.closeAtoms(block1_id)
        refPairs = [(0, 3), (1, 3), (2, 3)] # Also used to check PBC
        closePairs = []
        for iatom, ioblock, ioatom in closeList:
            closePairs.append( (iatom,ioatom) )
            
        #cell.writeXyz("close1.xyz", label=False)
            
        self.assertEqual(closePairs,
                        refPairs,
                         "Many contacts: {}".format(closePairs))
        
        # Too far for any contacts
        cell.delBlock(block2_id)
        block2 = cell.initBlock.copy()
        b2 = numpy.array([10,10,10], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        block2_id = cell.addBlock(block2)
        
        close = cell.closeAtoms(block1_id)
        self.assertEqual(close, None, "No contacts: ".format(close))
        
        # Now check across periodic boundary
        cell.delBlock(block2_id)
        block2 = cell.initBlock.copy()
        x = 2.2 + 2 * cell.A[0]
        y = 2.2 + 2 * cell.B[1]
        z = 2.2 + 2 * cell.C[2]
        
        b2 = numpy.array([x,y,z], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        block2_id = cell.addBlock(block2)
        
        #cell.writeXyz("close2.xyz", label=False)
        
        closeList =  cell.closeAtoms(block1_id)
        closePairs = []
        for iatom, ioblock, ioatom in closeList:
            closePairs.append( (iatom,ioatom) )

        self.assertEqual(closePairs, refPairs, "Periodic boundary: {}".format(closePairs))

#        cell.delBlock(block2_id)
#        block2 = cell.initBlock.copy()
#        b2 = numpy.array([29,1,1], dtype=numpy.float64 )
#        block2.translateCentroid( b2 )
#        block2_id = cell.addBlock(block2)
#        
#        cell.writeXyz("close.xyz", label=False)
#        
#        closeList =  cell.closeAtoms(block1_id)
#        closePairs = []
#        for iatom, ioblock, ioatom in closeList:
#            closePairs.append( (iatom,ioatom) )
#
#        self.assertEqual(closePairs, [(0, 2), (1, 2), (3, 0), (3, 2), (4, 0), (4, 2)], "Periodic boundary: ".format(closePairs))

    def testBonding(self):
        
        cell = Cell( )
        cell.fromXyz("canBond_lab.xyz")
        ib1 = cell.blocks.keys()[0]
        ok = cell.checkMove(ib1)
        #for b in cell.blocks.keys():
        #    print cell.checkMove(b)

    def testSurroundBoxes(self):
        """
        """
        cell = Cell( )
        cell.cellAxis( A=5, B=5, C=5 )
        # box size=1 - need to set manually as not reading in a block
        cell.maxAtomR = 0.5
        cell.atomMargin = 0.0
        
        cell.boxSize = ( cell.maxAtomR * 2 ) + cell.atomMargin
        cell.numBoxA = int(math.floor( cell.A[0] / cell.boxSize ) )
        cell.numBoxB = int(math.floor( cell.B[1] / cell.boxSize ) )
        cell.numBoxC = int(math.floor( cell.C[2] / cell.boxSize ) )
        
        s = [(2, 2, 2), (2, 2, 1), (2, 2, 3), (2, 1, 2), (2, 1, 1), (2, 1, 3), (2, 3, 2), (2, 3, 1), 
         (2, 3, 3), (1, 2, 2), (1, 2, 1), (1, 2, 3), (1, 1, 2), (1, 1, 1), (1, 1, 3), (1, 3, 2),
          (1, 3, 1), (1, 3, 3), (3, 2, 2), (3, 2, 1), (3, 2, 3), (3, 1, 2), (3, 1, 1), (3, 1, 3), 
          (3, 3, 2), (3, 3, 1), (3, 3, 3)]
        self.assertEqual(s, cell.surroundBoxes( (2,2,2) ), "in center")
        
        sb = cell.surroundBoxes( (0,0,0) )
        s = [ (0, 0, 0), (0, 0, 4), (0, 0, 1), (0, 4, 0), (0, 4, 4), (0, 4, 1), (0, 1, 0), (0, 1, 4), 
        (0, 1, 1), (4, 0, 0), (4, 0, 4), (4, 0, 1), (4, 4, 0), (4, 4, 4), (4, 4, 1), (4, 1, 0), (4, 1, 4),
         (4, 1, 1), (1, 0, 0), (1, 0, 4), (1, 0, 1), (1, 4, 0), (1, 4, 4), (1, 4, 1), (1, 1, 0), (1, 1, 4), (1, 1, 1)]
        self.assertEqual(s, sb , "periodic: {0}".format( sb ) )
        
    def XtestCloseDistance(self):
        """
        Test distance and close together
        """
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        cell.initCell("../ch4_typed.car",incell=False)
        block1 = cell.initBlock.copy()
        
        b1 = numpy.array([2,2,2], dtype=numpy.float64 )
        block1.translateCentroid( b1 )
        block1_id = cell.addBlock(block1)
        
        block2=cell.initBlock.copy()
        b2 = numpy.array([3,3,3], dtype=numpy.float64 )
        block2.translateCentroid( b2 )
        block2_id = cell.addBlock(block2)
        
        #cell.writeXyz("close1.xyz", label=False)
        
#        closeList =  cell.closeAtoms(block1_id)
#        for iatom, ioblock, ioatom in closeList:
#            b1c = block1.coords[iatom]
#            b2c = cell.blocks[ioblock].coords[ioatom]
#            distance = cell.distance(b1c,b2c)
#            print "{}: {} -> {}:{} = {}".format(iatom,b1c,ioatom,b2c,distance)
            
        # Distance measured with Avogadro so it MUST be right...
        refd = 0.673354948616
        distance = cell.distance(block1.coords[1], cell.blocks[block2_id].coords[3])
        self.assertAlmostEqual(refd,distance,12,"Closest atoms: {}".format(distance))

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
