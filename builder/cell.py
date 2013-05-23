'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import copy
import cPickle
import logging
import math
import os
import random
import sys
import unittest

# External modules
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
        
        # For time being origin always 0,0,0
        self.origin = numpy.array( [0,0,0], dtype=numpy.float64 )
        
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
        
        self.targetDensity = 10
        self.targetEndGroups = 100 # number of free endgroups left
        
        # Dict mapping box key to a list of tuples defining the atoms in that box
        self.box1={}
        
        # Dict mapping key to list of boxes surrounding the keyed box
        self.box3={}
        
        # the number of unit blocks - not the number of bonded blocks
        self.numBlocks = 0
        
        # max atom radius - used to calculate box size
        self.boxSize = None
        self.maxAtomR=None
        self.blockMass = None # the mass of an individual block
        
        # number of boxes in A,B,C axes - used for calculating PBC
        self.numBoxA = None
        self.numBoxB = None
        self.numBoxC = None
        
        # first block is kept separate as everything is built up from it
        self.initBlock = None
        
        # dictionary mapping id of the block to the block - can't use a list and indices
        # as we add and remove blocks and need to keep track of them
        self.blocks = {}
        
        self._setupLogging()
        
        self.fileCount=0
    
    def addBlock( self, block):
        """
        Add the block and put all atoms in their cells
        
        Args:
        block -- block to add
        
        Returns:
        index of the new block in the list
        
        jmht  - rather then the check here, we could create a copy of the _coords of the 
        block in the cell and then use these when checking distances - this way we wouldn't
        need to use PBC in the distance check, which would possibly be quicker at the expense of
        storing 2 sets of coordinates, which is negligible - test!
        """
        
        # The index of the new block
        idxBlock = id(block)
        
        # Add to the dict
        self.blocks[ idxBlock ] = block
        
        #print "nbox ",self.numBoxA,self.numBoxB,self.numBoxC
        #for idxCoord,coord in enumerate(block._coords):
        for idxCoord,coord in enumerate( block.iterCoord() ):
            
            #print "ADDING COORD ",coord
            
            # Periodic Boundaries
            x = coord[0] % self.A[0]
            y = coord[1] % self.B[1]
            z = coord[2] % self.C[2]
            
            a=int( math.floor( x / self.boxSize ) )
            b=int( math.floor( y / self.boxSize ) ) 
            c=int( math.floor( z / self.boxSize ) )
            
            key = (a,b,c)
            block.atomCell[ idxCoord ] = key
            try:
                self.box1[ key ].append( ( idxBlock, idxCoord ) )
            except KeyError:
                # Add to main list
                self.box1[ key ] = [ ( idxBlock, idxCoord ) ]
                # Map surrounding boxes
                self.box3[ key ] = self.surroundBoxes( key )
                
        return idxBlock
    

    def XbondBlock(self, bond):
        """ Bond the second block to the first and update the data structures"""
        
        iblock, iatom, ioblock, ioatom = bond
        block = self.blocks[iblock]
        lenbcoords = len(block._coords)
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
    
    def bondBlock(self, idxBlock, idxAtom, idxOblock, idxOatom ):
        """ Bond the second block to the first and update the data structures
        
        Previously we manipulted the cell data structures directly - now we take the simpler
        but cruder approach of deleting both blocks and then re-adding the joined one
        """
        
        self.logger.debug("bondBlock: {0}".format( (idxBlock, idxAtom, idxOblock, idxOatom) ) )
        
        block = self.blocks[ idxBlock ]
        oblock = self.blocks[ idxOblock ]
        
        self.delBlock( idxBlock )
        self.delBlock( idxOblock )
        
        block.addBlock( idxAtom, oblock, idxOatom )
        
        self.addBlock( block )
        
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

    def checkFinished(self):
        """
        See if we've finished according to our criteria
        """
        
        print "Got density: ",self.density()
        if self.density() < self.targetDensity:
            return True
        
        print "Got _endGroups: ",self._endGroups()
        if self._endGroups() < self.targetEndGroups:
            return True
        
        return False

    def checkMove(self,idxBlock):
        """
        See what happened with this move
        
        Return:
        True if the move succeeded
        """
        
        # Get a list of the close atoms
        close = self.closeAtoms(idxBlock)
        
        #print "GOT {} CLOSE ATOMS ".format(len(close))
        
        if not close:
            # Nothing to see so move along
            return True
        
        # convert bondAngle and bondMargin to angstroms
        bondAngleMargin = self.bondAngleMargin / util.RADIANS2DEGREES
        bondAngle = self.bondAngle / util.RADIANS2DEGREES
        
        # Get the block
        block = self.blocks[idxBlock]
        
        bonds = [] # list of possible bond atoms
        for ( idxAtom, idxOblock, idxOatom ) in close:
            
            symbol = block.atomSymbol( idxAtom )
            radius = block.atomRadius( idxAtom )
            coord = block.atomCoord( idxAtom )
            oblock = self.blocks[ idxOblock ]
            ocoord = oblock.atomCoord( idxOatom )
            
            #print "CHECKING  ATOMS {}:{}->{}:{} = {}".format(idxAtom,idxBlock, idxOatom,idxOblock,self.distance( coord, ocoord ) )
            
            # First see if both atoms are _endGroups
            #if idxAtom in block._endGroups and idxOatom in oblock._endGroups:
            if block.isEndGroup( idxAtom ) and oblock.isEndGroup( idxOatom ):
                # Checking for a bond
                # NB ASSUMPION FOR BOND LENGTH CHECK IS BOTH BLOCKS HAVE SAME ATOM TYPES
                osymbol = oblock.atomSymbol( idxOatom )
                
                bond_length = util.bondLength( symbol, osymbol )
                
                #print "CHECKING BOND ATOMS ",bond_length,self.distance( coord, ocoord )
                
                # THINK ABOUT BETTER SHORT BOND LENGTH CHECK
                if  bond_length - self.bondMargin < self.distance( coord, ocoord ) < bond_length + self.bondMargin:
                    
                    #print "Possible bond for ",idxAtom,idxOblock,idxOatom
                    # Possible bond so check the angle
                    angleAtom = block.angleAtom( idxAtom )
                    
                    #print "CHECKING ANGLE BETWEEN: {0} | {1} | {2}".format( contact, coord, ocoord )
                    angle = util.angle( angleAtom, coord, ocoord )
                    #print "{} < {} < {}".format( bondAngle-bondAngleMargin, angle, bondAngle+bondAngleMargin  )
                    
                    if ( bondAngle-bondAngleMargin < angle < bondAngle+bondAngleMargin ):
                        bonds.append( (idxBlock, idxAtom, idxOblock, idxOatom) )
                    else:
                        self.logger.debug( "Cannot bond due to angle: {}".format(angle  * util.RADIANS2DEGREES) )
                        return False
                        
                # Finished checking for bonds so move onto the next atoms
                continue
           
            # No bond so just check if the two atoms are close enough for a clash
            #oradius = oblock._atomRadii[idxOatom]
            oradius = oblock.atomRadius( idxOatom )
            
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
            idxBlock, idxAtom, idxOblock, idxOatom = bonds[0]
            self.bondBlock( idxBlock, idxAtom, idxOblock, idxOatom )
            block = self.blocks[ idxBlock ]
            #self.logger.debug(  "Block after bond: {0}".format( block ) )
            self.logger.info("Added bond: {}".format( bonds[0] ) )
        
        # Either got bonds or no clashes
        return True
    
    def _endGroups(self):
        """
        The number of free _endGroups in the cell
        """
        
        numEndGroups = 0
        for idxBlock in self.blocks.keys():
            block = self.blocks[ idxBlock ]
            numEndGroups += len( block._endGroups )
        
        return numEndGroups
        
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
#            for icoord,coord in enumerate(block._coords):
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
        Find all atoms that are close to the atoms in the given block.
        
        Args:
        iblock: index of the block in self.blocks
        
        Returns:
        a list of tuples: (thisAtomIndex, otherBlockIndex, otherAtomIndex) or None if there we no close contacts
        """
        
        contacts=[]
        
        #print "box1 ",self.box1
        
        block=self.blocks[iblock]
        #for idxCoord,coord in enumerate( block._coords ):
        for idxCoord,coord in enumerate( block.iterCoord() ):
            
            # Get the box this atom is in
            key = block.atomCell[ idxCoord ]
            #print "Close checking [{}] {}: {} : {}".format(key,iblock,idxCoord,coord)
            
            # Get a list of the boxes surrounding this one
            surrounding = self.box3[key]
            
            #For each box loop through all its atoms chekcking for clashes
            for i, sbox in enumerate(surrounding):
                
                #print "KEY ",i,sbox
                # For each box, get the list of the atoms as (block,coord) tuples
                # Check if we have a box with anything in it
                if not self.box1.has_key(sbox):
                    continue
                    
                for (idxOblock, idxOcoord) in self.box1[ sbox ]:
                    
                    # Check we are not checking ourself - need to check block index too!
                    if iblock == idxOblock:
                        continue
                    
                    oblock = self.blocks[idxOblock]
                    #ocoord = oblock._coords[idxOcoord]
                    ocoord = oblock.atomCoord( idxOcoord )
                    #print "AGAINST        [{}] {}: {} : {}".format(sbox,idxOblock, idxOcoord,ocoord)
                    #x = ocoord[0] % self.A[0]
                    #y = ocoord[1] % self.B[1]
                    #z = ocoord[2] % self.C[2]
                    #print "PBC: {}                         {}".format(self.distance( ocoord,coord ),[x,y,z] )
                    
                    if ( self.distance( ocoord,coord ) < self.boxSize ):
                        #print "CLOSE {}-{}:({}) and {}-{}:({}): {}".format( iblock,idxCoord,coord,idxOblock,idxOcoord,ocoord, self.distance( ocoord,coord ))
                        contacts.append( ( idxCoord, idxOblock, idxOcoord ) )
                        
        if len(contacts):
            return contacts
        else:
            return None

    def delBlock(self,blockId):
        """
        Remove the block with the given index from the cell
        """
        
        #print "Deleting block: {} from {}".format(blockId,self.blocks.keys())
        
        block =  self.blocks[ blockId ]
        
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
        
    def density(self):
        """
        The density of the cell
        """
        
        # each block has the same mass as we count unit blocks
        
        return ( self.numBlocks * self.blockMass ) / ( self.A[0] * self.B[1] * self.C[2] ) 

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
#            # anything - not sure if this quicker then saving the _coords & updating tho
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

    def dump(self,prefix="cell"):
        """Write out our current state"""
        
        self.fileCount+=1
        prefix=prefix+"_{0}".format(self.fileCount)
        self.writeXyz(prefix+".xyz",periodic=False)
        self.writeXyz(prefix+"_P.xyz",periodic=True)
        self.writePickle(prefix+".pkl")
        return
        
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
            block = buildingBlock.Block()
            block.fromLabelAndCoords( labels, coords )
            if not initBlock:
                self.setInitBlock( block.copy() )
                initBlock=True
            self.addBlock(block)
        
    def _readXyz(self, xyzFile ):
        """"Read in an xyz file containing a cell - cell axes is title line and 
        atoms are written out with their labels.
        This routine sets the axes and returns a list of the labels and coordinates
        """
        
        # Each block is a list of [label, _coords
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

    def _growBlock(self, growBlock, idxGrowAtom, idxStaticBlock, idxStaticAtom):
        """
        Position growBlock so it can bond to blockS, using the given _endGroups
        
        Arguments:

        
        We take responsibility for adding and removing the growBlock from the cell on 
        success or failure
        """
        
        #self.positionGrowBlock(growBlock, idxGrowAtom, idxStaticBlock, idxStaticAtom)
        staticBlock = self.blocks[ idxStaticBlock ]
        staticBlock.positionGrowBlock( idxStaticAtom, growBlock, idxGrowAtom )
        
        # Now add growBlock to the cell so we can check for clashes
        blockId = self.addBlock( growBlock )
        
        # Check it doesn't clash
        if self.checkMove( blockId ):
            return True
        
        # Didn't work so try rotating the growBlock about the bond to see if that lets it fit
        
        # Determine the axis and center to rotate about
        blockEndGroup = growBlock.atomCoord( idxGrowAtom )
        blockS = self.blocks[ idxStaticBlock ]
        blockSEndGroup = blockS.atomCoord( idxStaticAtom )
        axis = blockSEndGroup - blockEndGroup
        center = blockSEndGroup
        
        step = math.pi/18 # 10 degree increments
        for angle in util.frange( step, math.pi*2, step):
            
            #print "_growBlock rotating as clash: {}".format(angle*util.RADIANS2DEGREES)
            
            # remove the growBlock from the cell
            self.delBlock(blockId)
            
            # rotate by increment
            growBlock.rotate( axis, angle, center=center)
            
            # add it and check
            self.addBlock(growBlock)
            
            if self.checkMove( blockId ):
                self.logger.debug("_growBlock rotate worked")
                return True
        
        # remove the growBlock from the cell
        self.delBlock(blockId)
        
        return False
    
    def growNewBlocks(self, toGrow, maxTries=50 ):
        """
        Add toGrow new blocks to the cell based on the initBlock
        
        Args:
        toGrow: number of blocks to add
        maxTries: number of tries to add before we give up
        """
        
        self.logger.info( "Growing {0} new blocks".format( toGrow ) )
        
        added=0
        tries=0
        while added < toGrow:
            
            if tries > maxTries:
                self.logger.critical("growNewBlocks - exceeded maxtries when joining blocks!" )
                return False
            
            newblock = self.initBlock.copy()
            
            # Apply random rotation in 3 axes to randomise the orientation before we align
            newblock.randomRotate( origin=self.origin )
            
            ok = self.randomGrowBlock( newblock )
            
            if ok:
                self.numBlocks += 1
                added+=1
                tries=0
            else:
                tries+=1
                
        self.logger.info("After growNewBlocks numBlocks: {0} ({1})".format( len(self.blocks), self.numBlocks ) )
        
        return True

    def initCell(self,inputFile, incell=True):
        """
        Read in the first block from the input file and save it
        as the initBlock
        Also center the block in the cell
        """
        
        ib = buildingBlock.Block()
        ib.fromCarFile( inputFile )
        
        if incell:
            # Make sure it is positioned in the cell - if not keep moving it about randomly till
            # it fits
            while True:
                self.logger.debug( "moving initblock" )
                #print "{} | {} | {} | {}".format( ib.centroid()[0], ib.centroid()[1], ib.centroid()[2], ib.radius()  )
                if ib.centroid()[0] - ib.radius() < 0 or \
                   ib.centroid()[0] + ib.radius() > self.A[0] or \
                   ib.centroid()[1] - ib.radius() < 0 or \
                   ib.centroid()[1] + ib.radius() > self.B[1] or \
                   ib.centroid()[2] - ib.radius() < 0 or \
                   ib.centroid()[2] + ib.radius() > self.C[2]:
                    self.randomMoveBlock( ib, margin=ib.radius() )
                else:
                    # Its in!
                    break

        self.setInitBlock(ib)
        self.logger.debug( "initCell - block radius: {0}".format( ib.radius() ) )
        return
    
    def _joinBlocks(self):
        """
        See if we can join two randomly chosen blocks at random _endGroups
        """
        idxMoveBlock = self.randomBlockId()
        moveBlock = self.blocks[ idxMoveBlock ]
        
        # Copy the original block so we can replace it if the join fails
        blockCopy = moveBlock.copy()
        
        # Remove from cell so we don't check against itself and pick a different out
        self.delBlock( idxMoveBlock )
        
        idxStaticBlock = self.randomBlockId()
        staticBlock =  self.blocks[ idxStaticBlock ]
        
        # pick random endgroup and contact in target
        idxStaticBlockEG = staticBlock.randomEndGroup()
        
        # pick random endGroup on the block we are attaching
        idxMoveBlockEG = moveBlock.randomEndGroup()
        
        # now attach it
        ok = self._growBlock( moveBlock, idxMoveBlockEG, idxStaticBlock, idxStaticBlockEG )
        
        #mycell.randomRotateBlock( moveBlock )
        if not ok:
            self.addBlock( blockCopy )
        
        return ok
    
    def joinBlocks(self, toJoin, maxTries=50 ):
        """
        Try joining number of blocks together
        
        Args:
        toJoin: number of blocks to join
        maxTries: the maximum number of moves to try when joining
        """
        
        self.logger.info( "Joining {0} new blocks".format( toJoin ) )
        
        added=0
        tries=0
        while added < toJoin:
            
            if len ( self.blocks ) == 1:
                self.logger.info("joinBlocks - no more blocks to join!")
                return False
                
            if tries > maxTries:
                self.logger.critical("joinBlocks - exceeded maxtries when joining blocks!")
                return False
            
            # now attach it
            ok = self._joinBlocks()
            
            if ok:
                added+=1
                tries=0
            else:
                tries+=1
        
        self.logger.info("After joinBlocks numBlocks: {0} ({1})".format( len(self.blocks), self.numBlocks ) )
        
        return True
        
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
        selected _endGroups in each
        """
        
        # Select a random block to bond to
        idxStaticBlock = self.randomBlockId()
        staticBlock = self.blocks[ idxStaticBlock ]
        
        # pick random endgroup and contact in target
        idxStaticBlockEG = staticBlock.randomEndGroup()
        
        #self.logger.debug( "initBlock: {0}".format(self.initBlock))
        #self.logger.debug( "static block: {0}".format(staticBlock))
        #self.logger.debug( "block: {0}".format(block))
        
        # pick random endGroup on the block we are attaching
        idxBlockEG = block.randomEndGroup()
        
        self.logger.debug( "randomGrowblock calling _growBlock: {0}".format( (block, idxBlockEG, idxStaticBlock, idxStaticBlockEG) ) )
        # now attach it
        return self._growBlock( block, idxBlockEG, idxStaticBlock, idxStaticBlockEG )

    def randomMoveBlock(self, block, margin=None ):
        """Randomly move the given block
         If buffer is given, use this as a buffer from the edges of the cell
         when selecting the coord
        """
        
        # Get _coords of random point in the cell
        if margin:
            x = random.uniform(margin,self.A[0]-margin)
            y = random.uniform(margin,self.B[1]-margin)
            z = random.uniform(margin,self.C[2]-margin)
        else:
            x = random.uniform(0,self.A[0])
            y = random.uniform(0,self.B[1])
            z = random.uniform(0,self.C[2])
            
        coord = numpy.array([x,y,z], dtype=numpy.float64 )
        
        #print "Got random coord: {}".format(coord)
        
        # Move to origin, rotate there and then move to new coord
        # Use the cell axis definitions
        block.translateCentroid( self.origin )
        
        block.randomRotate( origin=self.origin, atOrigin=True )
        
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
        
        move_block.translateCentroid( coord )
        
        move_block.randomRotateBlock( origin=self.origin )
        return

    def seed( self, nblocks, inputFile, maxTries=1000 ):
        """ Seed a cell with nblocks based on the block that will be created
        from the input block
        """
        
        if self.A == None or self.B == None or self.C == None:
            raise RuntimeError,"Need to specify cell before seeding"
        
        self.initCell( inputFile )
        
        # If it's just one, add the init block - after a random move
        block = self.initBlock.copy()
        self.randomMoveBlock( block )
        self.addBlock( block )
        self.numBlocks += 1
        self.logger.debug("seed added first block: {0}".format( id(block) ) )
        if nblocks == 1:
            return
        
        # Loop through the nblocks adding the blocks to
        # the cell - nblocks-1 as we've already added the first
        for seedCount in range( nblocks-1 ):
            # Create new block
            newblock = self.initBlock.copy()
            tries = 0
            ok = False
            while not ok:
                # quit on maxTries
                if tries >= maxTries:
                    self.logger.critical("Exceeded maxtries when seeding")
                    sys.exit(1)
                
                # Move the block and rotate it
                self.randomMoveBlock(newblock)
                #print "RANDOM TO MOVE TO: {}".format( newblock.centroid() )
                
                #Add the block so we can check for clashes/bonds
                idxBlock = self.addBlock(newblock)
                
                # Test for Clashes with other molecules
                ok = self.checkMove( idxBlock )
                
                # Break out of try loop if no clashes
                if ok:
                    self.logger.debug("seed added block {0} after {1} tries.".format( seedCount+2, tries ) )
                    self.numBlocks += 1
                    break
                
                # Unsuccessful so remove the block from cell
                self.delBlock(idxBlock)
                
                # increment tries counter
                tries += 1
            
            # End Clash loop
        # End of loop to seed cell
        
        self.logger.info("After seed numBlocks: {0} ({1})".format( len(self.blocks), self.numBlocks ) )
        return
    # End seed
    
    def setInitBlock(self, block ):
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
        
        self.blockMass = block.mass()
        
        self.initBlock = block
        return
    
    def _setupLogging( self ):
        """
        Set up the various log files/console logging and return the logger
        """
        
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        
        # create file handler and set level to debug
        fl = logging.FileHandler("ambuild.log")
        fl.setLevel(logging.DEBUG)
        
        # create formatter for fl
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # add formatter to fl
        fl.setFormatter(formatter)
        
        # add fl to logger
        logger.addHandler(fl)
        
        # Now create console logger for outputting stuff
        # create file handler and set level to debug
        # Seems they changed the api in python 2.6->2.7
        try:
            cl = logging.StreamHandler(stream=sys.stdout)
        except TypeError:
            cl = logging.StreamHandler(strm=sys.stdout)

        cl.setLevel(logging.INFO)

        # create formatter for fl
        # Always add a blank line after every print
        formatter = logging.Formatter('%(message)s\n')

        # add formatter to fl
        cl.setFormatter(formatter)

        # add fl to logger
        logger.addHandler(cl)

        self.logger = logger

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
                self.writeXyz( os.path.splitext(filename)[0]+"_incell.xyz", periodic=True )
                self.writeXyz( filename+".lab", label=True )
                break
            
            #if True:
            if not step % 100:
                print "step {}".format(step)
                filename = util.newFilename(filename)
                self.writeXyz( filename )
                self.writeXyz( filename+".lab", label=True )
                self.writeXyz( os.path.splitext(filename)[0]+"_incell.xyz", periodic=True )
                
            istatic_block = None
            if stype == "block" or stype == "bond":
                imove_block, istatic_block = self.randomBlockId(count=2)
                static_block = self.blocks[istatic_block]    
            else:
                imove_block = self.randomBlockId()
                
            move_block = self.blocks[imove_block]
            
            # Copy the original coordinates so we can reject the move
            orig_coords = copy.deepcopy( move_block._coords )
            
            center = None
            if stype == "block":
                # Calculate how far to move
                circ = move_block.radius() + static_block.radius()
                radius = (circ/2) + self.atomMargin
                center = static_block.centroid()
            elif stype == "bond":
                staticEndGroupIndex = static_block.randomEndGroup()
                moveEndGroupIndex = move_block.randomEndGroup()
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
                    self.randomMoveBlock( move_block )
                
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
                    move_block._coords = copy.deepcopy(orig_coords)
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
      
    def writePickle( self, fileName=None ):
        """Pickle ourselves"""
        
        if not fileName:
            fileName="cell.pkl"
        
        # Need to close all open filehandles
        logging.shutdown()
        
        pfile = open( fileName, "w" )
        cPickle.dump(self,pfile)
        pfile.close()
        
        return
             
    def writeXyz(self, ofile, label=False, periodic=False ):
        """Write out the cell atoms to an xyz file
        If label is true we write out the atom label and block, otherwise the symbol
        """
        
        natoms=0
        xyz = ""
        for i,block in self.blocks.iteritems():
            for j, c in enumerate( block.iterCoord() ):
                if label:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( "{}_block#{}".format( block.atomLabel(j),
                                                                                                       c[0], c[1], c[2] ) )
                else:
                    if periodic:
                        # PBC
                        x = c[0] % self.A[0]
                        y = c[1] % self.B[1]
                        z = c[2] % self.C[2]
                        #xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.symbols[j], c[0], c[1], c[2] )
                        xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol( block.atomLabel(j) ), x, y, z )
                    else:
                        xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol( block.atomLabel(j) ),
                                                                                     c[0], c[1], c[2] )
                natoms += 1
        
        # Write out natoms and axes as title
        xyz = "{}\nAxes:{}:{}:{}\n".format(natoms, self.A[0], self.B[1], self.C[2] ) + xyz
        
        with open( ofile, 'w' ) as f:
            fpath = os.path.abspath(f.name)
            f.writelines( xyz )
            
        self.logger.info( "Wrote cell file: {0}".format(fpath) )
    
    def __str__(self):
        """
        """
        s = ""
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
        return d
    
    def __setstate__(self, d):
        """Called when we are unpickled """
        self.__dict__.update(d)
        self._setupLogging()
    
    
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
        
        paf = self.buildingBlock.Block()
        paf.fromCarFile("../PAF_bb_typed.car")
        return paf
        
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
        test_coord = cell.blocks[3]._coords[4]
        
        self.assertEqual( nblocks, len(cell.blocks), "Incorrect number of cell blocks at start: {}".format( len(cell.blocks) ))
        
        outfile = "./testCell.xyz"
        cell.writeXyz( outfile, label=True )
        
        newCell = Cell()
        
        newCell.fromXyz( outfile )
        self.assertEqual( nblocks, len(newCell.blocks), "Incorrect number of cell blocks after read: {}".format(len(newCell.blocks) ))
        
        self.assertTrue( numpy.allclose( test_coord, cell.blocks[3]._coords[4], rtol=1e-9, atol=1e-9 ),
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
        
    def XtestRotateBlock(self):
        """Test we can add a block correctly"""
        
        pass
        
        CELLA = 30
        CELLB = 30
        CELLC = 30
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        cell.seed( 1, "../PAF_bb_typed.car" )
        cell.writeXyz("1.xyz")
        
        block = cell.initBlock.copy()
        iblockEndGroup = block.randomEndGroup()
        blockEndGroup = block._coords[ iblockEndGroup ]
        
        iblockS = cell.blocks.keys()[0]
        blockS = cell.blocks[ iblockS ]
        iblockSEndGroup = blockS.randomEndGroup()
        blockSEndGroup = blockS._coords[ iblockSEndGroup ]
        
        cell.positionGrowBlock(block, iblockEndGroup, iblockS, iblockSEndGroup)
        
        bid = cell.addBlock(block)
        cell.writeXyz("2.xyz")
        cell.delBlock(bid)
        
        axis = blockSEndGroup - blockEndGroup
        center = blockSEndGroup
        block.rotate( axis, math.pi*2, center=center)
        
        bid = cell.addBlock(block)
        cell.writeXyz("3.xyz")
        
        
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
            for c in b._coords:
                if ( 0-radius > c[0] > CELLA[0]+ radius ) or \
                   ( 0-radius > c[1] > CELLB[1]+ radius ) or \
                   ( 0-radius > c[2] > CELLC[2]+ radius ):
                    
                    bad.append( b )
        
        self.assertEqual( 0, len(bad), "Got {} blocks outside cell: {}".format( len(bad), bad ) )
        
        #cell.writeXyz("seedTest.xyz")
        
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
#            b1c = block1._coords[iatom]
#            b2c = cell.blocks[ioblock]._coords[ioatom]
#            distance = cell.distance(b1c,b2c)
#            print "{}: {} -> {}:{} = {}".format(iatom,b1c,ioatom,b2c,distance)
            
        # Distance measured with Avogadro so it MUST be right...
        refd = 0.673354948616
        distance = cell.distance(block1._coords[1], cell.blocks[block2_id]._coords[3])
        self.assertAlmostEqual(refd,distance,12,"Closest atoms: {}".format(distance))

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
