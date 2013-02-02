'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import sys
import os
import copy
import random
import unittest

import numpy


class Cell():
    '''
    classdocs
    '''


    def __init__( self,  A, B, C ):
        '''
        Constructor
        '''
        
        # A, B and C are cell vectors - for the time being we assume they are orthogonal
        assert A[1] == 0 and A[2] == 0, "Cell vectors need to be orthogonal and in correct order"
        assert B[0] == 0 and B[2] == 0, "Cell vectors need to be orthogonal and in correct order"
        assert C[0] == 0 and C[1] == 0, "Cell vectors need to be orthogonal and in correct order"
        self.A = A
        self.B = B
        self.C = C
        
        self._blocks = []
        
    def blocks(self):
        return self._blocks
    
    def addBlock( self, block ):
        self._blocks.append( block )
    
    def getRandomBlockIndex(self):
        """Return the index of one of the blocks"""
        return random.randint( 0, len(self._blocks)-1 )
        
    def randomMove(self, block ):
        """Randomly move the given block
         Defintely needs more work on the rotation
        """
        
        # Get coords of random point in the cell
        x = random.uniform(0,self.A[0])
        y = random.uniform(0,self.B[1])
        z = random.uniform(0,self.C[2])
        position = numpy.array([x,y,z], dtype=numpy.float64 )
        
        #print "Got random position: {}".format(position)
        
        # Overklll - move to origin, rotate there and then move to new position
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
        
        
    def randomSmallMove(self, block, cog, prange=1.0 ):
        
        """Get a random position that moves the COG prange
        units from its current position, The cog argument is the
        original centre of geometry so that we sample around it"""
        
        
        # Calculate new position
        x = random.uniform(-prange/2,prange/2)
        y = random.uniform(-prange/2,prange/2)
        z = random.uniform(-prange/2,prange/2)
        xyz = numpy.array( [x,y,z], dtype=numpy.float64 )
        position = numpy.add( cog, xyz )
        
        # Overklll - move to origin, rotate there and then move to new position
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
        
    
    def seed( self, nblocks, firstBlock ):
        """ Seed a cell with nblocks based on firstBlock
        """
        
        MAXTRIES = 1000
        # Loop through the nblocks adding the blocks to
        # the cell
        for seedCount in range( nblocks ):
            # Create new block
            newblock = firstBlock.copy()
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
                
                # Test for Clashes with other molecules - always set clash
                # to False at start so it only becomes true if no clashes
                clash = False
                for block in self._blocks:
                    if block.clash(newblock):
                        clash = True
                        break
                
                # Break out if no clashes
                if clash == False:
                    break
                
                # increment tries counter
                tries += 1
            
            # End Seeding loop
            if (clash == False):
                #print "ADDED BLOCK AT POS: {}".format( newblock.centroid() )
                self.addBlock(newblock)
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
                
            iblock = self.getRandomBlockIndex()
            
            block = self._blocks[iblock]
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
             
            # Make a random move
            self.randomMove( block )
            
            # Remember the original cog of the block 
            cog = block.centroid()
            
            # Loop over all blocks to see if we can bond or if we clash
            
            #  list of a list containing the index of the block this block can bond to 
            # with a tuple of the index of the two atoms in this block and other 
            bonds = []
            clash = False
            for i in range( len( self.blocks() ) ):
                
                # skip the block we are using
                if i == iblock:
                    continue
                
                # Get the next block
                oblock = self._blocks[i]
                
                # First see if we are close enough to consider bonding
                if not block.close( oblock, margin = CLOSE_MARGIN ):
                    continue
                
                # See if we can bond
                bond = block.canBond( oblock, bondAngle=bondAngle )
                
                # No bond so move  to next block
                if not bond:
                    continue
                
                # Either a bond or a clash
                if bond == "clash":
                    # A clash so break out of the whole loop so we can reject this step
                    clash=True
                    break
                else:
                    # Got a bond, so add it to the list of possible bonds
                    print "Found bond between block {} and {}".format( iblock, i)
                    bonds.append( [i, bond] )

                #End Loop Over Blocks
            
            # See what happend
            if clash:
                # Reject this move and overwrite the changed block with the original one
                self._blocks[iblock] = orig_block
            elif len( bonds ):
                # Got some bonds so deal with them
                # Need to decrement index if we are removing blocks
                bcount=0
                for b,bond in bonds:
                    block.bond( self._blocks[b-bcount], bond )
                    # Remove the block from the list as it is now part of the other one
                    self._blocks.pop(b-bcount)
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
                
            iblock = self.getRandomBlockIndex()
            
            block = self._blocks[iblock]
            
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
            for i in range( len( self.blocks() ) ):
                
                # skip the block we are using
                if i == iblock:
                    continue
                
                # Get the next block
                oblock = self._blocks[i]
                
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
                        self.randomSmallMove( block, cog, prange=PRANGE )
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
                        self._blocks[i] = orig_block
                    else:
                        block.bond( oblock, bond )
                        print "Bonded block {} with block {}".format( iblock, i)
                        # Now delete block
                        removed.append(i)
                i+=1
                #End Loop Over Blocks
            
            #
            for r in removed:
                self._blocks.pop(r)
            #End step loop

             
    def write(self, ofile ):
        """Write out the cell atoms to an xyz file"""
        
        natoms=0
        xyz = ""
        for block in self._blocks:
            for i, c in enumerate( block.coords ):
                xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.labelToSymbol(block.labels[i]), c[0], c[1], c[2] )
                natoms += 1
        
        xyz = "{}\nCell Atoms\n".format(natoms) + xyz
        
        with open( ofile, 'w' ) as f:
            fpath = os.path.abspath(f.name)
            f.writelines( xyz )
            
        print "Wrote cell file: {0}".format(fpath)
    
    def __str__(self):
        """
        """
        s = ""
        for block in self._blocks:
            s+= str(block)
            
        return s
    
    
class TestCell(unittest.TestCase):
    
    # Import only here as only needed for testing
    import buildingBlock

    def makeCh4(self):
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
    
    def makePaf(self):
        """Return the PAF molecule for testing"""
        
        paf = self.buildingBlock.BuildingBlock()
        paf.fromCarFile("../PAF_bb_typed.car")
        return paf
        
    def testSeed(self):
        """Test we can seed correctly"""
        
        
        nblocks = 10
        CELLA = [ 10,  0,  0 ]
        CELLB = [ 0, 10,  0 ]
        CELLC = [ 0,  0, 10 ]
        
        ch4 = self.makeCh4()
        
        cell = Cell( CELLA, CELLB, CELLC )
        cell.seed ( nblocks, ch4 )
        
        self.assertEqual( nblocks, len(cell.blocks()), "Incorrect number of cell blocks" )
        
        # Check no block is outside the unit cell
        # We seed by the centroid so the edges could stick out by radius
        radius = ch4.radius()
        
        bad = []
        for b in cell.blocks():
            for c in b.coords:
                if ( 0-radius > c[0] > CELLA[0]+ radius ) or \
                   ( 0-radius > c[1] > CELLB[1]+ radius ) or \
                   ( 0-radius > c[2] > CELLC[2]+ radius ):
                    
                    bad.append( b )
        
        self.assertEqual( 0, len(bad), "Got {} blocks outside cell: {}".format( len(bad), bad ) )