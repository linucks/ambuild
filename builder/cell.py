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
        
        self._blocks = []
        
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
            
        
    def blocks(self):
        return self._blocks
    
    def addBlock( self, block ):
        self._blocks.append( block )
        
    def checkMove(self, block, iblock, closeMargin=None, bondAngle=None, blockList=None ):
        """
        See how this move went.
        Arguments
        block - block that we are checking the move for
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
            blockList = [ i for i in range( len(self.blocks()) ) ]
        
        bonds=[]
        #for i, oblock in enumerate( self.blocks() ):
        for i in blockList:
            
            print "Checking block: {}".format(i)
            # skip the block we are using
            if i == iblock:
                continue
            
            # Get the next block
            oblock = self._blocks[i]
            
            # First see if we are close enough to consider bonding
            # Don't need to do this as we have already checked that these blocks are close
            #if not block.close( oblock, margin = closeMargin ):
            #   continue
            
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
        for i,b in enumerate( self.blocks() ):
            # See if we are close enough
            dist = numpy.linalg.norm( centroid - b.centroid() )
            if dist < margin + b.radius():
                closeBlocks.append( i )
        
        return closeBlocks
        
    def fromXyz(self, xyzFile ):
        """ Read in an xyz file containing a cell and recreate the cell object"""
        
        blocks = self._readXyz(xyzFile)
        print blocks
        print blocks[0]
        for (label, coord) in blocks:
            block = buildingBlock.BuildingBlock()
            block.createFromLabelAndCoords( label, coord )
            self.addBlock(block)
        
    
    def _readXyz(self, xyzFile ):
        """"Read in an xyz file containing a cell - cell axes is title line and 
        atoms are written out with their labels.
        This routine sets the axes and returns a list of the labels and coordinates
        """
        
        # Each block is a list of [label, coords
        blocks = []
        
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
                    blocks.append( [] )
                
                label = labelf[0]+"_"+labelf[1]
                
                blocks[-1].append( [label, numpy.array(fields[1:4], dtype=numpy.float64)] )
                #labels.append(label) 
                #coords.append( numpy.array(fields[1:4], dtype=numpy.float64) )
        
        return blocks
    
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
            
            if len(self.blocks()) == 1:
                print "NO MORE BLOCKS TO BOND _ HOORAY!"
                return
            
            if not step % 100:
                print "Step: {}".format(step)
                filename = util.newFilename(filename)
                self.writeXyz( filename )
            
            # Get two blocks to sample round
            iblock, jblock = self.getRandomBlocks()
            block = self._blocks[iblock]
            oblock = self._blocks[jblock]
            print "Sampling block {} about {}".format(iblock,jblock)
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
            
            #jmht FIND OUT WHY THIS DOEN"ST WORK
            # Get a list of which blocks are close
            #closeBlocks = self.findClose( oblock, block, closeMargin=CLOSE_MARGIN, blockMargin=BLOCK_MARGIN )
            #print "Got close blocks: {}".format(closeBlocks)
            closeBlocks=None
            
            clash=0
            noclash=0
            gotBond=False
            for _ in range( nmoves ):
                
                self.randomMoveAroundBlock( oblock, block, margin=BLOCK_MARGIN )
                
                print "Move: {} about {}".format( oblock.centroid(), block.centroid() )
                
                # Loop over all blocks to see if we can bond or if we clash
                bonds = self.checkMove( block, iblock, closeMargin=CLOSE_MARGIN, bondAngle=bondAngle, blockList=closeBlocks )
                
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
                    block.bond( self._blocks[b-bcount], bond )
                    # Remove the block from the list as it is now part of the other one
                    self._blocks.pop(b-bcount)
                    bcount+=1
                
                # No need to loop anymore
                gotBond=True
                break
                            
            print "End of moves  clash/noclash: {}/{}".format(clash,noclash)
            
            if not gotBond:
                # If no bonds place the bond back at it's original position
                self._blocks[iblock] = orig_block
            
            #End move loop
        #End step loop
        
        print "END OF DIRECTED SHIMMY\nMade {} bonds and got {} clusters".format(nbonds,len(self.blocks()))
    #End directedShimmy
            

    def getRandomBlockIndex(self):
        """Return the index of one of the blocks"""
        return random.randint( 0, len(self._blocks)-1 )
    
    def getRandomBlocks(self):
        """Return two random block indices"""
        
        # Pick a block to move
        iblock = self.getRandomBlockIndex()
        
        # Find a different one to sample around
        jblock = iblock
        while jblock == iblock:
            jblock = self.getRandomBlockIndex()
        
        return (iblock, jblock)
        
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
        
    
    def seed( self, nblocks, firstBlock ):
        """ Seed a cell with nblocks based on firstBlock
        """
        
        if self.A == None or self.B == None or self.C == None:
            raise RuntimeError,"Need to specify cell before seeding"
        
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
                
            iblock = self.getRandomBlockIndex()
            block = self._blocks[iblock]
            
            # Copy the original coordinates so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker then saving the coords & updating tho
            orig_block = copy.deepcopy( block )
             
            # Make a random move
            self.randomMove( block )
            
            # Loop over all blocks to see if we can bond or if we clash
            bonds = self.checkMove( block, iblock, closeMargin=CLOSE_MARGIN, bondAngle=bondAngle )
            
            # See what happend
            if bonds:
                if bonds == "clash":
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

             
    def writeXyz(self, ofile, label=False ):
        """Write out the cell atoms to an xyz file
        If label is true we write out the atom label and block, otherwise the symbol
        """
        
        natoms=0
        xyz = ""
        for i,block in enumerate(self._blocks):
            for j, c in enumerate( block.coords ):
                if label:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( "{}_block#{}".format(block.labels[j], i), c[0], c[1], c[2] )
                else:
                    xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.labelToSymbol(block.labels[j]), c[0], c[1], c[2] )
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
        
        cell = Cell( )
        cell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
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
