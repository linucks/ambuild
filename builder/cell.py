'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import sys
import os
import copy
import random

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
        
        self.blocks = []
        
    def blocks(self):
        return self.blocks
    
    def addBlock( self, block ):
        self.blocks.append( block )
    
    def getRandomRotation(self):
        """Return a random rotation in radians"""
        
        x = random.uniform(0,self.A[0])
        y = random.uniform(0,self.B[1])
        z = random.uniform(0,self.C[2])
        axis = numpy.array( [ x, y, z ], dtype=numpy.float64 )
        angle = random.uniform( 0, 2*numpy.pi)
        return axis, angle
    
    def getRandomBlockIndex(self):
        """Return the index of one of the blocks"""
        return random.randint( 0, len(self.blocks)-1 )
        
        
    def randomMove(self, block ):
        """Randomly move the given block
        If given a prange parameter, only move by at most prange angstroms
        """
        
        # Get coords of random point in the cell
        x = random.uniform(0,self.A[0])
        y = random.uniform(0,self.B[1])
        z = random.uniform(0,self.C[2])
        position = numpy.array([x,y,z])
        
        # Translate COM of molecule to coords
        block.translateCenterOfGeometry( position )
        
        # Rotate by a random number
        axis, angle = self.getRandomRotation()
        
        block.rotate( axis, angle )
        
    def randomSmallMove(self, block, cog, prange=1.0 ):
        
        """Get a random position that moves the COG prange
        units from its current position, The cog argument is the
        original centre of geometry so that we sample around it"""
        
        x = random.uniform(-prange/2,prange/2)
        y = random.uniform(-prange/2,prange/2)
        z = random.uniform(-prange/2,prange/2)
        xyz = numpy.array( [x,y,z] )
        
        position = numpy.add( cog, xyz )
        
        block.translateCenterOfGeometry( position )
        
        # Rotate by a random number
        axis, angle = self.getRandomRotation()
        
        block.rotate( axis, angle )
    
    def seed( self, nblocks, firstBlock ):
        """ Seed a cell with nblocks based on firstBlock
        """
        
        MAXTRIES = 10000
        # Loop through the nblocks adding the blocks to
        # the cell
        for seedCount in range( nblocks ):
            # Create new block
            newblock = firstBlock.copy()
            tries = 0
            clash = True
            print "Adding block: "+str(seedCount)
            while clash:
                # quit on maxTries
                if tries >= MAXTRIES:
                    print "EXCEEDED MAXTRIES WHEN SEEDING"
                    sys.exit(1)
                
                # Move the block and rotate it
                self.randomMove(newblock)
                
                # Test for Clashes with other molecules - always set clash
                # to False at start so it only becomes true if no clashes
                clash = False
                for block in self.blocks:
                    if block.clash(newblock):
                        clash = True
                
                # Break out if no clashes
                if clash == False:
                    break
                
                # increment tries counter
                tries += 1
            
                # End Seeding loop
            if (clash == False):
                self.addBlock(newblock)
            else:
                print "ERROR ADDING BLOCK"
                sys.exit(1)
    
        # End of loop to seed cell
    
    # End seed
    

        
    def shimmy(self, nsteps ):
        """ Shuffle the molecules about making bonds where necessary for nsteps"""
        
        for step in range( nsteps ):
            
            # Number of sub-moves to attempt when the blocks are close
            minimoves = 20
            
            if not step % 20:
                print "step {}".format(step)
                
            iblock = self.getRandomBlockIndex()
            
            block = self.blocks[iblock]
            
            # Copy the origianl coordintes so we can reject the move
            # we copy the whole block so we don't need to recalculate
            # anything - not sure if this quicker the saving the coords & updating tho
            orig_block = copy.deepcopy( block )
             
            # Make a random move
            self.randomMove( block )
            
            # Remeber the original cog of the block 
            cog = block.centerOfGeometry()
            
            # Check all atoms and 
            removed = []
            for i in range( len(self.blocks) ):
                
                # skip the block we are using
                if i == iblock:
                    i+=1
                    continue
                
                oblock = self.blocks[i]
                
                # First see if we are close enough to consider bonding
                if not block.close( oblock, margin = 1.0 ):
                    i+=1
                    continue
                
                # See if we can bond and shimmy minimoves times
                # Get the original cog so we sample about it
                for j in range( minimoves ):
                    # Try to make a bond
                    bond = block.canBond( oblock )
                    if not bond or bond == "clash":
                        # Try a smaller move
                        print "minimove: {}->{}".format(j,bond)
                        self.randomSmallMove( block, cog, prange=1.5 )
                    else:
                        # Bonding worked so break out of loop
                        break
                
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
            
            for r in removed:
                self.blocks.pop(r)
                #r=r-1
                    

             
    def write(self, ofile ):
        """Write out the cell atoms to an xyz file"""
        
        natoms=0
        xyz = ""
        for block in self.blocks:
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
        for block in self.blocks:
            s+= str(block)
            
        return s