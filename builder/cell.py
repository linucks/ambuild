'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import sys
import os
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
        
    def getRandomPosition( self ) :
        '''
        Return a random position in the cell
        '''
        #
        x = random.uniform(0,self.A[0])
        y = random.uniform(0,self.B[1])
        z = random.uniform(0,self.C[2])
        
        return numpy.array([x,y,z])
    
    def getRandomRotationAngle(self):
        """Return a random rotation in radians"""
        return random.uniform( 0, 2*numpy.pi)
        
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
                    
                # Get coords of random point in the cell
                position = self.getRandomPosition()
                
                # Translate COM of molecule to coords
                newblock.translateCenterOfGeometry( position )
                
                # Rotate by a random number
                rotation = self.getRandomRotationAngle()
                
                newblock.rotate( numpy.array([1,0,0],dtype=numpy.float64), rotation )
                
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
    

        
    def shimmy(self, nsteps = 100 ):
        """ Shuffle the molecules about making bonds where necessary for nsteps"""
        
        for step in range( nsteps ):
            iblock = self.randomBlockIndex()
            block = self.blocks[iblock]
             
            # Make a random move
            self.randomMove( block )
             
             # Check all atoms and 
            for i, oblock in enumerate( self.blocks ):
                 
                # skip the block we are using
                if i == iblock:
                    continue
                
                # Try to make a bond
                bond = block.tryBond( oblock )
                
                if bond:
                    print "Bonded block {} with block {}\n".format( iblock, i) 

             
    def write(self, ofile ):
        """Write out the cell atoms to an xyz file"""
        
        natoms=0
        xyz = ""
        for block in self.blocks:
            for i, c in enumerate( block.coords ):
                xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( block.labels[i], c[0], c[1], c[2] )
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