'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import sys
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

    def getRandomPosition( self ) :
        '''
        Return a random position in the cell
        '''
        #
        x = random.uniform(0,self.A[0])
        y = random.uniform(0,self.B[1])
        z = random.uniform(0,self.C[2])
        
        return numpy.array([x,y,z])
    
    def getRandomRotation(self):
        """Return a random rotation in radians"""
        return random.uniform( 0, 2*numpy.pi)
        
    def blocks(self):
        return self.blocks
    
    def addBlock( self, block ):
        self.blocks.append( block )
        
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
                rotation = self.getRandomRotation()
                
                newblock.rotate(rotation)
                
                # Test for Clashes with other molecules - always set clash
                # to False at start so it only becomes true if no clashes
                clash = False
                for block in self.blocks():
                    if block.clashes(newblock):
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