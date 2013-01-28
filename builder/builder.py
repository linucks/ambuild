'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import buildingBlock
import cell

INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/PAF_bb_typed.car" 
INFILE = "/Users/jmht/1.xyz" 
INFILE = "/Users/jmht/Dropbox/Amorphousbuilder/PAF_bb_typed.car" 
#INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/builder/ch4.xyz"
OUTFILE1 = "/Users/jmht/cell_1.xyz"
OUTFILE2 = "/Users/jmht/cell_2.xyz"
nblocks = 50
CELLA = [ 50,  0,  0 ]
CELLB = [ 0, 50,  0 ]
CELLC = [ 0,  0, 50 ]
MAXTRIES =5000
STEPS = 1000
NMOVES=100


# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock( infile = INFILE )

# Create Cell and seed it with the blocks
cell = cell.Cell( CELLA, CELLB, CELLC )
cell.seed (nblocks, firstBlock )
#print cell
cell.write( OUTFILE1 )

cell.shimmy( STEPS, nsteps=STEPS, nmoves=NMOVES  )

cell.write( OUTFILE2 )
#if __name__ == '__main__':
#    pass