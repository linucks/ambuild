'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import buildingBlock
import cell

INFILE = "/Users/jmht/Dropbox/Amorphousbuilder/builder//PAF_bb_typed.car" 
INFILE = "/Users/jmht/Dropbox/Amorphousbuilder/builder/ch4.xyz"
OUTFILE = "/Users/jmht/cell.xyz"
nblocks = 100
CELLA = [ 50, 0,   0 ]
CELLB = [ 0,   50, 0 ]
CELLC = [ 0,   0,   50 ]
MAXTRIES = 10000



# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock( infile = INFILE )

# Create Cell and seed it with the blocks
cell = cell.Cell( CELLA, CELLB, CELLC )
cell.seed (nblocks, firstBlock )
#print cell
cell.write( OUTFILE )


#if __name__ == '__main__':
#    pass