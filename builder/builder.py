'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import buildingBlock
import cell

INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/pyrene_typed.car" 
INFILE = "/Users/jmht/Dropbox/Amorphousbuilder/pyrene_typed.car" 
INFILE = "/Users/jmht/1.xyz" 
#INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/builder/ch4.xyz"
OUTFILE1 = "/Users/abbietrewin/Dropbox/Amorphousbuilder/abbie_output/cell_1.xyz"
OUTFILE1 = "/Users/jmht/cell_1.xyz"
OUTFILE2 = "/Users/abbietrewin/Dropbox/Amorphousbuilder/abbie_output/cell_2.xyz"
OUTFILE2 = "/Users/jmht/cell_2.xyz"
nblocks = 1000 
CELLA = [ 50,  0,  0 ]
CELLB = [ 0, 50,  0 ]
CELLC = [ 0,  0, 50 ]
STEPS = 100
NMOVES=50


# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock( infile = INFILE )

# Create Cell and seed it with the blocks
cell = cell.Cell( CELLA, CELLB, CELLC )
cell.seed (nblocks, firstBlock )
#print cell
cell.write( OUTFILE1 )

cell.shimmy( nsteps=STEPS, nmoves=NMOVES  )

cell.write( OUTFILE2 )
#if __name__ == '__main__':
#    pass