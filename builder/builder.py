'''
Created on Jan 15, 2013

@author: abbietrewin



'''
#import numpy
#import sys


import buildingBlock
import cell


FILE = "/Users/jmht/Dropbox/Amorphousbuilder/builder//PAF_bb_typed.car" 
FILE = "/Users/jmht/Dropbox/Amorphousbuilder/builder/ch4.xyz" 
nblocks =50
CELLA = [ 100, 0,   0 ]
CELLB = [ 0,   100, 0 ]
CELLC = [ 0,   0,   100 ]
MAXTRIES = 10000



# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock()
#firstBlock.fromCarFile( CARFILE )
firstBlock.fromXyzFile( FILE )


# Create Cell and seed it with the blocks
cell = cell.Cell( CELLA, CELLB, CELLC )
cell.seed(nblocks, firstBlock)
print cell


if __name__ == '__main__':
    pass