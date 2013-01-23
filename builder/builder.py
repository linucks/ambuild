'''
Created on Jan 15, 2013

@author: abbietrewin



'''
#import numpy
import sys


import buildingBlock
import cell


CARFILE = "FOO" 
nblocks =50
CELLA = [ 100, 0,   0 ]
CELLB = [ 0,   100, 0 ]
CELLC = [ 0,   0,   100 ]
MAXTRIES = 10000



# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock()
firstBlock.fromCarFile( CARFILE )
com = firstBlock.centerOfMass()
com = firstBlock.centerOfMass()

# Create Cell and seed it with the blocks
cell = cell.Cell( CELLA, CELLB, CELLC )
cell.seed(nblocks, firstBlock)


if __name__ == '__main__':
    pass