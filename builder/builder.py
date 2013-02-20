'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import buildingBlock
import cell
import util

#1000
# 13/0
# 8/0
# 5/0
# 10/0
# 9/0
# 7/0
# 
# 9/0
# 6/0
# 8/0
# 18/0


INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/pyrene_typed.car" 
#INFILE = "/Users/jmht/Dropbox/Amorphousbuilder/pyrene_typed.car" 
#INFILE = "/Users/jmht/1.xyz" 
#INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/builder/ch4.xyz"
#OUTFILE1 = "/Users/abbietrewin/Dropbox/Amorphousbuilder/abbie_output/cell_1.xyz"
OUTFILE = "/Users/abbietrewin/AMBI/cell_1.xyz"
#OUTFILE2 = "/Users/abbietrewin/Dropbox/Amorphousbuilder/abbie_output/cell_2.xyz"
#OUTFILE2 = "/Users/jmht/cell_2.xyz"
nblocks = 400 
CELLA = [ 100,  0,  0 ]
CELLB = [ 0, 100,  0 ]
CELLC = [ 0,  0, 100 ]

STEPS = 50
NMOVES=100

BONDANGLE=180

# Create Cell and seed it with the blocks
new_cell = cell.Cell( CELLA, CELLB, CELLC )

# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock( infile = INFILE )

new_cell.seed (nblocks, firstBlock )
new_cell.write( OUTFILE )

# Loop through as many shimmy stages as required
while True:
    OUTFILE=util.newFilename(OUTFILE)
    #cell.shimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE  )
    new_cell.directedShimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE )
    new_cell.write( OUTFILE )
    response = raw_input('Do we have to do this _again_? (y/n)')
    if response.lower() == 'n':
        break
    
            
#if __name__ == '__main__':
#    pass
