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


INFILE = "../ch4_typed.car" 
INFILE = "../PAF_bb_typed.car" 
OUTFILE1 = "../afterSeed.xyz"
OUTFILE2 = "../afterShimmy.xyz"
nblocks = 11
CELLA = 100
CELLB = 100
CELLC = 100

STEPS = 5
MOVES=10

# Create Cell and seed it with the blocks
mycell = cell.Cell( bondMargin=0.2, blockMargin=2.0, atomMargin=0.5, bondAngle=180, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )

mycell.seed( nblocks, INFILE )
mycell.writeXyz( OUTFILE1 )
#mycell.shimmy(nsteps=STEPS,nmoves=MOVES)
#mycell.writeXyz( OUTFILE2 )

## Loop through as many shimmy stages as required
#while True:
#    OUTFILE=util.newFilename(OUTFILE)
#    #cell.shimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE  )
#    new_cell.directedShimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE )
#    new_cell.write( OUTFILE )
#    response = raw_input('Do we have to do this _again_? (y/n)')
#    if response.lower() == 'n':
#        break
    

