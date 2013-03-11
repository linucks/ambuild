'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import buildingBlock
import cell
import util


INFILE = "../ch4_typed.car" 
INFILE = "../PAF_bb_typed.car" 
OUTFILE1 = "../afterSeed.xyz"
OUTFILE2 = "../afterShimmy.xyz"
nblocks = 50
CELLA = 100
CELLB = 100
CELLC = 100

STEPS = 1000
MOVES=100

# Create Cell and seed it with the blocks
mycell = cell.Cell( bondMargin=0.2, blockMargin=2.0, atomMargin=0.5, bondAngle=180, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )

#import pdb
#pdb.set_trace()
mycell.seed( nblocks, INFILE )
mycell.writeXyz( OUTFILE1 )
mycell.shimmy(nsteps=STEPS,nmoves=MOVES)
mycell.writeXyz( OUTFILE2 )

## Loop through as many shimmy stages as required
#while True:
#    OUTFILE=util.newFilename(OUTFILE)
#    #cell.shimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE  )
#    new_cell.directedShimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE )
#    new_cell.write( OUTFILE )
#    response = raw_input('Do we have to do this _again_? (y/n)')
#    if response.lower() == 'n':
#        break
    

