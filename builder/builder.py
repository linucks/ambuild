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
OUTFILE2 = "../afterSeedLabel.xyz"
OUTFILE3 = "../afterShimmy.xyz"
nblocks = 15
CELLA = 30
CELLB = 30
CELLC = 30

STEPS = 10000
MOVES=100

# Create Cell and seed it with the blocks
mycell = cell.Cell( bondMargin=0.2, blockMargin=2.0, atomMargin=0.5, bondAngle=180, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )


#import pdb
#pdb.set_trace()
mycell.seed( nblocks, INFILE )
mycell.writeXyz( OUTFILE1 )
#mycell.writeXyz( OUTFILE2, label=True )
#mycell.fromXyz(OUTFILE2)
#mycell.shimmy(nsteps=STEPS)
mycell.newDirectedShimmy(nsteps=STEPS,nmoves=MOVES)
#mycell.writeXyz( OUTFILE3 )

## Loop through as many shimmy stages as required
#while True:
#    OUTFILE=util.newFilename(OUTFILE)
#    #cell.shimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE  )
#    new_cell.directedShimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE )
#    new_cell.write( OUTFILE )
#    response = raw_input('Do we have to do this _again_? (y/n)')
#    if response.lower() == 'n':
#        break
    

