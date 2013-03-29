'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import cell
import cPickle
import sys


INFILE = "../PAF_bb_typed.car" 
#INFILE = "../ch4_typed.car" 
#INFILE = "./afterSeed.bust.xyz"
#OUTFILE1 = "afterSeed.xyz"
#OUTFILE2 = "afterGrow.xyz"
#OUTFILE2a = "afterGrowP.xyz"
OUTFILE3 = "afterShimmy.xyz"
OUTFILE3a = "afterShimmyP.xyz"
# for ch4: 4,4,4,4
# 35,10,10,10
#
#for paf: 4, 15, 15, 15
CELLA = 40
CELLB = 40
CELLC = 40

STEPS = 1000
MOVES=100

# Create Cell and seed it with the blocks
mycell = cell.Cell( atomMargin=1.0, boxMargin=1.2, bondMargin=0.5, bondAngle=180, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )

#import pdb
#pdb.set_trace()
mycell.seed( 10, INFILE )
mycell.writeXyz( "afterSeed.xyz" )

ok = mycell.growNewBlocks(10, maxTries=50 )

mycell.writeXyz( "afterGrow.xyz" )
mycell.writeXyz( "afterGrowP.xyz", periodic=True )

# now trying growing the individual blocks
ok = mycell.joinBlocks( 5, maxTries=50 )

mycell.writeXyz( "afterJoin.xyz" )
mycell.writeXyz("afterJoinP.xyz", periodic=True )

#mycell.writeXyz("canBond.xyz", label=False)
#pfile = open("pickle1.pkl","w")
#cPickle.dump(mycell,pfile)
#pfile.close()
sys.exit(0)

#stype - 3 types: "bond", "block" or None (random)
mycell.shimmy(nsteps=STEPS, nmoves=MOVES, stype=None)
mycell.writeXyz( OUTFILE3 )
mycell.writeXyz( OUTFILE3a, periodic=True )



## Loop through as many shimmy stages as required
#while True:
#    OUTFILE=util.newFilename(OUTFILE)
#    #cell.shimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE  )
#    new_cell.directedShimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE )
#    new_cell.write( OUTFILE )
#    response = raw_input('Do we have to do this _again_? (y/n)')
#    if response.lower() == 'n':
#        break
    

