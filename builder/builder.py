'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import cell
import cPickle
import sys


INFILE = "/Users/abbietrewin/Dropbox/AmorphousBuilder/PAF_bb_typed.car" 
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


# Create Cell and seed it with the blocks
mycell = cell.Cell( atomMargin=1.0, boxMargin=1.2, bondMargin=0.5, bondAngle=180, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )

mycell.seed( 15, INFILE )
#mycell.checkFinished()
mycell.dump()

while True:
    
    for i in range(10):
        
        print "Growing new blocks in loop {0}".format(i)
        ok = mycell.growNewBlocks(15, maxTries=50 )
        mycell.dump()
        if not ok:
            print "FAILED TO GROW"
        #mycell.checkFinished()

        ok = mycell.joinBlocks( 10, maxTries=50 )
        print "Joining new blocks in loop {0}".format(i)
        mycell.dump()
        if not ok:
            print "FAILED TO JOIN"
        
    response = raw_input('Do we have to do this _again_? (y/n)')
    if response.lower() == 'n':
        break
    

