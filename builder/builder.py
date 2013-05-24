'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import cell
import cPickle
import sys


INFILE = "/Users/abbietrewin/Dropbox/AmorphousBuilder/pyrene_typed.car" 
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

seedCount=150
added = mycell.seed( seedCount, INFILE )
if added != seedCount:
    response = raw_input('Tried to seed with {0} but could only add {0}. Continue? (y/n)'.format(seedCount,added))
    if response.lower() == 'n':
        sys.exit(1)
    
#mycell.checkFinished()
mycell.dump()

while True:
    
    for i in range(40):
        
        print "Growing new blocks in loop {0}".format(i)
        ok = mycell.growNewBlocks(1, maxTries=50 )
        mycell.dump()
        if not ok:
            print "FAILED TO GROW"
        #mycell.checkFinished()

        ok = mycell.joinBlocks( 20, maxTries=5000 )
        print "Joining new blocks in loop {0}".format(i)
        mycell.dump()
        if not ok:
            print "FAILED TO JOIN"
        
    print "Got density: ",mycell.density()
    print "Got _endGroups: ",mycell._endGroups()
    response = raw_input('Do we have to do this _again_? (y/n)')
    if response.lower() == 'n':
        break
    

