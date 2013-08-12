'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import cell
import cPickle
import sys


#INFILE = "/Users/abbietrewin/Dropbox/AmorphousBuilder/pyrene_typed.car" 
INFILE = "../ch4_typed.car" 
INFILE = "../PAF_bb_typed.car" 
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
CELLA = 50
CELLB = 50
CELLC = 50


# Create Cell and seed it with the blocks
mycell = cell.Cell( atomMargin=0.1, boxMargin=0.5, bondMargin=0.5, bondAngle=180, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )

mycell.addInitBlock( filename="../PAF_bb_typed.car", fragmentType='A' )
mycell.addInitBlock( filename="../ch4_typed.car", fragmentType='B' )
mycell.addBondType( 'A-B' )
#mycell.addBondType( 'A-A' )

mycell.seed( 1, fragmentType='B' )
ok = mycell.growNewBlocks( 1, fragmentType='A', maxTries=500 )
ok = mycell.growNewBlocks( 1, fragmentType='A', maxTries=500 )
#ok = mycell.growNewBlocks( 1, fragmentType='A', maxTries=500 )
#ok = mycell.growNewBlocks( 1, fragmentType='B', maxTries=500 )
#ok = mycell.joinBlocks( 10, maxTries=500 )

mycell.dump()

#mycell.optimiseGeometry()

#mycell.dump()


sys.exit(0)

seedCount=1
added = mycell.seed( seedCount, INFILE )
print "Added {0} blocks".format(added)
if added != seedCount:
    response = raw_input('Tried to seed with {0} but could only add {1}. Continue? (y/n)'.format(seedCount,added))
    if response.lower() == 'n':
        sys.exit(1)

ok = mycell.growNewBlocks(100, maxTries=500 )
mycell.dump()
sys.exit()

#mycell.checkFinished()
mycell.dump()
#mycell.writeHoomdXml("hoomd.xml")

while True:
    
    for i in range(10):
        
        print "Growing new blocks in loop {0}".format(i)
        ok = mycell.growNewBlocks(1, maxTries=50 )
        mycell.dump()
        if not ok:
            print "FAILED TO GROW"
            sys.exit(1)
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