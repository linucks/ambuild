#!/usr/bin/env python
import sys
sys.path.append("/opt/ambuild/builder")
'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import cell
import cPickle
import sys


def cellFromPickle(pickleFile):
    with open(pickleFile) as f:
        myCell=cPickle.load(f)
    return myCell


def shake( mycell, tries=5 ):
    for i in range( tries ):
        
        print "SHAKING CELL TRY NUMBER ",i
        
        if not mycell.runMD( mdCycles=1000,
                      rCut=None,
                      T=0.1,
                      tau=0.5,
                      dt=0.005,
                      doDihedral=True,
                      doImproper=False,
                      quiet=False,
                      ):
            print "SHAKE RUNMD FAILED!"
            mycell.dump()
            sys.exit()
            
        if mycell.zipBlocks( bondMargin=5, bondAngleMargin=180 ) > 0:
            print "SHAKE ZIPBLOCKS WORKED SO OPTIMISING"
            if not mycell.optimiseGeometry( doDihedral=True, quiet=True ):
                mycell.dump()
                print "SHAKE ZIPBLOCKS OPTIMISATION FAILED!"
                sys.exit()
            
    return
 
#mycell = cellFromPickle("/opt/ambuild/work/CTF/CTF-2-interpen/jensNewCap/step_17.pkl")

#cell dimensions:
CELLA = CELLB = CELLC = 100
       
# Create Cell and seed it with the blocks
mycell = cell.Cell( atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15 )
mycell.cellAxis( CELLA, CELLB, CELLC )
       
#import the two fragment files if you have 2 different building blocks
#fragA="/opt/ambuild/work/building_blocks/PAF_bb_typed.car"
fragA="../node_2.car"
fragB="../linker_2.car"

mycell.addInitBlock( filename=fragA, fragmentType='A' )
mycell.addInitBlock( filename=fragB, fragmentType='B' )
mycell.addBondType( 'A-B' )

# Use center argument to place the first block at the center of the cell. Only for first seed.
mycell.seed( 4, fragmentType='B', center=True )

# Now loop adding blocks
for i in range( 1 ):
    
    toGrow=20
    maxTries=500
    if mycell.growBlocks( toGrow, fragmentType=None, maxTries=maxTries ) != toGrow:
        print "growBlocks failed after {0} tries".format( maxTries )
        mycell.dump()
        sys.exit()
        
    if not mycell.runMDAndOptimise( doDihedral=True,
                                    doImproper=False,
                                    mdCycles=1000,
                                    optCycles=100000,
                                    maxOptIter=100,
                                    T=1000.0,
                                    tau=0.5,
                                    dt=0.005,
                                    rCut=None,
                                    quiet=False,
                                    ):
        print "RUNMDANDOPTIMISE FAILED!"
        mycell.dump()
        sys.exit()
    
    #shake( mycell, tries=5 )
    mycell.dump()

sys.exit()

    
mycell.capBlocks( fragmentType='A', filename="../../cap_node_2.car" )
mycell.dump()
mycell.capBlocks( fragmentType='B', filename="../../cap_linker.car" )
mycell.dump()
if not mycell.optimiseGeometry(doDihedral=True, minCell=False):
    print "FINAL OPTIMISATION FAILED!"

mycell.dump()


