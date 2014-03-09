#!/usr/bin/env python
import sys
sys.path.append("/home/abbiet/ambuild.new/builder")
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

#mycell = cellFromPickle("/opt/ambuild/work/CTF/CTF-2-interpen/jensNewCap/step_17.pkl")


def run():

    #cell dimensions:
    CELLA = CELLB = CELLC = 400
       
    # Create Cell and seed it with the blocks
    mycell = cell.Cell( atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15 )
    mycell.cellAxis( CELLA, CELLB, CELLC )
       
    #import the two fragment files if you have 2 different building blocks
    #fragA="/opt/ambuild/work/building_blocks/PAF_bb_typed.car"
    fragA="../blocks/node_2.car"
    fragB="../blocks/linker_2.car"

    mycell.addInitBlock( filename=fragA, fragmentType='A' )
    mycell.addInitBlock( filename=fragB, fragmentType='B' )
    mycell.addBondType( 'A-B' )

    # Use center argument to place the first block at the center of the cell. Only for first seed.
    mycell.seed( 20, fragmentType='B', center=True )

    # Now loop adding blocks
    while True:
    
        toGrow=20
        if mycell.growBlocks( toGrow, fragmentType=None, maxTries=500 ) != toGrow:
            print "GROW FAILED!"
            mycell.dump()
            break

        print "GROW BLOCKS DUMP"
        mycell.dump()
    
        # Zip doesn't fail as such so don't need to check the result
        mycell.zipBlocks( bondMargin=5, bondAngleMargin=45 )
        print "ZIP BLOCKS DUMP"
        mycell.dump()
    
        if not mycell.runMDAndOptimise( doDihedral=True, quiet=False ):
            mycell.dump()
            print "OPTIMISATION FAILED!"
            sys.exit()
    
        mycell.dump()

    
    if False:
        mycell.capBlocks( fragmentType='A', filename="../blocks/cap_node_2.car" )
        mycell.dump()
        mycell.capBlocks( fragmentType='B', filename="../blocks/cap_linker.car" )
        mycell.dump()

    if not mycell.optimiseGeometry(doDihedral=True):
        print "FINAL OPTIMISATION FAILED!"
    mycell.dump()


if __name__ == "__main__":
    run()