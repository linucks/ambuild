#!/usr/bin/env python

# Annotated example amorphous builder script
import sys

# If you run ambi from a directory other than the builder directory, you will need to put the 
# path to the builder directory in builderHome variable below
builderHome="/home/abbiet/ambuild.new/builder"
sys.path.insert(0,builderHOme)

# This imports the builder cell module - this is the only module that should be required
import cell



# mycell.runMD(mdCycles=5000,dump=True)
##mycell = cellFromPickle("step_1.pkl")
#sys.exit()

#cell dimensions:
CELLA = CELLB = CELLC = 20
   
# Create Cell and seed it with the blocks
mycell = cell.Cell( atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15, doLog=False )
mycell.cellAxis( CELLA, CELLB, CELLC )
   
#import the two fragment files if you have 2 different building blocks
#fragA="/opt/ambuild/work/building_blocks/PAF_bb_typed.car"
fragA="../blocks/ch4.car"
fragB="../blocks/benzene.car"

mycell.addInitBlock( filename=fragA, fragmentType='A' )
mycell.addInitBlock( filename=fragB, fragmentType='B' )
mycell.addBondType( 'A:a-B:a' )

# Use center argument to place the first block at the center of the cell. Only for first seed.
mycell.seed( 10, fragmentType='B', center=True )
mycell.seed( 10, fragmentType='A', center=True )

mycell.growBlocks( 100)
mycell.dump()

#cProfile.run('mycell.zipBlocks( bondMargin=5, bondAngleMargin=20)', 'restats' )
#cProfile.run('mycell.growBlocks( 1000)', 'restats' )

sys.exit()

toGrow=10
if mycell.growBlocks( toGrow, libraryEndGroups=['B:a'], maxTries=500 ) != toGrow:
    print "GROW FAILED!"
mycell.dump()

added = mycell.zipBlocks( bondMargin=5, bondAngleMargin=45 )
if added > 0:
    mycell.dump()

sys.exit()

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
    mycell.zipBlocks( bondMargin=5, bondAngleMargin=10 )
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


def cellFromPickle(pickleFile):
    with open(pickleFile) as f:
        myCell=cPickle.load(f)
    return myCell



#if __name__ == "__main__":
#    run()