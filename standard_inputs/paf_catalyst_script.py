#!/usr/bin/env python

"""
1. addBondType('PAF:a-cat:a' ) - PAFs can bond to the free catalyst
2. addBondType('cat:a*-cat:a*' )  - catalysts with PAFs bonded to them can bond to each other
3. When PAF:a has bonded to cat:a, we change the bond type of cat:a  to cat:a*
4. When two catalysts with PAF bonded to them form a bond, we unbond the catalysts from the two PAF
   blocks and and bond the PAF blocks to each other

Questions
Do we do a multi-step operation - e.g. bond the catalysts together and then remove one and join the PAFS
 - or, as soon as two PAF-bonded cat groups can bond, we delete the two PAF/cat bonds and directly bond
 the PAF units. What happens to the catalysts? Presuamably we just leave them? 
 
 
 
 If the last bond was between two cat:a* groups:
     get the 4 fragments 2*cat & 2*paf
     unbond all 4, but keep track of the two paf endGroups
     bond the two paf endgroups
 
 
 MISC
 Fix debugging so it's easier to turn on and off in a standard way
 
 

"""

import sys
sys.path.append("/opt/ambuild/builder")

# Python imports
import cPickle
import csv

# Our imports
import buildingBlock
import cell
import util

#cell dimensions:
boxDim=[40,40,40]

#Create Cell and seed it with the blocks
mycell = cell.Cell(boxDim, atomMargin=0.1, bondMargin=0.5, bondAngleMargin=5, doLog=True )

#import the two fragment files if you have 2 different building blocks
#mycell.libraryAddFragment( filename='/opt/ambuild/blocks/PAF.car', fragmentType='PAF' )
mycell.libraryAddFragment( filename='/opt/ambuild/blocks/ch4.car', fragmentType='PAF' )
#mycell.libraryAddFragment( filename='pafh2.car', fragmentType='cat' )
mycell.libraryAddFragment( filename='/opt/ambuild/blocks/ch4.car', fragmentType='cat' )

mycell.addBondType( 'PAF:a-cat:a' )
mycell.addBondType( 'cat:a*-cat:a*' )

# Add a cat block and bond it to a PAF block
mycell.seed( 1, fragmentType='cat', center=True)
mycell.growBlocks(toGrow=1, cellEndGroups='cat:a', libraryEndGroups='PAF:a', maxTries=500)

# copy the block and get the two endGroups that we will use to position the two catalysts so they can bond
#newblock = self.getLibraryBlock(fragmentType=fragmentType) # Create new block
b1 = mycell.blocks.values()[0]

# Make copy of first block
b2 = b1.copy()

# Get the two cat* endGroups
endGroup1 = b1.freeEndGroups(endGroupTypes='cat:a*')[0]
endGroup2 = b2.freeEndGroups(endGroupTypes='cat:a*')[0]

# Position so they can bond
b1.positionGrowBlock(endGroup1, endGroup2)

# Put b2 in the cell
mycell.addBlock(b2)

mycell.checkMove(b2.id)
mycell.processBonds()
mycell.dump()

mycell.unbondCat()

mycell.dump()


# if mycell.unbondCat():
#     print "unbondCat2"

