#!/usr/bin/env python

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
mycell = cell.Cell(boxDim, atomMargin=0.1, bondMargin=0.5, bondAngleMargin=5, doLog=False )

#import the two fragment files if you have 2 different building blocks
mycell.libraryAddFragment( filename='PAF.car', fragmentType='PAF' )
mycell.libraryAddFragment( filename='pafh2.car', fragmentType='cat' )
mycell.addBondType( 'PAF:a-cat:a' )

def unbondCat(mycell):
    """Function to unbond a Ni-catalyst bonded to two PAF groups"""
    
    if not mycell.lastAdded: return
    
    lblock = mycell.blocks[mycell.lastAdded]
    # Last added bond is always last in the list
    bond = lblock._blockBonds[-1]
    
    # See if either of the blocks connected is the cat
    t1 = bond.endGroup1.type()
    t2 = bond.endGroup2.type()
    cfrag = None
    if t1 == 'cat:a':
        cfrag = bond.endGroup1.fragment
    elif t2 == 'cat:a':
        cfrag = bond.endGroup2.fragment
    else: return # Nothing to do
    
    # The block has multiple fragments, but we are only interested in this cat fragment
    # and if this has two bonds made to it.
    # For time being assume only two allowed bonds to cat
    endGroups = cfrag.endGroups()
    if not all([ eg.bonded for eg in endGroups ]): return
    assert len(endGroups) == 2, "Assumption is cat only has 2 endGroups!"
    eg1, eg2 = endGroups
    
    # Both are bonded so we now need to go through each bond in the block in turn and unbond the blocks
    
    # Get two bonds in this block that bond it to PAF and also the two paf endgroups
    b1 = None
    b2 = None
    paf1Eg = None
    paf2Eg = None
    for bond in lblock._blockBonds:
        #if bond.endGroup1 == eg1 or bond.endGroup2 == eg1 or bond.endGroup1 == eg2 or bond.endGroup2 == eg2:
        if bond.endGroup1 in [eg1, eg2]:
            assert not (b1 and b2), "More than 2 bonds linked to cat"
            if b1:
                b2 = bond
                paf2Eg = bond.endGroup2
            else:
                b1 = bond
                paf1Eg = bond.endGroup2
        elif bond.endGroup2 in [eg1, eg2]:
            assert not (b1 and b2), "More than 2 bonds linked to cat"
            if b1:
                b2 = bond
                paf2Eg = bond.endGroup1
            else:
                b1 = bond
                paf1Eg = bond.endGroup1
    
    assert b1 and b2,"Could not find both bonds"
    
    # Remove the block from the cell
    mycell.delBlock(lblock.id)
    
    # Break the two bonds
    paf1 = lblock.deleteBond(b1)
    #print "DELETED PAF1 ",paf1.id
    paf2 = lblock.deleteBond(b2)
    #print "DELETED PAF2 ",paf2.id
    
    # Add the unbonded blocks back to the cell
    mycell.addBlock(lblock)
    mycell.addBlock(paf1)
    mycell.addBlock(paf2)
    
    # We now need to bond the two PAF groups
    assert paf1Eg.free() and paf2Eg.free(),"PAF endgroups aren't free!"
    bond = buildingBlock.Bond(paf1Eg, paf2Eg)
    mycell.bondBlock(bond)
    
    return

mycell.seed( 1, fragmentType='cat', center=True)
mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups='PAF:a', maxTries=500)
unbondCat(mycell)
mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups='PAF:a', maxTries=500)
mycell.dump()
unbondCat(mycell)
mycell.dump()

