#!/usr/bin/env python
import os, sys
ambuild_home = "/opt/ambuild.git"
sys.path.insert(0, ambuild_home)

# This imports the builder cell module - this is the only module that should be required
from ambuild import ab_cell


# Create Cell and seed it with the blocks
cellDim=[100,100,100]
mycell = ab_cell.Cell(cellDim,atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15, paramsDir='/opt/paramsDir')
   
#import the two fragment files if you have 2 different building blocks
fragA = os.path.join(ambuild_home, "blocks/amine_typed.car")
fragB = os.path.join(ambuild_home, "blocks/triquin_typed.car")

mycell.libraryAddFragment( filename=fragA, fragmentType='amine' )
mycell.libraryAddFragment( filename=fragB, fragmentType='triquin' )
mycell.addBondType( 'amine:a-triquin:b' )

mycell.seed( 1, fragmentType='triquin', center=True )
mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups=None, maxTries=500)
mycell.dump()
