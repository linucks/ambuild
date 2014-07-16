#!/usr/bin/env python

import sys
sys.path.append("/opt/ambuild/builder")
import cell

# Create Cell and seed it with the blocks
cellDim=[100,100,100]
mycell = cell.Cell(cellDim,atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15)
   
#import the two fragment files if you have 2 different building blocks
fragA="/opt/ambuild/blocks/amine_typed.car"
fragB="/opt/ambuild/blocks/triquin_typed.car"

mycell.libraryAddFragment( filename=fragA, fragmentType='amine' )
mycell.libraryAddFragment( filename=fragB, fragmentType='triquin' )
mycell.addBondType( 'amine:a-triquin:b' )

mycell.seed( 1, fragmentType='triquin', center=True )
mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups=None, maxTries=500)
mycell.dump()
