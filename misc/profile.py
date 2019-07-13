#!/usr/bin/env python
"""
Script for profiling code.
See: https://github.com/jrfonseca/gprof2dot

To install dependencies on OSX
python -m pip install gprof2dot
brew install graphviz

Then:
gprof2dot -f pstats restats | dot -Tpng -o output.png
"""
import os
import sys
sys.path.insert(0, "../ambuild")
from ab_cell import Cell
from ab_util import cellFromPickle
from ab_paths import BLOCKS_DIR

def make_cell():
    boxDim = [20, 20, 20]
    mycell = Cell(boxDim)
    mycell.libraryAddFragment(filename=os.path.join(BLOCKS_DIR, 'ch4.car'), fragmentType='A')
    mycell.addBondType('A:a-A:a')
    
    mycell.seed(1, random=False, point=[12.5, 12.5, 12.5])
    mycell.seed(1, random=False, point=[12.5, 12.5, 17.5])
    mycell.seed(1, random=False, point=[12.5, 17.5, 12.5])
    mycell.seed(1, random=False, point=[12.5, 17.5, 17.5])
    mycell.seed(1, random=False, point=[17.5, 12.5, 12.5])
    mycell.seed(1, random=False, point=[17.5, 12.5, 17.5])
    mycell.seed(1, random=False, point=[17.5, 17.5, 12.5])
    mycell.seed(1, random=False, point=[17.5, 17.5, 17.5])
    mycell.growBlocks(8*30, random=False)
    mycell.dump()
    return mycell

#mycell = make_cell()
mycell = cellFromPickle('step_1.pkl')

import cProfile
cProfile.run('mycell.zipBlocks()', 'restats')
