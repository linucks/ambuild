"""
Created on 4 Mar 2018

@author: jmht
"""
import os
import unittest

#
import context

BLOCKS_DIR = context.ab_paths.BLOCKS_DIR
PARAMS_DIR = context.ab_paths.PARAMS_DIR
from context import ab_util, ab_cell

import numpy as np


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ch4Car = os.path.join(BLOCKS_DIR, "ch4.car")

    @unittest.skipUnless(
        ab_util.HOOMDVERSION and ab_util.HOOMDVERSION[0] == 2,
        "Need HOOMD-BLUE 2 to run",
    )
    def testReadSnapshot(self):
        import hoomd2

        boxDim = [20, 20, 20]
        mycell = ab_cell.Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType="A")
        mycell.seed(1, center=False)
        mycell.growBlocks(1)
        coords = []
        block = list(mycell.blocks.values())[0]
        for c in block.iterCoord():
            coords.append(c)

        # 1 cycle at low T so that nothing should move
        mycell.runMD(mdCycles=1, rigidBody=True, T=0.1, dump=False, dumpPeriod=1)

        newcoords = []
        block = list(mycell.blocks.values())[0]
        for c in block.iterCoord():
            newcoords.append(c)
        for before, after in zip(coords, newcoords):
            self.assertTrue(np.allclose(before, after))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
