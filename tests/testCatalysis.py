"""The tests for the cell class"""
import sys
import math
import os
import unittest

import numpy as np

import context
from context import ab_block
from context import ab_bond
from context import ab_cell
from context import ab_util

BLOCKS_DIR = context.BLOCKS_DIR
PARAMS_DIR = context.PARAMS_DIR
TESTDATA_DIR = context.TESTDATA_DIR


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ch4Car = os.path.join(BLOCKS_DIR, "ch4.car")
        cls.benzeneCar = os.path.join(BLOCKS_DIR, "benzene.car")
        cls.benzene2Car = os.path.join(BLOCKS_DIR, "benzene2.car")
        cls.pafCar = os.path.join(BLOCKS_DIR, "PAF_bb_typed.car")

        return

    @unittest.skipUnless(ab_util.HOOMDVERSION is not None, "Need HOOMD-BLUE to run")
    def testCat1Paf2(self):
        boxDim = [40, 40, 40]
        mycell = ab_cell.Cell(boxDim, paramsDir=PARAMS_DIR)
        paf = self.benzeneCar
        cat = self.benzene2Car
        mycell.libraryAddFragment(filename=paf, fragmentType="PAF")
        mycell.libraryAddFragment(filename=cat, fragmentType="cat", catalyst=True)
        mycell.addBondType("PAF:a-cat:a")
        mycell.addBondType("PAF:a-PAF:a")

        # Add PAF and grow so that we have multi-PAF blocks
        mycell.seed(2, fragmentType="PAF", point=[10.0, 10.0, 10.0], radius=5.0)
        mycell.growBlocks(
            toGrow=5, cellEndGroups="PAF:a", libraryEndGroups="PAF:a", maxTries=500
        )

        # Get the ids of the blocks
        paf1id, paf2id = mycell.blocks.keys()
        paf1 = mycell.blocks[paf1id]
        paf2 = mycell.blocks[paf2id]
        # Remove from the cell
        mycell.delBlock(paf1id)
        mycell.delBlock(paf2id)

        # Add catalyst
        mycell.seed(1, fragmentType="cat", center=True)
        # Now join the pafs to the cat
        cat = list(mycell.blocks.values())[0]
        paf1eg = paf1.freeEndGroups()[0]
        paf2eg = paf2.freeEndGroups()[0]
        cat1eg = cat.freeEndGroups()[0]
        cat.positionGrowBlock(cat1eg, paf1eg)
        bond = Bond(cat1eg, paf1eg)
        bond.engage()

        cat2eg = cat.freeEndGroups()[0]
        cat.positionGrowBlock(cat2eg, paf2eg)
        bond = ab_bond.Bond(cat2eg, paf2eg)
        bond.engage()
        mycell.newBonds = [bond]  # Hack to set newBonds

        self.assertEqual(len(mycell.blocks), 1)
        mycell.cat1Paf2(["PAF"], dt=0.00001, optCycles=10000)
        self.assertEqual(len(mycell.blocks), 2)
        return

    @unittest.skipUnless(ab_util.HOOMDVERSION is not None, "Need HOOMD-BLUE to run")
    def testCat2Paf2(self):
        """Given two catalysts bonded to each other, each with PAF blocks bonded, break the bond
        between the catalysts, move the PAFS from one catalysts to the other, and then join the PAFS
        on that catalyst with all the PAFS"""
        boxDim = [40, 40, 40]
        mycell = ab_cell.Cell(boxDim, paramsDir=PARAMS_DIR)
        paf = self.benzeneCar
        cat = self.benzene2Car
        mycell.libraryAddFragment(filename=paf, fragmentType="PAF")
        mycell.libraryAddFragment(
            filename=cat, fragmentType="cat", markBonded=True, catalyst=True
        )

        mycell.addBondType("PAF:a-PAF:a")
        mycell.addBondType("PAF:a-cat:a")
        mycell.addBondType("cat:a*-cat:a*")

        # Add a cat block and bond it to a PAF block
        mycell.seed(1, fragmentType="cat", center=True)
        mycell.growBlocks(
            toGrow=1, cellEndGroups="cat:a", libraryEndGroups="PAF:a", maxTries=500
        )
        # Add three PAF blocks to the PAF
        mycell.growBlocks(
            toGrow=3, cellEndGroups="PAF:a", libraryEndGroups="PAF:a", maxTries=500
        )

        # copy the block and get the two endGroups that we will use to position the two catalysts so they can bond
        # newblock = self.getLibraryBlock(fragmentType=fragmentType) # Create new block
        b1 = list(mycell.blocks.values())[0]

        # Make copy of first block
        b2 = b1.copy()

        # Get the two cat* endGroups
        endGroup1 = b1.freeEndGroups(endGroupTypes="cat:a*")[0]
        endGroup2 = b2.freeEndGroups(endGroupTypes="cat:a*")[0]

        # Position so they can bond
        b1.positionGrowBlock(endGroup1, endGroup2)

        # Put b2 in the cell and bond them together
        mycell.addBlock(b2)
        mycell.checkMove(b2.id)
        mycell.processBonds()
        # End setup

        self.assertEqual(len(mycell.blocks), 1)

        # Now see if we can split off the two cat blocks and join the two PAF blocks
        mycell.cat2Paf2(["PAF"], dt=0.00001, optCycles=10000)
        # mycell.dump()

        self.assertEqual(len(mycell.blocks), 3)
        return

    def testUnbonding(self):
        """Atom positions are correct on unbonding"""

        boxDim = [40, 40, 40]
        mycell = ab_cell.Cell(boxDim, paramsDir=PARAMS_DIR)
        dbb = self.benzeneCar
        cat = self.benzene2Car
        mycell.libraryAddFragment(filename=dbb, fragmentType="DBB")
        mycell.libraryAddFragment(
            filename=cat, fragmentType="cat", markBonded=False, catalyst=True
        )

        mycell.addBondType("cat:a-DBB:a")

        # Create block manually
        cat = ab_block.Block(
            filePath=cat,
            fragmentType="cat",
            catalyst=True,
            markBonded=True,
        )
        dbb1 = ab_block.Block(filePath=dbb, fragmentType="DBB", markBonded=True)
        dbb2 = ab_block.Block(filePath=dbb, fragmentType="DBB", markBonded=True)

        # Align bond along x-axis
        cat.alignAtoms(0, 1, [1, 0, 0])
        # cat in center of cell
        cat.translateCentroid([mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2])
        mycell.addBlock(cat)
        catEndGroup1 = cat.freeEndGroups()[0]
        catCapCoord1 = cat.coord(catEndGroup1.blockCapIdx)
        dbb1EndGroup = dbb1.freeEndGroups()[0]

        # Position dbb1
        cat.positionGrowBlock(catEndGroup1, dbb1EndGroup)
        mycell.addBlock(dbb1)

        catEndGroup2 = cat.freeEndGroups()[1]
        catCapCoord2 = cat.coord(catEndGroup2.blockCapIdx)
        dbb2EndGroup = dbb2.freeEndGroups()[0]

        # Position dbb2
        cat.positionGrowBlock(catEndGroup2, dbb2EndGroup)
        mycell.addBlock(dbb2)

        made = mycell.zipBlocks(bondMargin=0.5, bondAngleMargin=5)
        self.assertEqual(made, 2)

        b1 = list(mycell.blocks.values())[0]
        dbbf1 = b1.fragments[1]
        dbbf2 = b1.fragments[2]

        dbbf1.translate([1, 2, 3])
        dbbf2.translate([3, 1, 1])

        mycell.cat1Paf2(["TEB", "DBB"], optCycles=0)

        catEndGroup1 = cat.freeEndGroups()[0]
        catCapCoord1_post = cat.coord(catEndGroup1.blockCapIdx)
        catEndGroup2 = cat.freeEndGroups()[1]
        catCapCoord2_post = cat.coord(catEndGroup2.blockCapIdx)

        self.assertTrue(np.allclose(catCapCoord1, catCapCoord1_post, atol=1.0e-5))
        self.assertTrue(np.allclose(catCapCoord2, catCapCoord2_post, atol=1.0e-5))


if __name__ == "__main__":
    """
    Run the unit tests
    """
    unittest.main()
