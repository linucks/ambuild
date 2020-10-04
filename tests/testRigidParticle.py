"""
Created on 4 Mar 2018

@author: jmht
"""
import math
import os
import unittest
import numpy as np

import context
from context import ab_block
from context import ab_bond
from context import ab_cell
from context import ab_fragment
from context import ab_rigidparticle
from context import xyz_core
from context import xyz_util

BLOCKS_DIR = context.BLOCKS_DIR
PARAMS_DIR = context.PARAMS_DIR

# Not sure where best to do this
xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))


class Test(unittest.TestCase):
    def testConfig(self):
        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        ch4ca = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        ftype = "A"
        f1 = ab_fragment.fragmentFactory(ftype, ch4ca)
        rigidParticles = []
        for body in f1.bodies():
            rigidParticles.append(rigidParticleMgr.createParticle(body))

        pfx = rigidParticleMgr.PREFIX_STR
        configStr = {"1A0000": pfx + "AB", "2A0000": pfx + "AC", "0A0000": pfx + "AA"}
        self.assertEqual(rigidParticleMgr._configStr, configStr)
        self.assertEqual(
            len(rigidParticleMgr._configStr), len(rigidParticleMgr._positions)
        )

        b1 = list(f1.bodies())[0]
        quat_origin = np.array([1.0, 0.0, 0.0, 0.0])
        self.assertEqual(rigidParticleMgr.configStr(b1), pfx + "AA")

        # Rotate fragment and see if we get a different orientation
        axis = np.array([1, 0, 0])
        angle = math.pi / 3
        rotationMatrix = xyz_core.rotation_matrix(axis, angle)
        f1.rotate(rotationMatrix, f1._centerOfMass)
        b1 = list(f1.bodies())[0]
        ref_q = np.array([0.866, 0.5, 0.0, 0.0])
        self.assertTrue(
            np.allclose(rigidParticles[0].orientation, quat_origin, rtol=0.0001)
        )
        return

    def testManagerConfigStr(self):
        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        idx = 0
        cid = rigidParticleMgr.calcConfigStr(idx)
        pfx = rigidParticleMgr.PREFIX_STR
        self.assertTrue(cid, pfx + "AA")

        idx = 26 * 26 - 1
        cid = rigidParticleMgr.calcConfigStr(idx)
        self.assertTrue(cid, pfx + "ZZ")

        idx = idx + 1
        try:
            rigidParticleMgr.calcConfigStr(idx)
        except AssertionError:
            pass
        except Exception as e:
            self.fail("Unexpected exception: {}".format(e))

    def testOrientation1(self):
        """Two bonded fragments: check that rotation rotates one to the other"""

        # 2 blocks
        ch4_f = os.path.join(BLOCKS_DIR, "ch4.car")
        ch4_1 = ab_block.Block(filePath=ch4_f, fragmentType="A")
        ch4_2 = ab_block.Block(filePath=ch4_f, fragmentType="A")
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = ab_bond.Bond(eg1, eg2)
        bond.engage()

        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        rigidParticles = []
        for frag in ch4_1.fragments:
            for body in frag.bodies():
                rigidParticles.append(rigidParticleMgr.createParticle(body))

        # Check that the orientation quaternion can rotate the master coords onto the body
        for rp in rigidParticles:
            ref_coords = rigidParticleMgr._positions[rp.b_rigidConfigStr]
            rot_coords = xyz_core.rotate_quaternion(ref_coords, rp.orientation)
            self.assertTrue(np.allclose(rp.b_positions, rot_coords, atol=0.0001))

    def testOrientation2(self):
        """Two bonded fragments: check that rotation rotates one to the other"""
        #         class Cell(object):
        #             """Mock Cell object for testing"""
        #             def __init__(self, rigidParticleMgr):
        #                 self.rigidParticleMgr = rigidParticleMgr
        #
        #         rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        #         cell = Cell(rigidParticleMgr)

        boxWidth = 20
        boxDim = [boxWidth, boxWidth, boxWidth]
        mycell = ab_cell.Cell(boxDim, paramsDir=PARAMS_DIR)
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        mycell.libraryAddFragment(filename=ch4, fragmentType="A")
        mycell.addBondType("A:a-A:a")

        mycell.seed(1, random=False, center=True)
        mycell.growBlocks(1, random=False)
        d = mycell.cellData()
        # Check that the orientation quaternion can rotate the master coords onto the body
        for rp in d.rigidParticles:
            ref_coords = mycell.rigidParticleMgr._positions[rp.b_rigidConfigStr]
            rot_coords = xyz_core.rotate_quaternion(ref_coords, rp.orientation)
            self.assertTrue(np.allclose(rp.b_positions, rot_coords, atol=0.0001))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
