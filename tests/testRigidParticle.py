'''
Created on 4 Mar 2018

@author: jmht
'''
import math
import os
import unittest
import numpy as np

import context
BLOCKS_DIR = context.ab_paths.BLOCKS_DIR
PARAMS_DIR = context.ab_paths.PARAMS_DIR
from context import ab_block
from context import ab_bond
from context import ab_cell
from context import ab_fragment
from context import ab_rigidparticle
from context import xyz_core
from context import xyz_util

# Not sure where best to do this
xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))

class Test(unittest.TestCase):
    
    def testConfig(self):
        class Cell(object):
            """Mock Cell object for testing"""
            def __init__(self, rigidParticleMgr):
                self.rigidParticleMgr = rigidParticleMgr
        
        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        cell = Cell(rigidParticleMgr)
        
        ch4ca = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        ftype = 'A'
        f1 = ab_fragment.Fragment(filePath=ch4ca, fragmentType=ftype, cell=cell)
        
        configStr = {'1A0000': 'AB', '2A0000': 'AC', '0A0000': 'AA'}
        self.assertEqual(rigidParticleMgr._configStr, configStr)
        self.assertEqual(len(rigidParticleMgr._configStr), len(rigidParticleMgr._configs))
        
        b1 = list(f1.bodies())[0]
        quat_origin = np.array([1.0, 0.0, 0.0, 0.0])
        self.assertTrue(np.allclose(rigidParticleMgr.orientation(b1), quat_origin, rtol=0.0001))
        self.assertEqual(rigidParticleMgr.configStr(b1), "AA")
        
        # Rotate fragment and see if we get a different orientation
        axis = np.array([1, 0, 0])
        angle = math.pi / 3
        rotationMatrix = xyz_core.rotation_matrix(axis, angle)
        f1.rotate(rotationMatrix, f1._centerOfMass)
        b1 = list(f1.bodies())[0]
        ref_q = np.array([0.866, 0.5, 0.0, 0.0]) 
        self.assertTrue(np.allclose(rigidParticleMgr.orientation(b1), ref_q, atol=0.0001))
        
        return
    
    def testRigidParticle(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        ftype = 'A'
        f1 = ab_fragment.Fragment(filePath=ch4, fragmentType=ftype, cell=cell)
        rp = list(f1.bodies())[0].rigidParticle()
        quat_origin = np.array([1.0, 0.0, 0.0, 0.0])
        self.assertTrue(np.allclose(rp.orientation, quat_origin, rtol=0.0001))

    def testManager(self):
        class Cell(object):
            """Mock Cell object for testing"""
            def __init__(self, rigidParticleMgr):
                self.rigidParticleMgr = rigidParticleMgr
        
        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        cell = Cell(rigidParticleMgr)

        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        f1 = ab_fragment.Fragment(filePath=ch4, fragmentType='A', cell=cell)
        self.assertTrue(len(rigidParticleMgr._configs) == 1)
        for body in f1.bodies():
            self.assertTrue(body.rigidConfigStr in rigidParticleMgr._configs)
         
        f2 = f1.copy()
        self.assertTrue(len(rigidParticleMgr._configs) == 1)
 
        eg1 = f1.freeEndGroups()[0]
        eg2 = f2.freeEndGroups()[0]
        bond = ab_bond.Bond(eg1, eg2)
        bond.endGroup1.setBonded(bond)
        bond.endGroup2.setBonded(bond)
        self.assertTrue(len(rigidParticleMgr._configs) == 2)
 
        eg1 = f1.freeEndGroups()[-1]
        eg2 = f2.freeEndGroups()[0]
        bond = ab_bond.Bond(eg1, eg2)
        bond.endGroup1.setBonded(bond)
        bond.endGroup2.setBonded(bond)
        self.assertTrue(len(rigidParticleMgr._configs) == 4)
 
        f1 = ab_fragment.Fragment(filePath=ch4, fragmentType='B')
        self.assertTrue(len(rigidParticleMgr._configs) == 5)
         
    def testManagerConfigStr(self):
        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        idx = 0
        cid = rigidParticleMgr.calcConfigStr(idx)
        self.assertTrue(cid, 'AA')
         
        idx = 26 * 26 - 1
        cid = rigidParticleMgr.calcConfigStr(idx)
        self.assertTrue(cid, 'ZZ')
         
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
        ch4_1 = ab_block.Block(filePath=ch4_f, fragmentType='A')
        ch4_2 = ab_block.Block(filePath=ch4_f, fragmentType='A')
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
        mycell = ab_cell.Cell(boxDim)
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        mycell.libraryAddFragment(filename=ch4, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        mycell.seed(1, random=False, center=True)
        mycell.growBlocks(1, random=False)
        d = mycell.dataDict()
        # Check that the orientation quaternion can rotate the master coords onto the body
        for rp in d.rigidParticles:
            ref_coords = mycell.rigidParticleMgr._potitions[rp.b_rigidConfigStr]
            rot_coords = xyz_core.rotate_quaternion(ref_coords, rp.orientation)
            self.assertTrue(np.allclose(rp.b_positions, rot_coords, atol=0.0001))  

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
