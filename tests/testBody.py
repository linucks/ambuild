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
from context import ab_fragment
from context import ab_rigidparticle
from context import xyz_core
from context import xyz_util

# Not sure where best to do this
xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))

class Test(unittest.TestCase):
    
    def testBodyConfigStr(self):
        ch4ca = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        ftype = 'A'
        f1 = ab_fragment.Fragment(filePath=ch4ca, fragmentType=ftype)
        for i, b in enumerate(f1.bodies()):
            bstr = "{}{}{}".format(i, ftype, "0000")
            self.assertEqual(bstr, b.rigidConfigStr)
    
    def testOrientation1(self):
        
        class Cell(object):
            """Mock Cell object for testing"""
            def __init__(self, rigidParticleMgr):
                self.rigidParticleMgr = rigidParticleMgr
        
        rigidParticleMgr = ab_rigidparticle.RigidParticleManager()
        cell = Cell(rigidParticleMgr)
        
        ch4ca = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        ftype = 'A'
        f1 = ab_fragment.Fragment(filePath=ch4ca, fragmentType=ftype, cell=cell)
        b1 = list(f1.bodies())[0]
        
        tref = "{}{}{}".format(0, ftype, "0000")
        self.assertEqual(b1.rigidConfigStr, tref)
        tref = "AA"
        self.assertEqual(b1.rigidType, tref)
        self.assertTrue(np.allclose(b1.orientation, np.array([1.0, 0.0, 0.0, 0.0]), atol=0.0001))
        
    def testPrincipalMoments(self):
        """Test principal moments don't change with rotation"""
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = ab_fragment.Fragment(filePath=ch4, fragmentType='A')
 
        body = list(frag.bodies())[0]
        pi1 = body.principalMoments
        
        axis = [1, 0, 1]
        angle = math.pi/3
        rotationMatrix = xyz_core.rotation_matrix(axis, angle)
        center = body.centreOfMass
        frag.rotate(rotationMatrix, center)
        pi2 = body.principalMoments
        
        self.assertTrue(np.allclose(pi1, pi2))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
