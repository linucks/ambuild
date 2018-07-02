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
            self.assertEqual(bstr, b.rigidType)
    
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
        self.assertEqual(b1.rigidType, tref)
        self.assertTrue(np.allclose(b1.orientation, np.array([1.0, 0.0, 0.0, 0.0]), atol=0.0001))
        
    
    def XtestMomentOfInertia(self):
        Iref = np.array([[  1.58122888e+00,  -2.49643156e-17,   3.63000000e-07],
                        [ -2.49643156e-17,   1.58122879e+00,   1.36956935e-17],
                        [  3.63000000e-07,   1.36956935e-17,   1.58122800e+00]])

        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = ab_fragment.Fragment(filePath=ch4, fragmentType='A')
        if False:
                axis = [0,0,1]
                angle = math.pi/3
                rotationMatrix = xyz_core.rotation_matrix(axis, angle)
                center = np.array([ 0, 0, 0 ])
                frag.rotate(rotationMatrix, center)

        body = list(frag.bodies())[0]
        bodyFrame = True
        I = body.momentOfInertia(bodyFrame=bodyFrame)
        self.assertTrue(np.allclose(I, Iref))
        # https://github.com/pierrepo/principal_axes/blob/master/principal_axes.py
        #http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node67.html
        #https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
        
        def principal_axes(I):
            
            print("I ",I)
    
            # calculate principal moments of inertia
            e_values, e_vectors = np.linalg.eig(I)
            # warning eigen values are not necessary ordered!
            # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html
            print("(Unordered) eigen values:")
            print(e_values)
            print("(Unordered) eigen vectors:")
            print(e_vectors)
    
            #--------------------------------------------------------------------------
            # order eigen values (and eigen vectors)
            #
            # axis1 is the principal axis with the biggest eigen value (eval1)
            # axis2 is the principal axis with the second biggest eigen value (eval2)
            # axis3 is the principal axis with the smallest eigen value (eval3)
            #--------------------------------------------------------------------------
            order = np.argsort(e_values)
            eval3, eval2, eval1 = e_values[order]
            axis3, axis2, axis1 = e_vectors[:, order].transpose()
            print ("GOT AXES, ", axis3, axis2, axis1)
            print("Inertia axis are now ordered !")
            return (axis3, axis2, axis1)
        
        principal_axes(I)

        axis = [0,0,1]
        angle = math.pi/3
        rotationMatrix = xyz_core.rotation_matrix(axis, angle)
        center = np.array([ 0, 0, 0 ])
        frag.rotate(rotationMatrix, center)
        frag.translate([10.0,5.0, 3.0])
        body = list(frag.bodies())[0]

        I = body.momentOfInertia(bodyFrame=bodyFrame)
        principal_axes(I)
        
    def XtestMomentOfInertia2(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = ab_fragment.Fragment(filePath=ch4, fragmentType='A')
 
        body = list(frag.bodies())[0]
        com1 = body.momentOfInertia()
        
        axis = [1, 0, 1]
        angle = math.pi/3
        rotationMatrix = xyz_core.rotation_matrix(axis, angle)
        center = body.centreOfMass()
        frag.rotate(rotationMatrix, center)
        com2 = body.momentOfInertia()
        print(com1)
        print(com2)
        for x, y in zip(com1, com2):
            for a, b in zip(x, y):
                print(abs(a - b))
        print body.rigidType()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
