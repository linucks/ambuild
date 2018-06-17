'''
Created on 4 Mar 2018

@author: jmht
'''
import math
import os
import unittest

from paths import BLOCKS_DIR, PARAMS_DIR
import fragment
import numpy as np
import util

# Not sure where best to do this
util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))

class Test(unittest.TestCase):
    
    def testBodyConfig(self):
        ch4ca = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        f1 = fragment.Fragment(filePath=ch4ca, fragmentType='A')
        self.assertEqual(len(list(f1.bodies())), 3)
        self.assertEqual(list(f1.bodies())[-1].rigidType(), "AAA2")
    
    def testMomentOfInertia(self):
        Iref = np.array([[  1.58122888e+00,  -2.49643156e-17,   3.63000000e-07],
                        [ -2.49643156e-17,   1.58122879e+00,   1.36956935e-17],
                        [  3.63000000e-07,   1.36956935e-17,   1.58122800e+00]])

        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = fragment.Fragment(filePath=ch4, fragmentType='A')
        if False:
                axis = [0,0,1]
                angle = math.pi/3
                rotationMatrix = util.rotation_matrix(axis, angle)
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
        rotationMatrix = util.rotation_matrix(axis, angle)
        center = np.array([ 0, 0, 0 ])
        frag.rotate(rotationMatrix, center)
        frag.translate([10.0,5.0, 3.0])
        body = list(frag.bodies())[0]

        I = body.momentOfInertia(bodyFrame=bodyFrame)
        principal_axes(I)
        
    def testMomentOfInertia2(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = fragment.Fragment(filePath=ch4, fragmentType='A')
 
        body = list(frag.bodies())[0]
        com1 = body.momentOfInertia()
        
        axis = [1, 0, 1]
        angle = math.pi/3
        rotationMatrix = util.rotation_matrix(axis, angle)
        center = body.centreOfMass()
        frag.rotate(rotationMatrix, center)
        
        com2 = body.momentOfInertia()
        
        print com1
        print com2
        for x, y in zip(com1, com2):
            for a, b in zip(x, y):
                print abs(a - b)
        print body.rigidType()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
