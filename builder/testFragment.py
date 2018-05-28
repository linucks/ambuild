'''
Created on 4 Mar 2018

@author: jmht
'''
import os
import unittest

from paths import AMBUILD_DIR, BLOCKS_DIR, PARAMS_DIR
import fragment
import numpy as np
import util


# Not sure where best to do this
util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))

class Test(unittest.TestCase):

    def XtestMomentOfInertia(self):
        Iref = np.array([[  1.58122888e+00,  -2.49643156e-17,   3.63000000e-07],
                        [ -2.49643156e-17,   1.58122879e+00,   1.36956935e-17],
                        [  3.63000000e-07,   1.36956935e-17,   1.58122800e+00]])

        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = fragment.Fragment(filePath=ch4, fragmentType='A')
        body = list(frag.bodies())[0]

        coords = body.coords(dim=None, center=False, bodyspace=True)
        I = body.momentOfInertia(coords)
        self.assertTrue(np.allclose(I, Iref))
        # https://github.com/pierrepo/principal_axes/blob/master/principal_axes.py
        #http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node67.html
        #https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/

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
        print ("GOT ACXES, ", axis3, axis2, axis1)
        print("Inertia axis are now ordered !")

    def testMomentOfInertia2(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = fragment.Fragment(filePath=ch4, fragmentType='A')
        body = list(frag.bodies())[0]
        dim = np.array([10.0, 10.0, 10.0])
        coords, images = body.coords(dim=dim, center=True, bodyspace=True)
        I = body.momentOfInertia(coords)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
