'''
Created on 4 Mar 2018

@author: jmht
'''
import os
import unittest

from paths import AMBUILD_DIR, BLOCKS_DIR, PARAMS_DIR
import fragment
import numpy
import util


# Not sure where best to do this
util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))

class Test(unittest.TestCase):

    def testMomentOfInertia(self):
        Iref = numpy.array([[  1.58122888e+00,  -2.49643156e-17,   3.63000000e-07],
                            [ -2.49643156e-17,   1.58122879e+00,   1.36956935e-17],
                            [  3.63000000e-07,   1.36956935e-17,   1.58122800e+00]])
        
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        frag = fragment.Fragment(filePath=ch4, fragmentType='A')
        for body in frag.bodies():
            I = body.momentInertia(dim=None, center=False, bodyspace=True)
            self.assertTrue(numpy.allclose(I, Iref))
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()