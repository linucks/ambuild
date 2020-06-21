"""
Created on 4 Mar 2018

@author: jmht
"""
import os
import unittest

import context
from context import ab_fragment
from context import xyz_util

BLOCKS_DIR = context.BLOCKS_DIR
PARAMS_DIR = context.PARAMS_DIR

xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))


class Test(unittest.TestCase):
    def testFragmentConfigStr(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        ftype = "A"
        f1 = ab_fragment.Fragment(filePath=ch4, fragmentType=ftype)
        cstr = ftype + "0000"
        self.assertEqual(cstr, f1.configStr)

        for eg in f1.endGroups():
            eg.bonded = True
        f1.update()
        cstr = ftype + "1111"
        self.assertEqual(cstr, f1.configStr)

    def testMultipleCapAtoms(self):
        h3car = os.path.join(BLOCKS_DIR, "H3.car")
        with self.assertRaises(RuntimeError):
            ab_fragment.Fragment(filePath=h3car, fragmentType="A")


if __name__ == "__main__":
    unittest.main()
