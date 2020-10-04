"""
Created on 4 Mar 2018

@author: jmht
"""
import os
import unittest

import context
from context import ab_fragment
from context import ab_util
from context import xyz_util

BLOCKS_DIR = context.BLOCKS_DIR
PARAMS_DIR = context.PARAMS_DIR

xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))


class Test(unittest.TestCase):

    def testFragmentType(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        ftype = 'A'
        f1 = ab_fragment.fragmentFactory(ftype, ch4)
        self.assertEqual(ftype, f1.fragmentType)

    def testDisallowedFragmentType(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        ftype = 'Fragment'
        with self.assertRaises(AssertionError):
            f1 = ab_fragment.fragmentFactory(ftype, ch4)

    def testPickleRoundtrip(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        ftype = 'A'
        f1 = ab_fragment.fragmentFactory(ftype, ch4)
        f2 = f1.copy()
        labels = ['A', 'A', 'A', 'A', 'A']
        f2._labels = labels
        fname = 'f1.pkl.gz'
        ab_util.pickleObj(f2, fname)
        f3 = ab_util.unpickleObj(fname)
        self.assertEqual(f3._labels, labels)
        os.unlink(fname)

    def testConfigStr(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        ftype = 'A'
        f1 = ab_fragment.fragmentFactory(ftype, ch4)
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
            ab_fragment.fragmentFactory('A', h3car)


if __name__ == "__main__":
    unittest.main()
