'''
Created on 4 Mar 2018

@author: jmht
'''
import os
import unittest

import context
BLOCKS_DIR = context.paths.BLOCKS_DIR
PARAMS_DIR = context.paths.PARAMS_DIR
from context import ab_bond
from context import ab_fragment
from context import xyz_util


# Not sure where best to do this
xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))

class Test(unittest.TestCase):

    def setUp(self):
        ab_fragment.configManager.reset()
        return
    
    def testFragmentConfig(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        f1 = ab_fragment.Fragment(filePath=ch4, fragmentType='A')
        self.assertTrue(f1.fragmentType in ab_fragment.configManager.configs)
        self.assertTrue(len(ab_fragment.configManager.configs) == 1)
        
        f2 = f1.copy()
        self.assertTrue(len(ab_fragment.configManager.configs) == 1)

        eg1 = f1.freeEndGroups()[0]
        eg2 = f2.freeEndGroups()[0]
        bond = ab_bond.Bond(eg1, eg2)
        bond.endGroup1.setBonded(bond)
        bond.endGroup2.setBonded(bond)
        self.assertTrue(len(ab_fragment.configManager.configs[f1.fragmentType]) == 2)

        eg1 = f1.freeEndGroups()[-1]
        eg2 = f2.freeEndGroups()[0]
        bond = ab_bond.Bond(eg1, eg2)
        bond.endGroup1.setBonded(bond)
        bond.endGroup2.setBonded(bond)
        self.assertTrue(len(ab_fragment.configManager.configs[f1.fragmentType]) == 4)

        f1 = ab_fragment.Fragment(filePath=ch4, fragmentType='B')
        self.assertTrue(len(ab_fragment.configManager.configs) == 2)
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
