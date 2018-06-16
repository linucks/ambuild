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
    
    def testFragmentConfig(self):
        ch4 = os.path.join(BLOCKS_DIR, "ch4.car")
        f1 = fragment.Fragment(filePath=ch4, fragmentType='A')
        self.assertTrue(f1.fragmentType in fragment.configManager.configs)
        self.assertTrue(len(fragment.configManager.configs) == 1)
        
        f2 = f1.copy()
        self.assertTrue(len(fragment.configManager.configs) == 1)

        eg1 = f1.freeEndGroups()[0]
        eg2 = f2.freeEndGroups()[0]
        import buildingBlock
        bond = buildingBlock.Bond(eg1, eg2)
        bond.endGroup1.setBonded(bond)
        bond.endGroup2.setBonded(bond)
        self.assertTrue(len(fragment.configManager.configs[f1.fragmentType]) == 2)

        eg1 = f1.freeEndGroups()[-1]
        eg2 = f2.freeEndGroups()[0]
        bond = buildingBlock.Bond(eg1, eg2)
        bond.endGroup1.setBonded(bond)
        bond.endGroup2.setBonded(bond)
        self.assertTrue(len(fragment.configManager.configs[f1.fragmentType]) == 4)

        f1 = fragment.Fragment(filePath=ch4, fragmentType='B')
        self.assertTrue(len(fragment.configManager.configs) == 2)
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
