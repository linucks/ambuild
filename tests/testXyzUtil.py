# python imports
import os
import unittest

# external imports
import numpy as np

import context
PARAMS_DIR = context.ab_paths.PARAMS_DIR
from context import xyz_util

class Test(unittest.TestCase):
    
    def setUp(self):
        #logging.basicConfig(level=logging.DEBUG)
        xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))
        return
    
    def testLabel2Symbol(self):
        label = 'hc'
        self.assertEqual(xyz_util.label2symbol(label), 'H')
        
    
    def testCalcBonds(self):
        coords = np.array([[ 0.      ,  0.      ,  0.      ],
                           [ 0.      ,  0.      ,  1.089   ],
                           [ 1.026719,  0.      , -0.363   ],
                           [-0.51336 , -0.889165, -0.363   ],
                           [-0.51336 ,  0.889165, -0.363   ]])
        symbols = ['c', 'h', 'h', 'h', 'h']
        dim = None
        maxAtomRadius = 0.70380574116999994
        bondMargin = 0.2
        boxMargin = 1.0
        bonds = xyz_util.calcBonds(coords,
                               symbols,
                               dim=dim,
                               maxAtomRadius=maxAtomRadius,
                               bondMargin=bondMargin,
                               boxMargin=boxMargin)
    
        ref_bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
        self.assertEqual(bonds, ref_bonds)

    def testCloseAtoms(self):
        coords = np.array([[ 0.      ,  0.      ,  0.      ],
                           [ 0.      ,  0.      ,  1.089   ],
                           [ 1.026719,  0.      , -0.363   ],
                           [-0.51336 , -0.889165, -0.363   ],
                           [-0.51336 ,  0.889165, -0.363   ]])
        symbols = ['c', 'h', 'h', 'h', 'h']
        dim = None
        close = xyz_util.closeAtoms(coords,
                               symbols,
                               dim=dim)
    
        ref_close = [(0, 2), (0, 4), (0, 1), (0, 3), (1, 2), (1, 4), (1, 3), (2, 4), (2, 3), (3, 4)]
        # order of atoms doesn't 'matter
        self.assertEqual(set(close), set(ref_close))


if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()