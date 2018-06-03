# python imports
import math
import os
import unittest

# external imports
import numpy as np

# our imports
import util
from paths import AMBUILD_DIR, PARAMS_DIR

class Test(unittest.TestCase):

    def setUp(self):
        #logging.basicConfig(level=logging.DEBUG)
        util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))
        return
    
    def testDistance(self):
        v1 = np.array([0, 0, 0])
        v2 = np.array([2.5, 2.5, 2.5])
        v3 = np.array([10, 10, 10])
        v4 = np.array([7.5, 7.5, 7.5])
        
        # Test without cell
        ref = [17.32050808, 8.66025404]
        result = util.distance([v1,v2], [v3,v4], dim=None)
        self.assertTrue(np.allclose(ref, result),
                         msg="Incorrect with no cell: {0}".format(result))
        
        # Test with cell
        cell = np.array([10.0, 10.0, 10.0])
        ref = [0, 8.66025404]
        result = util.distance([v1,v2], [v3,v4], dim=cell)
        self.assertTrue(np.allclose(ref, result),
                         msg="Incorrect with cell: {0}".format(result))
        
        # Test with only some conditions having PBC
        cell = np.array([10.0, 10.0, 10.0])
        ref = [10.0, 8.66025404]
        result = util.distance([v1,v2], [v3,v4], dim=cell, pbc=[False,True,True])
        self.assertTrue(np.allclose(ref,result ),
                         msg="Incorrect with partial cell: {0}".format(result))
        return
    
    def testLabel2Symbol(self):
        label = 'hc'
        self.assertEqual(util.label2symbol(label), 'H')
        
    
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
        bonds = util.calcBonds(coords,
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
        close = util.closeAtoms(coords,
                               symbols,
                               dim=dim)
    
        ref_close = [(0, 2), (0, 4), (0, 1), (0, 3), (1, 2), (1, 4), (1, 3), (2, 4), (2, 3), (3, 4)]
        # order of atoms doesn't 'matter
        self.assertEqual(set(close), set(ref_close))

    @unittest.skipUnless(util.PYTHONFLAVOUR == 2, "pkl file created with Python2")
    def testCellFromPickle(self):
        """Get old pickle file from Pierre"""
        pickleFile = os.path.join(AMBUILD_DIR,"tests","oldversion.pkl")
        mycell = util.cellFromPickle(pickleFile)
        self.assertEqual(len(mycell.blocks), 20, "Incorrect number of blocks: {0}".format(len(mycell.blocks)))
        
        # Just check we can build onto the cell as it demonstrates all the values are ok
        toGrow = 2
        grown = mycell.growBlocks(toGrow)
        self.assertEqual(toGrow,grown,"Failed to grow blocks after unpickling")
        os.unlink(mycell.logfile)
        os.unlink(mycell.logcsv)
        return

    @unittest.skipUnless(util.PYTHONFLAVOUR == 2, "pkl file created with Python2")
    def testDumpDLPOLY(self):
        """Get old pickle file from Pierre"""
        pickleFile = os.path.join(AMBUILD_DIR,"tests","oldversion.pkl")
        rigidBody = True
        skipDihedrals = True
        util.dumpDLPOLY(pickleFile, rigidBody=rigidBody, skipDihedrals=skipDihedrals)
        # For now just make sure we write something out...
        for fname in  ['CONFIG', 'CONTROL', 'FIELD']:
            self.assertTrue(os.path.isfile(fname))
            os.unlink(fname)
        return

    def testVectorAngle(self):
        """Test we can measure angles"""
        v1 = np.array([0, 0, 0])
        v2 = np.array([1, 0, 0])
        theta = util.vectorAngle(v1, v2)
        self.assertEqual(180,math.degrees(theta))
        return
    
    def testVecDiff(self):
        v1 = np.array([0, 0, 0])
        v2 = np.array([2.5, 2.5, 2.5])
        v3 = np.array([10, 10, 10])
        v4 = np.array([7.5, 7.5, 7.5])
        v5 = np.array([27.5, 27.5, 27.5])
        
        dim = np.array([10.0, 10.0, 10.0])
        
        # Test without cell
        ref = [[ -10.0, -10.0,  -10.0], [-5.0, - 5.0, -5.0]]
        result = util.vecDiff([v1,v2], [v3,v4], dim=dim, pbc=[False,False,False])
        self.assertTrue(np.allclose(ref, result),
                         msg="Incorrect with no cell: {0}".format(result))
         
        # Test with single vectors in cell
        ref = [[ 0.0, 0.0,  0.0]]
        result = util.vecDiff(v1, v3, dim=dim, pbc=[True,True,True])
        self.assertTrue(np.allclose(ref, result),
                         msg="Incorrect with cell single: {0}".format(result))
         
        # Test with cell
        ref = [[ -0.0, -0.0,  -0.0], [5.0, 5.0, 5.0]]
        result = util.vecDiff([v1,v2], [v3,v4], dim=dim, pbc=[True,True,True])
        self.assertTrue(np.allclose(ref, result),
                         msg="Incorrect with cell: {0}".format(result))
         
        # Test with only some conditions having PBC
        ref = [[ -10.0, -0.0,  -0.0], [-5.0, 5.0, 5.0]]
        result = util.vecDiff([v1,v2], [v3,v4], dim=dim, pbc=[False,True,True])
        self.assertTrue(np.allclose(ref, result),
                         msg="Incorrect with partial cell: {0}".format(result))
        
        # Test with only some conditions having PBC across multiple cells
        ref = [-25.0, 5.0, 5.0]
        result = util.vecDiff(v2, v5, dim=dim, pbc=[False,True,True])
        self.assertTrue(np.allclose(ref, result), msg="Incorrect with partial cell: {0}".format(result))
        return
    
    def testWrapCoord3_1(self):
        dim = np.array([10,20,30])
        c1 = np.array([101.0,202.0,303.0 ])
        center=True
        c1u, image = util.wrapCoord3(c1, dim, center=center)
        cref = np.array([-4.0,-8.0,-12.0])
        iref = np.array([10,10,10],dtype=np.int)
        self.assertTrue(np.allclose(c1u,cref))
        self.assertTrue(np.array_equal(image,iref))
        
    def testWrapCoord3_2(self):
        dim = np.array([10,20,30])
        c1 = np.array([-101.0,-202.0,-303.0 ])
        center=True
        c1u, image = util.wrapCoord3(c1, dim, center=center)
        cref = np.array([ 4.,8.,12.])
        iref = np.array([-11,-11,-11],dtype=np.int)
        self.assertTrue(np.allclose(c1u,cref))
        self.assertTrue(np.array_equal(image,iref))
        
    def testWrapCoord3_list(self):
        dim = np.array([10,20,30])
        c1 = np.array([[101.0,202.0,303.0 ], [-101.0,-202.0,-303.0 ]])
        center=True
        c1u, image = util.wrapCoord3(c1, dim, center=center)
        cref = np.array([[-4.0,-8.0,-12.0], [ 4.,8.,12.]])
        iref = np.array([[10,10,10], [-11,-11,-11]], dtype=np.int)
        self.assertTrue(np.allclose(c1u,cref))
        self.assertTrue(np.array_equal(image,iref))    

    def testUnWrapCoord3_1(self):
        cin = np.array([-4.0,-8.0,-12.0])
        idxin = np.array([10,10,10],dtype=np.int)
        dim = np.array([10,20,30])
        coord = util.unWrapCoord3(cin, idxin, dim, centered=True)
        self.assertTrue(np.allclose(coord,np.array([101.0,202.0,303.0 ])))
        
if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()