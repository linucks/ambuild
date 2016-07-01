# python imports
import math
import os
import unittest

# external imports
import numpy

# our imports
import util
from paths import AMBUILD_DIR

class Test(unittest.TestCase):

    def testDistance(self):
        v1 = numpy.array([0, 0, 0])
        v2 = numpy.array([2.5, 2.5, 2.5])
        v3 = numpy.array([10, 10, 10])
        v4 = numpy.array([7.5, 7.5, 7.5])
        
        # Test without cell
        ref = [17.32050808, 8.66025404]
        result = util.distance([v1,v2], [v3,v4], dim=None)
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with no cell: {0}".format(result))
        
        # Test with cell
        cell = numpy.array([10.0, 10.0, 10.0])
        ref = [0, 8.66025404]
        result = util.distance([v1,v2], [v3,v4], dim=cell)
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with cell: {0}".format(result))
        
        # Test with only some conditions having PBC
        cell = numpy.array([10.0, 10.0, 10.0])
        ref = [10.0, 8.66025404]
        result = util.distance([v1,v2], [v3,v4], dim=cell, pbc=[False,True,True])
        self.assertTrue(numpy.allclose(ref,result ),
                         msg="Incorrect with partial cell: {0}".format(result))
        return
    
    def testCellFromPickle(self):
        """Get old pick file from Pierre"""
        pickleFile = os.path.join(AMBUILD_DIR,"tests","oldversion.pkl")
        mycell = util.cellFromPickle(pickleFile)
        self.assertEqual(len(mycell.blocks), 20, "Incorrect number of blocks: {0}".format(len(mycell.blocks)))
        
        # Just check we can build onto the cell as it demonstrates all the values are ok
        mycell.growBlocks(2)
        return

    def testVectorAngle(self):
        """Test we can measure angles"""
        v1 = numpy.array([0, 0, 0])
        v2 = numpy.array([1, 0, 0])
        theta = util.vectorAngle(v1, v2)
        self.assertEqual(180,math.degrees(theta))
        return
    
    def testVecDiff(self):
        v1 = numpy.array([0, 0, 0])
        v2 = numpy.array([2.5, 2.5, 2.5])
        v3 = numpy.array([10, 10, 10])
        v4 = numpy.array([7.5, 7.5, 7.5])
        v5 = numpy.array([27.5, 27.5, 27.5])
        
        dim = numpy.array([10.0, 10.0, 10.0])
        
        # Test without cell
        ref = [[ -10.0, -10.0,  -10.0], [-5.0, - 5.0, -5.0]]
        result = util.vecDiff([v1,v2], [v3,v4], dim=dim, pbc=[False,False,False])
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with no cell: {0}".format(result))
         
        # Test with single vectors in cell
        ref = [[ 0.0, 0.0,  0.0]]
        result = util.vecDiff(v1, v3, dim=dim, pbc=[True,True,True])
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with cell single: {0}".format(result))
         
        # Test with cell
        ref = [[ -0.0, -0.0,  -0.0], [5.0, 5.0, 5.0]]
        result = util.vecDiff([v1,v2], [v3,v4], dim=dim, pbc=[True,True,True])
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with cell: {0}".format(result))
         
        # Test with only some conditions having PBC
        ref = [[ -10.0, -0.0,  -0.0], [-5.0, 5.0, 5.0]]
        result = util.vecDiff([v1,v2], [v3,v4], dim=dim, pbc=[False,True,True])
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with partial cell: {0}".format(result))
        
        # Test with only some conditions having PBC across multiple cells
        ref = [-25.0, 5.0, 5.0]
        result = util.vecDiff(v2, v5, dim=dim, pbc=[False,True,True])
        self.assertTrue(numpy.allclose(ref, result),
                         msg="Incorrect with partial cell: {0}".format(result))
        
        
        return

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()