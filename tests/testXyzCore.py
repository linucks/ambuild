# python imports
import math
import unittest

# external imports
import numpy as np

# our imports
from context import xyz_core


class Test(unittest.TestCase):
    def testDistance(self):
        v1 = np.array([0, 0, 0])
        v2 = np.array([2.5, 2.5, 2.5])
        v3 = np.array([10, 10, 10])
        v4 = np.array([7.5, 7.5, 7.5])

        # Test without cell
        ref = [17.32050808, 8.66025404]
        result = xyz_core.distance([v1, v2], [v3, v4], dim=None)
        self.assertTrue(
            np.allclose(ref, result), msg="Incorrect with no cell: {0}".format(result)
        )

        # Test with cell
        cell = np.array([10.0, 10.0, 10.0])
        ref = [0, 8.66025404]
        result = xyz_core.distance([v1, v2], [v3, v4], dim=cell)
        self.assertTrue(
            np.allclose(ref, result), msg="Incorrect with cell: {0}".format(result)
        )

        # Test with only some conditions having PBC
        cell = np.array([10.0, 10.0, 10.0])
        ref = [10.0, 8.66025404]
        result = xyz_core.distance(
            [v1, v2], [v3, v4], dim=cell, pbc=[False, True, True]
        )
        self.assertTrue(
            np.allclose(ref, result),
            msg="Incorrect with partial cell: {0}".format(result),
        )
        return

    def testRigidRotateSelf(self):
        p1 = np.array(
            [
                [7.30084273e-02, 0.00000000e00, -6.70417149e-08],
                [-2.89991573e-01, 0.00000000e00, -1.02671907e00],
                [-2.89991573e-01, -8.89165000e-01, 5.13359933e-01],
                [-2.89991573e-01, 8.89165000e-01, 5.13359933e-01],
            ]
        )
        ref = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        M = xyz_core.rigid_rotate(p1, p1)
        self.assertTrue(np.allclose(ref, M, rtol=0.0001))

        # Rotate and make sure all overlap
        p_1 = np.dot(p1, M.T)
        self.assertTrue(np.allclose(p1, p_1, rtol=0.0001))
        return

    def testRigidRotate(self):
        """Both coordinates on centre of mass"""
        p1 = np.array(
            [
                [7.30084273e-02, 0.00000000e00, -6.70417149e-08],
                [-2.89991573e-01, 0.00000000e00, -1.02671907e00],
                [-2.89991573e-01, -8.89165000e-01, 5.13359933e-01],
                [-2.89991573e-01, 8.89165000e-01, 5.13359933e-01],
            ]
        )

        p2 = np.array(
            [
                [-7.30084273e-02, 0.00000000e00, 6.70417145e-08],
                [2.89991573e-01, 0.00000000e00, 1.02671907e00],
                [2.89991573e-01, -8.89165000e-01, -5.13359933e-01],
                [2.89991573e-01, 8.89165000e-01, -5.13359933e-01],
            ]
        )
        m2 = np.array([12.0107, 1.00794, 1.00794, 1.00794])

        ref = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]
        M = xyz_core.rigid_rotate(p1, p2)
        self.assertTrue(np.allclose(ref, M, rtol=0.0001))

        # Rotate p1 to p2 and make sure all overlap
        p_2 = np.dot(p1, M.T)
        self.assertTrue(np.allclose(p2, p_2, rtol=0.0001))
        return

    def testRigidRotate2(self):
        """Both coordinates on centre of mass"""
        p1 = np.array(
            [
                [7.30084273e-02, 0.00000000e00, -6.70417149e-08],
                [-2.89991573e-01, 0.00000000e00, -1.02671907e00],
                [-2.89991573e-01, -8.89165000e-01, 5.13359933e-01],
                [-2.89991573e-01, 8.89165000e-01, 5.13359933e-01],
            ]
        )
        M1 = xyz_core.rotation_matrix(np.array([0, 1, 0]), np.pi / 3)
        p2 = np.dot(p1, M1.T)
        M2 = xyz_core.rigid_rotate(p1, p2)
        # Rotate p1 to p2 and make sure all overlap
        p_2 = np.dot(p1, M2.T)
        self.assertTrue(np.allclose(p2, p_2, rtol=0.0001))
        return

    def testRotMat2Quat(self):
        M = np.matrix([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        Q = xyz_core.quaternion_from_matrix(M)
        ref = np.array([1.0, 0.0, 0.0, 0.0])
        self.assertTrue(np.allclose(Q, ref, rtol=0.0001))
        return

    def testWrapCoord3_1(self):
        dim = np.array([10, 20, 30])
        c1 = np.array([101.0, 202.0, 303.0])
        center = True
        c1u, image = xyz_core.wrapCoord3(c1, dim, center=center)
        cref = np.array([-4.0, -8.0, -12.0])
        iref = np.array([10, 10, 10], dtype=np.int)
        self.assertTrue(np.allclose(c1u, cref))
        self.assertTrue(np.array_equal(image, iref))

    def testWrapCoord3_2(self):
        dim = np.array([10, 20, 30])
        c1 = np.array([-101.0, -202.0, -303.0])
        center = True
        c1u, image = xyz_core.wrapCoord3(c1, dim, center=center)
        cref = np.array([4.0, 8.0, 12.0])
        iref = np.array([-11, -11, -11], dtype=np.int)
        self.assertTrue(np.allclose(c1u, cref))
        self.assertTrue(np.array_equal(image, iref))

    def testWrapCoord3_list(self):
        dim = np.array([10, 20, 30])
        c1 = np.array([[101.0, 202.0, 303.0], [-101.0, -202.0, -303.0]])
        center = True
        c1u, image = xyz_core.wrapCoord3(c1, dim, center=center)
        cref = np.array([[-4.0, -8.0, -12.0], [4.0, 8.0, 12.0]])
        iref = np.array([[10, 10, 10], [-11, -11, -11]], dtype=np.int)
        self.assertTrue(np.allclose(c1u, cref))
        self.assertTrue(np.array_equal(image, iref))

    def testUnWrapCoord3_1(self):
        cin = np.array([-4.0, -8.0, -12.0])
        idxin = np.array([10, 10, 10], dtype=np.int)
        dim = np.array([10, 20, 30])
        coord = xyz_core.unWrapCoord3(cin, idxin, dim, centered=True)
        self.assertTrue(np.allclose(coord, np.array([101.0, 202.0, 303.0])))

    def testVectorAngle(self):
        """Test we can measure angles"""
        v1 = np.array([0, 0, 0])
        v2 = np.array([1, 0, 0])
        theta = xyz_core.vectorAngle(v1, v2)
        self.assertEqual(180, math.degrees(theta))
        return

    def testVecDiff(self):
        v1 = np.array([0, 0, 0])
        v2 = np.array([2.5, 2.5, 2.5])
        v3 = np.array([10, 10, 10])
        v4 = np.array([7.5, 7.5, 7.5])
        v5 = np.array([27.5, 27.5, 27.5])

        dim = np.array([10.0, 10.0, 10.0])

        # Test without cell
        ref = [[-10.0, -10.0, -10.0], [-5.0, -5.0, -5.0]]
        result = xyz_core.vecDiff(
            [v1, v2], [v3, v4], dim=dim, pbc=[False, False, False]
        )
        self.assertTrue(
            np.allclose(ref, result), msg="Incorrect with no cell: {0}".format(result)
        )

        # Test with single vectors in cell
        ref = [[0.0, 0.0, 0.0]]
        result = xyz_core.vecDiff(v1, v3, dim=dim, pbc=[True, True, True])
        self.assertTrue(
            np.allclose(ref, result),
            msg="Incorrect with cell single: {0}".format(result),
        )

        # Test with cell
        ref = [[-0.0, -0.0, -0.0], [5.0, 5.0, 5.0]]
        result = xyz_core.vecDiff([v1, v2], [v3, v4], dim=dim, pbc=[True, True, True])
        self.assertTrue(
            np.allclose(ref, result), msg="Incorrect with cell: {0}".format(result)
        )

        # Test with only some conditions having PBC
        ref = [[-10.0, -0.0, -0.0], [-5.0, 5.0, 5.0]]
        result = xyz_core.vecDiff([v1, v2], [v3, v4], dim=dim, pbc=[False, True, True])
        self.assertTrue(
            np.allclose(ref, result),
            msg="Incorrect with partial cell: {0}".format(result),
        )

        # Test with only some conditions having PBC across multiple cells
        ref = [-25.0, 5.0, 5.0]
        result = xyz_core.vecDiff(v2, v5, dim=dim, pbc=[False, True, True])
        self.assertTrue(
            np.allclose(ref, result),
            msg="Incorrect with partial cell: {0}".format(result),
        )
        return


if __name__ == "__main__":
    """
    Run the unit tests
    """
    unittest.main()
