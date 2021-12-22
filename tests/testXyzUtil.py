# python imports
import os
import unittest

# external imports
import numpy as np

from context import xyz_core
from context import xyz_util
from context import PARAMS_DIR


class Test(unittest.TestCase):
    def setUp(self):
        # logging.basicConfig(level=logging.DEBUG)
        xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))
        return

    def testLabel2Symbol(self):
        label = "hc"
        self.assertEqual(xyz_util.label2symbol(label), "H")

    def testCalcBonds(self):
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.089],
                [1.026719, 0.0, -0.363],
                [-0.51336, -0.889165, -0.363],
                [-0.51336, 0.889165, -0.363],
            ]
        )
        symbols = ["c", "h", "h", "h", "h"]
        dim = None
        maxAtomRadius = 0.70380574116999994
        bondMargin = 0.2
        boxMargin = 1.0
        bonds = xyz_util.calcBonds(
            coords,
            symbols,
            dim=dim,
            maxAtomRadius=maxAtomRadius,
            bondMargin=bondMargin,
            boxMargin=boxMargin,
        )

        ref_bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
        self.assertEqual(bonds, ref_bonds)

    def testCloseAtoms(self):
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.089],
                [1.026719, 0.0, -0.363],
                [-0.51336, -0.889165, -0.363],
                [-0.51336, 0.889165, -0.363],
            ]
        )
        symbols = ["c", "h", "h", "h", "h"]
        dim = None
        close = xyz_util.closeAtoms(coords, symbols, dim=dim)

        ref_close = [
            (0, 2),
            (0, 4),
            (0, 1),
            (0, 3),
            (1, 2),
            (1, 4),
            (1, 3),
            (2, 4),
            (2, 3),
            (3, 4),
        ]
        # order of atoms doesn't 'matter
        self.assertEqual(set(close), set(ref_close))

    def testTrilateration(self):
        """
        Test Trilateration to determine atom position

                 3[0,1,0]
        2[-1,0,0] 1[0,0,0] 4[1,0,0]
                 t[0,-1,0]
        """

        p1 = np.array([0.0, 0.0, 0.0])
        p2 = np.array([-1.0, 0.0, 0.0])
        p3 = np.array([0.0, 1.0, 0.0])
        p4 = np.array([1.0, 0.0, 0.0])
        points = [p1, p2, p3, p4]

        t = np.array([0.0, -1.0, 0.0])
        p1t = xyz_core.distance(t, p1)
        p2t = xyz_core.distance(t, p2)
        p3t = xyz_core.distance(t, p3)
        p4t = xyz_core.distance(t, p4)
        distances = [p1t, p2t, p3t, p4t]

        position = xyz_util.trilaterate3D(distances, points)
        self.assertTrue(np.allclose(position, t))

        #
        # Real example with Benzene atoms
        #
        p1 = np.array([12.5632, 6.03464, -1.74782])
        p2 = np.array([13.86152, 4.39379, 0.51081])
        p3 = np.array([14.74447, 6.33331, -2.51304])
        p4 = np.array([11.74412, 4.14417, 1.20077])
        points = [p1, p2, p3, p4]

        t = np.array([12.12731, 6.58007, -2.50777])
        p1t = xyz_core.distance(t, p1)
        p2t = xyz_core.distance(t, p2)
        p3t = xyz_core.distance(t, p3)
        p4t = xyz_core.distance(t, p4)
        distances = [p1t, p2t, p3t, p4t]

        position = xyz_util.trilaterate3D(distances, points)
        self.assertTrue(np.allclose(position, t))


if __name__ == "__main__":
    """
    Run the unit tests
    """
    unittest.main()
