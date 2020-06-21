# python imports
import os
import unittest

# our imports
from context import ab_util
from context import PARAMS_DIR, TESTDATA_DIR


class Test(unittest.TestCase):
    @unittest.skipUnless(ab_util.PYTHONFLAVOUR == 2, "pkl file created with Python2")
    def testCellFromPickle(self):
        """Get old pickle file from Pierre"""
        pickleFile = os.path.join(TESTDATA_DIR, "oldversion.pkl")
        mycell = ab_util.cellFromPickle(pickleFile)
        self.assertEqual(
            len(mycell.blocks),
            20,
            "Incorrect number of blocks: {0}".format(len(mycell.blocks)),
        )
        # Just check we can build onto the cell as it demonstrates all the values are ok
        toGrow = 2
        grown = mycell.growBlocks(toGrow)
        self.assertEqual(toGrow, grown, "Failed to grow blocks after unpickling")
        os.unlink(mycell.logfile)
        os.unlink(mycell.logcsv)
        return

    @unittest.skipUnless(ab_util.PYTHONFLAVOUR == 2, "pkl file created with Python2")
    def testCellFromPickleGzip(self):
        pickleFile = os.path.join(TESTDATA_DIR, "new.pkl.gz")
        mycell = ab_util.cellFromPickle(pickleFile)
        self.assertEqual(
            len(mycell.blocks),
            3,
            "Incorrect number of blocks: {0}".format(len(mycell.blocks)),
        )
        # Just check we can build onto the cell as it demonstrates all the values are ok
        toGrow = 2
        grown = mycell.growBlocks(toGrow)
        self.assertEqual(toGrow, grown, "Failed to grow blocks after unpickling")
        os.unlink(mycell.logfile)
        os.unlink(mycell.logcsv)
        return

    @unittest.skipUnless(ab_util.PYTHONFLAVOUR == 2, "pkl file created with Python2")
    def testDumpDLPOLY(self):
        """Get old pickle file from Pierre"""
        pickleFile = os.path.join(TESTDATA_DIR, "oldversion.pkl")
        rigidBody = True
        skipDihedrals = True
        ab_util.dumpDLPOLY(
            pickleFile,
            rigidBody=rigidBody,
            skipDihedrals=skipDihedrals,
            paramsDir=PARAMS_DIR,
        )
        # For now just make sure we write something out...
        for fname in ["CONFIG", "CONTROL", "FIELD"]:
            self.assertTrue(os.path.isfile(fname))
            os.unlink(fname)
        return


if __name__ == "__main__":
    """
    Run the unit tests
    """
    unittest.main()
