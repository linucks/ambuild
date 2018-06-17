import logging
import math
import os
import unittest

import numpy

import context
AMBUILD_DIR = context.paths.AMBUILD_DIR
BLOCKS_DIR = context.paths.BLOCKS_DIR
PARAMS_DIR = context.paths.PARAMS_DIR
Bond = context.ab_bond.Bond
Block = context.buildingBlock.Block
from context import fragment
from context import xyz_util


logger = logging.getLogger(__name__)

class Test(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.INFO)

        self.ch4Xyz = os.path.join(BLOCKS_DIR, "ch4.xyz")
        self.ch4Car = os.path.join(BLOCKS_DIR, "ch4.car")
        self.cx4Car = os.path.join(BLOCKS_DIR, "cx4.car")
        self.ch4_1Car = os.path.join(BLOCKS_DIR, "ch4_1.car")
        # self.pafCar = os.path.join( BLOCKS_DIR, "PAF_bb_typed.car" )
        self.benzeneCar = os.path.join(BLOCKS_DIR, "benzene.car")
        self.benzene2Car = os.path.join(BLOCKS_DIR, "benzene2.car")
        self.ch4Ca2Car = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        
        xyz_util.setModuleBondLength(os.path.join(PARAMS_DIR, "bond_params.csv"))
        fragment.configManager.reset()
        return

    def catBlocks(self, blocks, filename):

        symbols = []
        coords = []
        for b in blocks:
            for i, c in enumerate(b.iterCoord()):
                symbols.append(b._symbol(i))
                coords.append(c)

        with open(filename, "w") as f:
            fpath = os.path.abspath(f.name)
            f.write("{}\n".format(len(coords)))
            f.write("id={}\n".format(str(id(self))))
            for i, c in enumerate(coords):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format(symbols[ i ], c[0], c[1], c[2]))

        logger.info("Wrote file: {0}".format(fpath))
        return

    def testBodies(self):
        b1 = Block(filePath=self.ch4Ca2Car, fragmentType='A')
        b2 = b1.copy()
        eg1 = b1.freeEndGroups()[0]
        eg2 = b2.freeEndGroups()[0]
        b1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()
        ref = [0, 0, 0, 0, 1, 2, 4, 4, 4, 4, 5, 6]
        self.assertEqual(b1._bodies, ref)

        return

    def testCH4(self):
        """Test the creation of a CH4 molecule"""
        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        endGroups = [ 0, 0, 0, 0 ]
        self.assertEqual(endGroups, [ e.blockEndGroupIdx for e in ch4.freeEndGroups() ])
        angleAtoms = [ 1, 2, 3, 4 ]
        self.assertEqual(angleAtoms, [ e.blockCapIdx for e in ch4.freeEndGroups() ])
        return

    def testCX4(self):
        """Test the creation of a CX4 molecule"""

        cx4_1 = Block(filePath=self.cx4Car, fragmentType='A')

        self.assertEqual([ 0, 0, 0, 0 ], [ e.blockEndGroupIdx for e in cx4_1.freeEndGroups() ])
        self.assertEqual([ 1, 2, 3, 4, ], [ e.blockCapIdx for e in cx4_1.freeEndGroups() ])

        cx4_2 = Block(filePath=self.cx4Car, fragmentType='A')

        eg1 = cx4_1.freeEndGroups()[0]
        eg2 = cx4_2.freeEndGroups()[0]


        cx4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()
        return

    def testCH4_Fragmentbond(self):
        """First pass"""

        ch4 = Block(filePath=self.ch4Car, fragmentType='A')

        self.assertEqual(len(ch4.fragments[0]._bonds), 4)
        return

    def testAnglesAndDihedrals(self):
        ch4_1 = Block(filePath=self.benzeneCar, fragmentType='A')

        ad = ch4_1.anglesAndDihedrals()
        #
        # UNCHECKED JUST HERE SO I CAN SPOT IF ANYTHING CHANGES
        ref = ([(0, 1, 2), (0, 1, 6), (0, 5, 4), (0, 5, 11), (1, 0, 5), (1, 0, 8), (1, 2, 3), (1, 2, 10), (2, 1, 6), (2, 3, 4), (2, 3, 7), (3, 2, 10), (3, 4, 5), (3, 4, 9), (4, 3, 7), (4, 5, 11), (5, 0, 8), (5, 4, 9)], [(0, 1, 2, 3), (0, 1, 2, 10), (0, 5, 4, 3), (0, 5, 4, 9), (1, 0, 5, 4), (1, 0, 5, 11), (1, 2, 3, 4), (1, 2, 3, 7), (2, 1, 0, 5), (2, 1, 0, 8), (2, 3, 4, 5), (2, 3, 4, 9), (3, 2, 1, 6), (3, 4, 5, 11), (4, 3, 2, 10), (4, 5, 0, 8), (5, 0, 1, 6), (5, 4, 3, 7), (6, 1, 0, 8), (6, 1, 2, 10), (7, 3, 2, 10), (7, 3, 4, 9), (8, 0, 5, 11), (9, 4, 5, 11)], [(0, 8, 1, 5), (1, 0, 2, 6), (2, 1, 10, 3), (3, 2, 4, 7), (4, 9, 3, 5), (5, 0, 11, 4)])

        self.assertEqual(ad, ref, "untested angles and dihedrals")
        return

    def testBond1(self):
        """First pass"""

        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]

        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()

        self.assertEqual([0, 0, 0, 4, 4, 4], [ eg.blockEndGroupIdx for eg in ch4_1.freeEndGroups() ])

        self.assertEqual(len(ch4_1._blockBonds), 1)
        self.assertEqual([ (0, 4) ], [ (b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx) for b in ch4_1._blockBonds ])

        # Check block Bonds
        self.assertEqual([ (0, 4) ], ch4_1.blockBonds())

        # Check all bonds
        ref_bonds = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7), (0, 4)]
        # order irrelevant
        self.assertEqual(ref_bonds, ch4_1.bonds())

        return

    def testBondSelf(self):
        """Silly test as bonds aren't feasible"""

        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')

        eg1 = ch4_1.freeEndGroups()[1]
        eg2 = ch4_2.freeEndGroups()[2]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_1.freeEndGroups()[-1]
        bond = Bond(eg1, eg2)
        bond.engage()
        return

    def testDeleteBondSimple(self):
        """Bfoo"""
        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        ch4_1.deleteBond(bond)
        return

    def testDeleteBondCircular(self):
        """Bfoo"""
        def egFromF(block, f1):
            # Need to find endGroups that match the fragments at either end
            for eg in block.freeEndGroups():
                if eg.fragment == f1:
                    return eg
            assert False

        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        f1 = ch4_1.fragments[0]
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')
        f2 = ch4_2.fragments[0]
        ch4_3 = Block(filePath=self.ch4Car, fragmentType='A')
        f3 = ch4_3.fragments[0]
        ch4_4 = Block(filePath=self.ch4Car, fragmentType='A')
        f4 = ch4_4.fragments[0]
        ch4_5 = Block(filePath=self.ch4Car, fragmentType='A')
        f5 = ch4_5.fragments[0]

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()

        eg1 = egFromF(ch4_1, f2)
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bondM = Bond(eg1, eg2)
        bondM.engage()

        eg1 = egFromF(ch4_1, f3)
        eg2 = ch4_4.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()

        eg1 = egFromF(ch4_1, f4)
        eg2 = ch4_5.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()

        # Now just bond into a loop
        eg1 = egFromF(ch4_1, f1)
        eg2 = egFromF(ch4_1, f5)
        bond = Bond(eg1, eg2)
        bond.engage()

        #ch4_1.writeCml("foo2.cml")
        self.assertFalse(ch4_1.deleteBond(bondM))
        #ch4_1.writeCml("foo3.cml")
        return

    def testDeleteBondSplit(self):
        """Create a central block with 4 attached blocks and the split one of the bonds so we get 2 blocks"""

        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_3 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_4 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_5 = Block(filePath=self.ch4Car, fragmentType='A')

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_4.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_5.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()

        #ch4_1.writeCml("foo1.cml")
        self.assertTrue(bool(ch4_1.deleteBond(bond)))
        #ch4_1.writeCml("foo2.cml")
        return
    
    def testDeleteFragment(self):
        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        f1 = ch4_1.fragments[0]
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_3 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_4 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_5 = Block(filePath=self.ch4Car, fragmentType='A')

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_4.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_5.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()
        
#         ch4_1.writeXyz("foo.xyz")
        blocks = ch4_1.deleteFragment(f1)
        self.assertEqual(len(blocks), 4, "Not enough blocks returned")
#         for i, b in enumerate(blocks):
#             b.writeXyz("foo_{0}.xyz".format(i))
        
        return

    def XtestAlignBlocks(self):
        """Test we can align two _blocks correctly"""

        blockS = Block(filePath=self.benzeneCar, fragmentType='A')
        block = blockS.copy()

        block.translateCentroid([ 3, 4 , 5 ])
        block.randomRotate()

        # Get the atoms that define things
        eg1 = blockS.freeEndGroups()[0]
        idxSatom = eg1.blockEndGroupIdx
        idxAatom = eg1.blockCapAtomIdx
        blockSEndGroup = blockS._coord(idxSatom)
        blockSangleAtom = blockS._coord(idxAatom)


        # idxAtom = 7
        # idxAatom2 = 1
        # blockAngleAtom = block._coord( idxAatom2 )
        eg2 = block.freeEndGroups()[0]
        blockAngleAtom = eg1.blockCapAtomIdx

        # we want to align along blockSangleAtom -> blockSEndGroup
        refVector = blockSEndGroup - blockSangleAtom

        # Position block so contact is at origin
        block.translate(-blockAngleAtom)
        block.alignAtoms(idxAatom2, idxAtom, refVector)

        # Check the relevant atoms are in the right place
        blockEndGroup = block._coord(idxAtom)
        blockAngleAtom = block._coord(idxAatom2)

        newVector = blockEndGroup - blockAngleAtom

        # Normalise two vectors so we can compare them
        newNorm = newVector / numpy.linalg.norm(newVector)
        refNorm = refVector / numpy.linalg.norm(refVector)

        # Slack tolerances - need to work out why...
        self.assertTrue(numpy.allclose(newNorm, refNorm),
                         msg="End Group incorrectly positioned: {0} | {1}".format(newNorm, refNorm))
        return

    def testAlignAtoms(self):
        block = Block(filePath=self.benzeneCar, fragmentType='A')

        # Check atoms are not aligned along axis
        c1Idx = 2
        c2Idx = 5
        c1 = block.coord(c1Idx)
        c2 = block.coord(c2Idx)

        # self.assertTrue( numpy.allclose( c1-c2 , [ 3.0559,  -0.36295,  0.07825], atol=1E-7  ), "before" )
        self.assertTrue(numpy.allclose(c1 - c2 , [ 3.0559, -0.36295, 0.07825]), "before")

        # Align along z-axis
        block.alignAtoms(c1Idx, c2Idx, [ 0, 0, 1 ])

        # check it worked
        c1 = block.coord(c1Idx)
        c2 = block.coord(c2Idx)
        z = numpy.array([  0.0, 0.0, -3.07837304 ])

        self.assertTrue(numpy.allclose(c1 - c2 , z), "after")
        return

    def testCentroid(self):
        """
        Test calculation of Center of Geometry
        """

        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        # We need to move so that all values are not the same or allclose will succeed with a single number
        ch4.translate(numpy.array([1.0,2.0,3.0]))
        cog = ch4.centroid()
        correct = numpy.array([  1.000000, 2.000000, 3.000000 ])
        self.assertTrue(numpy.allclose(correct, cog, rtol=1e-9, atol=1e-6),
                         msg="testCentroid incorrect: {0} -> {1}.".format(correct, cog))

        return

    def testCenterOfMass(self):
        """
        Test calculation of Center of Mass
        """
        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4.translate(numpy.array([1.0,2.0,3.0]))
        com = ch4.centerOfMass()
        correct = numpy.array([  1.000000, 2.000000, 3.000000 ])
        self.assertTrue(numpy.allclose(correct, com, rtol=1e-6, atol=1e-6),
                         msg="testCenterOfMass incorrect COM: {0}".format(com))
        return

    def testDihedrals(self):
        """foo"""

        ch4_1 = Block(filePath=self.benzeneCar, fragmentType='A')
        ch4_2 = ch4_1.copy()

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()

        # Check just across bonds
        ref = [ (1, 0, 11, 16),
                (1, 0, 11, 12),
                (5, 0, 11, 16),
                (5, 0, 11, 12) ]
        dihedrals = ch4_1.dihedrals(eg1.endGroupIdx(), eg2.endGroupIdx(), bondOnly=True)
        self.assertEqual(dihedrals, ref, "across bond: {} {}".format(ref, dihedrals))

        # Now all dihedrals
        ref = [(11, 0, 1, 2),
               (11, 0, 1, 6),
               (11, 0, 5, 10),
               (11, 0, 5, 4),
               (0, 11, 16, 21),
               (0, 11, 16, 15),
               (0, 11, 12, 17),
               (0, 11, 12, 13),
               (1, 0, 11, 16),
               (1, 0, 11, 12),
               (5, 0, 11, 16),
               (5, 0, 11, 12)]
        dihedrals = ch4_1.dihedrals(eg1.endGroupIdx(), eg2.endGroupIdx(), bondOnly=False)
        self.assertEqual(dihedrals, ref, "all dihedrals:\n{}\n{}".format(ref, dihedrals))

        return
    
    def testFreeEndGroups(self):
        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = ch4_1.copy()
        b1 = Block(filePath=self.benzeneCar, fragmentType='B')

        # create a chain of ch4 - c6h6 - ch4
        eg1 = b1.freeEndGroups()[0]
        eg2 = ch4_1.freeEndGroups()[0]
        b1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()

        eg1s = b1.freeEndGroups()
        self.assertEqual(len(eg1s), 4)
        eg = eg1s[0]
        ref_idxs = [0, 11, 2, 12]
        idxs = [eg.fragmentEndGroupIdx, eg.blockEndGroupIdx, eg.fragmentCapIdx, eg.blockCapIdx]
        self.assertEqual(idxs, ref_idxs)
        
        eg = eg1s[-1]
        ref_idxs = [3, 3, 7, 7]
        idxs = [eg.fragmentEndGroupIdx, eg.blockEndGroupIdx, eg.fragmentCapIdx, eg.blockCapIdx]
        self.assertEqual(idxs, ref_idxs)
        
        eg2s = ch4_2.freeEndGroups()
        self.assertEqual(len(eg2s), 4)
        eg = eg2s[0]
        ref_idxs = [0, 0, 1, 1]
        idxs = [eg.fragmentEndGroupIdx, eg.blockEndGroupIdx, eg.fragmentCapIdx, eg.blockCapIdx]
        self.assertEqual(idxs, ref_idxs)
        
        eg = eg2s[-1]
        ref_idxs = [0, 0, 4, 4]
        idxs = [eg.fragmentEndGroupIdx, eg.blockEndGroupIdx, eg.fragmentCapIdx, eg.blockCapIdx]       
        self.assertEqual(idxs, ref_idxs)
    
    def testMaxBond(self):
        carfile = os.path.join(BLOCKS_DIR, "DCX.car")
        f = fragment.Fragment(filePath=carfile, fragmentType='A')
        f.setMaxBond('A:CH', 1)
        block1 = Block(initFragment=f)
        
        block2 = Block(filePath=self.ch4Car, fragmentType='B')
        
        # Check we have 6 free endGroups at the start
        self.assertEqual(len(block1.freeEndGroups()),6)
        self.assertEqual(len(block2.freeEndGroups()),4)
        
        # Create a bond to the maxBonded type
        eg1 = block1.freeEndGroups(endGroupTypes='A:CH')[0]
        eg2 = block2.freeEndGroups()[0]
        bond = Bond(eg1, eg2)
        bond.engage()
        
        # Now see if the other three are blocked - we need to include the three from the second block
        self.assertEqual(len(block1.freeEndGroups()),5)
        return
    
    def testBondingFunction(self):
        carfile = os.path.join(BLOCKS_DIR, "benzene6.car")
        f = fragment.Fragment(filePath=carfile, fragmentType='A')
        def x(endGroup):
            fragment = endGroup.fragment
            egt = endGroup.type()
            for eg in fragment.endGroups():
                if egt == 'A:a':
                    if not eg.bonded and eg.type() in ['A:c', 'A:e']:
                        eg.blocked = True
                elif egt == 'A:b':
                    if not eg.bonded and eg.type() in ['A:d', 'A:f']:
                        eg.blocked = True
                elif egt == 'A:c':
                    if not eg.bonded and eg.type() in ['A:a', 'A:e']:
                        eg.blocked = True
                elif egt == 'A:d':
                    if not eg.bonded and eg.type() in ['A:b', 'A:f']:
                        eg.blocked = True
                elif egt == 'A:e':
                    if not eg.bonded and eg.type() in ['A:a', 'A:c']:
                        eg.blocked = True
                elif egt == 'A:f':
                    if not eg.bonded and eg.type() in ['A:b', 'A:d']:
                        eg.blocked = True
            return
        
        f.onbondFunction = x
        block1 = Block(initFragment=f)
        
        block2 = Block(filePath=self.ch4Car, fragmentType='B')
        
        # Check we have 6 free endGroups at the start
        self.assertEqual(len(block1.freeEndGroups()),6)
        self.assertEqual(len(block2.freeEndGroups()),4)
        
        # Create a bond to the first endGroup
        eg1 = block1.freeEndGroups(endGroupTypes='A:a')[0]
        eg2 = block2.freeEndGroups()[0]
        bond = Bond(eg1, eg2)
        bond.engage()
        
        # Now see if the other two are blocked - we need to include the three from the second block
        self.assertEqual(len(block1.freeEndGroups()),6)
        return

    def testMultiEndGroups(self):
        """Test we can move correctly"""
        # Try with no settings
        f = fragment.Fragment(filePath=self.ch4_1Car, fragmentType='A')
        m1 = Block(initFragment=f)
        m2 = m1.copy()
        eg1 = m1.freeEndGroups()[0]
        eg2 = m2.freeEndGroups()[0]
        m1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()
        self.assertEqual(6, len(m1.freeEndGroups()))

        # Try with specifying bond
        f = fragment.Fragment(filePath=self.ch4_1Car, fragmentType='A')
        m1 = Block(initFragment=f)
        f.setMaxBond('A:a', 1)
        m2 = m1.copy()
        eg1 = m1.freeEndGroups()[0]
        eg2 = m2.freeEndGroups()[0]
        m1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        bond.engage()
        self.assertEqual(4, len(m1.freeEndGroups()))
        return

    def testMove(self):
        """Test we can move correctly"""

        paf = Block(filePath=self.benzeneCar, fragmentType='A')
        m = paf.copy()
        m.translate(numpy.array([5, 5, 5]))
        c = m.centroid()
        paf.translateCentroid(c)
        p = paf.centroid()

        self.assertTrue(numpy.allclose(p, c, rtol=1e-9, atol=1e-9), "simple move")
        return

    def testPositionGrowBlock(self):

        blockS = Block(filePath=self.benzeneCar, fragmentType='A')

        growBlock = blockS.copy()

        growBlock.translateCentroid([ 3, 4, 5 ])
        growBlock.randomRotate()


        # Get the atoms that define things
        endGroup1 = blockS.freeEndGroups()[ 0 ]
        endGroup2 = growBlock.freeEndGroups()[ 1 ]

        # Get position to check
        newPos = blockS.newBondPosition(endGroup1, growBlock.symbol(endGroup2.blockEndGroupIdx))

        # Position the block
        blockS.positionGrowBlock(endGroup1, endGroup2)

        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock.coord(endGroup2.blockEndGroupIdx)
        self.assertTrue(numpy.allclose(newPos, endGroupCoord, rtol=1e-9, atol=1e-7),
                         msg="testCenterOfMass incorrect COM.")

        return

    def testPositionGrowBlock2(self):

        staticBlock = Block(filePath=self.benzeneCar, fragmentType='A')

        growBlock = staticBlock.copy()

        growBlock.translateCentroid([ 3, 4, 5 ])
        growBlock.randomRotate()

        # Get the atoms that define things
        endGroup1 = staticBlock.freeEndGroups()[ 0 ]
        endGroup2 = growBlock.freeEndGroups()[ 0 ]

        # Get position to check
        newPos = staticBlock.newBondPosition(endGroup1, growBlock.symbol(endGroup2.blockEndGroupIdx))

        # staticBlock._symbols.append( 'N' )
        # staticBlock._coords.append( newPos )
        # staticBlock.writeXyz("FOO.xyz")


        # Position the block
        # staticBlock.XXpositionGrowBlock( endGroup1, growBlock, endGroup2 )
        staticBlock.positionGrowBlock(endGroup1, endGroup2)

        # self.catBlocks( [staticBlock, growBlock ], "both.xyz")

        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock.coord(endGroup2.blockEndGroupIdx)
        self.assertTrue(numpy.allclose(newPos, endGroupCoord, rtol=1e-9, atol=1e-7),
                         msg="testCenterOfMass incorrect COM.")

        return

    def testPositionDihedral(self):

        staticBlock = Block(filePath=self.benzeneCar, fragmentType='A')

        growBlock = staticBlock.copy()

        growBlock.translateCentroid([ 3, 4, 5 ])
        growBlock.randomRotate()

        # Get the atoms that define things
        endGroup1 = staticBlock.freeEndGroups()[ 0 ]
        endGroup2 = growBlock.freeEndGroups()[ 0 ]

        # Get position to check
        staticBlock.newBondPosition(endGroup1,
                                    growBlock.symbol(endGroup2.blockEndGroupIdx))
        # Position the block
        staticBlock.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi / 2)

        # Hacky - just use one of the coords I checked manually
        hcheck = numpy.array([11.98409351860, 8.826721156800, -1.833703434310])
        endGroupCoord = growBlock.coord(11)
        self.assertTrue(numpy.allclose(hcheck, endGroupCoord, rtol=1e-9, atol=1e-7),
                         msg="testCenterOfMass incorrect COM.")

        # self.catBlocks( [staticBlock, growBlock ], "both2.xyz")
        return

    def testRadius(self):
        """
        Test calculation of the radius
        """

        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        r = ch4.blockRadius()
        self.assertAlmostEqual(r, 1.792806, 6, "Incorrect radius: {}".format(str(r)))
        return

    def XtestMaxAtomRadius(self):
        """
        Test calculation of the radius
        """

        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        r = ch4.maxAtomRadius()
        # jmht - check...- old was: 1.78900031214
        self.assertAlmostEqual(r, 0.70380574117, 7, "Incorrect radius: {}".format(str(r)))

    def testRotate(self):
        """
        Test the rotation
        """

        ch4 = Block(filePath=self.ch4Car, fragmentType='A')

        array1 = numpy.array([ -0.51336 , 0.889165, -0.363 ])
        self.assertTrue(numpy.array_equal(ch4.coord(4), array1),
                         msg="testRotate arrays before rotation incorrect.")

        axis = numpy.array([1, 2, 3])
        angle = 2
        ch4.rotate(axis, angle)

        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue(numpy.allclose(ch4.coord(4), array2, rtol=1e-9, atol=1e-8),
                         msg="testRotate arrays after rotation incorrect.")

        # Check rotation by 360
        axis = numpy.array([1, 2, 3])
        angle = numpy.pi * 2
        ch4.rotate(axis, angle)

        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue(numpy.allclose(ch4.coord(4), array2, rtol=1e-9, atol=1e-8),
                         msg="testRotate arrays after rotation incorrect.")

        return

    def testSplitFragment(self):
        """
        Test the rotation
        """

        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = ch4_1.copy()
        b1 = Block(filePath=self.benzeneCar, fragmentType='B')

        # create a chain of ch4 - c6h6 - ch4
        eg1 = b1.freeEndGroups()[0]
        eg2 = ch4_1.freeEndGroups()[0]
        b1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()

        eg1 = b1.freeEndGroups()[-1]
        eg2 = ch4_2.freeEndGroups()[0]
        b1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()

        # b1.writeCml("foo1.cml")
        # b1.writeXyz("foo.xyz")

        coords, symbols, bonds = b1.dataByFragment('A')
        cmlFilename = "test.cml"
        xyz_util.writeCml(cmlFilename,
                          coords,
                          symbols,
                          bonds=bonds,
                          prettyPrint=True)

        with open(cmlFilename) as f:
            test = f.readlines()

        with open(os.path.join(AMBUILD_DIR, 'tests', 'test_data', 'testSplitFragment.cml')) as f:
            ref = f.readlines()

        self.assertEqual(test, ref, "cml compare")
        os.unlink(cmlFilename)

        return

    def testWriteCml(self):
        """foo"""
        ch4_1 = Block(filePath=self.benzeneCar, fragmentType='A')
        ch4_2 = ch4_1.copy()

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        bond.engage()

        # Write out the cml and see if it matches what we've saved
        fname = "test.cml"
        ch4_1.writeCml(fname)
        with open(fname) as f:
            test = f.readlines()

        with open(os.path.join(AMBUILD_DIR, "tests", 'test_data', "benzeneBond.cml")) as f:
            ref = f.readlines()
        self.assertEqual(test, ref, "cml compare")
        os.unlink(fname)
        return
if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()

