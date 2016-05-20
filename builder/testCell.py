"""The tests for the cell class"""

import cPickle
import math
import os
import unittest

import numpy

from cell import Cell
import buildingBlock
import opt
from paths import AMBUILD_DIR, BLOCKS_DIR
import util

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.cx4Car = os.path.join(BLOCKS_DIR, "cx4.car")
        cls.ch4Car = os.path.join(BLOCKS_DIR, "ch4.car")
        cls.ch4_1Car = os.path.join(BLOCKS_DIR, "ch4_1.car")
        cls.capLinker = os.path.join(BLOCKS_DIR, "cap_linker.car")
        cls.benzeneCar = os.path.join(BLOCKS_DIR, "benzene.car")
        cls.benzene2Car = os.path.join(BLOCKS_DIR, "benzene2.car")
        cls.pafCar = os.path.join(BLOCKS_DIR, "PAF_bb_typed.car")
        cls.dcxCar = os.path.join(BLOCKS_DIR, "DCX.car")
        cls.ch4Ca2Car = os.path.join(BLOCKS_DIR, "ch4Ca2.car")
        cls.amineCar = os.path.join(BLOCKS_DIR, "amine_typed.car")
        cls.triquinCar = os.path.join(BLOCKS_DIR, "triquin_typed.car")
        cls.graphiteCar = os.path.join(BLOCKS_DIR, "2_graphite_cont.car")

        return

    def createTestCell(self):
        # Cell dimensions need to be: L > 2*(r_cut+r_buff) and L < 3*(r_cut+r_buff)
        # From code looks like default r_buff is 0.4 and our default r_cut is 5.0
        boxDim = [20, 20, 20]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.benzene2Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')
        mycell.seed(5, fragmentType='A')
        mycell.growBlocks(8)
        return mycell
    
    def clashes(self, mycell, minDist=1.0, pbc=[True, True, True]):
        """Return True of if any atoms are < minDist apart"""
        # coords=[ c for b in cell.blocks.itervalues() for i, c in enumerate(b.iterCoord())]
        symbols = []
        coords = []
        for b in mycell.blocks.itervalues():
            for i, c in enumerate(b.iterCoord()):
                coords.append(c)
                symbols.append(b.symbol(i))
        dim = numpy.array([mycell.dim[0], mycell.dim[1], mycell.dim[2]])
        close = util.closeAtoms(coords, symbols, dim=dim, boxMargin=1.0)
        v1 = []
        v2 = []
        for idxAtom1, idxAtom2 in close:
            v1.append(coords[idxAtom1])
            v2.append(coords[idxAtom2])
        distances = util.distance(v1, v2, dim=dim, pbc=pbc)
        return any(map(lambda x: x < minDist, distances))

    def testCX4(self):
        """First pass"""

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(fragmentType='A', filename=self.cx4Car)
        mycell.addBondType('A:a-A:a')
        mycell.seed(1)
        mycell.growBlocks(1)

        return

#     def XXXtimeCheck(self):
#         """NOT A TEST JUST CODE TO TIME CHECKMOVE"""
#
#         def cellFromPickle(pickleFile):
#             with open(pickleFile) as f:
#                 myCell=cPickle.load(f)
#             return myCell
#
#         mycell = cellFromPickle("step_1.pkl")
#         # Get the last block
#         idxb = mycell.blocks.keys()[-1]
#         b = mycell.blocks[ idxb ]
#         be = b.freeEndGroups()[0]
#
#         i =  mycell.getInitBlock()
#         ie = i.freeEndGroups()[0]
#
#         def run():
#             global b, be, i, ie
#             for _ in xrange( 1000):
#                 b.positionGrowBlock( be, ie, dihedral=None )
#                 blockId = mycell.addBlock( i )
#                 mycell.checkMove( blockId )
#                 mycell.delBlock(blockId)
#
#         cProfile.run('run()','restats')
#         return

    def testAmbody(self):

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.ch4Ca2Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')
        added = mycell.seed(3)
        self.assertEqual(added, 3)

        return

    def testBond(self):
        """First pass"""

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)

        fragmentType = 'A'
        mycell.libraryAddFragment(fragmentType=fragmentType, filename=self.benzeneCar)
        mycell.addBondType('A:a-A:a')

        centralBlock = mycell.getInitBlock(fragmentType=fragmentType)

        block1 = mycell.getInitBlock(fragmentType=fragmentType)
        block2 = mycell.getInitBlock(fragmentType=fragmentType)

        centralBlock.positionGrowBlock(centralBlock.freeEndGroups()[ 0 ], block1.freeEndGroups()[ 0 ])
        centralBlock.positionGrowBlock(centralBlock.freeEndGroups()[ 1 ], block2.freeEndGroups()[ 0 ])

        # Now add growBlock to the cell so we can check for clashes
        mycell.addBlock(block1)
        mycell.addBlock(block2)

        # Now add the central block - it will have 2 bonds to make
        centralId = mycell.addBlock(centralBlock)

        self.assertTrue(mycell.checkMove(centralId), "checkMove failed!")
        self.assertEqual(mycell.processBonds(), 2)

        return

    def testBlockTypes(self):
        """Test we can add a block correctly"""


        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(fragmentType='A', filename=self.ch4Car)
        mycell.libraryAddFragment(fragmentType='B', filename=self.ch4Car)
        mycell.libraryAddFragment(fragmentType='C', filename=self.ch4Car)
        mycell.libraryAddFragment(fragmentType='D', filename=self.ch4Car)

        mycell.addBondType('A:a-B:a')
        mycell.addBondType('A:a-C:a')
        mycell.addBondType('A:a-D:a')

        mycell.seed(1, fragmentType='A')
        toGrow = 3
        grew = mycell.growBlocks(toGrow, endGroupType=None, maxTries=5)
        self.assertEqual(toGrow, grew, "testBlockTypes not ok")
        return

    def XtestCapBlocks(self):
        mycell = self.createTestCell()
        mycell.capBlocks(fragmentType='A', filename=self.capLinker)
        # mycell.dump()
        mycell.blocks[ mycell.blocks.keys()[0] ]
        return

    def testCellIO(self):
        """Check we can write out and then read in a cell
        """

        # Remember a coordinate for checking
        mycell = self.createTestCell()
        test_coord = mycell.blocks[ mycell.blocks.keys()[0] ].coord(4)

        outfile = "./testCellIO.pkl"
        mycell.writePickle(outfile)
        with open(outfile) as f:
            newCell = cPickle.load(f)

        self.assertTrue(numpy.allclose(test_coord, mycell.blocks[ mycell.blocks.keys()[0] ].coord(4),
                                         rtol=1e-9, atol=1e-9),
                         msg="Incorrect testCoordinate of cell.")

        # mycell.growBlocks( 5 )
        os.unlink(outfile)

        return

    def testCloseAtoms(self):

        mycell = Cell([30, 30, 30], atomMargin=0.1, bondMargin=0.1, bondAngleMargin=15)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')
        block1 = mycell.getInitBlock('A')

        block1.translateCentroid([ 1, 1, 1 ])
        block1_id = mycell.addBlock(block1)
        block2 = mycell.getInitBlock('A')
        block2.translateCentroid([ 2.2, 2.2, 2.2 ])
        block2_id = mycell.addBlock(block2)
 
        closeList, wallClash = mycell.closeAtoms(block1_id)
         
        self.assertFalse(wallClash,"Got clash with wall")
         
        # Updated refPairs as now have
        refPairs = [ (0, 3), (1, 3), (2, 3) ]  # Also used to check PBC
        closePairs = []
        for iatom, _, ioatom, _ in closeList:
            closePairs.append((iatom, ioatom))
 
        # mycell.writeXyz("close1.xyz", label=False)
 
        self.assertEqual(closePairs,
                        refPairs,
                         "Many contacts: {0}".format(closePairs))
 
        # Too far for any contacts
        mycell.delBlock(block2_id)
        block2 = mycell.getInitBlock('A')
        block2.translateCentroid([10, 10, 10])
        block2_id = mycell.addBlock(block2)
 
        close, wallClash = mycell.closeAtoms(block1_id)
        self.assertFalse(wallClash, "Wallclash")
        self.assertFalse(close, "No contacts: ".format(close))
 
        # Now check across periodic boundary
        mycell.delBlock(block2_id)
        block2 = mycell.getInitBlock('A')
        x = 2.2 + 2 * mycell.dim[0]
        y = 2.2 + 2 * mycell.dim[1]
        z = 2.2 + 2 * mycell.dim[2]

        block2.translateCentroid([x, y, z])
        block2_id = mycell.addBlock(block2)

        # mycell.writeXyz("close2.xyz", label=False)

        closeList, wallClash = mycell.closeAtoms(block1_id)
        self.assertFalse(wallClash, "Wallclash")
        closePairs = []
        for iatom, block, ioatom, distance in closeList:
            closePairs.append((iatom, ioatom))

        self.assertEqual(closePairs, refPairs, "Periodic boundary: {}".format(closePairs))

        return

    def testCloseAtoms2(self):

        mycell = Cell([2.1, 2.1, 2.1], atomMargin=0.1, bondMargin=0.1, bondAngleMargin=15)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')
        block1 = mycell.getInitBlock('A')

        natoms = block1.numAtoms()

        # block radius is 1.8
        block1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
        block1Idx = mycell.addBlock(block1)

        # Alone in cell but in center
        close, wallClash = mycell.closeAtoms(block1Idx)
        self.assertFalse(close)
        self.assertFalse(wallClash)

        # Add second block overlapping first but in other image
        block2 = mycell.getInitBlock('A')
        block2.translateCentroid([ mycell.dim[0] / 2 + mycell.dim[0],
                                   mycell.dim[1] / 2 + mycell.dim[1],
                                   mycell.dim[2] / 2 + mycell.dim[2] ])
        block2Idx = mycell.addBlock(block2)

        # Make sure every atom overlaps with ever other
        close, wallClash = mycell.closeAtoms(block1Idx)
        self.assertFalse(wallClash)
        self.assertEquals(natoms * natoms, len(close))
        
        # See we have enough clashing atoms - NOT CHECKED THIS NUMBER
        close, wallClash = mycell.closeAtoms(block2Idx)
        self.assertFalse(wallClash)
        self.assertEquals(natoms * natoms, len(close))

        return

    def testCloseDistance(self):
        """
        Test distance and close together
        """
        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)
        # mycell.initCell("../ch4_typed.car",incell=False)
        # block1 = mycell.initBlock.copy()
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')
        block1 = mycell.getInitBlock('A')

        b1 = numpy.array([2, 2, 2], dtype=numpy.float64)
        block1.translateCentroid(b1)
        block1_id = mycell.addBlock(block1)

        # block2=mycell.initBlock.copy()
        block2 = mycell.getInitBlock('A')
        b2 = numpy.array([3, 3, 3], dtype=numpy.float64)
        block2.translateCentroid(b2)
        block2_id = mycell.addBlock(block2)

        # mycell.writeXyz("close1.xyz", label=False)

#        closeList =  mycell.closeAtoms(block1_id)
#        for iatom, ioblock, ioatom in closeList:
#            b1c = block1._coords[iatom]
#            b2c = mycell.blocks[ioblock]._coords[ioatom]
#            distance = mycell.distance(b1c,b2c)
#            print "{}: {} -> {}:{} = {}".format(iatom,b1c,ioatom,b2c,distance)

        # Distance measured with Avogadro so it MUST be right...
        refd = 0.673354948616
        distance = mycell.distance(block1.coord(1), mycell.blocks[ block2_id ].coord(3))
        self.assertAlmostEqual(refd, distance, 12, "Closest atoms: {}".format(distance))
        return
    
    def testDeleteBlocksIndices(self):
        
        boxDim = [50, 50, 50]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')

        # See if we can delete a range of blocks
        toSeed = 10
        seeded = mycell.seed(toSeed)
        self.assertEqual(seeded,toSeed,"Failed to seed all blocks")
        
        toDelete = [0,3,5]
        
        bidxs = [ blockId for i, blockId in enumerate(mycell.blocks.iterkeys()) if i in toDelete]
        mycell.deleteBlocksIndices(toDelete,save=True)
        self.assertEqual(len(mycell.blocks), toSeed-len(toDelete))
        
        for idx in bidxs:
            self.assertNotIn(idx, mycell.blocks.keys(), "Block was left in list")
            
        mycell.restoreBlocks()
        self.assertEqual(len(mycell.blocks), toSeed)
        return

    def testDeleteBlocksType(self):
        # Cell dimensions need to be: L > 2*(r_cut+r_buff) and L < 3*(r_cut+r_buff)
        # From code looks like default r_buff is 0.4 and our default r_cut is 5.0
        boxDim = [50, 50, 50]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.benzene2Car, fragmentType='A')
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='B')
        mycell.addBondType('A:a-A:a')
        mycell.addBondType('A:a-B:a')

        # See if we can delete single blocks
        toSeed = 10
        mycell.seed(toSeed, fragmentType='A')
        toDelete = toSeed - 1
        deleted = mycell.deleteBlocksType('A', numBlocks=toDelete)
        self.assertEqual(deleted, toDelete)
        self.assertEqual(mycell.numBlocks(), toSeed - deleted)

        # Now try deleting multiple blocks of a single type
        mycell.growBlocks(8, cellEndGroups=['A:a'], libraryEndGroups=['A:a'])
        deleted = mycell.deleteBlocksType('A', numBlocks=1)
        self.assertEqual(deleted, 0)
        deleted = mycell.deleteBlocksType('A', numBlocks=1, multiple=True)
        self.assertEqual(deleted, 1)
        
        return
    
    def testDeleteBondType(self):
        boxDim = [50, 50, 50]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.benzene2Car, fragmentType='A')
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='B')
        mycell.addBondType('A:a-A:a')
        mycell.addBondType('A:a-B:a')
        self.assertEqual(len(mycell.bondTypes),2)
        mycell.seed(5, fragmentType='A')
        toGrow = 1
        added = mycell.growBlocks(toGrow, cellEndGroups=['A:a'], libraryEndGroups=['B:a'])
        self.assertEqual(added, toGrow, "Failed to grow in testDeleteBondType")
        mycell.deleteBondType('A:a-B:a')
        self.assertEqual(len(mycell.bondTypes),1)
        added = mycell.growBlocks(toGrow, cellEndGroups=['A:a'], libraryEndGroups=['B:a'])
        self.assertEqual(added, 0, "Added disallowed bond in testDeleteBondType")
        return
        
    def testDeleteFragment1(self):
        boxDim = [20, 20, 20]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # See if we can delete single blocks
        mycell.seed(1, fragmentType='A', center=True)
        f1 = mycell.blocks.values()[0].fragments[0]
        mycell.deleteFragment(f1)
        self.assertEqual(len(mycell.blocks),0)
        return 
    
    def testDeleteFragment2(self):
        boxDim = [20, 20, 20]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # See if we can delete single blocks
        mycell.seed(1, fragmentType='A', center=True)
        f1 = mycell.blocks.values()[0].fragments[0]
        mycell.growBlocks(1, cellEndGroups=['A:a'], libraryEndGroups=['A:a'])
        mycell.deleteFragment(f1)
        return 
    
    def testDeleteFragment3(self):
        boxDim = [20, 20, 20]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # See if we can delete single blocks
        mycell.seed(1, fragmentType='A', center=True)
        f1 = mycell.blocks.values()[0].fragments[0]
        mycell.growBlocks(20, cellEndGroups=['A:a'], libraryEndGroups=['A:a'])
        mycell.deleteFragment(f1)
        #mycell.writeCml("foo2.cml")
        return 

    def testDLPOLY(self):
        """testDLPOLY"""
        boxDim = [100.0, 100.0, 100.0]
        mycell = Cell(boxDim)

        ch4Car = self.ch4Car
        # mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b2 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b3 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')

        # b1 in center
        b1.translateCentroid([ 25, 25, 25 ])
        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]
        b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)

        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock(bond)
        mycell.addBlock(b1)

        b3.translateCentroid([ 25, 25, 20 ])
        mycell.addBlock(b3)

        d = opt.DLPOLY()

        # data = mycell.dataDict(periodic=True, center=True, rigidBody=True)
        d.writeCONTROL()
        d.writeFIELDandCONFIG(mycell)

        os.unlink('CONTROL')
        os.unlink('CONFIG')
        os.unlink('FIELD')

        return

    def testDistance(self):
        """Test the distance under periodic boundary conditions"""

        CELLA = CELLB = CELLC = 10
        boxDim = [CELLA, CELLB, CELLC]
        mycell = Cell(boxDim)

        v1 = [ 2.46803012, 1.67131881, 1.96745421]
        v2 = [ 1.07988345, 0.10567109, 1.64897769]

        nv1 = numpy.array(v1)
        nv2 = numpy.array(v2)

        dc1 = mycell.distance(nv1, nv2)
        dn = numpy.linalg.norm(nv2 - nv1)
        self.assertEqual(dc1, dn, "Distance within cell:{} | {}".format(dc1, dn))

        x = v2[0] + 2 * CELLA
        y = v2[1] + 2 * CELLB
        z = v2[2] + 2 * CELLC
        nv2 = numpy.array([ x, y, z ])
        dc2 = mycell.distance(nv1, nv2)
        self.assertAlmostEqual(dc1, dc2, 11, "Distance across multiple cells +ve: {} | {}".format(dc1, dc2))

        x = v2[0] - 2 * CELLA
        y = v2[1] - 2 * CELLB
        z = v2[2] - 2 * CELLC
        nv2 = numpy.array([ x, y, z ])
        dc3 = mycell.distance(nv1, nv2)
        self.assertAlmostEqual(dc1, dc3, 11, "Distance across multiple cells -ve: {} | {}".format(dc1, dc3))

        v1 = numpy.array([ 0.0, 0.0, 0.0 ])
        v2 = numpy.array([ 0.0, 0.0, 8.0 ])
        dc = mycell.distance(v1, v2)
        self.assertEqual(dc, 2.0, "Distance across boundary cell:{}".format(dc))

        return

    def testDihedral(self):

        CELLDIM = 30
        boxDim = [CELLDIM, CELLDIM, CELLDIM]
        mycell = Cell(boxDim)

        p1 = numpy.array([ 0.0, 0.0, 0.0 ])
        p2 = numpy.array([ 10.0, 0.0, 0.0 ])
        p3 = numpy.array([ 10.0, 10.0, 0.0 ])
        p4 = numpy.array([ 20.0, 10.0, 10.0 ])

        ref = util.dihedral(p1, p2, p3, p4)
        
        self.assertEqual(ref, mycell.dihedral(p1, p2, p3, p4))

        # Move by a full cell along x-axis - result should be the same
        p3 = numpy.array([ 10.0 + CELLDIM, 10.0, 0.0 ])
        p4 = numpy.array([ 20.0 + CELLDIM, 10.0, 10.0 ])

        self.assertEqual(ref, mycell.dihedral(p1, p2, p3, p4))

        return

    def testEndGroupTypes(self):
        """Test we can add a block correctly"""

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(fragmentType='A', filename=self.ch4Car)
        mycell.libraryAddFragment(fragmentType='B', filename=self.ch4Car)
        mycell.libraryAddFragment(fragmentType='C', filename=self.ch4Car)
        mycell.libraryAddFragment(fragmentType='D', filename=self.ch4Car)

        # Everything can bond to A (apart from A itself), but nothing can bond to anything else
        mycell.addBondType('A:a-B:a')
        mycell.addBondType('A:a-C:a')
        mycell.addBondType('A:a-D:a')

        mycell.seed(1, fragmentType='A')
        mycell.seed(1, fragmentType='B')
        mycell.seed(1, fragmentType='C')

        banned = ['B:a', 'C:a', 'D:a']
        for _ in xrange(5):
            cellEndGroup, libraryEndGroup = mycell.libraryEndGroupPair()
            if cellEndGroup.type() != 'A:a':
                self.assertNotIn(libraryEndGroup.type(),
                                  banned,
                                  "cell: {0} library: {1}".format(cellEndGroup.type(), libraryEndGroup.type())
                                  )

        # Check it works if we give a cellEndGroup to bond
        cellEndGroup, libraryEndGroup = mycell.libraryEndGroupPair(cellEndGroups='B:a')

        # Check it fails if we ask for a block not in the cell
        self.assertRaises(RuntimeError, mycell.libraryEndGroupPair, {'cellEndGroups', 'D:a'})

        # Check it works if we give it a list of libraryBlocks
        cellEndGroup, libraryEndGroup = mycell.libraryEndGroupPair(libraryEndGroups='B:a')
        self.assertEqual(libraryEndGroup.type(), 'B:a')

        return

    def testFragMaxEnergy(self):
        """

        run calc
        find highest energy fragment => need to link from a compute/tag to a block/fragment

        - need list of block computes
        - need list of fragment computes

        each groupname is of the form X:Y where X is block index in cell and Y is fragment in index in the block

        The list of groups isn't needed, just the labels

        have a list of labels as long as there are groups - same as number of computes

        find highest energy compute - the index links back to the list of labels, and from the label we get the block/framgent

        block compute could be the first in the list

        each block has tags
        each fragment has tags

        tags = [
                 [
                   [s,e],[s,e],[s,e]
                   ]

        createGroups



        """

        boxDim = [100.0, 100.0, 100.0]
        mycell = Cell(boxDim)

        ch4Car = self.ch4Car
        # mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b2 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b3 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')


        # b1 in center
        b1.translateCentroid([ 25, 25, 25 ])
        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]
        b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)

        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock(bond)
        mycell.addBlock(b1)

        b3.translateCentroid([ 25, 25, 20 ])
        mycell.addBlock(b3)

        # mycell.writeXyz("foo.xyz")

        # Calculate the energy
        mycell.fragMaxEnergy()

        return

    def testGrowBlocks(self):
        """Test we can add blocks correctly"""
        boxDim = [30, 30, 30]
        mycell = Cell(boxDim, doLog=False)
        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        # mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        added = mycell.seed(1, center=True, random=False)
        self.assertEqual(added, 1, 'seed')
        natoms = mycell.numAtoms()
        nblocks = 10
        added = mycell.growBlocks(nblocks, endGroupType=None, maxTries=1, random=False)
  
        # mycell.writeCml("foo.cml", periodic=True, pruneBonds=False)
        self.assertEqual(added, nblocks, "growBlocks did not return ok")
        self.assertEqual(1, len(mycell.blocks), "Growing blocks found {0} blocks".format(len(mycell.blocks)))
  
        natoms2 = mycell.numAtoms()
        nblocks += 1
        # Need to subtract cap atoms
        self.assertEqual(natoms2, (natoms * nblocks) - (nblocks - 1) * 2)
        block = mycell.blocks.values()[0]
        self.assertTrue(numpy.allclose(block.centroid(), [ 12.91963557,  17.39975016,  11.65120767]))
        self.assertFalse(self.clashes(mycell))
        return
    
    def testGetBox(self):
        mycell = Cell([10, 10, 10])   
        boxSize = 1
        mycell.boxSize = boxSize
        mycell.numBoxes = [int(math.ceil(mycell.dim[0] / boxSize)),
                         int(math.ceil(mycell.dim[1] / boxSize)),
                         int(math.ceil(mycell.dim[2] / boxSize)) ]
        
        p1 = numpy.array([0, 0, 0])
        self.assertEqual(mycell._getBox(p1), (0, 0, 0))
        
        p1 = numpy.array([5.5, 5.5, 5.5])
        self.assertEqual(mycell._getBox(p1), (5, 5, 5))
        
        p1 = numpy.array([-5.5, -5.5, -5.5])
        self.assertEqual(mycell._getBox(p1), (4, 4, 4))
        
        p1 = numpy.array([10.5, 10.5, 10.5])
        self.assertEqual(mycell._getBox(p1), (0, 0, 0))
        
        p1 = numpy.array([100.5, 100.5, 100.5])
        self.assertEqual(mycell._getBox(p1), (0, 0, 0))
        
        p1 = numpy.array([-100.5, -100.5, -100.5])
        self.assertEqual(mycell._getBox(p1), (9, 9, 9))
        
        p1 = numpy.array([-3, 0, 0])
        self.assertEqual(mycell._getBox(p1), (7, 0, 0))
        
        return

    def testGrowLimited(self):
        """Test we can add blocks correctly"""
        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4_1Car, fragmentType='A')
        mycell.addBondType('A:a-A:b')
        mycell.setMaxBond('A:a', 1)

        mycell.seed(1, center=True, random=False)
        mycell.growBlocks(10, cellEndGroups=None, maxTries=10, random=False)
        block = mycell.blocks.values()[0]
        self.assertTrue(numpy.allclose(block.centroid(), [ 14.64787487, 15.43129673, 16.00520584]))
        self.assertFalse(self.clashes(mycell))
        return

    def testGrowLimited2(self):
        """Test we can add blocks correctly"""
        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.dcxCar, fragmentType='A')
        mycell.addBondType('A:CH-A:CCl')
        mycell.setMaxBond('A:CH', 1)
        mycell.seed(1, center=True, random=False)
        mycell.growBlocks(10, cellEndGroups=['A:CH'], random=False)
        block = mycell.blocks.values()[0]
        self.assertTrue(numpy.allclose(block.centroid(), [  6.17554446, 4.10903953, 17.51814496]))
        self.assertFalse(self.clashes(mycell))
        return

    def testGrowBlocks2(self):
        """Test we can add blocks correctly"""

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        block1 = mycell.getInitBlock('A')

        # Position block so that it's aligned along x-axis
        # - use two opposing C-atoms 0 & 3
        block1.alignAtoms(0, 3, [ 1, 0, 0 ])

        # Now put it at the center of the cell
        block1.translateCentroid([mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2])

        # Add to the cell
        mycell.addBlock(block1)

        # Try adding 6 blocks - only 5 will fit
        nblocks = 6
        added = mycell.growBlocks(nblocks, endGroupType=None, maxTries=5, random=False)

        self.assertEqual(added, 5, "growBlocks did not add 5 blocks")
        self.assertEqual(1, len(mycell.blocks), "Growing blocks found {0} blocks".format(len(mycell.blocks)))
        self.assertFalse(self.clashes(mycell))
        self.assertTrue(numpy.allclose(block1.centroid(), [ 12.69119085, 14.93774993, 14.96541503]))
        return

    def testGrowBlocksDihedral(self):
        """Test we can add blocks correctly"""

        boxDim = [10, 10, 10]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')
        mycell.seed(1, center=True)

        nblocks = 2
        dihedral = 90
        mycell.growBlocks(nblocks, cellEndGroups=None, libraryEndGroups=None, dihedral=dihedral, maxTries=5)

        # Check angle between two specified dihedrals is correct
        block1 = mycell.blocks[ mycell.blocks.keys()[0] ]
        bond1 = block1._blockBonds[0]

        p1 = block1.coord(bond1.endGroup1.dihedralIdx())
        p2 = block1.coord(bond1.endGroup1.endGroupIdx())
        p3 = block1.coord(bond1.endGroup2.endGroupIdx())
        p4 = block1.coord(bond1.endGroup2.dihedralIdx())

        self.assertAlmostEqual(math.degrees(util.dihedral(p1, p2, p3, p4)), dihedral)
        self.assertFalse(self.clashes(mycell))

        return

    def testGrowBlocksUw(self):
        """Test we can add blocks correctly"""

        boxDim = [10, 10, 10]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.amineCar, fragmentType='amine')
        mycell.libraryAddFragment(filename=self.triquinCar, fragmentType='triquin')
        mycell.addBondType('amine:a-triquin:b')

        mycell.seed(1, fragmentType='triquin', center=True)
        mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups=['amine:a'], maxTries=1)
        self.assertFalse(self.clashes(mycell))

        return
    
    def testIntersectedCells(self):
        
        mycell = Cell([10, 10, 10])   
        boxSize = 1
        mycell.boxSize = boxSize
        mycell.numBoxes = [int(math.ceil(mycell.dim[0] / boxSize)),
                         int(math.ceil(mycell.dim[1] / boxSize)),
                         int(math.ceil(mycell.dim[2] / boxSize)) ]
         
#         p1 = numpy.array([0.5, 0.5, 0.5])
#         p2 = numpy.array([8.5, 8.5, 8.5])
#         ref = [(0, 0, 0), (0, 0, 9), (0, 0, 8), (0, 9, 8), (0, 8, 8), (9, 8, 8), (8, 8, 8)]
#         self.assertEqual(mycell._intersectedCells(p1, p2), ref)
#          
#         p1 = numpy.array([0.5, 0.5, 0.5])
#         p2 = numpy.array([0.5, 8.5, 8.5])
#         ref = [(0, 0, 0), (0, 0, 9), (0, 0, 8), (0, 9, 8), (0, 8, 8)]
#         self.assertEqual(mycell._intersectedCells(p1, p2), ref)
#          
#         p1 = numpy.array([0.5, 0.5, 0.5])
#         p2 = numpy.array([0.6, 8.5, 8.5])
#         ref = [(0, 0, 0), (0, 0, 9), (0, 0, 8), (0, 9, 8), (0, 8, 8)]
#         self.assertEqual(mycell._intersectedCells(p1, p2), ref)
#          
#         p1 = numpy.array([0.5, 0.5, 0.5])
#         p2 = numpy.array([0.5, 0.5, 8.5])
#         ref = [(0, 0, 0), (0, 0, 9), (0, 0, 8)]
#         self.assertEqual(mycell._intersectedCells(p1, p2), ref)
#  
#         p1 = numpy.array([-2.5, -2.5, -2.5])
#         p2 = numpy.array([8.5, 8.5, 8.5])
#         ref = [(7, 7, 7), (7, 7, 8), (7, 8, 8), (8, 8, 8)]
#         self.assertEqual(mycell._intersectedCells(p1, p2), ref)
#      
#         p1 = numpy.array([5., 0., 0.])
#         p2 = numpy.array([-3., 0., 0.])
#         self.assertEqual(mycell._intersectedCells(p1, p2), [(5, 0, 0), (6, 0, 0), (7, 0, 0)])
         
        p1 = numpy.array([9., 9., 9.])
        p2 = numpy.array([1., 1., 1.])
        self.assertEqual(mycell._intersectedCells(p1, p2), [(9, 9, 9), (9, 9, 0), (9, 0, 0),
                                                          (0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)])


        mycell = Cell([25, 25, 25])   
        mycell.boxSize = 1.91761148234
        mycell.numBoxes = [14, 14, 14]

        p1 = numpy.array([ -3.46529620e+00, 7.99752596e-03, 1.86788673e-03])
        p2 = numpy.array([  3.45583501e+00, 7.99752596e-03, 1.86788673e-03])
        self.assertEqual(mycell._intersectedCells(p1, p2), [(11, 0, 0), (12, 0, 0), (13, 0, 0),
                                                          (0, 0, 0), (1, 0, 0)])

        return

    def testJoinBlocks(self):
        """Test we can add a block correctly"""

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)

        # mycell.libraryAddFragment(filename=self.pafCar, fragmentType='A')
        mycell.libraryAddFragment(filename=self.benzene2Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        toAdd = 5
        added = mycell.seed(toAdd, center=True)
        self.assertEqual(added, toAdd, 'seed')
        nc = 0
        for block in mycell.blocks.itervalues():
            nc += block.numAtoms()

        toJoin = 4
        added = mycell.joinBlocks(toJoin, cellEndGroups=None, maxTries=1)
        self.assertEqual(added, toJoin, "joinBlocks did join enough")
        self.assertEqual(1, len(mycell.blocks), "joinBlocks found {0} blocks".format(len(mycell.blocks)))

        nc2 = 0
        for block in mycell.blocks.itervalues():
            nc2 += block.numAtoms()

        # Need to subtract cap atoms
        self.assertEqual(nc - (toJoin * 2), nc2, "Growing blocks found {0} coords".format(nc2))
        self.assertFalse(self.clashes(mycell))

        return

    def testOptimiseGeometryRigid(self):
        """
        """
        # self.testCell.writeCml("foo.cml")
        mycell = self.createTestCell()
        mycell.optimiseGeometry(rigidBody=True,
                                        doDihedral=True,
                                        quiet=True,
                                        rCut=5.0,
                                        optCycles=1000000,
                                        dt=0.005,
                                        Etol=1e-5,
                                         )
        self.assertFalse(self.clashes(mycell))
        os.unlink("hoomdOpt.xml")
        return

    def testOptimiseGeometryAll(self):
        """
        """
        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.amineCar, fragmentType='amine')
        mycell.libraryAddFragment(filename=self.triquinCar, fragmentType='triquin')
        mycell.addBondType('amine:a-triquin:b')

        mycell.seed(1, fragmentType='triquin', center=True)
        mycell.growBlocks(toGrow=2, cellEndGroups=None, libraryEndGroups=['amine:a'], maxTries=1)
#         ok=mycell.optimiseGeometry(rigidBody=False,
#                                 doDihedral=True,
#                                 optCycles=10000,
#                                 dump=True,
#                                 quiet=False,
#                                 dt=0.005,
#                                 Nmin=5,
#                                 alpha_start=0.1,
#                                 ftol=1e-2,
#                                 Etol=1e-4,
#                                 finc=1.1,
#                                 fdec=0.5,
#                                  )
        ok = mycell.optimiseGeometry(rigidBody=False,
                                doDihedral=True,
                                optCycles=10000,
                                dump=False,
                                quiet=False
                                 )

        self.assertTrue(ok, "All Optimisation failed!")
        self.assertFalse(self.clashes(mycell))
        # os.unlink("hoomdOpt.xml")
        return


    def testRunMDAll(self):
        """
        """

        boxDim = [30, 30, 30]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.amineCar, fragmentType='amine')
        mycell.libraryAddFragment(filename=self.triquinCar, fragmentType='triquin')
        mycell.addBondType('amine:a-triquin:b')

        mycell.seed(1, fragmentType='triquin', center=True)
        mycell.growBlocks(toGrow=1, cellEndGroups=None, libraryEndGroups=['amine:a'], maxTries=1)
        # self.testCell.writeCml("foo.cml")
        mycell.runMD(rigidBody=False,
                     doDihedral=True,
                     rCut=5.0,
                     mdCycles=1000,
                     T=1.0,
                     tau=0.5,
                     dt=0.0005)
        self.assertFalse(self.clashes(mycell))
        # os.unlink("hoomdOpt.xml")
        return

    def XtestOptimiseGeometryMinCell(self):
        """
        """

        # Create a large cell and populate it with a block
        CELLA = CELLB = CELLC = 100
        mycell = Cell()
        mycell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        mycell.seed(1)
        mycell.growBlocks(3, 'A')

        # Get the block and put it in the centre of the cell
        block = mycell.blocks[ mycell.blocks.keys()[0] ]
        block.translateCentroid([ CELLA / 2, CELLB / 2, CELLC / 2 ])

        # mycell.dump()
        # Center of mass shouldn't have changed
        com = 0.0
        for _, block in mycell.blocks.iteritems():
            com += block.centerOfMass()

        mycell.optimiseGeometry(minCell=True, optAttempts=1)

        comA = 0.0
        for _, block in mycell.blocks.iteritems():
            comA += block.centerOfMass()

        self.assertTrue(numpy.allclose(com, comA))

        return

    def testOptimiseGeometryDihedral(self):
        """
        """
        mycell = self.createTestCell()
        mycell.optimiseGeometry(doDihedral=True, quiet=True)
        self.assertFalse(self.clashes(mycell))
        os.unlink("hoomdOpt.xml")
        return
    
    def testOptimiseGeometryStatic(self):
        """
        """

        mycell = Cell(filePath=self.graphiteCar)
        mycell.libraryAddFragment(filename=self.amineCar, fragmentType='amine')
        mycell.libraryAddFragment(filename=self.triquinCar, fragmentType='triquin')
        mycell.addBondType('amine:a-triquin:b')

        mycell.seed(3, fragmentType='triquin', center=True)
        mycell.growBlocks(toGrow=2, cellEndGroups=None, libraryEndGroups=['amine:a'], maxTries=10)
        
        # mycell.dump()

        ok = mycell.optimiseGeometry(rigidBody=False,
                                doDihedral=True,
                                optCycles=10000,
                                dump=False,
                                quiet=False
                                 )
        self.assertFalse(self.clashes(mycell))
        # mycell.dump()

        return

    def testRunMD(self):
        """
        """
        mycell = self.createTestCell()
        mycell.runMD(doDihedral=True,
                             quiet=True,
                             rCut=5.0,
                             mdCycles=100,
                             T=1.0,
                             tau=0.5,
                             dt=0.0005,
                         )
        self.assertFalse(self.clashes(mycell))
        os.unlink("hoomdMD.xml")
        return
    
    def testRunMdNpt(self):
        """
        """
        mycell = self.createTestCell()
        mycell.runMD(doDihedral=True,
                     quiet=False,
                     rCut=5.0,
                     mdCycles=100,
                     T=1.0,
                     tau=0.5,
                     P=1.0,
                     tauP=0.5,
                     dt=0.0005,
                     dump=True,
                     dumpPeriod=20,
                     integrator='npt',
                    )
        self.assertFalse(self.clashes(mycell))
        os.unlink("hoomdMD.xml")
        #os.unlink("runmd.dcd")
        return

    def testRunMDAndOptimise(self):
        """
        """
        mycell = self.createTestCell()
        mycell.runMDAndOptimise(doDihedral=True, quiet=True)
        self.assertFalse(self.clashes(mycell))
        os.unlink("hoomdMDOpt.xml")
        return

    def testPeriodic(self):

        import hoomdblue

        mycell = self.createTestCell()

        # Grab coords
        coords = []
        for block in mycell.blocks.itervalues():
            for i, coord in enumerate(block.iterCoord()):
                coords.append(coord)

        # Wrap them
        wcoords = []
        images = []
        for c in coords:
            x, xi = util.wrapCoord(c[0], mycell.dim[0], center=False)
            y, yi = util.wrapCoord(c[1], mycell.dim[1], center=False)
            z, zi = util.wrapCoord(c[2], mycell.dim[2], center=False)
            wcoords.append([ x, y, z ])
            images.append([ xi, yi, zi ])

        # Now umwrap them
        for i, c in enumerate(wcoords):
            x = util.unWrapCoord(c[0], images[i][0], mycell.dim[0], centered=False)
            y = util.unWrapCoord(c[1], images[i][1], mycell.dim[1], centered=False)
            z = util.unWrapCoord(c[2], images[i][2], mycell.dim[2], centered=False)

            self.assertAlmostEqual(x,
                                    coords[ i ][ 0 ],
                                    places=9,
                                    msg="{0} {1}".format(x, coords[ i ][ 0 ]))
            self.assertAlmostEqual(y,
                                    coords[ i ][ 1 ],
                                    places=9,
                                    msg="{0} {1}".format(y, coords[ i ][ 1 ]))
            self.assertAlmostEqual(z,
                                    coords[ i ][ 2 ],
                                    places=9,
                                    msg="{0} {1}".format(z, coords[ i ][ 2 ]))

        # Now wrap them with centering
        wcoords = []
        images = []
        for c in coords:
            x, xi = util.wrapCoord(c[0], mycell.dim[0], center=True)
            y, yi = util.wrapCoord(c[1], mycell.dim[1], center=True)
            z, zi = util.wrapCoord(c[2], mycell.dim[2], center=True)
            wcoords.append([ x, y, z ])
            images.append([ xi, yi, zi ])

        # Now umwrap them
        for i, c in enumerate(wcoords):
            x = util.unWrapCoord(c[0], images[i][0], mycell.dim[0], centered=True)
            y = util.unWrapCoord(c[1], images[i][1], mycell.dim[1], centered=True)
            z = util.unWrapCoord(c[2], images[i][2], mycell.dim[2], centered=True)

            self.assertAlmostEqual(x,
                                    coords[ i ][ 0 ],
                                    places=9,
                                    msg="{0} {1}".format(x, coords[ i ][ 0 ]))

            self.assertAlmostEqual(y,
                                    coords[ i ][ 1 ],
                                    places=9,
                                    msg="{0} {1}".format(y, coords[ i ][ 1 ]))
            self.assertAlmostEqual(z,
                                    coords[ i ][ 2 ],
                                    places=9,
                                    msg="{0} {1}".format(z, coords[ i ][ 2 ]))

        # Now test with HOOMD-Blue
        filename = "periodicTest.xml"
        data = mycell.dataDict(periodic=True, center=True, rigidBody=True)
        o = opt.HoomdOptimiser()
        o.writeXml(data,
                   xmlFilename=filename,
                   rigidBody=True,
                   doDihedral=True,
                   doImproper=False,
                   doCharges=True,
                   )

        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        system = hoomdblue.init.read_xml(filename=filename)

        wcoords = []
        for i, p in enumerate(system.particles):
            x, y, z = p.position
            ix, iy, iz = p.image

            x = util.unWrapCoord(x, ix, mycell.dim[0], centered=True)
            y = util.unWrapCoord(y, iy, mycell.dim[1], centered=True)
            z = util.unWrapCoord(z, iz, mycell.dim[2], centered=True)

            self.assertAlmostEqual(x,
                                    coords[ i ][ 0 ],
                                    places=6,
                                    msg="{0} {1}".format(x, coords[ i ][ 0 ]))

            self.assertAlmostEqual(y,
                                    coords[ i ][ 1 ],
                                    places=6,
                                    msg="{0} {1}".format(y, coords[ i ][ 1 ]))
            self.assertAlmostEqual(z,
                                    coords[ i ][ 2 ],
                                    places=6,
                                    msg="{0} {1}".format(z, coords[ i ][ 2 ]))
        os.unlink(filename)
        return
        
    def testSeed(self):
        """Test we can seed correctly"""

        boxDim = [50, 50, 50]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        # mycell.libraryAddFragment( filename=self.ch4Car, fragmentType='A' )
        mycell.addBondType('A:a-A:a')

        nblocks = 10
        added = mycell.seed(nblocks)

        self.assertEqual(nblocks, added, "Incorrect number of cell blocks")
        self.assertEqual(nblocks, len(mycell.blocks), "Incorrect number of cell blocks")

        # Check no block is outside the unit cell
        # We seed by the centroid so the edges could stick out by radius
        # radius = ch4.radius()
        
        bad = []
        for i, b in mycell.blocks.iteritems():
            radius = b.blockRadius()
            for c in b.iterCoord():
                if not 0 - radius < c[0] < mycell.dim[0] + radius  or \
                   not 0 - radius < c[1] < mycell.dim[1] + radius  or \
                   not 0 - radius < c[2] < mycell.dim[2] + radius :

                    bad.append(c)

        self.assertEqual(0, len(bad), "Got {} atoms outside cell: {}".format(len(bad), bad))
        self.assertFalse(self.clashes(mycell))

        # mycell.writeXyz("seedTest.xyz")

        return
    
    def testSetBoxSize(self):

        boxDim = [10, 10, 10]
        mycell = Cell(boxDim)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        mycell.seed(5)
        mycell.growBlocks(5)
        mycell.setBoxSize([20, 20, 20])
        mycell.growBlocks(5)

        return
    
    def testSolvent(self):
        """Create a block as solvent and then ensure there are no clashes
        - uses same setup as ZipClash test
        """

        boxDim = [25, 25, 25]
        mycell = Cell(boxDim, doLog=True)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='B', solvent=True)
        mycell.addBondType('A:a-A:a')

        # Create blocks manually
        b1 = mycell.getInitBlock(fragmentType='A')

        # Align bond along x-axis
        b1.alignAtoms(0, 3, [ 0, 0, 1 ])
        b1.alignAtoms(0, 3, [ 0, 1, 0 ])
        b1.alignAtoms(0, 3, [ 1, 0, 0 ])

        b3 = b1.copy()
        
        # b1 in center of cell -5 on x-axis
        b1.translateCentroid([ (mycell.dim[0] / 2) - 5, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
        
        # Get b2 from library and set up so that it's in the middle
        b2 = mycell.getInitBlock(fragmentType='B')
        
        b2.alignAtoms(0, 3, [ 0, 0, 1 ])
        b2.alignAtoms(0, 3, [ 0, 1, 0 ])
        b2.alignAtoms(0, 3, [ 1, 0, 0 ])

        # b2 in center
        b2.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # rotate so ring facing bond axis
        b2.rotate([0, 1, 0], math.pi / 2, center=b2.centroid())
        
        # b3 + 5 on x
        b3.translateCentroid([ (mycell.dim[0] / 2) + 5, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # Add blocks either side to cell and make sure we can bond
        mycell.addBlock(b1)
        mycell.addBlock(b3)

        # Add the pesky middle block and see if it now works
        mycell.addBlock(b2)
        
        made = mycell.zipBlocks(bondMargin=6, bondAngleMargin=1, clashCheck=True, clashDist=1.8)
        
        self.assertEqual(made,1)

        pass

    def testSurroundBoxes(self):
        """
        """
        boxDim = [5, 5, 5]
        mycell = Cell(boxDim)
        # box size=1 - need to set manually as not reading in a block
        mycell.maxAtomRadius = 0.5
        mycell.atomMargin = 0.0

        mycell.boxSize = (mycell.maxAtomRadius * 2) + mycell.atomMargin
        mycell.numBoxes = [int(math.ceil(mycell.dim[0] / mycell.boxSize)),
                         int(math.ceil(mycell.dim[1] / mycell.boxSize)),
                         int(math.ceil(mycell.dim[2] / mycell.boxSize)) ]

        s = [(2, 2, 2), (2, 2, 1), (2, 2, 3), (2, 1, 2), (2, 1, 1), (2, 1, 3), (2, 3, 2), (2, 3, 1),
         (2, 3, 3), (1, 2, 2), (1, 2, 1), (1, 2, 3), (1, 1, 2), (1, 1, 1), (1, 1, 3), (1, 3, 2),
          (1, 3, 1), (1, 3, 3), (3, 2, 2), (3, 2, 1), (3, 2, 3), (3, 1, 2), (3, 1, 1), (3, 1, 3),
          (3, 3, 2), (3, 3, 1), (3, 3, 3)]
        
        self.assertEqual(sorted(s), sorted(mycell.haloCells((2, 2, 2))), "in center")

        sb = sorted(mycell.haloCells((0, 0, 0)))
        s = sorted([ (0, 0, 0), (0, 0, 4), (0, 0, 1), (0, 4, 0), (0, 4, 4), (0, 4, 1), (0, 1, 0), (0, 1, 4),
        (0, 1, 1), (4, 0, 0), (4, 0, 4), (4, 0, 1), (4, 4, 0), (4, 4, 4), (4, 4, 1), (4, 1, 0), (4, 1, 4),
         (4, 1, 1), (1, 0, 0), (1, 0, 4), (1, 0, 1), (1, 4, 0), (1, 4, 4), (1, 4, 1), (1, 1, 0), (1, 1, 4), (1, 1, 1)])
        self.assertEqual(s, sb , "periodic: {0}".format(sb))
        return
    
        
    def testWall(self):
        """Test that we can implement a wall correctly"""
        
        boxDim = [10, 10, 10]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment( filename=self.ch4Car, fragmentType='A' )
        #mycell.addBondType('A:a-A:a')
        
        # Create the wall
        mycell.setWall(XOY=True, XOZ=False, YOZ=False)

        # first add a block in the center
        nblocks = 1
        added = mycell.seed(nblocks, center=True)
        self.assertEqual(added, nblocks, "Failed to add central block")
        
        # Now add a block against the XOY / a wall and confirm it fails
        nblocks = 1
        added = mycell.seed(nblocks, point=[0.0, 5.0, 5.0 ], maxTries=1)
        self.assertEqual(added, 0, "Added block at wall")
        
        #mycell.dump()
        
        return

    def testZipBlocks(self):

        boxDim = [12.0, 12.0, 12.0]
        mycell = Cell(boxDim)

        ch4Car = self.ch4Car
        # mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b2 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b3 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b4 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')

        # b1 in center of cell
        b1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
        mycell.addBlock(b1)
        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Position b2 - all blocks will be positioned around this one
        b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)
        mycell.addBlock(b2)

        # Position b3
        endGroup1 = b2.freeEndGroups()[ 1 ]
        endGroup2 = b3.freeEndGroups()[ 0 ]
        b2.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)
        mycell.addBlock(b3)

        # Position b4
        endGroup1 = b2.freeEndGroups()[ 2 ]
        endGroup2 = b4.freeEndGroups()[ 0 ]
        b2.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)
        mycell.addBlock(b4)

        made = mycell.zipBlocks(bondMargin=0.5, bondAngleMargin=0.5)

        self.assertEqual(made, 3)
        self.assertFalse(self.clashes(mycell))

        return

    def testZipBlocks2(self):

        boxDim = [12.0, 12.0, 12.0]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        b2 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')

        # Position block so that it's aligned along x-axis
        # - use two opposing C-atoms 0 & 3
        b1.alignAtoms(0, 3, [ 1, 0, 0 ])

        b1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Position the block
        b1.positionGrowBlock(endGroup1, endGroup2)

        # and bond
        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock(bond)
        mycell.addBlock(b1)
        
        made = mycell.zipBlocks(bondMargin=5, bondAngleMargin=5)

        self.assertEqual(made, 1)
        self.assertFalse(self.clashes(mycell))
        return

    def testZipBlocks3(self):
        """Bonding with margin"""

        boxDim = [12.0, 12.0, 12.0]
        mycell = Cell(boxDim)

        ch4Car = self.ch4Car
        mycell.libraryAddFragment(filename=ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create block manually
        b1 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')
        b2 = buildingBlock.Block(filePath=ch4Car, fragmentType='A')

        # Align bond along x-axis
        b1.alignAtoms(0, 1, [ 1, 0, 0 ])

        # b1 in center of cell
        b1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
        mycell.addBlock(b1)

        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Position b2
        b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)
        mycell.addBlock(b2)

        # Rotate on y-axis so that it's slightly off
        b2.rotateT([0, 1, 0], math.radians(15))

        made = mycell.zipBlocks(bondMargin=0.5, bondAngleMargin=5)
        self.assertEqual(made, 0)

        made = mycell.zipBlocks(bondMargin=0.5, bondAngleMargin=16)
        self.assertEqual(made, 1)
        self.assertFalse(self.clashes(mycell))

        return


    def testZipClash1(self):
        """Test no clashes"""

        boxDim = [25, 25, 25]
        mycell = Cell(boxDim, doLog=True)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create blocks manually
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        # b2 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )
        # b3 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )

        # Align bond along x-axis
        b1.alignAtoms(0, 3, [ 0, 0, 1 ])
        b1.alignAtoms(0, 3, [ 0, 1, 0 ])
        b1.alignAtoms(0, 3, [ 1, 0, 0 ])

        b2 = b1.copy()
        b3 = b1.copy()
        
        # b1 in center of cell -5 on x-axis
        b1.translateCentroid([ (mycell.dim[0] / 2) - 5, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # b2 in center
        b2.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # rotate so ring facing bond axis
        b2.rotate([0, 1, 0], math.pi / 2, center=b2.centroid())
        
        # b3 + 5 on x
        b3.translateCentroid([ (mycell.dim[0] / 2) + 5, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # Add blocks either side to cell and make sure we can bond
        mycell.addBlock(b1)
        mycell.addBlock(b3)

        made = mycell.zipBlocks(bondMargin=6, bondAngleMargin=1, clashCheck=True)
        self.assertEqual(made, 1)
        self.assertFalse(self.clashes(mycell))
        return
    
    def testZipClash2(self):
        """Test clash when bond is through a benzene ring"""

        boxDim = [25, 25, 25]
        mycell = Cell(boxDim, doLog=True)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create blocks manually
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        # b2 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )
        # b3 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )

        # Align bond along x-axis
        b1.alignAtoms(0, 3, [ 0, 0, 1 ])
        b1.alignAtoms(0, 3, [ 0, 1, 0 ])
        b1.alignAtoms(0, 3, [ 1, 0, 0 ])

        b2 = b1.copy()
        b3 = b1.copy()
        
        # b1 in center of cell -5 on x-axis
        b1.translateCentroid([ (mycell.dim[0] / 2) - 5, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # b2 in center
        b2.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # rotate so ring facing bond axis
        b2.rotate([0, 1, 0], math.pi / 2, center=b2.centroid())
        
        # b3 + 5 on x
        b3.translateCentroid([ (mycell.dim[0] / 2) + 5, mycell.dim[1] / 2, mycell.dim[2] / 2 ])

        # Add blocks either side to cell and make sure we can bond
        mycell.addBlock(b1)
        mycell.addBlock(b3)

        # Add the pesky middle block and see if it now works
        mycell.addBlock(b2)
        # print "ATOM BLOCK ",sorted(b2.atomCell)
                
        made = mycell.zipBlocks(bondMargin=6, bondAngleMargin=1, clashCheck=True, clashDist=1.8)
        self.assertEqual(made, 0)
        self.assertFalse(self.clashes(mycell))
        
        return    

    def testZipClashPBC1(self):
        """Test clash for bond across PBC"""

        boxDim = [25, 25, 25]
        mycell = Cell(boxDim, doLog=True)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create blocks manually
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        # b2 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )
        # b3 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )

        # Align bond along x-axis
        b1.alignAtoms(0, 3, [ 0, 0, 1 ])
        b1.alignAtoms(0, 3, [ 0, 1, 0 ])
        b1.alignAtoms(0, 3, [ 1, 0, 0 ])

        b2 = b1.copy()
        b3 = b1.copy()

        # b1 in center of cell -5 on x-axis
        b1.translateCentroid([ -5, 0, 0 ])

        # b2 on cell bounudary
        b2.translateCentroid([ 0, 0, 0 ])

        # rotate so ring facing bond axis
        b2.rotate([0, 1, 0], math.pi / 2, center=b2.centroid())

        # b3 + 5 on x
        b3.translateCentroid([ 5, 0, 0 ])

        # Add blocks either side to cell and make sure we can bond
        mycell.addBlock(b1)
        mycell.addBlock(b3)
        
        # Add the pesky middle block and see if it now works
        mycell.addBlock(b2)
        # mycell.dump()
        # print "ATOM BLOCK ",sorted(b2.atomCell)

        made = mycell.zipBlocks(bondMargin=6, bondAngleMargin=1, clashCheck=True, clashDist=1.8)
        # mycell.dump()
        self.assertEqual(made, 0)

        return
    
    def testZipClashPBC2(self):
        """Test clash for bond across PBC"""

        boxDim = [25, 25, 25]
        mycell = Cell(boxDim, doLog=True)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create blocks manually
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        # b2 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )
        # b3 = buildingBlock.Block( filePath=self.benzeneCar, fragmentType='A' )

        # Align bond along x-axis
        b1.alignAtoms(0, 3, [ 0, 0, 1 ])
        b1.alignAtoms(0, 3, [ 0, 1, 0 ])
        b1.alignAtoms(0, 3, [ 1, 0, 0 ])

        b2 = b1.copy()
        b3 = b1.copy()

        # b1 in center of cell -5 on x-axis
        b1.translateCentroid([ -5, 0, 0 ])

        # b2 on cell bounudary
        b2.translateCentroid([ 0, 0, 0 ])

        # rotate so ring facing bond axis
        b2.rotate([0, 1, 0], math.pi / 2, center=b2.centroid())

        # b3 + 5 on x
        b3.translateCentroid([ 5, 0, 0 ])

        # Add blocks either side to cell and make sure we can bond
        mycell.addBlock(b1)
        mycell.addBlock(b3)
        
        # Add the pesky middle block and see if it now works
        mycell.addBlock(b2)
        # mycell.dump()
        # print "ATOM BLOCK ",sorted(b2.atomCell)

        made = mycell.zipBlocks(bondMargin=6, bondAngleMargin=1, clashCheck=False)
        # mycell.dump()
        self.assertEqual(made, 1)

        return

    def testWriteCml(self):
        """
        write out cml
        """

        boxDim = [100, 100, 100]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='B')
        mycell.addBondType('B:a-B:a')

        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='C')
        mycell.addBondType('C:a-C:a')

        # Create block manually - this is so we have reproducible results
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        b2 = b1.copy()
        b3 = b1.copy()
        b4 = b1.copy()

        # b1 in center of cell
        b1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
        mycell.addBlock(b1)

        for b in [ b2, b3, b4]:
            endGroup1 = b1.freeEndGroups()[ 0 ]
            endGroup2 = b.freeEndGroups()[ 0 ]

            # Add b2
            b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)

            # bond them
            bond = buildingBlock.Bond(endGroup1, endGroup2)
            b1.bondBlock(bond)

        # Add another block that's not overlapping with the first
        b1 = buildingBlock.Block(filePath=self.ch4Car, fragmentType='B')
        b2 = b1.copy()
        b3 = b1.copy()
        b4 = b1.copy()
        b1.translateCentroid([ 10, 10, 10])
        mycell.addBlock(b1)

        for b in [ b2, b3, b4]:
            endGroup1 = b1.freeEndGroups()[ 0 ]
            endGroup2 = b.freeEndGroups()[ 0 ]

            # Add b2
            b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)

            # bond them
            bond = buildingBlock.Bond(endGroup1, endGroup2)
            b1.bondBlock(bond)

        # Add another block that's not overlapping with the others
        b1 = buildingBlock.Block(filePath=self.ch4Car, fragmentType='C')
        b2 = b1.copy()
        b3 = b1.copy()
        b4 = b1.copy()
        b1.translateCentroid([ 90, 90, 90])
        mycell.addBlock(b1)

        for b in [ b2, b3, b4]:
            endGroup1 = b1.freeEndGroups()[ 0 ]
            endGroup2 = b.freeEndGroups()[ 0 ]

            # Add b2
            b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)

            # bond them
            bond = buildingBlock.Bond(endGroup1, endGroup2)
            b1.bondBlock(bond)


        fname = "test.cml"
        mycell.writeCml(fname, periodic=False, rigidBody=True)
        # Test is same as reference
        with open(fname) as f:
            test = f.readlines()
        with open(os.path.join(AMBUILD_DIR, "tests", "testCellRigid.cml")) as f:
            ref = f.readlines()

        self.assertEqual(test, ref, "cml compare rigid")

        fname = "test.cml"
        mycell.writeCml(fname, periodic=False, rigidBody=False)
        # Test is same as reference
        with open(fname) as f:
            test = f.readlines()
        with open(os.path.join(AMBUILD_DIR, "tests", "testCellAll.cml")) as f:
            ref = f.readlines()

        self.assertEqual(test, ref, "cml compare all")
        os.unlink(fname)
        return


    def testWriteHoomdblue(self):
        """
        write out hoomdblue xml
        """
        import hoomdblue

        boxDim = [20, 20, 20]
        mycell = Cell(boxDim)

        mycell.libraryAddFragment(filename=self.benzeneCar, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        # Create block manually - this is so we have reproducible results
        b1 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')
        b2 = buildingBlock.Block(filePath=self.benzeneCar, fragmentType='A')

        # b1 in center of cell
        b1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
        endGroup1 = b1.freeEndGroups()[ 0 ]
        endGroup2 = b2.freeEndGroups()[ 0 ]

        # Bond b2 to it
        b1.positionGrowBlock(endGroup1, endGroup2, dihedral=math.pi)
        bond = buildingBlock.Bond(endGroup1, endGroup2)
        b1.bondBlock(bond)
        mycell.addBlock(b1)

        initcoords = []
        for block in mycell.blocks.itervalues():
            for c in block.iterCoord():
                initcoords.append(c)

        xmlFilename = "testWriteHoomdblue.xml"
        data = mycell.dataDict(periodic=True, center=True, rigidBody=True)
        o = opt.HoomdOptimiser()
        o.writeXml(data,
                   xmlFilename=xmlFilename,
                   rigidBody=True,
                   doDihedral=True,
                   doImproper=False,
                   doCharges=True,
                   )

        # Test what we've written out matches the reference file
        with open(xmlFilename) as f:
            test = f.readlines()
        with open(os.path.join(AMBUILD_DIR, "tests", xmlFilename)) as f:
            ref = f.readlines()
        self.assertEqual(test, ref, "xml compare")

        # Init the sytem from the file
        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        system = hoomdblue.init.read_xml(filename=xmlFilename)

        # Read it back in to make sure we get the same values
        mycell.fromHoomdblueSystem(system)
        finalcoords = []
        for block in mycell.blocks.itervalues():
            for c in block.iterCoord():
                finalcoords.append(c)

        self.assertTrue(all(map(lambda x : numpy.allclose(x[0], x[1]), zip(initcoords, finalcoords))),
                         "coords don't match")

        # ok = opt.HoomdOptimiser().optimiseGeometry( xmlFilename,doDihedral=doDihedral)

        os.unlink(xmlFilename)

        return

    def XtestWriteHoomdblueMinCell(self):
        """
        write out hoomdblue xml
        """

        import hoomdblue

        # Create a large cell and populate it with a block
        CELLA = CELLB = CELLC = 100
        mycell = Cell()
        mycell.cellAxis(A=CELLA, B=CELLB, C=CELLC)
        mycell.libraryAddFragment(filename=self.ch4Car, fragmentType='A')
        mycell.addBondType('A:a-A:a')

        mycell.seed(1)
        mycell.growBlocks(3, 'A')

        # Get the block and put it in the centre of the cell
        block = mycell.blocks[ mycell.blocks.keys()[0] ]
        block.translateCentroid([ CELLA / 2, CELLB / 2, CELLC / 2 ])

        initcoords = []
        for block in mycell.blocks.itervalues():
            initcoords += block._coords

        filename = "testWriteHoomdblue.xml"
        natoms = mycell.writeHoomdXml(xmlFilename=filename)

        # Init the sytem from the file
        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        system = hoomdblue.init.read_xml(filename=filename)

        # Read it back in to make sure we get the same values
        mycell.fromHoomdblueSystem(system, minCell=True)

        finalcoords = []
        for block in mycell.blocks.itervalues():
            finalcoords += block._coords

        self.assertEqual(initcoords, finalcoords, "coords don't match")

        os.unlink(filename)

        return

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
