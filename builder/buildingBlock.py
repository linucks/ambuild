'''
Created on Jan 15, 2013

@author: abbietrewin

Things to look at:
http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co

http://en.wikipedia.org/wiki/Periodic_boundary_conditions

http://mail.scipy.org/pipermail/scipy-dev/2012-March/017177.html
https://groups.google.com/forum/?fromgroups=#!topic/scipy-user/P6k8LEo30ws
https://github.com/patvarilly/periodic_kdtree
'''

# Python imports
import copy
import itertools
import os
import math
import random
import unittest

# external imports
import numpy

# local imports
import fragment
import util

class Bond(object):
    """An object to hold all info on a bond
    """
    def __init__(self, endGroup1, endGroup2):
        self.endGroup1 = endGroup1
        self.endGroup2 = endGroup2
        return
    def __str__(self):
        """List the data attributes of this object"""
        s = "Bond {0}: {1}:{2}:{3}-{4} -> {5}:{6}:{7}-{8}".format(id(self),
                                                   self.endGroup1.block().id, id(self.endGroup1.fragment),
                                                   id(self.endGroup1), self.endGroup1.blockEndGroupIdx,
                                                   self.endGroup2.block().id, id(self.endGroup2.fragment),
                                                   id(self.endGroup2), self.endGroup2.blockEndGroupIdx)
        return s

class Block(object):
    '''
    Structure is:
    * a list of fragments
    * each fragment my contain endGroups
    * each endGroup handles the change to atoms required on bonding
    * fragments may be bonded to each other with bonds
    * bonds just point to two endGroups that are bonded to each other
    '''

    def __init__(self, filePath=None, fragmentType=None, initFragment=None,):
        '''
        Constructor
        '''

        # Need to change so cannot create block withough fragmentType
        if filePath:
            assert os.path.isfile(filePath) and fragmentType
            initFragment = fragment.Fragment(filePath, fragmentType)

        # List of the fragments contained in this one
        self.fragments = []
        if initFragment: self.fragments.append(initFragment )

        # List of bond objects between blocks
        self._blockBonds = []

        # List of tuples of atoms that are bonded
        self._bonds = []

        # List of the bonds within fragments as a tuple (fragmentType, bond)
        self._bondsByFragmentType = []

        # List of which atom is bonded to which
        self._bondedToAtom = []

        # list of tuples of ( idFrag, idxData )
        self._dataMap = []
        self._bodies = []  # List of which body in the block each atom belongs too - frags can contain multiple bodies

        # The list of atoms that are endGroups and their corresponding angleAtoms
        self._freeEndGroups = {}
        self._numFreeEndGroups = 0
        self._endGroupType2EndGroups = {}

        # Flag to indicate if block has changed and parameters (such as centerOfMass)
        # need to be changed
        self._changed = True

        # Below need to be updated when we move etc
        self._centroid = numpy.zeros(3)
        self._centerOfMass = numpy.zeros(3)
        self._maxAtomRadius = -1
        self._radius = None
        self._blockMass = 0
        self.id = id(self)

        if self.fragments: self._update()
        return

    def alignAtoms(self, atom1Idx, atom2Idx, refVector):
        """Move molecule so two atoms are aligned along refVector"""

        atom1 = self.coord(atom1Idx)
        atom2 = self.coord(atom2Idx)
        if isinstance(refVector, list):  # Should check if numpy array
            refVector = numpy.array(refVector, dtype=numpy.float64)

        return self.alignVector(atom1, atom2, refVector)

    def alignVector(self, pos1, pos2, refVector):
        """
        Align this block, so that the axis defined by the two atoms is aligned with
        the refVector

        pos1 is the coordinate of atom1 and pos2 the coordinate of atom2
        """

        # Move so that pos1 is at origin so the vector of pos2 can be
        # aligned with the refVector

        # Shift block so angleAtom at center,
        self.translate(-pos1)

        # Check neither is zero
        if numpy.array_equal(refVector, [0, 0, 0]) or numpy.array_equal(pos2, [0, 0, 0]):
            raise RuntimeError, "alignBlock - one of the vectors is zero!\nrefVector: {0} endGroup: {1}".format(refVector, pos2)

        # Makes no sense if they are already aligned
        if numpy.array_equal(pos2 / numpy.linalg.norm(pos2), refVector / numpy.linalg.norm(refVector)):
            print "alignBlock - block already aligned along vector. May not be a problem, but you should know..."
            self.translate(pos1)  # NEW - put it back
            return

        # print "alignBlock BEFORE: {0} | {1}".format( endGroup, refVector )

        # Calculate normalised cross product to find an axis orthogonal
        # to both that we can rotate about
        cross = numpy.cross(refVector, pos2)

        if numpy.array_equal(cross, [0, 0, 0]):
            # They must be already aligned but anti-parallel, so we flip
            print "alignBlock - vectors are anti-parallel so flipping"
            self.flip(refVector)
            self.translate(pos1)  # NEW - put it back
        else:

            # Find angle
            angle = util.vectorAngle(refVector, pos2)

            # Normalised cross to rotate about
            ncross = cross / numpy.linalg.norm(cross)

            # Rotate
            self.rotate(ncross, angle)

        # Now shift back
        self.translate(pos1)

        # print "alignBlock AFTER: {0} | {1}".format( self._coord( idxAtom ), refVector )
        return

    def atomBonded1(self, idxAtom):
        """Return the indices of all atoms directly bonded to idxAtom"""
        # We return a copy or else modifying the list changes the actual bond list
        return copy.copy(self._bondedToAtom[ idxAtom ])

    def atomBonded2(self, idxAtom):
        """Return the indices of all atoms bonded by <= 2 bonds to idxAtom"""
        bonded = copy.copy(self.atomBonded1(idxAtom))
        for a1 in list(bonded):  # Copy to list so we're not changing while looping thru it
            bonded.update(self.atomBonded1(a1))
        return bonded

    def atomBonded3(self, idxAtom):
        """Return the indices of all atoms bonded by <= 3 bonds to idxAtom"""
        bonded = copy.copy(self.atomBonded2(idxAtom))
        for a1 in list(bonded):  # Copy to list so we're not changing while looping thru it
            bonded.update(self.atomBonded1(a1))
        return bonded

    def body(self, idxAtom):
        return self._bodies[idxAtom]

    def charge(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag.charge(idxData)

    def coord(self, idxAtom, coord=None):
        """Get and set coordinate in external indices"""
        frag, idxData = self._dataMap[idxAtom]
        if coord is not None:
            if isinstance(coord, list):
                coord = numpy.array(coord)
            frag.coord(idxData, coord)
        else:
            return frag.coord(idxData)

    def fragment(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag

    def fragmentType(self, idxAtom):
        """The type of the fragment that this atom belongs to."""
        frag, idxData = self._dataMap[idxAtom]
        return frag.fragmentType

    def label(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag.label(idxData)

    def mass(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag.mass(idxData)

    def radius(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag.radius(idxData)

    def symbol(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag.symbol(idxData)

    def type(self, idxAtom):
        frag, idxData = self._dataMap[idxAtom]
        return frag.type(idxData)

    def anglesAndDihedrals(self):
        """
        Borrowed from openMM
        """
        # Make a list of all unique angles
        uniqueAngles = set()
        for atom1, atom2 in self._bonds:
            for atom in self._bondedToAtom[atom1]:
                if atom != atom2:
                    if atom < atom2:
                        uniqueAngles.add((atom, atom1, atom2))
                    else:
                        uniqueAngles.add((atom2, atom1, atom))
            for atom in self._bondedToAtom[atom2]:
                if atom != atom1:
                    if atom > atom1:
                        uniqueAngles.add((atom1, atom2, atom))
                    else:
                        uniqueAngles.add((atom, atom2, atom1))

        # Sort and reindex
        angles = sorted(list(uniqueAngles))

        # Make a list of all unique proper torsions

        uniquePropers = set()
        for angle in angles:
            for atom in self._bondedToAtom[angle[0]]:
                if atom != angle[1]:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in self._bondedToAtom[angle[2]]:
                if atom != angle[1]:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))

        propers = sorted(list(uniquePropers))

        # Make a list of all unique improper torsions
        impropers = []
        # for atom in range(len(bondedToAtom)):
        for atom in range(len(self._dataMap)):
            bondedTo = self._bondedToAtom[atom]
            if len(bondedTo) > 2:
                for subset in itertools.combinations(bondedTo, 3):
                    impropers.append((atom, subset[0], subset[1], subset[2]))
        return angles, propers, impropers

    def atomEndGroups(self, idxAtom):
        """Return a list of the endGroup objects for this atom
        Index in external coordinates
        """
        try:
            return self._freeEndGroups[idxAtom]
        except KeyError:
            return False

    def bonds(self):
        """All bonds in external indices"""
        return self._bonds

    def blockBonds(self):
        """External indices"""
        return [ (b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx) for b in self._blockBonds ]

    def bondBlock(self, bond):
        """ Add newBlock to this one
        """
        assert bond.endGroup1 != bond.endGroup2
        assert bond.endGroup1.block() == self
        assert bond.endGroup1.free()
        assert bond.endGroup2.free()
        assert not bond.endGroup1.isBonded()
        assert not bond.endGroup2.isBonded()
        assert not bond.endGroup1.saturated()
        assert not bond.endGroup2.saturated()
        assert not bond.endGroup1.fragment == bond.endGroup2.fragment

        # Mark both endGroups as used
        bond.endGroup1.setBonded()
        bond.endGroup2.setBonded()

        # Tried optimising this by passing in the bond to update and only updating those fragments/
        # endGroups that had changed but it rapidly got rather complicated so we keep to a simple
        # update and add the data for the new block here
        # Append fragments and bonds of other block to this one
        if bond.endGroup1.block() != bond.endGroup2.block():
            self.fragments += bond.endGroup2.block().fragments
            self._blockBonds += bond.endGroup2.block()._blockBonds

        # add the new bond
        self._blockBonds.append(bond)

        return self._update()
    
    def blockRadius(self):
        if self._changed: self._calcProperties()
        return self._radius
    
    def _calcCenters(self):
        sumG = numpy.zeros(3)
        sumM = numpy.zeros(3)
        totalMass = 0.0
        for f in self.fragments:
            mass = f.totalMass()
            totalMass += mass
            sumG += f.centroid()
            sumM += mass * f.centroid()

        self._centroid = sumG / len(self.fragments)
        self._centerOfMass = sumM / totalMass
        return

    def _calcRadius(self):
        """
        Calculate a simple size metric of the block so that we can screen for whether
        two blocks are within touching distance
        Assumes centroid already calculated
        """

        distances = []
        self._maxAtomRadius = 0.0
        for f in self.fragments:
            for coord in f.iterCoord():
                distances.append(util.distance(self._centroid, coord))
                self._maxAtomRadius = max(f.maxAtomRadius(), self._maxAtomRadius)

        assert distances
        imax = numpy.argmax(distances)
        dist = distances[ imax ]

        # Add on the radius of the largest atom
        self._radius = dist + self.maxAtomRadius()
        return

    def _calcProperties(self):
        self._calcCenters()
        self._calcRadius()
        self._changed = False
        return

    def centerOfMass(self):
        """
        Return or calculate the center of mass for this building block.
        """
        #print "centerOfMass changed ", self._changed
        if self._changed:
            self._calcProperties()
        return self._centerOfMass

    def centroid(self):
        """
        Return or calculate the center of geometry for the building block.
        """
        if self._changed: self._calcProperties()
        return self._centroid

    def copy(self):
        """Return a copy of ourselves."""
        new = copy.deepcopy(self)
        new.id = id(new)
        return new

    def dataByFragment(self, fragmentType):
        """Return the data for a specific fragmentType within the block"""

        coords = []
        symbols = []

        atomCount = 0
        fbondRen = {}
        for i in range(len(self._dataMap)):
            if  self.fragmentType(i) == fragmentType:
                fbondRen[i] = atomCount
                coords.append(self.coord(i))
                symbols.append(self.symbol(i))
                atomCount += 1

        bonds = [ (fbondRen[b1], fbondRen[b2]) \
                 for ftype, (b1, b2) in self._bondsByFragmentType if ftype == fragmentType ]

        return coords, symbols, bonds

    def deleteBond(self, bond):
        # See if breaking the bond separates the block into two separate blocks

        # Create dictionary of which fragments are bonded to which fragments
        bondedToFragment = {}
        for b in self._blockBonds:
            if b.endGroup1.fragment not in bondedToFragment:
                bondedToFragment[b.endGroup1.fragment] = set()
            if b.endGroup2.fragment not in bondedToFragment:
                bondedToFragment[b.endGroup2.fragment] = set()
            if b != bond:
#             if (b.endGroup1.blockEndGroupIdx==idxAtom1 and b.endGroup2.blockEndGroupIdx==idxAtom2) or \
#                 (b.endGroup1.blockEndGroupIdx==idxAtom2 and b.endGroup2.blockEndGroupIdx==idxAtom1):
#                 bond=b
                bondedToFragment[b.endGroup1.fragment].add(b.endGroup2.fragment)
                bondedToFragment[b.endGroup2.fragment].add(b.endGroup1.fragment)

        # Take the two fragments on either side of the bond
        f1 = bond.endGroup1.fragment
        f2 = bond.endGroup2.fragment

        # Sets of which fragments can be reached from each fragment
        f1set = set()
        f2set = set()

        def addFragments(startf, f1set):
            # Trundle through the bond topology for f1 and set True for all fragments we reach
            for f in bondedToFragment[startf]:
                if f not in f1set and f != startf:
                    f1set.add(f)
                    f1set.update(addFragments(f, f1set))
            return f1set

        f1set = addFragments(f1, f1set)
        f2set = addFragments(f2, f2set)

        if bool(f1set.intersection(f2set)):
            # Fragments in common with both, so just delete the bond
            #print "Block remains contiguous"
            # Need to unmask the fragment atoms
            bond.endGroup1.unBond()
            bond.endGroup2.unBond()

            # Now delete the bond from the block
            self._blockBonds.remove(bond)
            self._update()
            return None
        
        # Breaking the bond splits the block in two, so we separate the two fragments, keeep the largest
        # for ourselves and return the new block. We set f1 and f1set to be the biggest
        if f2set > f1set:
            tmp = f1set
            f1set = f2set
            f2set = tmp
            tmp = f1
            f1 = f2
            f2 = tmp
            
        # How to partition the bonds?
        f1bonds = set()
        f2bonds = set()
        for b in self._blockBonds:
            if b == bond: continue
            if b.endGroup1.fragment in f1set and b.endGroup2.fragment in f1set:
                f1bonds.add(b)
            elif b.endGroup1.fragment in f2set and b.endGroup2.fragment in f2set:
                f2bonds.add(b)
            else: raise RuntimeError,"Bond not in set: {0}".format(b)
            
        bond.endGroup1.unBond()
        bond.endGroup2.unBond()
        
        # Create a new block with the smaller fragments
        newBlock = Block()
        newBlock.fragments = [f2]
        if len(f2set): newBlock.fragments += list(f2set)
        if len(f2bonds): newBlock._blockBonds = list(f2bonds)
        newBlock._update()
        
        # Update our list of bonds 
        self.fragments = [f1]
        if len(f1set): self.fragments += list(f1set)
        if len(f1bonds): self._blockBonds = list(f1bonds)
        self._update()

        return newBlock
    
    def deleteFragment(self, frag):
        if frag not in self.fragments:
            raise RuntimeError,"Cannot find fragment {0} in block".format(frag.fragmentType)
        
        if len(self.fragments) == 1:
            # This is the only fragment in the block, so we need to delete ourselves!
            # Return False to indicate this is the case
            return False
            
        # Get the list of bonds that this fragment makes
        dbonds = [ b for b in self._blockBonds if b.endGroup1.fragment == frag and b.endGroup1.fragment == frag ]
        
        assert len(dbonds),"Fragment is not involved in any bonds!"
        
        # Loop through deleting all bonds and returning any blocks created
        # The fragment we seek to delete might end up in one of the other blocks if it happens to be bigger then us
        # so we need to check for this and remove the fragment from that block
        blocks = []
        removed=False
        for bond in dbonds:
            block = self.deleteBond(bond)
            if block:
                if frag in block.fragments:
                    if len(block.fragments) > 1:
                        block.fragments.remove(frag)
                        block._update()
                    else:
                        block = None
                    removed=True
                if block: blocks.append(block)
            
        # Now we might need to delete the fragment from ourself and update
        if not removed:
            if len(self.fragments) > 1:
                self.fragments.remove(frag)
                self._update()
            else:
                raise RuntimeError,"FIX FOR WHEN THERE IS ONLY ONE BLOCK"
        
        return blocks

    def dihedrals(self, atom1Idx, atom2Idx, bondOnly=False):
        """Return a list of all the dihedrals around these two bonded atoms
        input & output in internacl coodrdinates

        This needs more work to check for when things are looped back and bonded to each other
        """

        # Create a list of lists of all the atoms that each endGroup is bonded to - excluding the other endGroup
        # Add all dihedrals on each side of the bond - both bond atoms plus 1, 2 connected to the endGrup
        # Add all dihedrals across the bond - both endGroups plus 1 either side

        # Create list of what's bonded to atom1 - we exclude anything that loops back on itself
        atom1Bonded = {}
        for a1 in self.atomBonded1(atom1Idx):
            if a1 == atom2Idx:
                continue
            atom1Bonded[ a1 ] = []
            for a2 in self.atomBonded1(a1):
                if a2 == atom1Idx:
                    continue
                assert not a2 == atom2Idx, "Dihedral atom loops back onto bond!"
                atom1Bonded[ a1 ].append(a2)

        # Create list of what's bonded to atom2 - we exclude anything that loops back on itself
        atom2Bonded = {}
        for a1 in self.atomBonded1(atom2Idx):
            if a1 == atom1Idx:
                continue
            atom2Bonded[ a1 ] = []
            for a2 in self.atomBonded1(a1):
                if a2 == atom2Idx:
                    continue
                assert not a2 == atom1Idx, "Dihedral atom loops back onto bond!"
                atom2Bonded[ a1 ].append(a2)

        dindices = []
        if not bondOnly:
            # Add all dihedrals on endGroup1's side of the bond
            for a3 in atom1Bonded:
                for a4 in atom1Bonded[ a3 ]:
                    dindices.append((atom2Idx, atom1Idx, a3, a4))

            # Add all dihedrals on endGroup2's side of the bond
            for a3 in atom2Bonded:
                for a4 in atom2Bonded[ a3 ]:
                    dindices.append((atom1Idx, atom2Idx, a3, a4))

        # Now add dihedrals across the bond
        for a1 in atom1Bonded:
            for a2 in atom2Bonded:
                dindices.append((a1, atom1Idx, atom2Idx, a2))

        return dindices

    def flip(self, fvector):
        """Rotate perpendicular to fvector so we  facing the opposite way along the fvector
        """

        # Find vector perpendicular to the bond axis
        # Dot product needs to be 0
        # xu + yv + zw = 0 - set u and v to 1, so w = (x + y)/z
        # vector is 1, 1, w
        w = -1.0 * (fvector[0] + fvector[1]) / fvector[2]
        orth = numpy.array([1.0, 1.0, w])

        # Find axis that we can rotate about
        rotAxis = numpy.cross(fvector, orth)

        # Rotate by 180
        self.rotate(rotAxis, numpy.pi)

        return
    
    def endGroupConfig(self, fragmentType):
        """See endGroupConfig in cell"""
        c = []
        for f in self.fragments:
            if f.fragmentType == fragmentType:
                # calculate number of bonded endGroups
                c.append(len(f.endGroups()) - len(f.freeEndGroups()))
        if len(c):
            return ":".join([str(n) for n in c])
        else:
            return None

    def fragmentBonds(self):
        """Return a list of all the internal fragment bonds for this block
        This excludes the bonds between fragments
        """
        # loop through all fragments and get a list of which atoms are internally bonded
        # Then go through the data map and map these to the 'external' atom indices
        fbonds = []
        for fragment in self.fragments:
            for b1, b2 in fragment.bonds():
                fbonds.append((b1 + fragment.blockIdx, b2 + fragment.blockIdx))
        return fbonds

    def fragmentTypeDict(self):
        """A dictionary with the number of the different types of fragment we contain"""
        return self._fragmentTypeDict

    def freeEndGroups(self):
        # http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
        return  [item for sublist in self._freeEndGroups.values() for item in sublist]

    def freeEndGroupsFromTypes(self, endGroupTypes):
        # Make sure we have a list to check against
        if isinstance(endGroupTypes, str):
            endGroupTypes = [ endGroupTypes ]
        endGroups = []
        for t in endGroupTypes:
            endGroups += self._endGroupType2EndGroups[ t ]
        return endGroups

    def freeEndGroupTypes(self):
        """Return a list of the fragmentTypes for the available endGroups"""
        return self._endGroupType2EndGroups.keys()

    def isEndGroup(self, idxAtom):
        """Return True if this atom is a free endGroup
        """
        # No need to do conversion as atomEndGroups is external interface
        if self.atomEndGroups(idxAtom):
            return True
        return False

    def iterCoord(self):
        """Generator to return the coordinates"""
        for i in range(len(self._dataMap)):
            yield self.coord(i)
        return

    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius

    def newBondPosition(self, endGroup, symbol):
        """Return the position where a bond to an atom of type 'symbol'
        would be placed if bonding to the target endgroup
         I'm sure this algorithm is clunky in the extreme...
        """

        targetEndGroup = self.coord(endGroup.blockEndGroupIdx)
        targetSymbol = self.symbol(endGroup.blockEndGroupIdx)
        targetCapAtom = self.coord(endGroup.blockCapIdx)

        # Get the bond length between these two atoms
        bondLength = util.bondLength(targetSymbol, symbol)

        # Find unit vector pointing from targetAngleAtom to targetEndGroup

        # vector from targetEndgroup to targetCapAtom
        # v1 = targetEndGroup - targetCapAtom
        v1 = targetCapAtom - targetEndGroup

        # Now get unit vector
        uv = v1 / numpy.linalg.norm(v1)

        # Multiply unit vector by bond length to get the component to add on
        newPosition = targetEndGroup + (uv * bondLength)

        return newPosition

    def numFreeEndGroups(self):
        return self._numFreeEndGroups

    def numAtoms(self):
        """Number of atoms visible externally"""
        return len(self._dataMap)

    def blockMass(self):
        return self._blockMass

    def positionGrowBlock(self, endGroup, growEndGroup, dihedral=None):
        """
        Position growBlock so it can bond to us
        """

        # The vector we want to align along is the vector from the endGroup
        # to the capAtom
        endGroupAtom = self.coord(endGroup.blockEndGroupIdx)
        capAtom = self.coord(endGroup.blockCapIdx)
        refVector = endGroupAtom - capAtom

        growBlock = growEndGroup.block()

        # get the coord where the next block should bond
        # symbol of endGroup tells us the sort of bond we are making which determines
        # the bond length
        symbol = growBlock.symbol(growEndGroup.blockEndGroupIdx)
        bondPos = self.newBondPosition(endGroup, symbol)
        # print "got bondPos for {0}: {1}".format( symbol, bondPos )

        # Align along the staticBlock bond
        growBlock.alignAtoms(growEndGroup.blockEndGroupIdx,
                              growEndGroup.blockCapIdx,
                              refVector)

        # Now need to place the endGroup at the bond coord to do this now we just add the bondPos
        growBlock.translate(bondPos)

        # We need to rotate to adhere to the specified dihedral angle
        if dihedral is not None:
            assert 0 <= dihedral < math.pi * 2
            # Get current angle
            current = util.dihedral(self.coord(endGroup.blockDihedralIdx),
                                     self.coord(endGroup.blockEndGroupIdx),
                                     growBlock.coord(growEndGroup.blockEndGroupIdx),
                                     growBlock.coord(growEndGroup.blockDihedralIdx))

            assert endGroup.blockDihedralIdx != -1 and growEndGroup.blockDihedralIdx != -1, \
            "Need to have specified dihedrals as 3rd column in ambi file first!"

            # Find how much we need to rotate by
            angle = dihedral - current
            if angle != 0:
                # Need to rotate so get the axis to rotate about
                axis = self.coord(endGroup.blockEndGroupIdx) - growBlock.coord(growEndGroup.blockEndGroupIdx)
                growBlock.rotate(axis, angle, center=self.coord(endGroup.blockEndGroupIdx))

        return

    def randomEndGroup(self, endGroupTypes=None):
        """Return a random free endGroup in the block"""

        if endGroupTypes == None:
            # We pick a random endGroup
            endGroup = random.choice(self.freeEndGroups())
        else:
            # Make sure we have a list to check against
            if isinstance(endGroupTypes, str):
                endGroupTypes = [ endGroupTypes ]

            # See if any in common
            common = frozenset(self.freeEndGroupTypes()).intersection(frozenset(endGroupTypes))
            if not bool(common):
                raise RuntimeError, "Cannot find {0} in available types {1}".format(endGroupTypes, self.freeEndGroupTypes())

            # We can definitely return something so pick a random fragment type and get a random endGroup
            ftype = random.choice(list(common))
            endGroup = random.choice(self.freeEndGroupsFromTypes(endGroupTypes=[ ftype ]))

        return endGroup

    def randomRotate(self, origin=[ 0, 0, 0 ], atOrigin=False):
        """Randomly rotate a block.

         Args:
         atOrigin -- flag to indicate if the block is already positioned at the origin
        """

        if not atOrigin:
            position = self.centroid()
            self.translateCentroid(origin)


        xAxis = [ 1, 0, 0 ]
        yAxis = [ 0, 1, 0 ]
        zAxis = [ 0, 0, 1 ]

        angle = random.uniform(0, 2 * numpy.pi)
        self.rotate(xAxis, angle)

        angle = random.uniform(0, 2 * numpy.pi)
        self.rotate(yAxis, angle)

        angle = random.uniform(0, 2 * numpy.pi)
        self.rotate(zAxis, angle)

        if not atOrigin:
            self.translateCentroid(position)

    def randomBlock(self):
        """Return a random block"""
        return random.choice(self.allBlocks)

    def rotate(self, axis, angle, center=None):

        if center is None:
            center = numpy.array([ 0, 0, 0 ])

        rotationMatrix = util.rotation_matrix(axis, angle)

        for f in self.fragments:
            f.rotate(rotationMatrix, center)
        return

    def rotateT(self, axis, angle, center=None):
        """Rotation with translation to center"""

        position = self.centroid()
        origin = numpy.array([ 0, 0, 0 ])
        self.translateCentroid(origin)

        rotationMatrix = util.rotation_matrix(axis, angle)
        for f in self.fragments:
            f.rotate(rotationMatrix, origin)

        self.translateCentroid(position)
        return

    def translate(self, tvector):
        """ translate the molecule by the given vector"""

        # CHANGE SO WE CHECK IF IS A NUMPY ARRAY
        if isinstance(tvector, list):
            tvector = numpy.array(tvector)

        # Loop through each fragment and translate each in turn
        for f in self.fragments:
            f.translate(tvector)

        self._changed = True
        return

    def translateCentroid(self, position):
        """Translate the molecule so the center of geometry moves
        to the given coord"""
        self.translate(position - self.centroid())
        return

    def _update(self):
        """Set the list of _endGroups & update data for new block
        """
        #
        # Now build up the dataMap listing where each fragment starts in the block and linking the
        # overall block atom index to the fragment and fragment atom index
        #
        self._dataMap = []
        self._bodies = []
        self._blockMass = 0
        self._fragmentTypeDict = {}
        bodyCount = -1
        lastBody = 0
        for fragment in self.fragments:

            # Set the block
            fragment.block = self

            # Increment body count for each fragment
            bodyCount += 1

            # Count the number of each type of fragment in the block (see Analyse)
            t = fragment.fragmentType
            if t not in self._fragmentTypeDict:
                self._fragmentTypeDict[ t ] = 1
            else:
                self._fragmentTypeDict[ t ] += 1

            fragment.blockIdx = len(self._dataMap)  # Mark where the data starts in the block
            for i in xrange(fragment.numAtoms()):
                self._dataMap.append((fragment, i))

                # Bring up the bodies
                b = fragment.body(i)
                if b != lastBody:
                    bodyCount += 1
                    lastBody = b
                self._bodies.append(bodyCount)
                self._blockMass += fragment.mass(i)
        #
        # Have dataMap so now update the endGroup information
        self._numFreeEndGroups = 0
        self._freeEndGroups = {}
        self._endGroupType2EndGroups = {}
        for fragment in self.fragments:
            for i, endGroup in enumerate(fragment.endGroups()):
                assert endGroup.fragment == fragment
                # assert id(endGroup) == id( fragment._endGroups[ i ] ) # no longer valid as we only return free ones
                # Set the block index - we sort out the others after we've done bonding
                endGroup.updateEndGroupIndex()
                if endGroup.free():
                    # Add to the list of all free endGroups
                    try:
                        self._freeEndGroups[ endGroup.blockEndGroupIdx ].append(endGroup)
                    except KeyError:
                        self._freeEndGroups[ endGroup.blockEndGroupIdx ] = [ endGroup ]
                    self._numFreeEndGroups += 1
                    # Now add to the type list
                    if endGroup.type() not in self._endGroupType2EndGroups:
                        self._endGroupType2EndGroups[ endGroup.type() ] = []
                    self._endGroupType2EndGroups[ endGroup.type() ].append(endGroup)

        # Now need to create the list of all bonds throughout the block
        self._bonds = []
        self._bondsByFragmentType = []
        # First all bonds within the fragments
        for fragment in self.fragments:
            for b1, b2 in fragment.bonds():
                # Convert to block indices
                b1 = b1 + fragment.blockIdx
                b2 = b2 + fragment.blockIdx
                self._bonds.append((b1, b2))
                self._bondsByFragmentType.append((fragment.fragmentType, (b1, b2)))

        # Then all bonds between fragments
        cap2EndGroup = {}
        for b in self._blockBonds:
            self._bonds.append((b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx))
            # Need to map cap atoms to their bonded counterparts so we can look these up when we
            # fix the endGroup indices - we map the fragment, fragmentIndex to the corresponding block index
            # This is somewhat untidy as we use the internal fragment index here - which really should be hidden
            cap2EndGroup[ (b.endGroup1.fragment, b.endGroup1.fragmentCapIdx) ] = b.endGroup2.blockEndGroupIdx
            cap2EndGroup[ (b.endGroup2.fragment, b.endGroup2.fragmentCapIdx) ] = b.endGroup1.blockEndGroupIdx

        # Now create the list of which atoms are bonded to which
        self._bondedToAtom = []
        for i in range(len(self._dataMap)):
            self._bondedToAtom.append(set())
        for b1, b2 in self._bonds:
            self._bondedToAtom[b1].add(b2)
            self._bondedToAtom[b2].add(b1)

        # Finally update the ancillary blockIndices for the endGroups - we need the bonding to have been done
        # as some of the atoms will now be defined by atoms in other fragments
        for fragment in self.fragments:
            for endGroup in fragment.endGroups():
                endGroup.updateAncillaryIndices(cap2EndGroup)

        # Recalculate the data for this new block
        self._calcProperties()

        return

    def writeCml(self, cmlFilename, cell=None):
        atomTypes = []
        coords = []
        symbols = []
        for i, coord in enumerate(self.iterCoord()):
            coords.append(coord)
            symbols.append(self.symbol(i))
            atomTypes.append(self.type(i))

        cmlFilename = util.writeCml(cmlFilename,
                                    coords,
                                    symbols,
                                    bonds=self.bonds(),
                                    atomTypes=atomTypes,
                                    cell=cell,
                                    pruneBonds=False)

        print "wrote block CML file ", cmlFilename
        return

    def writeXyz(self, name, cell=None):
        symbols = []
        coords = []
        for i, c in enumerate(self.iterCoord()):
            coords.append(c)
            symbols.append(self.symbol(i))
        fpath = util.writeXyz(name, coords, symbols, cell=cell)
        print "Wrote block xyz file: {0}".format(fpath)
        return

    def __str__(self):
        """
        Return a string representation of the block
        """

        mystr = "Block: {0}\n".format(self.id)

        mystr += "Num fragments: {0}\n".format(len(self.fragments))
        mystr += "endGroups: {0}\n".format([ str(e) for e in self.freeEndGroups() ])
        # mystr += "_bondedCapAtoms: {0}\n".format( self._bondedCapIdxs )
        mystr += "Block bonds: {0}\n".format(self._blockBonds)
        mystr += "bonds: {0}\n".format(self._bonds)
        mystr += "bonded: {0}\n".format(self._bondedToAtom)

        for i, c in enumerate(self.iterCoord()):
            # mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self._labels[i], c[0], c[1], c[2])
            mystr += "{0}  {1:5} {2:0< 15} {3:0< 15} {4:0< 15} \n".format(i, self.label(i), c[0], c[1], c[2])

        return mystr

class TestBlock(unittest.TestCase):

    def setUp(self):

        thisd = os.path.abspath(os.path.dirname(__file__))
        paths = thisd.split(os.sep)
        self.ambuildDir = os.sep.join(paths[ :-1 ])

        self.ch4Xyz = os.path.join(self.ambuildDir, "blocks", "ch4.xyz")
        self.ch4Car = os.path.join(self.ambuildDir, "blocks", "ch4.car")
        self.cx4Car = os.path.join(self.ambuildDir, "blocks", "cx4.car")
        self.ch4_1Car = os.path.join(self.ambuildDir, "blocks", "ch4_1.car")
        # self.pafCar = os.path.join( self.ambuildDir, "blocks", "PAF_bb_typed.car" )
        self.benzeneCar = os.path.join(self.ambuildDir, "blocks", "benzene.car")
        self.benzene2Car = os.path.join(self.ambuildDir, "blocks", "benzene2.car")
        self.ch4Ca2Car = os.path.join(self.ambuildDir, "blocks", "ch4Ca2.car")

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

        print "Wrote file: {0}".format(fpath)

        return

    def testBodies(self):

        b1 = Block(filePath=self.ch4Ca2Car, fragmentType='A')
        b2 = b1.copy()

        eg1 = b1.freeEndGroups()[0]
        eg2 = b2.freeEndGroups()[0]
        b1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        b1.bondBlock(bond)

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
        cx4_1.bondBlock(bond)
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
        ch4_1.bondBlock(bond)

        self.assertEqual([0, 0, 0, 4, 4, 4], [ eg.blockEndGroupIdx for eg in ch4_1.freeEndGroups() ])

        self.assertEqual(len(ch4_1._blockBonds), 1)
        self.assertEqual([ (0, 4) ], [ (b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx) for b in ch4_1._blockBonds ])

        # Check block Bonds
        self.assertEqual([ (0, 4) ], ch4_1.blockBonds())

        # Check all bonds
        self.assertEqual([(0, 1), (0, 3), (0, 2), (4, 5), (4, 7), (4, 6), (0, 4)], ch4_1.bonds())

#         # print
#         for fragment in ch4_1.fragments:
#             for i, eg in enumerate( fragment.endGroups() ):
#                 print eg


        return

    def testBondSelf(self):
        """Silly test as bonds aren't feasible"""

        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')

        eg1 = ch4_1.freeEndGroups()[1]
        eg2 = ch4_2.freeEndGroups()[2]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_1.freeEndGroups()[-1]
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        return

    def testDeleteBondSimple(self):
        """Bfoo"""
        ch4_1 = Block(filePath=self.ch4Car, fragmentType='A')
        ch4_2 = Block(filePath=self.ch4Car, fragmentType='A')

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)
        
        #ch4_1.writeCml("foo1.cml")
        x = ch4_1.deleteBond(bond)
        ch4_1.writeCml("foo2.cml")
        x.writeCml("foo3.cml")
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
        ch4_1.bondBlock(bond)

        eg1 = egFromF(ch4_1, f2)
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bondM = Bond(eg1, eg2)
        ch4_1.bondBlock(bondM)

        eg1 = egFromF(ch4_1, f3)
        eg2 = ch4_4.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        eg1 = egFromF(ch4_1, f4)
        eg2 = ch4_5.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        #ch4_1.writeCml("foo1.cml")

        # Now just bond into a loop
        eg1 = egFromF(ch4_1, f1)
        eg2 = egFromF(ch4_1, f5)
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        #ch4_1.writeCml("foo2.cml")
        self.assertFalse(ch4_1.deleteBond(bondM))
        #ch4_1.writeCml("foo3.cml")
        return

    def testDeleteBondSplit(self):
        """Create a central block with 4 attached blocks and the split one of the bonds so we get 2 blocks"""

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
        ch4_1.bondBlock(bond)
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_4.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_5.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        #ch4_1.writeCml("foo1.cml")
        self.assertTrue(bool(ch4_1.deleteBond(bond)))
        #ch4_1.writeCml("foo2.cml")
        return
    
    def testDeleteFragment(self):
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
        ch4_1.bondBlock(bond)
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_4.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_5.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)
        
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
        bcopy = block.copy()

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

        correct = numpy.array([  0.000000, 0.000000, 0.000000 ])
        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        cog = ch4.centroid()
        self.assertTrue(numpy.allclose(correct, cog, rtol=1e-9, atol=1e-6),
                         msg="testCentroid incorrect COM.".format(cog))

        return

    def testCenterOfMass(self):
        """
        Test calculation of Center of Mass
        """

        correct = numpy.array([  0.000000, 0.000000, 0.000000 ])
        ch4 = Block(filePath=self.ch4Car, fragmentType='A')
        com = ch4.centerOfMass()
        self.assertTrue(numpy.allclose(correct, com, rtol=1e-6, atol=1e-6),
                         msg="testCenterOfMass incorrect COM: {0}".format(com))
        return

    def testDihedrals(self):
        """foo"""

        ch4_1 = Block(filePath=self.benzeneCar, fragmentType='A')
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2)
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

#         # Check just across bonds
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
        m1.bondBlock(bond)

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
        m1.bondBlock(bond)

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
        newPos = staticBlock.newBondPosition(endGroup1,
                                              growBlock.symbol(endGroup2.blockEndGroupIdx),
                                               )

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
        # jmht - check...- old was: 1.78900031214
        # or maybe: 1.45942438719
        self.assertAlmostEqual(r, 1.79280605406, 7, "Incorrect radius: {}".format(str(r)))
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
        b1.bondBlock(bond)

        eg1 = b1.freeEndGroups()[-1]
        eg2 = ch4_2.freeEndGroups()[0]
        b1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        b1.bondBlock(bond)

        # b1.writeCml("foo1.cml")
        # b1.writeXyz("foo.xyz")

        coords, symbols, bonds = b1.dataByFragment('A')
        cmlFilename = "foo.cml"
        util.writeCml(cmlFilename,
                      coords,
                      symbols,
                      bonds=bonds)

        with open(cmlFilename) as f:
            test = f.readlines()

        with open(os.path.join(self.ambuildDir, "tests", "testSplitFragment.cml")) as f:
            ref = f.readlines()

        self.assertEqual(test, ref, "cml compare")
        os.unlink(cmlFilename)

        return

    def testWriteCml(self):
        """foo"""
        ch4_1 = Block(filePath=self.benzeneCar, fragmentType='A')
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()

        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock(eg1, eg2, dihedral=math.radians(180))
        bond = Bond(eg1, eg2)
        ch4_1.bondBlock(bond)

        # Write out the cml and see if it matches what we've saved
        fname = "test.cml"
        ch4_1.writeCml(fname)
        with open(fname) as f:
            test = f.readlines()

        with open(os.path.join(self.ambuildDir, "tests", "benzeneBond.cml")) as f:
            ref = f.readlines()

        self.assertEqual(test, ref, "cml compare")
        os.unlink(fname)

        return


if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()

