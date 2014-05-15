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
import collections
import copy
import os
import math
import random
import sys
import types
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
        s = "Bond {0}: {1}:{2}:{3}-{4} -> {5}:{6}:{7}-{8}".format( id(self),
                                                   self.endGroup1.block().id,id(self.endGroup1.fragment),
                                                   id(self.endGroup1),self.endGroup1.blockEndGroupIdx,
                                                   self.endGroup2.block().id,id(self.endGroup2.fragment),
                                                   id(self.endGroup2),self.endGroup2.blockEndGroupIdx )
        return s

class Block(object):
    '''
    foo
    '''

    def __init__( self, filePath=None, fragmentType=None, initFragment=None, ):
        '''
        Constructor
        '''
        
        # Need to change so cannot create block withough fragmentType
        if filePath:
            assert os.path.isfile( filePath ) and fragmentType
            initFragment = fragment.Fragment( filePath, fragmentType )
        
        # List of the fragments contained in this one
        self._fragments = [ initFragment ]
        
        # List of bond objects between blocks
        self._blockBonds = []
        
        # List of tuples of atoms that are bonded
        self._bonds = []
        
        # List of which atom is bonded to which
        self._bondedToAtom = []
        
        # list of tuples of ( idFrag, idxData )
        self._dataMap = []
        self._int2ext = {} # Maps internal block index to external with masked atoms remvoed
        self._ext2int = {} # reverse lookup
        
        # The list of atoms that are endGroups and their corresponding angleAtoms
        self._freeEndGroups          = {}
        self._numFreeEndGroups       = 0
        self._endGroupType2EndGroups = {}
        
        # Flag to indicate if block has changed and parameters (such as centerOfMass)
        # need to be changed
        self._changed = True
        
        # Below need to be updated when we move etc
        self._centroid = numpy.zeros( 3 )
        self._centerOfMass = numpy.zeros( 3 )
        self._maxAtomRadius = -1
        self._radius = None
        self._blockMass = 0
        self.id = id(self) 
        
        return self._update()
    
    def alignAtoms(self, atom1Idx, atom2Idx, refVector ):
        """Move molecule so two atoms are aligned along refVector"""
        
        atom1 = self._coord( atom1Idx )
        atom2 = self._coord( atom2Idx )
        if isinstance( refVector, list ): # Should check if numpy array
            refVector = numpy.array( refVector, dtype=numpy.float64 )
            
        return self.alignVector( atom1, atom2, refVector )
    
    def alignVector(self, pos1, pos2, refVector ):
        """
        Align this block, so that the axis defined by the two atoms is aligned with
        the refVector
        
        pos1 is the coordinate of atom1 and pos2 the coordinate of atom2
        """
        
        # Move so that pos1 is at origin so the vector of pos2 can be
        # aligned with the refVector
        
        # Shift block so angleAtom at center, 
        self.translate( -pos1 )
        
        # Check neither is zero
        if numpy.array_equal( refVector, [0,0,0] ) or numpy.array_equal( pos2, [0,0,0] ):
            raise RuntimeError, "alignBlock - one of the vectors is zero!\nrefVector: {0} endGroup: {1}".format( refVector, pos2 )
        
        # Makes no sense if they are already aligned
        if numpy.array_equal( pos2/numpy.linalg.norm(pos2), refVector/numpy.linalg.norm(refVector) ):
            print "alignBlock - block already aligned along vector. May not be a problem, but you should know..."
            self.translate( pos1 ) # NEW - put it back
            return

        #print "alignBlock BEFORE: {0} | {1}".format( endGroup, refVector )

        # Calculate normalised cross product to find an axis orthogonal 
        # to both that we can rotate about
        cross = numpy.cross( refVector, pos2 )
        
        if numpy.array_equal( cross, [0,0,0] ):
            # They must be already aligned but anti-parallel, so we flip
            print "alignBlock - vectors are anti-parallel so flipping"
            self.flip( refVector )
            self.translate( pos1 ) #NEW - put it back
        else:
            
            # Find angle
            angle = util.vectorAngle( refVector, pos2 )
            
            # Normalised cross to rotate about
            ncross = cross / numpy.linalg.norm( cross )
            
            # Rotate
            self.rotate( ncross, angle )
        
        # Now shift back 
        self.translate( pos1 )
        
        #print "alignBlock AFTER: {0} | {1}".format( self._coord( idxAtom ), refVector )
        return

    def atomBonded1(self, idxAtom):
        return [ self._int2ext[a] for a in self._atomBonded1(self._ext2int[idxAtom]) ]
        
    def _atomBonded1(self, idxAtom):
        """Return the indices of all atoms directly bonded to idxAtom"""
        return self._bondedToAtom[ idxAtom ]
#         # Get the atoms bonded to the atom in fragment space
#         fragment, idxData = self._dataMap[ idxAtom ]
#         fbonded = fragment.bonded( idxData )
#         
#         # Updated to block indices and fix and bonded cap atoms
#         bonded = []
#         for fatom in fbonded:
#             atom = fatom + fragment._blockIdx # Convert index from fragment to block
#             
#             # See if the atom is a bonded Cap - if so we return the index of the corresponding
#             # endGroup in the bond
#             try:
#                 atom = self.cap2EndGroup[ atom ]
#             except KeyError:
#                 pass
#             
#             bonded.append( atom )
#         
#         return bonded

    def atomBonded2(self, idxAtom):
        return [ self._int2ext[a] for a in self._atomBonded2(self._ext2int[idxAtom]) ]

    def _atomBonded2(self, idxAtom ):
        """Return the indices of all atoms bonded by <= 2 bonds to idxAtom"""
        bonded = copy.copy( self.atomBonded1( idxAtom ) )
        for a1 in list(bonded): # Convert to list so we're not changing while looping thru it
            bonded.update( self.atomBonded1( a1 ) )
        return bonded

    def atomBonded3(self, idxAtom):
        return [ self._int2ext[a] for a in self._atomBonded3(self._ext2int[idxAtom]) ]

    def _atomBonded3(self, idxAtom ):
        """Return the indices of all atoms bonded by <= 3 bonds to idxAtom"""
        bonded = copy.copy( self.atomBonded1( idxAtom ) )
        for a1 in list(bonded): # Convert to list so we're not changing while looping thru it
            for a2 in self.atomBonded1( a1 ):
                bonded.update( [ a1, a2 ] + self.atomBonded1( a2 ) )
        return bonded

    def body(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.bodies[idxData]
 
    def _body(self, idxAtom):
        frag, idxData = self._dataMap[ idxAtom ]
        return frag.bodies[idxData]
    
    def charge(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.charges[idxData]

    def _charge(self, idxAtom):
        frag, idxData = self._dataMap[ idxAtom ]
        return frag.charges[idxData]
    
    def coord(self, idxAtom, coord=None):
        """Get and set coordinate in external indices"""
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        if coord is not None:
            if isinstance(coord,list):
                coord=numpy.array(coord)
            frag.coords[idxData] = coord
        else:
            return frag.coords[idxData]

    def _coord(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return frag.coords[idxData]
    
    def label(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.labels[idxData]

    def _label(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag.labels[idxData]

    def mass(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.masses[idxData]

    def _mass(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag.masses[idxData]
    
    def radius(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.radii[idxData]

    def _radius(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag.radii[idxData]
    
    def symbol(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.symbols[idxData]

    def _symbol(self, idxAtom):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag.symbols[idxData]
    
    def type(self, idxAtom):
        frag, idxData = self._dataMap[ self._ext2int[ idxAtom ] ]
        return frag.types[idxData]

    def _type(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag.types[idxData]

    def atomEndGroups(self,idxAtom):
        """Return a list of the endGroup objects for this atom
        Index in external coordinates
        """
        try:
            return self._freeEndGroups[ self._ext2int[idxAtom] ]
        except KeyError,e:
            return False
    
    def atomFragType(self, idxAtom ):
        """The type of the fragment that this atom belongs to."""
        frag, idxData = self._dataMap[ idxAtom ]
        return frag.fragmentType

    def bonds(self):
        """All bonds in external indices"""
        return [ (self._int2ext[b1], self._int2ext[b2]) for b1, b2 in self._bonds ]

    def blockBonds(self):
        """External indices"""
        return [ (self._int2ext[b.endGroup1.blockEndGroupIdx], self._int2ext[b.endGroup2.blockEndGroupIdx]) for b in self._blockBonds ]
    
    def bondBlock( self, bond ):
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
            self._fragments += bond.endGroup2.block()._fragments
            self._blockBonds += bond.endGroup2.block()._blockBonds
        
        # add the new bond
        self._blockBonds.append( bond )
        
        return self._update()
    
    def _calcCenters(self):

        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        totalMass = 0.0
        for f in self._fragments:
            
            mass = f.totalMass()
            totalMass += mass
            sumG += f.centroid()
            sumM += mass * f.centroid()
        
        self._centroid = sumG / len( self._fragments )
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
        for f in self._fragments:
            for c in f.coords:
                distances.append( util.distance( self._centroid, c ) )
                self._maxAtomRadius = max( f.maxAtomRadius(), self._maxAtomRadius )
            
        imax = numpy.argmax( distances )
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
        
        print "centerOfMass changed ",self._changed
        if self._changed:
            self._calcProperties()
        return self._centerOfMass

    def centroid(self):
        """
        Return or calculate the center of geometry for the building block.
        """
        
        if self._changed:
            self._calcProperties()
        return self._centroid

    def copy( self ):
        """Return a copy of ourselves."""
        new = copy.deepcopy(self)
        new.id=id( new )
        return new

    def dihedrals(self, idxAtom1,idxAtom2, bondOnly=False):
        """Return a list of all the dihedrals around these two bonded atoms
        input and output is in external coordinates"""
        # Not sure if the lambda thing entirely kosha...
        return [ tuple( map( lambda x: self._int2ext[x], di ) ) for di in self._dihedrals( self._ext2int[idxAtom1],
                                  self._ext2int[idxAtom2],
                                  bondOnly=bondOnly ) ]
    
    def _dihedrals(self, atom1Idx, atom2Idx, bondOnly=False):
        """Return a list of all the dihedrals around these two bonded atoms
        input & output in internacl coodrdinates
        
        This needs more work to check for when things are looped back and bonded to each other
        """
        
        # Create a list of lists of all the atoms that each endGroup is bonded to - excluding the other endGroup
        # Add all dihedrals on each side of the bond - both bond atoms plus 1, 2 connected to the endGrup
        # Add all dihedrals across the bond - both endGroups plus 1 either side
        
        # Create list of what's bonded to atom1 - we exclude anything that loops back on itself
        atom1Bonded = {}
        for a1 in self._atomBonded1( atom1Idx ):
            if a1 == atom2Idx:
                continue
            atom1Bonded[ a1 ] = []
            for a2 in self._atomBonded1( a1 ):
                if a2 == atom1Idx:
                    continue
                assert not a2 == atom2Idx,"Dihedral atom loops back onto bond!"
                atom1Bonded[ a1 ].append( a2 )
        
        # Create list of what's bonded to atom2 - we exclude anything that loops back on itself
        atom2Bonded = {}
        for a1 in self._atomBonded1( atom2Idx ):
            if a1 == atom1Idx:
                continue
            atom2Bonded[ a1 ] = []
            for a2 in self._atomBonded1( a1 ):
                if a2 == atom2Idx:
                    continue
                assert not a2 == atom1Idx,"Dihedral atom loops back onto bond!"
                atom2Bonded[ a1 ].append( a2 )
        
        dindices = []
        if not bondOnly:
            # Add all dihedrals on endGroup1's side of the bond
            for a3 in atom1Bonded:
                for a4 in atom1Bonded[ a3 ]:
                    dindices.append( ( atom2Idx, atom1Idx, a3, a4 ) )
                        
            # Add all dihedrals on endGroup2's side of the bond
            for a3 in atom2Bonded:
                for a4 in atom2Bonded[ a3 ]:
                    dindices.append( ( atom1Idx, atom2Idx, a3, a4 ) )
        
        # Now add dihedrals across the bond
        for a1 in atom1Bonded:
            for a2 in atom2Bonded:
                dindices.append( ( a1, atom1Idx, atom2Idx, a2 ) )

        return dindices

    def flip( self, fvector ):
        """Rotate perpendicular to fvector so we  facing the opposite way along the fvector
        """
        
        # Find vector perpendicular to the bond axis
        # Dot product needs to be 0
        # xu + yv + zw = 0 - set u and v to 1, so w = (x + y)/z
        # vector is 1, 1, w
        w =  -1.0 * ( fvector[0] + fvector[1] ) / fvector[2]
        orth = numpy.array( [1.0, 1.0, w] )
        
        # Find axis that we can rotate about
        rotAxis = numpy.cross( fvector, orth )
        
        # Rotate by 180
        self.rotate( rotAxis, numpy.pi )
        
        return

    def fragmentBonds(self):
        """Return a list of all the internal fragment bonds for this block
        This excludes the bonds between fragments
        """
        # loop through all fragments and get a list of which atoms are internally bonded
        # Then go through the data map and map these to the 'external' atom indices
        fbonds = []
        for fragment in self._fragments:
            for b in fragment.bonds():
                fbonds.append( b )
        return fbonds
    
    def fragmentTypeDict(self):
        """A dictionary with the number of the different types of fragment we contain"""
        return self._fragmentTypeDict

    def freeEndGroups(self):
        # http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
        return  [item for sublist in self._freeEndGroups.values() for item in sublist]

    def freeEndGroupsFromTypes(self,endGroupTypes):
        # Make sure we have a list to check against
        if isinstance( endGroupTypes, str ):
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
        for i in range( len(self._ext2int) ):
            yield self.coord( i )
    
    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius

    def newBondPosition(self, endGroup, symbol ):
        """Return the position where a bond to an atom of type 'symbol'
        would be placed if bonding to the target endgroup
         I'm sure this algorithm is clunky in the extreme...
        """
        
        targetEndGroup  = self._coord( endGroup.blockEndGroupIdx )
        targetSymbol   = self._symbol( endGroup.blockEndGroupIdx )
        targetCapAtom  = self._coord( endGroup.blockCapIdx )
        
        # Get the bond length between these two atoms
        bondLength = util.bondLength( targetSymbol, symbol )
        
        # Find unit vector pointing from targetAngleAtom to targetEndGroup
        
        # vector from targetEndgroup to targetCapAtom
        #v1 = targetEndGroup - targetCapAtom
        v1 =  targetCapAtom - targetEndGroup
        
        # Now get unit vector
        uv = v1 / numpy.linalg.norm( v1 )
        
        # Multiply unit vector by bond length to get the component to add on
        newPosition = targetEndGroup + ( uv * bondLength )
        
        return newPosition

    def numFreeEndGroups(self):
        return self._numFreeEndGroups
    
    def numAtoms(self):
        """Number of atoms visible externally"""
        return len(self._ext2int)
    
    def blockMass(self):
        return self._blockMass

    def positionGrowBlock( self, endGroup, growEndGroup, dihedral=None ):
        """
        Position growBlock so it can bond to us
        """

        # The vector we want to align along is the vector from the endGroup
        # to the capAtom
        endGroupAtom = self._coord( endGroup.blockEndGroupIdx )
        capAtom      = self._coord( endGroup.blockCapIdx )
        refVector    = endGroupAtom - capAtom
        
        growBlock = growEndGroup.block()
        
        # get the coord where the next block should bond
        # symbol of endGroup tells us the sort of bond we are making which determines
        # the bond length
        symbol = growBlock._symbol( growEndGroup.blockEndGroupIdx )
        bondPos = self.newBondPosition( endGroup, symbol )
        #print "got bondPos for {0}: {1}".format( symbol, bondPos )
        
        # Align along the staticBlock bond
        growBlock.alignAtoms( growEndGroup.blockEndGroupIdx,
                              growEndGroup.blockCapIdx,
                              refVector )
        
        # Now need to place the endGroup at the bond coord to do this now we just add the bondPos
        growBlock.translate( bondPos )
        
        # We need to rotate to adhere to the specified dihedral angle
        if dihedral is not None:
            # Get current angle
            current = util.dihedral( self._coord( endGroup.blockDihedralIdx ),
                                     self._coord( endGroup.blockEndGroupIdx ),
                                     growBlock._coord( growEndGroup.blockEndGroupIdx ),
                                     growBlock._coord( growEndGroup.blockDihedralIdx ) )
            
            assert endGroup.blockDihedralIdx != -1 and growEndGroup.blockDihedralIdx != -1, \
            "Need to have specified dihedrals as 3rd column in ambi file first!"
            
            # Find how much we need to rotate by
            angle = dihedral - current 
            if angle != 0:
                # Need to rotate so get the axis to rotate about
                axis = self._coord( endGroup.blockEndGroupIdx ) - growBlock._coord( growEndGroup.blockEndGroupIdx )
                growBlock.rotate(axis, angle, center=self._coord( endGroup.blockEndGroupIdx ) )

        return

    def randomEndGroup( self, endGroupTypes=None ):
        """Return a random free endGroup in the block"""
        
        if endGroupTypes == None:
            # We pick a random endGroup
            endGroup = random.choice( self.freeEndGroups() )
        else:
            # Make sure we have a list to check against
            if isinstance( endGroupTypes, str ):
                endGroupTypes = [ endGroupTypes ]
            
            # See if any in common
            common = frozenset( self.freeEndGroupTypes() ).intersection( frozenset( endGroupTypes ) )
            if not bool( common ):
                raise RuntimeError,"Cannot find {0} in available types {1}".format( endGroupTypes, self.freeEndGroupTypes() )
            
            # We can definitely return something so pick a random fragment type and get a random endGroup
            ftype = random.choice( list( common ) )
            endGroup = random.choice( self.freeEndGroupsFromTypes( endGroupTypes=[ ftype ] ) )
            
        return endGroup
    
    def blockRadius(self):
        if self._changed:
            self._calcProperties()
        return self._radius

    def randomRotate( self, origin=[ 0,0,0 ], atOrigin=False ):
        """Randomly rotate a block.
        
         Args:
         atOrigin -- flag to indicate if the block is already positioned at the origin
        """
        
        if not atOrigin:
            position = self.centroid()
            self.translateCentroid( origin )
            
        
        xAxis = [ 1, 0, 0 ]
        yAxis = [ 0, 1, 0 ]
        zAxis = [ 0, 0, 1 ]
        
        angle = random.uniform( 0, 2*numpy.pi)
        self.rotate( xAxis, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        self.rotate( yAxis, angle )
        
        angle = random.uniform( 0, 2*numpy.pi)
        self.rotate( zAxis, angle )
        
        if not atOrigin:
            self.translateCentroid( position )

    def randomBlock(self):
        """Return a random block"""
        return random.choice( self.allBlocks )
    
    def rotate( self, axis, angle, center=None ):
        
        if center is None:
            center = numpy.array( [ 0, 0, 0 ] )
        
        rotationMatrix = util.rotation_matrix( axis, angle )

        for f in self._fragments:
            f.rotate( rotationMatrix, center )
            
        return
    
    def rotateT( self, axis, angle, center=None ):
        """Rotation with translation to center"""

        position = self.centroid()
        origin = numpy.array( [ 0, 0, 0 ] )
        self.translateCentroid( origin )
        
        rotationMatrix = util.rotation_matrix( axis, angle )
        for f in self._fragments:
            f.rotate( rotationMatrix, origin )
        
        self.translateCentroid( position )
        return
    
    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # CHANGE SO WE CHECK IF IS A NUMPY ARRAY
        if isinstance( tvector, list ):
            tvector = numpy.array( tvector )
        
        # Loop through each fragment and translate each in turn
        for f in self._fragments:
            f.translate( tvector )
        
        self._changed = True
        
        return

    def translateCentroid( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given coord"""
        self.translate( position - self.centroid() )
        return

    def _update(self):
        """Set the list of _endGroups & update data for new block
        """
        #
        # Now build up the dataMap listing where each fragment starts in the block and linking the
        # overall block atom index to the fragment and fragment atom index
        #
        self._dataMap = []
        self._int2ext = collections.OrderedDict()
        self._ext2int = collections.OrderedDict()
        self._blockMass = 0
        self._fragmentTypeDict = {}
        icount=0
        ecount=0
        for fragment in self._fragments:
            
            # Set the block
            fragment.block = self
            
            # Count the number of each type of fragment in the block (see Analyse)
            t = fragment.fragmentType
            if t not in self._fragmentTypeDict:
                self._fragmentTypeDict[ t ] = 1
            else:
                self._fragmentTypeDict[ t ] += 1
                
            fragment._blockIdx = len(self._dataMap) # Mark where the data starts in the block
            for i in xrange( fragment.numAtoms() ):
                self._dataMap.append( ( fragment, i ) )
                # Sort out indexing for masked/available atoms
                if not fragment.isMasked(i):
                    self._blockMass += fragment.masses[i]
                    self._ext2int[ecount] = icount
                    self._int2ext[icount] = ecount
                    ecount += 1
                icount += 1
        #
        # Have dataMap so now update the endGroup information
        self._numFreeEndGroups       = 0 
        self._freeEndGroups          = {}
        self._endGroupType2EndGroups = {}
        for fragment in self._fragments:
            for i, endGroup in enumerate( fragment.endGroups() ):
                assert endGroup.fragment == fragment
                #assert id(endGroup) == id( fragment._endGroups[ i ] ) # no longer valid as we only return free ones
                # Set the block index - we sort out the others after we've done bonding
                endGroup.updateEndGroupIndex()
                if endGroup.free():
                    # Add to the list of all free endGroups
                    try:
                        self._freeEndGroups[ endGroup.blockEndGroupIdx ].append( endGroup )
                    except KeyError:
                        self._freeEndGroups[ endGroup.blockEndGroupIdx ] = [ endGroup ]
                    self._numFreeEndGroups += 1
                    # Now add to the type list
                    if endGroup.type() not in self._endGroupType2EndGroups:
                        self._endGroupType2EndGroups[ endGroup.type() ] = []
                    self._endGroupType2EndGroups[ endGroup.type() ].append( endGroup )
        
        # Now need to create the list of all bonds throughout the block
        self._bonds = []
        # First all bonds within the fragments
        for fragment in self._fragments:
            for b in fragment.bonds():
                self._bonds.append( b )
        
        # Then all bonds between fragments
        cap2EndGroup = {}
        for b in self._blockBonds:
            self._bonds.append( (b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx) )
            # Need to map cap atoms to their bonded counterparts so we can look these up when we 
            # fix the endGroup indices - we map the fragment, fragmentIndex to the corresponding block index
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
        for fragment in self._fragments:
            for endGroup in fragment.endGroups():
                endGroup.updateAncillaryIndices(cap2EndGroup)
        
        # Recalculate the data for this new block
        self._calcProperties()
        
        return

    def writeCml(self,cmlFilename):

        atomTypes = []
        coords    = []
        symbols   = []
        count     = 0
        
        for i, coord in enumerate(self.iterCoord()):
            coords.append(coord)
            symbols.append(self.symbol(i))
            atomTypes.append(self.type(i))
        
        cell = None
        cmlFilename = util.writeCml(cmlFilename,
                                    coords,
                                    symbols,
                                    bonds=self.bonds(),
                                    atomTypes=atomTypes,
                                    cell=cell,
                                    pruneBonds=False)
        
        print "wrote CML file ",cmlFilename

        return
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
        
        with open(name,"w") as f:
            
            coords = []
            count=0
            for frag, i in self._dataMap:
                c = frag.coords[i]
                s = frag.symbols[i]
                
                #print count, id(frag), s, frag._masked[i], c[0], c[1], c[2]
                #if not frag._masked[i]:
                coords.append("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( s, c[0], c[1], c[2]) )
                count += 1
            
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( len(coords) ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
            f.writelines( coords )
            
        
        print "Wrote file: {0}".format(fpath)
        
        return
    def XwriteXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
        
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( self.numAtoms() ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
            
            for i, c in enumerate( self.iterCoord() ):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( self._symbol( i ), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
        return
        
    def __str__(self):
        """
        Return a string representation of the block
        """
        
        mystr = "Block: {0}\n".format(self.__repr__())
        
        mystr += "Num fragments: {0}\n".format( len( self._fragments) )
        mystr += "endGroups idxs: {0}\n".format( [ e.blockEndGroupIdx for e in self.freeEndGroups() ]  )
        mystr += "capAtoms: {0}\n".format( [ e.blockCapIdx for e in self.freeEndGroups() ] )
        #mystr += "_bondedCapAtoms: {0}\n".format( self._bondedCapIdxs )
        mystr += "bonds: {0}\n".format( self._blockBonds )
        
        for i,c in enumerate( self.iterCoord() ):
            #mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self._labels[i], c[0], c[1], c[2])
            mystr += "{0}  {1:5} {2:0< 15} {3:0< 15} {4:0< 15} \n".format( i, self._label( i ), c[0], c[1], c[2])

        return mystr
    
class TestBlock(unittest.TestCase):
    
    def setUp(self):
        
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ambuildDir = os.sep.join( paths[ : -1 ] )
        
        self.ch4Xyz = os.path.join( self.ambuildDir, "blocks", "ch4.xyz" )
        self.ch4Car = os.path.join( self.ambuildDir, "blocks", "ch4.car" )
        self.cx4Car = os.path.join( self.ambuildDir, "blocks", "cx4.car" )
        #self.pafCar = os.path.join( self.ambuildDir, "blocks", "PAF_bb_typed.car" )
        self.benzeneCar = os.path.join( self.ambuildDir, "blocks", "benzene.car" )
        self.benzene2Car = os.path.join( self.ambuildDir, "blocks", "benzene2.car" )
        
        return
    
    def catBlocks(self, blocks, filename ):
        
        symbols = []
        coords = []
        for b in blocks:
            for i, c in enumerate( b.iterCoord() ):
                symbols.append( b._symbol( i ) )
                coords.append( c )
                
        with open(filename,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( len(coords) ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
            for i, c in enumerate( coords ):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( symbols[ i ], c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
        return
    
    def testCH4(self):
        """Test the creation of a CH4 molecule"""
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        endGroups = [ 0, 0, 0, 0 ]
        self.assertEqual( endGroups, [ e.blockEndGroupIdx for e in ch4.freeEndGroups() ] )
        
        angleAtoms = [ 1, 2, 3, 4 ]
        self.assertEqual( angleAtoms, [ e.blockCapIdx for e in ch4.freeEndGroups() ] )
        
        return
    
    def testCX4(self):
        """Test the creation of a CX4 molecule"""
        
        cx4_1 = Block( filePath=self.cx4Car, fragmentType='A' )
        
        self.assertEqual( [ 0, 0, 0, 0 ], [ e.blockEndGroupIdx for e in cx4_1.freeEndGroups() ] )
        self.assertEqual( [ 1, 2, 3, 4, ], [ e.blockCapIdx for e in cx4_1.freeEndGroups() ] )
        
        cx4_2 = Block( filePath=self.cx4Car, fragmentType='A' )
        
        eg1 = cx4_1.freeEndGroups()[0]
        eg2 = cx4_2.freeEndGroups()[0]
        

        cx4_1.positionGrowBlock( eg1, eg2 )
        bond = Bond(eg1,eg2)
        cx4_1.bondBlock( bond )
        return

    def testCH4_Fragmentbond(self):
        """First pass"""
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        self.assertEqual( len( ch4._fragments[0]._bonds ), 4 )
        return
    
    def testBond1(self):
        """First pass"""
        
        ch4_1 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_2 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        
        ch4_1.positionGrowBlock( eg1, eg2 )
        bond = Bond(eg1,eg2)
        ch4_1.bondBlock( bond )
        
        self.assertEqual( [0, 0, 0, 5, 5, 5], [ eg.blockEndGroupIdx for eg in ch4_1.freeEndGroups() ] )
        #self.assertEqual( [1, 6], ch4_1._bondedCapIdxs )
        
        self.assertEqual( len(ch4_1._blockBonds), 1 )
        self.assertEqual( [ (0, 5) ], [ (b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx) for b in ch4_1._blockBonds ] )
        
        # Check external indices
        self.assertEqual( [ (0, 4) ], ch4_1.blockBonds() )
        
        # Check all bonds
        self.assertEqual(  [(0, 2), (0, 4), (0, 3), (5, 7), (5, 9), (5, 8), (0, 5)], ch4_1._bonds)
        
        # Check external indices
        self.assertEqual( [(0, 1), (0, 3), (0, 2), (4, 5), (4, 7), (4, 6), (0, 4)], ch4_1.bonds() )
        
        return
    
    def testBondSelf(self):
        """Silly test as bonds aren't feasible"""
        
        ch4_1 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_2 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        eg1 = ch4_1.freeEndGroups()[1]
        eg2 = ch4_2.freeEndGroups()[2]
        ch4_1.positionGrowBlock( eg1, eg2 )
        bond = Bond(eg1,eg2)
        ch4_1.bondBlock( bond )
    
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_1.freeEndGroups()[-1]
        bond = Bond(eg1,eg2)
        ch4_1.bondBlock( bond )
        
        return
    
    def testBondSelf2(self):
        """Bfoo"""
        
        ch4_1 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_2 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_3 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_4 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock( eg1, eg2, dihedral=180 )
        bond = Bond(eg1,eg2)
        ch4_1.bondBlock( bond )
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_3.freeEndGroups()[0]
        ch4_1.positionGrowBlock( eg1, eg2 )
        bond = Bond(eg1,eg2)
        
        ch4_1.bondBlock( bond )
        
        #print ch4_1.bonds()
        #print ch4_1.blockBonds()
        #print ch4_1._bondedToAtom
        return

    def XtestAlignBlocks(self):
        """Test we can align two _blocks correctly"""
    
        blockS = Block( filePath=self.benzeneCar, fragmentType='A'  )
        block = blockS.copy()
        
        block.translateCentroid( [ 3, 4 ,5 ] )
        block.randomRotate()
        
        # Get the atoms that define things
        eg1 = blockS.freeEndGroups()[0]
        idxSatom = eg1.blockEndGroupIdx
        idxAatom = eg1.blockCapAtomIdx
        blockSEndGroup = blockS._coord( idxSatom )
        blockSangleAtom = blockS._coord( idxAatom )
        
        
        #idxAtom = 7
        #idxAatom2 = 1
        #blockAngleAtom = block._coord( idxAatom2 )
        eg2 = block.freeEndGroups()[0]
        blockAngleAtom = eg1.blockCapAtomIdx
        
        # we want to align along blockSangleAtom -> blockSEndGroup
        refVector = blockSEndGroup - blockSangleAtom
        
        # Position block so contact is at origin
        block.translate( -blockAngleAtom )
        block.alignAtoms( idxAatom2, idxAtom, refVector )
        
        # Check the relevant atoms are in the right place
        blockEndGroup  = block._coord( idxAtom )
        blockAngleAtom = block._coord( idxAatom2 )
        
        newVector = blockEndGroup - blockAngleAtom
        
        # Normalise two vectors so we can compare them
        newNorm = newVector / numpy.linalg.norm(newVector)
        refNorm = refVector / numpy.linalg.norm(refVector)
        
        # Slack tolerances - need to work out why...
        self.assertTrue( numpy.allclose( newNorm, refNorm),
                         msg="End Group incorrectly positioned: {0} | {1}".format(newNorm, refNorm )  )
        return

    def testCentroid(self):
        """
        Test calculation of Center of Geometry
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        ch4 = Block( filePath=self.ch4Car, fragmentType='A'  )
        cog = ch4.centroid()
        self.assertTrue( numpy.allclose( correct, cog, rtol=1e-9, atol=1e-6 ),
                         msg="testCentroid incorrect COM.".format( cog ))
        
        return

    def testCenterOfMass(self):
        """
        Test calculation of Center of Mass
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        ch4 = Block( filePath=self.ch4Car, fragmentType='A'  )
        com = ch4.centerOfMass()
        self.assertTrue( numpy.allclose( correct, com, rtol=1e-6, atol=1e-6 ),
                         msg="testCenterOfMass incorrect COM: {0}".format( com ) )
        return

    def testDihedrals(self):
        """foo"""

        ch4_1 = Block( filePath=self.benzeneCar, fragmentType='A' )
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock( eg1, eg2 )
        bond = Bond(eg1,eg2)
        ch4_1.bondBlock( bond )
        
        # Check just across bonds
        ref = [ (1,0,12,17),
                (1,0,12,13),
                (5,0,12,17),
                (5,0,12,13) ]
        dihedrals = ch4_1._dihedrals(eg1.blockEndGroupIdx, eg2.blockEndGroupIdx, bondOnly=True)
        self.assertEqual(dihedrals,ref,"internal indices")
        
        ref = [ (1,0,11,16),
                (1,0,11,12),
                (5,0,11,16),
                (5,0,11,12) ]
        dihedrals = ch4_1.dihedrals(eg1.endGroupIdx(), eg2.endGroupIdx(), bondOnly=True)
        self.assertEqual(dihedrals,ref,"external indices")
        
        # Now all dihedrals
        ref = [(12, 0, 1, 2),
               (12, 0, 1, 6),
               (12, 0, 5, 11),
               (12, 0, 5, 4),
               (0, 12, 17, 16),
               (0, 12, 17, 23),
               (0, 12, 13, 18),
               (0, 12, 13, 14),
               (1, 0, 12, 17),
               (1, 0, 12, 13),
               (5, 0, 12, 17),
               (5, 0, 12, 13)]
        dihedrals = ch4_1._dihedrals(eg1.blockEndGroupIdx, eg2.blockEndGroupIdx, bondOnly=False)
        self.assertEqual(dihedrals,ref,"internal indices")
        
        ref = [(11, 0, 1, 2),
               (11, 0, 1, 6),
               (11, 0, 5, 10),
               (11, 0, 5, 4),
               (0, 11, 16, 15),
               (0, 11, 16, 21),
               (0, 11, 12, 17),
               (0, 11, 12, 13),
               (1, 0, 11, 16),
               (1, 0, 11, 12),
               (5, 0, 11, 16),
               (5, 0, 11, 12)]
        dihedrals = ch4_1.dihedrals(eg1.endGroupIdx(), eg2.endGroupIdx(), bondOnly=False)
        self.assertEqual(dihedrals,ref,"external indices")
        
        return

    def testMove(self):
        """Test we can move correctly"""
        
        paf = Block( filePath=self.benzeneCar, fragmentType='A'  )
        m = paf.copy()
        m.translate( numpy.array( [5,5,5] ) )
        c = m.centroid()
        paf.translateCentroid( c )
        p = paf.centroid()
        
        self.assertTrue( numpy.allclose( p, c, rtol=1e-9, atol=1e-9 ), "simple move")
        return
    
    def testAlignAtoms(self):
        
        block = Block( filePath=self.benzeneCar, fragmentType='A'  )
        bcopy = block.copy()
        
        # Check atoms are not aligned along axis
        c1Idx=2
        c2Idx=5
        c1 = block._coord( c1Idx )
        c2 = block._coord( c2Idx )
        
        #self.assertTrue( numpy.allclose( c1-c2 , [ 3.0559,  -0.36295,  0.07825], atol=1E-7  ), "before" )
        self.assertTrue( numpy.allclose( c1-c2 , [ 3.0559,  -0.36295,  0.07825] ), "before" )
        
        # Align along z-axis
        block.alignAtoms( c1Idx, c2Idx, [ 0, 0, 1 ] )
        
        # check it worked
        c1 = block._coord( c1Idx )
        c2 = block._coord( c2Idx )
        z = numpy.array([  0.0,  0.0, -3.07837304 ])
        
        self.assertTrue( numpy.allclose( c1-c2 , z ), "after" )

        return
    
    def testPositionGrowBlock(self):
        
        blockS = Block(filePath= self.benzeneCar, fragmentType='A'  )
        
        growBlock = blockS.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        
        # Get the atoms that define things
        endGroup1 = blockS.freeEndGroups()[ 0 ]
        endGroup2 = growBlock.freeEndGroups()[ 1 ]
        
        # Get position to check
        newPos = blockS.newBondPosition( endGroup1, growBlock._symbol( endGroup2.blockEndGroupIdx ) )
        
        # Position the block
        blockS.positionGrowBlock( endGroup1, endGroup2 )
        
        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock._coord( endGroup2.blockEndGroupIdx ) 
        self.assertTrue( numpy.allclose( newPos, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        return
    
    def testPositionGrowBlock2(self):
        
        staticBlock = Block( filePath=self.benzeneCar, fragmentType='A'  )
        
        growBlock = staticBlock.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        # Get the atoms that define things
        endGroup1 = staticBlock.freeEndGroups()[ 0 ]
        endGroup2 = growBlock.freeEndGroups()[ 0 ]
        
        # Get position to check
        newPos = staticBlock.newBondPosition( endGroup1, growBlock._symbol( endGroup2.blockEndGroupIdx ) )
        
        #staticBlock._symbols.append( 'N' )
        #staticBlock._coords.append( newPos )
        #staticBlock.writeXyz("FOO.xyz")
        
        
        # Position the block
        #staticBlock.XXpositionGrowBlock( endGroup1, growBlock, endGroup2 )
        staticBlock.positionGrowBlock( endGroup1, endGroup2 )
        
        #self.catBlocks( [staticBlock, growBlock ], "both.xyz")
        
        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock._coord( endGroup2.blockEndGroupIdx ) 
        self.assertTrue( numpy.allclose( newPos, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        return
    
    def testPositionDihedral(self):

        staticBlock = Block( filePath=self.benzeneCar, fragmentType='A' )
        
        growBlock = staticBlock.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        # Get the atoms that define things
        endGroup1 = staticBlock.freeEndGroups()[ 0 ]
        endGroup2 = growBlock.freeEndGroups()[ 0 ]
        
        # Get position to check
        newPos = staticBlock.newBondPosition( endGroup1,
                                              growBlock._symbol( endGroup2.blockEndGroupIdx ),
                                               )
        
        # Position the block
        staticBlock.positionGrowBlock( endGroup1, endGroup2, dihedral=math.pi/2 )
        
        # Hacky - just use one of the coords I checked manually
        hcheck = numpy.array( [11.98409351860, 8.826721156800, -1.833703434310] )
        endGroupCoord = growBlock._coord( 11 ) 
        self.assertTrue( numpy.allclose( hcheck, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        #self.catBlocks( [staticBlock, growBlock ], "both2.xyz")
        return

    def testRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A'  )
        r = ch4.blockRadius()
        #jmht - check...- old was: 1.78900031214
        # or maybe: 1.45942438719
        self.assertAlmostEqual(r, 1.79280605406, 7, "Incorrect radius: {}".format(str(r)) )
        return
        
    def XtestMaxAtomRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A' )
        r = ch4.maxAtomRadius()
        #jmht - check...- old was: 1.78900031214
        self.assertAlmostEqual(r, 0.70380574117, 7, "Incorrect radius: {}".format(str(r)) )
        
    def testRotate(self):
        """
        Test the rotation
        """
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A'  )
        
        array1 = numpy.array( [ -0.51336 ,  0.889165, -0.363 ] )
        self.assertTrue( numpy.array_equal( ch4._coord(4), array1 ),
                         msg="testRotate arrays before rotation incorrect.")
        
        axis = numpy.array([1,2,3])
        angle = 2
        ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose( ch4._coord(4), array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")
        
        # Check rotation by 360
        axis = numpy.array([1,2,3])
        angle = numpy.pi*2
        ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose( ch4._coord(4), array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")
        
        return
    
    def testWriteCml(self):
        """foo"""
        ch4_1 = Block( filePath=self.benzeneCar, fragmentType='A' )
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()
        
        eg1 = ch4_1.freeEndGroups()[0]
        eg2 = ch4_2.freeEndGroups()[0]
        ch4_1.positionGrowBlock( eg1, eg2 )
        bond = Bond(eg1,eg2)
        ch4_1.bondBlock( bond )
        
        # Write out the cml and see if it matches what we've saved
        fname="test.cml"
        ch4_1.writeCml(fname)
        with open(fname) as f:
            test = f.readlines()
        
        with open(os.path.join( self.ambuildDir, "tests", "benzeneBond.cml" )) as f:
            ref = f.readlines()
        
        self.assertEqual(test,ref,"cml compare")
        os.unlink(fname)
        
        return


if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
        
