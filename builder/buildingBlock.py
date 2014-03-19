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
    def __init__(self):
        
        self.block1         = None
        self.idxBlock1      = None
        self.endGroup1      = None
        
        self.block2         = None
        self.idxBlock2      = None
        self.endGroup2      = None
        
        return
    
    def __str__(self):
        """List the data attributes of this object"""
        
        s = "Bond {0}: {1}-{2} -> {3}-{4}".format( id(self), self.idxBlock1, self.endGroup1.blockEndGroupIdx,
                                                   self.idxBlock2, self.endGroup1.blockEndGroupIdx )
        return s

class EndGroup(object):
    
    def __init__(self):
        
        self.free = True
        
        self.fragment             = None
        
        self.fragmentEndGroupIdx  = None
        self.blockEndGroupIdx     = None
        
        self.fragmentCapIdx       = None
        self.blockCapIdx          = None
        
        self.fragmentDihedralIdx  = -1
        self.blockDihedralIdx     = -1
        
        self.fragmentUwIdx       = -1
        self.blockUwIdx          = -1
        
        return
    
    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
                
        t = []
        for k in sorted(me):
            t.append( str( ( k, me[k] ) ) )
                                   
        return "{0} : {1}".format(self.__repr__(), ",".join( t ) )
        

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
        
        # List of all the direct bonds to this atom
        self._bonds = []
        
        # Links cap atoms in bonds to the opposite endGroup
        self._cap2EndGroup = None
        
        # List of the different fragmentTypes contained in this block
        self._ftype2endGroup = []
        
        # list of tuples of ( idFrag, idxData )
        self._dataMap = []
        
        # The list of atoms that are endGroups and their corresponding angleAtoms
        self._endGroups = []
        self._freeEndGroupIdxs = []
        self._ignoreAtom = []
        
        # the number of free endGroups
        self._numFeeEndGroups = 0
        
        # Flag to indicate if block has changed and parameters (such as centerOfMass)
        # need to be changed
        self._changed = True
        
        # Below need to be updated when we move etc
        self._centroid = numpy.zeros( 3 )
        self._centerOfMass = numpy.zeros( 3 )
        self._maxAtomRadius = -1
        self._radius = None
        self._mass = 0
        self._numAtoms = 0

        self._blockId = id(self)
        
        return self._update()
    
    def alignAtoms(self, atom1Idx, atom2Idx, refVector ):
        """Move molecule so two atoms are aligned along refVector"""
        
        atom1 = self.atomCoord( atom1Idx )
        atom2 = self.atomCoord( atom2Idx )
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
        
        #print "alignBlock AFTER: {0} | {1}".format( self.atomCoord( idxAtom ), refVector )
        return

    def atomFragType(self, idxAtom ):
        """The type of the fragment that this atom belongs to."""
        frag, idxData = self._dataMap[ idxAtom ]
        return frag._fragmentType

    def atomBonded(self, idxAtom):
        
        # Make sure the atom isn't one that has been removed from the cell
        assert not self.invisibleAtom(idxAtom)
        
        # Get the atoms bonded to the atom in fragment space
        fragment, idxData = self._dataMap[ idxAtom ]
        fbonded = fragment.bonded( idxData )
        
        # Updated to block indices and fix and bonded cap atoms
        bonded = []
        for fatom in fbonded:
            
            atom = fatom + fragment._blockIdx # Convert index from fragment to block
            
            # See if the atom is a bonded Cap - if so we return the index of the corresponding
            # endGroup in the bond
            try:
                atom = self.cap2EndGroup[ atom ]
            except KeyError:
                # If it's not a bonded cap but is invisible (e.g.uwAtoms) we skip it
                if self.invisibleAtom( atom ):
                    continue
            
            bonded.append( atom )
        
        return bonded
    
    def atomCharge(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return frag._charges[ idxData ]
    
    def atomCoord(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return frag._coords[ idxData]
    
    def atomLabel(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag._labels[ idxData ]
        
    def atomType(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag._atomTypes[ idxData ]
    
    def atomMass(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag._masses[ idxData ]
    
    def atomRadius(self, idxAtom ):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag._atomRadii[ idxData ]
    
    def atomSymbol(self, idxAtom):
        frag, idxData = self._dataMap[ idxAtom ]
        return  frag._symbols[ idxData ]
    
    def bondBlock( self, bond ):
        """ Add newBlock to this one
        """
        
        assert bond.block1 == self
        
        # Mark both endGroups as used
        bond.endGroup1.free = False
        bond.endGroup2.free = False
        
        # Add all the other blocks fragments to this one
        if bond.block1 != bond.block2:
            self._fragments += bond.block2._fragments
        
            # Add all the bonds between the fragments in the other block to this one
            self._bonds += bond.block2._bonds
        
        # Add the bond we are making to list of bonds
        self._bonds.append( bond )
        
        return self._update()
    
    def bonds(self):
        return self._bonds
    
    def _calcCenters(self):

        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        totalMass = 0.0
        for f in self._fragments:
            
            mass = f.mass()
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
            for c in f._coords:
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
        """
        Create a copy of ourselves and return it.
        Hit problem when passing in the distance function from the cell as the deepcopy
        also had a copy of the cell
        """
        new = copy.deepcopy(self)
        new._blockId = id( new )
        return new

    def dihedrals(self, atom1Idx, atom2Idx):
        """Return a list of all the dihedrals around these two bonded atoms
        
        This needs more work to check for when things are looped back and bonded to each other
        """
        
        
        # Create a list of lists of all the atoms that each endGroup is bonded to - excluding the other endGroup
        # Add all dihedrals on each side of the bond - both bond atoms plus 1, 2 connected to the endGrup
        # Add all dihedrals across the bond - both endGroups plus 1 either side
        
        # Check are bonded (remove when sure of code)
        #assert atom2Idx in self.atomBonded( atom1Idx ) and atom1Idx in self.atomBonded( atom2Idx )
        
        # Create list of what's bonded to atom1 - we exclude anything that loops back on itself
        atom1Bonded = {}
        for a1 in self.atomBonded( atom1Idx ):
            if a1 == atom2Idx:
                continue
            atom1Bonded[ a1 ] = []
            for a2 in self.atomBonded( a1 ):
                if a2 == atom1Idx:
                    continue
                assert not a2 == atom2Idx,"Dihedral atom loops back onto bond!"
                atom1Bonded[ a1 ].append( a2 )
                
        # Create list of what's bonded to atom2 - we exclude anything that loops back on itself
        atom2Bonded = {}
        for a1 in self.atomBonded( atom2Idx ):
            if a1 == atom1Idx:
                continue
            atom2Bonded[ a1 ] = []
            for a2 in self.atomBonded( a1 ):
                if a2 == atom2Idx:
                    continue
                assert not a2 == atom1Idx,"Dihedral atom loops back onto bond!"
                atom2Bonded[ a1 ].append( a2 )
                
        dindices = []
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
        """Return a list of all the internal fragment bonds for this block"""
        
        # loop through all fragments and get a list of which atoms are internally bonded
        # Then go through the data map and map these to the 'external' atom indices
        fbonds = []
        for fragment in self._fragments:
            for b in fragment._bonds:
                fbonds.append( ( fragment._blockIdx + b[0], fragment._blockIdx + b[1] ) )
        return fbonds
    
    def fragmentTypeDict(self):
        """A dictionary with the number of the different types of fragment we contain"""
        return self._fragmentTypeDict


    def getEndGroups(self, idxAtom=None, fragmentTypes=None, checkFree=True ):
        """return a list of endGroups
        With no arguments, return all endGroups, given the index of an atom that is an endGroup, 
        return all the endGroup objects associated with it. If fragmentType is set only return
        endGroups of that fragment type"""
        
        if idxAtom is None and fragmentTypes is None and not checkFree:
            return self._endGroups
        
        if idxAtom is not None and not self.isEndGroup( idxAtom ):
            raise RuntimeError,"{0} is not an endGroup!".format( idxAtom )
        
        # Think about how to better optimise this so we don't have to loop each time
        # FIX TO USE ftype2endGroup!!!
        endGroups = []
        for endGroup in self._endGroups:
            
            if checkFree and not endGroup.free:
                continue
            
            if idxAtom is not None and endGroup.blockEndGroupIdx != idxAtom:
                continue
            
            if fragmentTypes is not None and endGroup.fragment.type() not in fragmentTypes:
                continue
            
            endGroups.append( endGroup )
        
        #if not len( endGroups ):
        #    raise RuntimeError,"getEndGroups found no endGroups for ftype {0} idxAtom {1}!".format( fragmentTypes, idxAtom  )
        
        return endGroups
    
    def getEndGroupTypes(self):
        """Return a list of the fragmentTypes for the available endGroups"""
        return self._ftype2endGroup.keys()
    
    def id(self):
        return self._blockId

    def ignoreAtom(self, idxAtom ):
        return self._ignoreAtom[ idxAtom ]
    
    def invisibleAtom(self, idxAtom ):
        return self.atomSymbol( idxAtom ).lower() == 'x' or self.ignoreAtom(idxAtom)

    def isEndGroup(self, idxAtom):
        """Return True if this atom is a free endGroup
        Currently this is the only place _freeEndGroupIdxs is used - 
        we could use getEndGroups with the checkFree argument to obviate
        the need for _freeEndGroupIdxs but for the time being I think the added
        complexity is probably warranted for the extra speed (not measured)
        """
        return idxAtom in self._freeEndGroupIdxs
    
    def iterCoord(self):
        """Generator to return the coordinates"""
        for idx in range( len(self._dataMap) ):
            yield self.atomCoord( idx )
    
    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius

    def newBondPosition(self, endGroup, symbol ):
        """Return the position where a bond to an atom of type 'symbol'
        would be placed if bonding to the target endgroup
         I'm sure this algorithm is clunky in the extreme...
        """
        
        targetEndGroup  = self.atomCoord( endGroup.blockEndGroupIdx )
        targetSymbol    = self.atomSymbol( endGroup.blockEndGroupIdx )
        targetCapAtom = self.atomCoord( endGroup.blockCapIdx )
        
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
        return len( self._freeEndGroupIdxs )
    
    def numAllAtoms(self):
        """Number of atoms including bonded caps"""
        return len(self._dataMap)
    
    def numAtoms(self):
        """Number of atoms excluding bonded caps"""
        return self._numAtoms
    
    def mass(self):
        return self._mass

    def positionGrowBlock( self, endGroup, growBlock, growEndGroup, dihedral=None ):
        """
        Position growBlock so it can bond to us
        """

        # The vector we want to align along is the vector from the endGroup
        # to the capAtom
        endGroupAtom = self.atomCoord( endGroup.blockEndGroupIdx )
        capAtom      = self.atomCoord( endGroup.blockCapIdx )
        refVector    =  endGroupAtom - capAtom
        
        # get the coord where the next block should bond
        # symbol of endGroup tells us the sort of bond we are making which determines
        # the bond length
        symbol = growBlock.atomSymbol( growEndGroup.blockEndGroupIdx )
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
            current = util.dihedral( self.atomCoord( endGroup.blockDihedralIdx ),
                                     self.atomCoord( endGroup.blockEndGroupIdx ),
                                     growBlock.atomCoord( growEndGroup.blockEndGroupIdx ),
                                     growBlock.atomCoord( growEndGroup.blockDihedralIdx ) )
            
            assert endGroup.blockDihedralIdx != -1 and growEndGroup.blockDihedralIdx != -1, \
            "Need to have specified dihedrals as 3rd column in ambi file first!"
            
            # Find how much we need to rotate by
            angle = dihedral - current 
            if angle != 0:
                # Need to rotate so get the axis to rotate about
                axis = self.atomCoord( endGroup.blockEndGroupIdx ) - growBlock.atomCoord( growEndGroup.blockEndGroupIdx )
                growBlock.rotate(axis, angle, center=self.atomCoord( endGroup.blockEndGroupIdx ) )

        return

    def randomEndGroup( self, fragmentTypes=None ):
        """Randomly select an atom that is an endGroup """
        
        if fragmentTypes == None:
            # We pick a random endGroup index
            #i = random.randint( 0, len( self._endGroups ) - 1 )
            endGroup = random.choice( self.getEndGroups( checkFree=True ) )
        else:
            # Make sure we have a list to check against
            if isinstance( fragmentTypes, str ):
                fragmentTypes = [ fragmentTypes ]
            
            # See if any in common
            common = frozenset( self.getEndGroupTypes() ).intersection( frozenset( fragmentTypes ) )
            if not bool( common ):
                return False
            
            # We can definitely return something so pick a random fragment type and get a random endGroup
            ftype = random.choice( list( common ) )
            endGroup = random.choice( self.getEndGroups( fragmentTypes=[ ftype ] ) )
            
        return endGroup
    
    def radius(self):
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
    
    def setCoord(self, idxAtom, coord ):
        frag, idxData = self._dataMap[ idxAtom ]
        frag._coords[ idxData] = coord
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
        to the given coord
        """
        self.translate( position - self.centroid() )
        
        return

    def _update( self ):
        """Set the list of _endGroups & update data for new block"""
        
        # Now build up the dataMap
        self._dataMap = []
        self._mass = 0
        self._fragmentTypeDict = {}
        count=0
        for fragment in self._fragments:
            t = fragment.type()
            if t not in self._fragmentTypeDict:
                self._fragmentTypeDict[ t ] = 1
            else:
                self._fragmentTypeDict[ t ] += 1
                
            fragment._blockIdx = count # Mark where the data starts in the block
            for j in range( len( fragment._coords ) ):
                self._dataMap.append( ( fragment, j ) )
                self._mass += fragment._masses[ j ]
                count += 1
        
        # Have dataMap so now update the endGroup information
        self._endGroups = []
        self._freeEndGroupIdxs = []
        self._ignoreAtom = [ False ] * len(self._dataMap) # Use bool array so we can just check an index and don't need to search
        # through the array each time. Could use indexes and then trigger an exception - not sure which is best yet
        self._ftype2endGroup = {}
        self._numAtoms = len(self._dataMap) # Get total number of atoms and subtract # endGroups
        for fragment in self._fragments:
            for i, endGroup in enumerate( fragment._endGroups ):
                
                assert endGroup.fragment == fragment
                assert id(endGroup) == id( fragment._endGroups[ i ] )
                
                # Update the block-wide indices for the endGroups
                endGroup.blockEndGroupIdx  = fragment._blockIdx + endGroup.fragmentEndGroupIdx
                endGroup.blockCapIdx = fragment._blockIdx + endGroup.fragmentCapIdx
                # -1 means no dihedral or uw atom set
                if endGroup.fragmentDihedralIdx != -1:
                    endGroup.blockDihedralIdx = fragment._blockIdx + endGroup.fragmentDihedralIdx
                if endGroup.fragmentUwIdx != -1:
                    endGroup.blockUwIdx = fragment._blockIdx + endGroup.fragmentUwIdx
                self._endGroups.append( endGroup )
                
                if endGroup.free:
                    self._freeEndGroupIdxs.append( endGroup.blockEndGroupIdx )
                    if not endGroup.fragment.type() in self._ftype2endGroup:
                        self._ftype2endGroup[ endGroup.fragment.type() ] = []
                    self._ftype2endGroup[ endGroup.fragment.type() ].append( endGroup )
                else:
                    # If the endGroup is involved in a bond the capAtom becomes invisible
                    self._ignoreAtom[ endGroup.blockCapIdx ] = True
                    
                    # Removed cap atom so remove from list of atoms and also adjust the mass
                    self._numAtoms -= 1
                    self._mass -= fragment._masses[ endGroup.fragmentCapIdx ]
                    # Uw - unwanted atoms are atoms that need to be removed when we make a bond - 
                    # is a bit of a hack at the moment
                    if endGroup.blockUwIdx != -1:
                        self._ignoreAtom[ endGroup.blockUwIdx ] = True
                        self._numAtoms -= 1
                        self._mass -= fragment._masses[ endGroup.fragmentCapIdx ]

        # Need to map bonded cap Atoms to the opposite endGroups
        self.cap2EndGroup = {}
        for bond in self._bonds:
            self.cap2EndGroup[ bond.endGroup1.blockCapIdx ] = bond.endGroup2.blockEndGroupIdx
            self.cap2EndGroup[ bond.endGroup2.blockCapIdx ] = bond.endGroup1.blockEndGroupIdx
        
        # Recalculate the data for this new block
        self._calcProperties()
        
        return
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
        
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( self.numAtoms() ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
                             
            for i, c in enumerate( self.iterCoord() ):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( self.atomSymbol( i ), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
        return
        
    def __str__(self):
        """
        Return a string representation of the block
        """
        
        mystr = "Block: {0}\n".format(self.__repr__())
        
        mystr += "Num fragments: {0}\n".format( len( self._fragments) )
        mystr += "endGroups idxs: {0}\n".format( [ e.blockEndGroupIdx for e in self._endGroups ]  )
        mystr += "_freeEndGroupIdxs: {0}\n".format( self._freeEndGroupIdxs )
        mystr += "capAtoms: {0}\n".format( [ e.blockCapIdx for e in self._endGroups ] )
        #mystr += "_bondedCapAtoms: {0}\n".format( self._bondedCapIdxs )
        mystr += "bonds: {0}\n".format( self._bonds )
        
        for i,c in enumerate( self.iterCoord() ):
            #mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self._labels[i], c[0], c[1], c[2])
            mystr += "{0}  {1:5} {2:0< 15} {3:0< 15} {4:0< 15} \n".format( i, self.atomLabel( i ), c[0], c[1], c[2])

        return mystr
    
class TestBlock(unittest.TestCase):
    
    def setUp(self):
        
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ambuildDir = os.sep.join( paths[ : -1 ] )
        
        self.ch4Xyz = os.path.join( self.ambuildDir, "blocks", "ch4.xyz" )
        self.ch4Car = os.path.join( self.ambuildDir, "blocks", "ch4.car" )
        self.cx4Car = os.path.join( self.ambuildDir, "blocks", "cx4.car" )
        self.pafCar = os.path.join( self.ambuildDir, "blocks", "PAF_bb_typed.car" )
        self.benzeneCar = os.path.join( self.ambuildDir, "blocks", "benzene.car" )
        self.benzene2Car = os.path.join( self.ambuildDir, "blocks", "benzene2.car" )
        
        return
    
    def catBlocks(self, blocks, filename ):
        
        symbols = []
        coords = []
        for b in blocks:
            for i, c in enumerate( b.iterCoord() ):
                symbols.append( b.atomSymbol( i ) )
                coords.append( c )
                
        with open(filename,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( len(coords) ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
            for i, c in enumerate( coords ):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( symbols[ i ], c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
        return
    
    def bondBlocks(self, block1, block2, endGroupIdx1, endGroupIdx2):
        
        endGroup1 = block1._endGroups[ endGroupIdx1 ]
        endGroup2 = block2._endGroups[ endGroupIdx2 ]
        
        block1.positionGrowBlock( endGroup1, block2, endGroup2 )
        
        bond = Bond()
        bond.block1          = block1
        bond.endGroup1       = endGroup1
        
        bond.block2          = block2
        bond.endGroup2       = endGroup2
        
        block1.bondBlock( bond )
        
        return
    
    def testCH4(self):
        """Test the creation of a CH4 molecule"""
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        #endGroups = [ 1, 2, 3, 4 ]
        endGroups = [ 0, 0, 0, 0 ]
        self.assertEqual( endGroups, [ e.blockEndGroupIdx for e in ch4._endGroups ] )
        
        #angleAtoms = [ 0, 0, 0, 0 ]
        angleAtoms = [ 1, 2, 3, 4 ]
        self.assertEqual( angleAtoms, [ e.blockCapIdx for e in ch4._endGroups ] )
        
        return
    
    def testCX4(self):
        """Test the creation of a CX4 molecule"""
        
        cx4_1 = Block( filePath=self.cx4Car, fragmentType='A' )
        
        self.assertEqual( [ 0, 0, 0, 0 ], [ e.blockEndGroupIdx for e in cx4_1._endGroups ] )
        self.assertEqual( [ 1, 2, 3, 4, ], [ e.blockCapIdx for e in cx4_1._endGroups ] )
        
        cx4_2 = Block( filePath=self.cx4Car, fragmentType='A'  )
        
        idxEg1=0
        idxEg2=0
        self.bondBlocks( cx4_1, cx4_2, idxEg1, idxEg2 )
        
        return

    def testCH4_Fragmentbond(self):
        """First pass"""
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        self.assertEqual( len( ch4._fragments[0]._bonds ), 4 )
        return
    
    def testCH4_bond(self):
        """First pass"""
        
        ch4_1 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_2 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        idxEg1=0 # 0 -> 1
        idxEg2=0 # 0 -> 1
        
        self.bondBlocks( ch4_1, ch4_2, idxEg1, idxEg2 )
        
        #ch4_1.writeXyz("testBond.xyz")
        
        self.assertEqual( [0, 0, 0, 5, 5, 5], ch4_1._freeEndGroupIdxs )
        #self.assertEqual( [1, 6], ch4_1._bondedCapIdxs )
        
        self.assertEqual( len(ch4_1._bonds), 1 )
        
        bonds = [ ( b.endGroup1.blockEndGroupIdx, b.endGroup2.blockEndGroupIdx ) for b in ch4_1._bonds ]
        self.assertEqual( [ (0, 5) ], bonds )
        
        return
    
    def testCH4_bond_self(self):
        """First pass"""
        
        
        ch4_1 = Block( filePath=self.ch4Car, fragmentType='A' )
        ch4_2 = Block( filePath=self.ch4Car, fragmentType='A' )
        
        idxEg1=1
        idxEg2=2
        self.bondBlocks( ch4_1, ch4_2, idxEg1, idxEg2 )
        
        idxEg1=0
        idxEg2=5
        self.bondBlocks( ch4_1, ch4_1, idxEg1, idxEg2 )
         
        return

    def testAlignBlocks(self):
        """Test we can align two _blocks correctly"""
    
        blockS = Block( filePath=self.pafCar, fragmentType='A'  )
        block = blockS.copy()
        
        block.translateCentroid( [ 3, 4 ,5 ] )
        block.randomRotate()
        
        # Get the atoms that define things
        idxSatom = 17
        idxAatom = 3
        blockSEndGroup = blockS.atomCoord( idxSatom )
        blockSangleAtom = blockS.atomCoord( idxAatom )
        
        idxAtom = 7
        idxAatom2 = 1
        blockAngleAtom = block.atomCoord( idxAatom2 )
        
        # we want to align along blockSangleAtom -> blockSEndGroup
        refVector = blockSEndGroup - blockSangleAtom
        
        # Position block so contact is at origin
        block.translate( -blockAngleAtom )
        block.alignAtoms( idxAatom2, idxAtom, refVector )
        
        # Check the relevant atoms are in the right place
        blockEndGroup  = block.atomCoord( idxAtom )
        blockAngleAtom = block.atomCoord( idxAatom2 )
        
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

    def testMove(self):
        """Test we can move correctly"""
        
        paf = Block( filePath=self.pafCar, fragmentType='A'  )
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
        c1 = block.atomCoord( c1Idx )
        c2 = block.atomCoord( c2Idx )
        
        #self.assertTrue( numpy.allclose( c1-c2 , [ 3.0559,  -0.36295,  0.07825], atol=1E-7  ), "before" )
        self.assertTrue( numpy.allclose( c1-c2 , [ 3.0559,  -0.36295,  0.07825] ), "before" )
        
        # Align along z-axis
        block.alignAtoms( c1Idx, c2Idx, [ 0, 0, 1 ] )
        
        # check it worked
        c1 = block.atomCoord( c1Idx )
        c2 = block.atomCoord( c2Idx )
        z = numpy.array([  0.0,  0.0, -3.07837304 ])
        
        self.assertTrue( numpy.allclose( c1-c2 , z ), "after" )

        return
    
    def testPositionGrowBlock(self):
        
        blockS = Block(filePath= self.pafCar, fragmentType='A'  )
        
        growBlock = blockS.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        
        # Get the atoms that define things
        endGroup1 = blockS._endGroups[ 0 ]
        endGroup2 = growBlock._endGroups[ 1 ]
        
        # Get position to check
        newPos = blockS.newBondPosition( endGroup1, growBlock.atomSymbol( endGroup2.blockEndGroupIdx ) )
        
        # Position the block
        blockS.positionGrowBlock( endGroup1, growBlock, endGroup2 )
        
        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock.atomCoord( endGroup2.blockEndGroupIdx ) 
        self.assertTrue( numpy.allclose( newPos, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        return
    
    def testPositionGrowBlock2(self):
        
        staticBlock = Block( filePath=self.benzeneCar, fragmentType='A'  )
        
        growBlock = staticBlock.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        # Get the atoms that define things
        endGroup1 = staticBlock._endGroups[ 0 ]
        endGroup2 = growBlock._endGroups[ 0 ]
        
        # Get position to check
        newPos = staticBlock.newBondPosition( endGroup1, growBlock.atomSymbol( endGroup2.blockEndGroupIdx ) )
        
        #staticBlock._symbols.append( 'N' )
        #staticBlock._coords.append( newPos )
        #staticBlock.writeXyz("FOO.xyz")
        
        
        # Position the block
        #staticBlock.XXpositionGrowBlock( endGroup1, growBlock, endGroup2 )
        staticBlock.positionGrowBlock( endGroup1, growBlock, endGroup2 )
        
        #self.catBlocks( [staticBlock, growBlock ], "both.xyz")
        
        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock.atomCoord( endGroup2.blockEndGroupIdx ) 
        self.assertTrue( numpy.allclose( newPos, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        return
    
    def testPositionDihedral(self):

        staticBlock = Block( filePath=self.benzeneCar, fragmentType='A'  )
        
        growBlock = staticBlock.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        # Get the atoms that define things
        endGroup1 = staticBlock._endGroups[ 0 ]
        endGroup2 = growBlock._endGroups[ 0 ]
        
        # Get position to check
        newPos = staticBlock.newBondPosition( endGroup1,
                                              growBlock.atomSymbol( endGroup2.blockEndGroupIdx ),
                                               )
        
        # Position the block
        staticBlock.positionGrowBlock( endGroup1, growBlock, endGroup2, dihedral=math.pi/2 )
        
        # Hacky - just use one of the coords I checked manually
        hcheck = numpy.array( [11.98409351860, 8.826721156800, -1.833703434310] )
        endGroupCoord = growBlock.atomCoord( 11 ) 
        self.assertTrue( numpy.allclose( hcheck, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        #self.catBlocks( [staticBlock, growBlock ], "both2.xyz")
        return

    def testRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = Block( filePath=self.ch4Car, fragmentType='A'  )
        r = ch4.radius()
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
        self.assertTrue( numpy.array_equal( ch4.atomCoord(4), array1 ),
                         msg="testRotate arrays before rotation incorrect.")
        
        axis = numpy.array([1,2,3])
        angle = 2
        ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose( ch4.atomCoord(4), array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")
        
        # Check rotation by 360
        axis = numpy.array([1,2,3])
        angle = numpy.pi*2
        ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose( ch4.atomCoord(4), array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")
        
        return

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
        
