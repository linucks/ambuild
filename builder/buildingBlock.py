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
    def __init__(self):
        
        self.block1 = None
        self.frag1 = None
        self.frag1atomIdx = None # fragment index 
        self.frag1angleIdx = None
        self.atom1Idx = None # index in whole block
        self.angle1Idx = None
        
        self.block2 = None
        self.frag2 = None
        self.frag2atomIdx = None
        self.frag2angleIdx = None
        self.atom2Idx = None
        self.angle2Idx = None
        
    def __str__(self):
        """
        Return a string representation of the bond
        """
        
        mystr = "Bond: {0} [\n".format(self.__repr__())
        
        mystr += "block1: {0}\n".format( id(self.block1) )
        mystr += "frag1: {0}\n".format( self.frag1 )
        mystr += "atomIdx1: {0}\n".format( self.frag1atomIdx )
        mystr += "block2: {0}\n".format( id(self.block2) )
        mystr += "frag2: {0}\n".format( self.frag2 )
        mystr += "atomIdx2: {0} ]\n".format( self.frag2atomIdx )
        
        return mystr

class Block(object):
    '''
    classdocs
    
    external API:
    atomSymbol
    atomRadius
    atomCoord
    atomLabel
    centroid
    iterCoord
    newBondPosition
    positionGrowBlock
    radius
    randomEndGroup
    randomRotate
    randomRotateBlock
    rotate
    translateCentroid
    
    '''

    def __init__( self, filename=None, fragmentType=None ):
        '''
        Constructor
        '''
        
        # The root (original) fragment
        self._rootFragment = fragment.Fragment( filename=filename, fragmentType=fragmentType )
        
        # List of all the direct bonds to this atom
        self._myBonds = []
        
        # The lists below here are constructed dynamically on bonding
        
        # List of the fragments contained in this one
        self._fragments = []
        
        # List of the different fragmentTypes contained in this block
        self._fragmentTypes = []
        
        # list of tuples of ( idFrag, idxData )
        self._dataMap = []
        
        # The list of all bonds in this block and all subblocks
        self._bonds = []
        
        # The list of atoms that are endGroups and their corresponding angleAtoms
        self._endGroups = []
        self._angleAtoms = []
   
        # Flag to indicate if block has changed and parameters (such as centerOfMass)
        # need to be changed
        self._changed = True
        
        # Below need to be updated when we move etc
        self._centroid = 0
        self._centerOfMass = 0
        self._maxAtomRadius = 0
        self._mass = 0
        
        self.atomCell = None # The cell that this atom is in (see cell.py)
        
        self._update()

    def alignBond(self, idxAtom, refVector ):
        """
        Align this block, so that the bond defined by idxAtom is aligned with
        the refVector
        
        This assumes that the block has already been positioned with the contact atom at the origin
        """
        
        endGroup = self.atomCoord( idxAtom )
        
        # Check neither is zero
        if numpy.array_equal( refVector, [0,0,0] ) or numpy.array_equal( endGroup, [0,0,0] ):
            raise RuntimeError, "alignBlock - one of the vectors is zero!\nrefVector: {0} endGroup: {1}".format( refVector, endGroup )
        
        # Makes no sense if they are already aligned
        if numpy.array_equal( endGroup/numpy.linalg.norm(endGroup), refVector/numpy.linalg.norm(refVector) ):
            print "alignBlock - block already aligned along vector. May not be a problem, but you should know..."
            return

        #print "alignBlock BEFORE: {0} | {1}".format( endGroup, refVector )

        # Calculate normalised cross product to find an axis orthogonal 
        # to both that we can rotate about
        cross = numpy.cross( refVector, endGroup )
        
        if numpy.array_equal( cross, [0,0,0] ):
            # They must be already aligned but anti-parallel, so we flip
            print "alignBlock - vectors are anti-parallel so flipping"
            self.flip( refVector )
        else:
            
            # Find angle
            angle = util.vectorAngle( refVector, endGroup )
            
            # Normalised cross to rotate about
            ncross = cross / numpy.linalg.norm( cross )
            
            # Rotate
            self.rotate( ncross, angle )
            
        #print "alignBlock AFTER: {0} | {1}".format( self.atomCoord( idxAtom ), refVector )
        return
    
    def angleAtom( self, idxAtom ):
        """Return the coordinate for the angleAtom of the given atom"""
        
        #return self._angleAtoms[ self._endGroups.index( idxAtom ) ]
        idxAngleAtom = self._angleAtoms[ self._endGroups.index( idxAtom ) ]
        return self.atomCoord( idxAngleAtom )
    
    def atomFragType(self, idxAtom ):
        """The type of the fragment that this atom belongs to."""
        idxFrag, idxData = self._dataMap[ idxAtom ]
        return self._fragments[ idxFrag ]._fragmentType
    
    def atomCoord(self, idxAtom ):
        idxFrag, idxData = self._dataMap[ idxAtom ]
        return self._fragments[ idxFrag ]._coords[ idxData ]
    
    def atomLabel(self, idxAtom ):
        idxFrag, idxData = self._dataMap[ idxAtom ]
        return self._fragments[ idxFrag ]._labels[ idxData ]
    
    def atomMass(self, idxAtom ):
        idxFrag, idxData = self._dataMap[ idxAtom ]
        return self._fragments[ idxFrag ]._masses[ idxData ]
    
    def atomRadius(self, idxAtom ):
        idxFrag, idxData = self._dataMap[ idxAtom ]
        return self._fragments[ idxFrag ]._atomRadii[ idxData ]
    
    def atomSymbol(self, idxAtom):
        idxFrag, idxData = self._dataMap[ idxAtom ]
        return self._fragments[ idxFrag ]._symbols[ idxData ]
        
    def bondBlock( self, idxAtom1, addBlock, idxAtom2 ):
        """ Add newBlock to this one
        """
        
        # Given a block and an atom index, we need to get the fragment
        
        idxFrag1, idxAtomFrag1 = self._dataMap[ idxAtom1 ]
        frag1 = self._fragments[ idxFrag1 ]
        
        idxFrag2, idxAtomFrag2 = addBlock._dataMap[ idxAtom2 ]
        frag2 = addBlock._fragments[ idxFrag2 ]
        
        bond = Bond()
        bond.block1 = self
        bond.frag1 = frag1
        bond.frag1atomIdx = idxAtomFrag1
        bond.frag1angleIdx = frag1.angleAtom( idxAtomFrag1 )
        
        bond.block2 = addBlock
        bond.frag2 = frag2
        bond.frag2atomIdx = idxAtomFrag2
        bond.frag2angleIdx = frag2.angleAtom( idxAtomFrag2 )
        
        self._myBonds.append( bond )
         
        self._update()
        
        return

    def _calcCenters(self): 
        """ Calc center for the combined block
        """
        
        self._centroid = 0
        self._centerOfMass = 0
        for frag in self._fragments:
            frag._calcCenters()
            self._centroid += frag._centroid
            self._centerOfMass += frag._centerOfMass
        return
        
    def _calcRadius(self):
        """
        Calculate a simple size metric of the block so that we can screen for whether
        two _blocks are within touching distance
        
        First try a simple approach with a loop just to get a feel for things
        - Find the largest distance between any atom and the center of geometry
        - Get the covalent radius of that atom
        - return that distance + radius + buffer
        
        Should move to use scipy as detailed here:
        http://stackoverflow.com/questions/6430091/efficient-distance-calculation-between-n-points-and-a-reference-in-numpy-scipy
        """
        
        cog = self.centroid()
        
        distances = []
        self._maxAtomRadius = 0
        for frag in self._fragments:
            if frag.maxAtomRadius() > self._maxAtomRadius:
                self._maxAtomRadius = frag.maxAtomRadius()
            for coord in frag._coords:
                #distances.append( numpy.linalg.norm(coord-cog) )
                distances.append( util.distance(cog, coord) )
            
        imax = numpy.argmax( distances )
        dist = distances[ imax ]
        
        # Add on the radius of the largest atom
        atomR = self.maxAtomRadius()
        
        # Set radius
        self._radius = dist + atomR
        return
        
    def centroid(self):
        """
        Return or calculate the center of geometry for the building block.
        """
        
        if self._changed:
            self._calcCenters()
            self._changed = False
        
        return self._centroid
    
        
    def centerOfMass(self):
        """
        Return or calculate the center of mass for this building block.
        """
        if self._changed:
            self._calcCenters()
            self._changed = False
        
        return self._centerOfMass

    def coords( self ):
        """Return a list of the coordinates for this block"""
        coords = []
        for frag in self._fragments:
            coords += frag._coords
        return coords

    def copy( self ):
        """
        Create a copy of ourselves and return it.
        Hit problem when passing in the distance function from the cell as the deepcopy
        also had a copy of the cell
        """
        return copy.deepcopy(self)

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

    def _getBonds(self):
        """Return all bonds for this block"""
        self._bonds = []
        for bond in self._myBonds:
            
            # First add the bonds to this block
            self._bonds.append( bond )
            
            # Then add the bonds in the bonded blocks - excluding those to ourselves
            if bond.block2 != self:
                self._bonds += bond.block2._getBonds()
            
        return self._bonds
        
    def hasFragmentType(self, fragmentType ):
        """Return True if this atom is an endGroup - doesn't check if free"""
        return fragmentType in self._fragmentTypes
    
    def isEndGroup(self, idxAtom):
        """Return True if this atom is an endGroup - doesn't check if free"""
        return idxAtom in self._endGroups

    def iterCoord(self):
        """Generator to return the coordinates"""
        for idx in range( len(self._dataMap) ):
            yield self.atomCoord( idx )

    def labels( self ):
        """Return a list of the coordinates for this block"""
        labels = []
        for frag in self._fragments:
            labels += frag._labels
        return labels

    def maxAtomRadius(self):
        """Return the maxium atom radius
        """
        
        if self._maxAtomRadius != 0:
            return self._maxAtomRadius
        
        raise RuntimeError,"maxAtomRadius not calculated!"
    
    def mass(self):
        """Total mass of the molecule
        """
        if not self._mass:
            mass = 0.0
            for frag in self._fragments:
                for m in frag._masses:
                    mass += m
            self._mass = mass
            
        return self._mass

    def newBondPosition(self, idxAtom, symbol ):
        """Return the position where a bond to an atom of type 'symbol'
        would be placed if bonding to the target endgroup
         I'm sure this algorithm is clunky in the extreme...
        """
        
        targetEndGroup = self.atomCoord( idxAtom )
        targetAngleAtom = self.angleAtom( idxAtom )
        targetSymbol = self.atomSymbol( idxAtom )
        
        # Get the bond length between these two atoms
        bondLength = util.bondLength( targetSymbol, symbol )
        
        # vector from target EndGroup contact to Endgroup
        bondVec =  targetEndGroup - targetAngleAtom 
        
        # Now extend it by Bondlength
        vLength = numpy.linalg.norm( bondVec )
        ratio = ( vLength + bondLength ) / vLength
        newVec = bondVec * ratio
        diff = newVec-bondVec
        newPosition = targetEndGroup + diff
        
        return newPosition

    def positionGrowBlock( self, idxAtom, growBlock, idxGrowAtom ):
        """
        Position growBlock so it can bond to us, using the given _endGroups
        
        Arguments:
        idxAtom: the index of the endGroup atom to use for the bond
        growBlock: the block we are positioning
        idxGrowAtom: the index of the endGroup to use for the bond
        """

        # The vector we want to align along is the vector from the angleAtom
        # to the endGroup
        endGroup = self.atomCoord( idxAtom )
        angleAtom = self.angleAtom( idxAtom )
        refVector =  endGroup - angleAtom
        
        # Get the angleAtom for the growBlock
        growBlockAngleAtom = growBlock.angleAtom( idxGrowAtom )
        
        # get the coord where the next block should bond
        # symbol of endGroup tells us the sort of bond we are making which determines
        # the bond length
        symbol = growBlock.atomSymbol( idxGrowAtom )
        bondPos = self.newBondPosition( idxAtom, symbol )
        #print "got bondPos: {}".format( bondPos )
        
        # Shift block so angleAtom at center, so the vector of the endGroup can be
        # aligned with the refVector
        growBlock.translate( -growBlockAngleAtom )
        
        # Align along the staticBlock bond
        growBlock.alignBond( idxGrowAtom, refVector )
        
        # Now turn the second block around
        # - need to get the new vectors as these will have changed due to the moves
        growBlockEG = growBlock.atomCoord( idxGrowAtom )
        
        # Flip by 180 along bond axis
        growBlock.flip( growBlockEG )
        
        # Now need to place the endGroup at the bond coord - at the moment
        # the contact atom is on the origin so we need to subtract the vector to the endGroup
        # from the translation vector
        #jmht FIX!
        growBlock.translate( bondPos + growBlockEG )
        
        return

    def radius(self):
        
        if self._changed:
            self._calcCenters()
            self._calcRadius()
            self._changed = False
        
        return self._radius

    def randomEndGroup( self, fragmentTypes=None ):
        """Randomly select at atom that is an endGroup - return index in global array"""
        
        MAXCOUNT=50 # to make sure we don't loop forever...
        count=0
        while count < MAXCOUNT:
            count += 1
            if count > MAXCOUNT:
                return False
            
            endGroupIdx = random.choice( self._endGroups )
            if fragmentTypes:
                if not self.atomFragType( endGroupIdx ) in fragmentTypes:
                    continue
                
            return endGroupIdx

    def randomRotate( self, origin=[0,0,0], atOrigin=False ):
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
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
        if center==None:
            center = numpy.array([0,0,0])
        
        rmat = util.rotation_matrix( axis, angle )
        
        # loop through all _blocks and change all _coords
        # Need to check that the center bit is the best way of doing this-
        # am almost certainly doing more ops than needed
        for frag in self._fragments:
            for i,coord in enumerate( frag._coords ):
                #self._coords[i] = numpy.dot( rmat, coord )
                coord = coord - center
                coord = numpy.dot( rmat, coord )
                frag._coords[i] = coord + center
        
        self._changed = True

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # Make sure its a numpy array - is this needed?
        if isinstance(tvector,list):
            tvector = numpy.array( tvector )

        for frag in self._fragments:
            # Use len as we don't need to return the _coords
            for i in range( len (frag._coords ) ):
                frag._coords[i] += tvector
        
        self._changed = True
    
    def translateCentroid( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given coord
        """
        self.translate( position - self.centroid() )
        
    def _update( self  ):
        """Set the list of _endGroups & update data for new block"""
        
        # get list of all fragments
        # recursively loop across all local bonds
        self._bonds = self._getBonds()
        
        # From the list of all bonds get a list of all blocks - excluding this one
        blocks = []
        for bond in self._bonds:
            if bond.block2 not in blocks and bond.block2 != self:
                blocks.append( bond.block2 )
        
        # Now get a list of all fragments
        self._fragments = [ self._rootFragment ]
        
        # - all fragments or just root?
        for block in blocks:
            self._fragments.append( block._rootFragment )
            
        # Now build up the dataMap
        self._dataMap = []
        self._endGroups = []
        self._angleAtoms = []
        self._fragmentTypes = []
        count=0
        for i, fragment in enumerate( self._fragments ):
            
            maxr = fragment.maxAtomRadius()
            if maxr > self._maxAtomRadius:
                self._maxAtomRadius = maxr
            
            if fragment._fragmentType not in self._fragmentTypes:
                self._fragmentTypes.append( fragment._fragmentType )
                
            for j in range( len(fragment._coords) ):
                
                self._dataMap.append( (i, j) )
                
                # Add endGroups - need to exclude those used in bonds
                if fragment.isEndGroup( j ):
                    
                    self._endGroups.append( count )
                    
                    # Get index of the angleAtom in the fragment
                    k = fragment.angleAtom( j )
                    # Now add that on starting from where the atoms for this fragment start
                    self._angleAtoms.append( count-j+k )
                
                count += 1
        
        # Have dataMap so now set bond indices
        for bond in self._bonds:
            
            idxFrag1 = self._fragments.index( bond.frag1 )
            bond.atom1Idx = self._dataMap.index( (idxFrag1, bond.frag1atomIdx) )
            bond.angle1Idx = self._dataMap.index( (idxFrag1, bond.frag1angleIdx) )
            
            idxFrag2 = self._fragments.index( bond.frag2 )
            bond.atom2Idx = self._dataMap.index( (idxFrag2, bond.frag2atomIdx) )
            bond.angle2Idx = self._dataMap.index( (idxFrag2, bond.frag2angleIdx) )
            
            # Remove these atoms from our endGroups
            i = self._endGroups.index( bond.atom1Idx )
            self._endGroups.pop( i )
            i = self._endGroups.index( bond.atom2Idx )
            self._endGroups.pop( i )
            
        
        # Fill atomCell list
        self.atomCell = [None]*len( self._dataMap )
        
        # Recalculate the data for this new block
        self._calcCenters()
        self._calcRadius()
        
        return
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
            
        
        coords = self.coords()
        labels = self.labels()
        
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( len(coords) ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
                             
            for i,c in enumerate(coords):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol(labels[i]), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
    def __str__(self):
        """
        Return a string representation of the block
        """
        
        mystr = "Block: {0}\n".format(self.__repr__())
        
        mystr += "Num fragments: {0}\n".format( len( self._fragments) )
        mystr += "endGroups: {0}\n".format( self._endGroups )
        mystr += "angleAtoms: {0}\n".format( self._angleAtoms )
        mystr += "bonds: {0}\n".format( self._bonds )
        
#         mystr = ""
#         mystr += "BlockID: {}\n".format(id(self))
#         
#         coords = []
#         labels = []
#         for block in self.allBlocks:
#             coords += block._coords
#             labels += block._labels
#         mystr += "{}\n".format( len(coords) )
#             
#         for i,c in enumerate( coords ):
#             #mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self._labels[i], c[0], c[1], c[2])
#             mystr += "{0:5} {1:0< 15} {2:0< 15} {3:0< 15} \n".format( labels[i], c[0], c[1], c[2])
#             
#         #mystr += "radius: {}\n".format( self.radius() )
#         mystr += "COM: {}\n".format( self._centerOfMass )
#         mystr += "COG: {}\n".format( self._centroid )
#         mystr += "_endGroups: {}\n".format( self._endGroups )
#         mystr += "_myEndGroups: {}\n".format( self._myEndGroups )
#         mystr += "_myAngleAtoms: {}\n".format( self._myAngleAtoms )
#         mystr += "_blocks: {}\n".format( self._blocks )
#         mystr += "allBlocks: {}\n".format( self.allBlocks )

        return mystr
    
class TestBlock(unittest.TestCase):
    
    
    def writeCombined(self, block1, block2, filename ):
        """Write an xyz file with the coordinates of both blocks"""
        
        block1coords = block1.coords()
        block1labels = block1.labels()
    
        block2coords = block2.coords()
        block2labels = block2.labels()
        
        natoms = len( block1coords ) + len( block2coords )
        
        f = open( filename,'w')
        f.write("{0}\n".format( natoms ) )
        f.write("Combined coords\n")
        for i, c in enumerate( block1coords ):
            f.write( "{0}    {1}    {2}    {3}\n".format( util.label2symbol(block1labels[i]), c[0], c[1], c[2] ) )
        for i, c in enumerate( block2coords ):
            f.write( "{0}    {1}    {2}    {3}\n".format( util.label2symbol(block2labels[i]), c[0], c[1], c[2] ) )
            
        f.close()
        return
    
    def testCH4(self):
        """Test the creation of a CH4 molecule"""
        
        infile="../ch4.xyz"
        ch4 = Block( infile )
        
        endGroups = [ 1, 2, 3, 4 ]
        self.assertEqual( endGroups, ch4._endGroups )
        angleAtoms = [ 0, 0, 0, 0, ]
        self.assertEqual( angleAtoms, ch4._angleAtoms )
        
        return

    def testCH4_bond(self):
        """First pass"""
        
        infile="../ch4.xyz"
        ch4_1 = Block( infile )
        
        ch4_2 = Block( infile )
        
        idxAtom1=2
        idxAtom2=3
        ch4_1.bondBlock( idxAtom1, ch4_2, idxAtom2 )
        
        endGroups = [ 1, 3, 4, 6, 7, 9 ]
        self.assertEqual( endGroups, ch4_1._endGroups )
        
        angleAtoms = [0, 0, 0, 0, 5, 5, 5, 5]
        self.assertEqual( angleAtoms, ch4_1._angleAtoms )
        
        ref_bonds = [ (2, 8) ]
        bonds = [ ( b.atom1Idx, b.atom2Idx ) for b in ch4_1._bonds ]
        self.assertEqual( ref_bonds, bonds )
        
        return
    
    def testCH4_bond2(self):
        """First pass"""
        
        infile="../ch4.xyz"
        ch4_1 = Block( infile )
        ch4_2 = Block( infile )
        ch4_3 = Block( infile )
        
        idxAtom1=2
        idxAtom2=3
        ch4_1.bondBlock( idxAtom1, ch4_2, idxAtom2 )
        
        idxAtom1=9
        idxAtom2=1
        ch4_1.bondBlock( idxAtom1, ch4_3, idxAtom2 )
        
        endGroups = [ 1, 3, 4, 6, 7, 12, 13, 14 ]
        self.assertEqual( endGroups, ch4_1._endGroups )
        
        angleAtoms = [0, 0, 0, 0, 5, 5, 5, 5, 10, 10, 10, 10]
        self.assertEqual( angleAtoms, ch4_1._angleAtoms )
        
        ref_bonds = [ (2, 8), (9, 11) ]
        bonds = [ ( b.atom1Idx, b.atom2Idx ) for b in ch4_1._bonds ]
        self.assertEqual( ref_bonds, bonds )
        
        return

    def testCH4_bond_self(self):
        """First pass"""
        
        
        infile="../ch4.xyz"
        ch4_1 = Block( infile )
        ch4_2 = Block( infile )
        
        idxAtom1=1
        idxAtom2=1
        ch4_1.bondBlock( idxAtom1, ch4_2, idxAtom2 )
        
        # now bond to self
        idxAtom1=2
        idxAtom2=9
        ch4_1.bondBlock( idxAtom1, ch4_1, idxAtom2 )
        
        endGroups = [3, 4, 7, 8]
        self.assertEqual( endGroups, ch4_1._endGroups )
         
        angleAtoms = [0, 0, 0, 0, 5, 5, 5, 5]
        self.assertEqual( angleAtoms, ch4_1._angleAtoms )
         
        ref_bonds = [ (1, 6), (2, 9) ]
        bonds = [ ( b.atom1Idx, b.atom2Idx ) for b in ch4_1._bonds ]
        self.assertEqual( ref_bonds, bonds )
         
        return

    def testAlignBlocks(self):
        """Test we can align two _blocks correctly"""
    
        blockS = Block( "../PAF_bb_typed.car" )
        block = blockS.copy()
        
        block.translateCentroid( [ 3, 4 ,5 ] )
        block.randomRotate()
        
        # Get the atoms that define things
        idxSatom = 17
        blockSEndGroup = blockS.atomCoord( idxSatom )
        blockSangleAtom = blockS.angleAtom( idxSatom )
        
        idxAtom = 7
        blockAngleAtom = block.angleAtom( idxAtom )
        
        # we want to align along block1Contact -> block1EndGroup
        refVector = blockSEndGroup - blockSangleAtom
        
        # Position block so contact is at origin
        block.translate( -blockAngleAtom )

        block.alignBond( idxAtom, refVector )
        
        # Check the relevant atoms are in the right place
        blockEndGroup  = block.atomCoord( idxAtom )
        blockAngleAtom = block.angleAtom( idxAtom )
        
        newVector = blockEndGroup - blockAngleAtom
        
        # Normalise two vectors so we can compare them
        newNorm = newVector / numpy.linalg.norm(newVector)
        refNorm = refVector / numpy.linalg.norm(refVector)
        
        # Slack tolerances - need to work out why...
        self.assertTrue( numpy.allclose( newNorm, refNorm),
                         msg="End Group incorrectly positioned: {0} | {1}".format(newNorm, refNorm )  )
        return

    def testCenterOfGeometry(self):
        """
        Test calculation of Center of Geometry
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        ch4 = Block("../ch4.xyz")
        cog = ch4.centroid()
        self.assertTrue( numpy.allclose( correct, cog, rtol=1e-9, atol=1e-6 ),
                         msg="testCenterOfGeometry incorrect COM.")

    def testCenterOfMass(self):
        """
        Test calculation of Center of Mass
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        ch4 = Block("../ch4.xyz")
        com = ch4.centerOfMass()
        self.assertTrue( numpy.allclose( correct, com, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        return

    def testMove(self):
        """Test we can move correctly"""
        
        paf = Block( "../PAF_bb_typed.car" )
        m = paf.copy()
        m.translate( numpy.array( [5,5,5] ) )
        c = m.centroid()
        paf.translateCentroid( c )
        p = paf.centroid()
        
        self.assertTrue( numpy.allclose( p, c, rtol=1e-9, atol=1e-9 ), "simple move")
        return
    
    def testPositionGrowBlock(self):
        
        blockS = Block( "../PAF_bb_typed.car" )
        growBlock = blockS.copy()
        
        growBlock.translateCentroid( [ 3,4,5 ] )
        growBlock.randomRotate()
        
        #self.writeCombined(blockS, growBlock, "b4.xyz")
        
        # Get the atoms that define things
        idxBlockAtom = 7
        idxGrowBlockAtom = 12
        
        # Get position to check
        newPos = blockS.newBondPosition( idxBlockAtom, growBlock.atomSymbol( idxGrowBlockAtom ) )
        
        # Position the block
        blockS.positionGrowBlock( idxBlockAtom, growBlock, idxGrowBlockAtom )
        
        # After move, the endGroup of the growBlock should be at newPos
        endGroupCoord = growBlock.atomCoord( idxGrowBlockAtom )
        self.assertTrue( numpy.allclose( newPos, endGroupCoord, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")
        
        #self.writeCombined(blockS, growBlock, "after.xyz")

        return

    def testRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = Block("../ch4.xyz")
        r = ch4.radius()
        #jmht - check...- old was: 1.78900031214
        # or maybe: 1.45942438719
        self.assertAlmostEqual(r, 1.79280605406, 7, "Incorrect radius: {}".format(str(r)) )
        return
        
    def XtestMaxAtomRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = Block("../ch4.xyz")
        r = ch4.maxAtomRadius()
        #jmht - check...- old was: 1.78900031214
        self.assertAlmostEqual(r, 0.70380574117, 7, "Incorrect radius: {}".format(str(r)) )
        
    def testRotate(self):
        """
        Test the rotation
        """
        
        ch4 = Block("../ch4.xyz")
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
        