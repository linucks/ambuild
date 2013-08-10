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
import random
import unittest

# external imports
import numpy

# local imports
import util


class Fragment(object):
    '''
    classdocs
    '''

    def __init__( self, filename=None, fragmentType=None ):
        '''
        Constructor
        '''
        
        # the coordinates for this block
        self._coords = []
        
        # ordered array of _labels
        self._labels = []
        
        # ordered array of _symbols (in upper case)
        self._symbols = []
        
        # ordered array of _masses
        self._masses = []
        
        # orderd array of atom radii
        self._atomRadii = []
        
        # List of the cell (3-tuple) to which each atom belongs
        self.atomCell = []
        
        # A list of block, endGroup indices of the atoms that are _endGroups
        self._endGroups = []
        
        # The type this fragment is (for bonding checks)
        self._fragmentType = fragmentType
        
        # Holds the center of mass for the combined block
        self._centerOfMass = numpy.zeros( 3 )
        
        # Holds the centroid of the combined block
        self._centroid = numpy.zeros( 3 )
        
        # The radius of the combined block assuming it is a circle centered on the centroid
        self._radius = 0
        
        # maximum Atom Radius
        self._maxAtomRadius = 0
        
        # Total mass of the combined block
        self._mass = None
        
        # Flag to indicate if block has changed and parameters (such as centerOfMass)
        # need to be changed
        self._changed = True
        
        if filename:
            if filename.endswith(".car"):
                self.fromCarFile( filename )
            elif filename.endswith(".xyz"):
                self.fromXyzFile( filename )
            else:
                raise RuntimeError("Unrecognised file suffix: {}".format(filename) )

    def angleAtom(self, endGroup ):
        """Return the index of the angleAtom for this endGroup"""
        if not self.isEndGroup( endGroup ):
            raise RuntimeError,"Index {0} is not an endGroup".format( endGroup )
        i = self._endGroups.index( endGroup )
        return self._angleAtoms[ i ]
        
    def createFromArgs(self, coords, labels, endGroups, angleAtoms ):
        """ Create from given arguments
        """
        # could check if numpy array here
        self._coords = coords
        self._labels = labels
        self._endGroups = endGroups
        self._angleAtoms = angleAtoms
        
        self.fillData()
        
        self.update()
        
    def fromLabelAndCoords(self, labels, coords):
        """ Given an array of _labels and another of _coords, create a block
        This requires determining the end groups and contacts from the label
        """
        
        # array of indexes
        endGroups = []
        egAaLabel = []
        
        # For tracking the mapping of the label used to mark the angle atom
        # to the true index of the atom
        aaLabel2index = {}
        
        for i, label in enumerate(labels):
            
            # End groups and the atoms which define their bond angles 
            # are of the form XX_EN for endgroups and XX_AN for the defining atoms
            # where N can be any number, but linking the two atoms together 
            # The indexing of the atoms starts from 1!!!!
            if "_" in label:
                _, ident = label.split("_")
                atype=ident[0]
                anum=int(ident[1:])
                
                if atype.upper() == "E":
                    # An endGroup
                    if i in endGroups:
                        raise RuntimeError,"Duplicate endGroup key for: {0} - {1}".format( i, endGroups )
                    endGroups.append(i)
                    egAaLabel.append( anum )
                    
                elif atype.upper() == "A":
                    if aaLabel2index.has_key(anum):
                        raise RuntimeError,"Duplicate angleAtom key for: {0} - {1}".format(anum, aaLabel2index)
                    aaLabel2index[ anum ] = i
                else:
                    raise RuntimeError,"Got a duff label! - {}".format(label)
        
#        #Now match up the endgroups with their angle atoms
        angleAtoms = []
        for label in egAaLabel:
            angleAtoms.append( aaLabel2index[ label ]  )
        
        #self.createFromArgs(_coords, _labels, endGroups, endGroupContacts)
        self.createFromArgs( coords, labels, endGroups, angleAtoms )
  
    def _calcCenters(self):
        """Calculate the center of mass and geometry for this fragment
        """
        
        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        totalMass = 0.0
        for i, coord in enumerate( self._coords ):
            mass = self._masses[i]
            totalMass += mass
            sumG += coord
            sumM += mass * coord
        
        self._centroid = sumG / (i+1)
        self._centerOfMass = sumM / totalMass
        
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
        for coord in self._coords:
            #distances.append( numpy.linalg.norm(coord-cog) )
            distances.append( util.distance(cog, coord) )
            
        imax = numpy.argmax( distances )
        dist = distances[ imax ]
        
        # Add on the radius of the largest atom
        atomR = self.maxAtomRadius()
        
        # Set radius
        self._radius = dist + atomR
        
        
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

    def fillData(self):
        """ Fill the data arrays from the label """
        
        if len(self._labels) != len(self._coords):
            raise RuntimeError("fillData needs _labels filled!")
        
        # now in -update
        # Fill atomCell list
        #self.atomCell = [None]*len(self._coords)
        
        symbol_types=[]
        for label in self._labels:
            
            # Symbols
            symbol = util.label2symbol( label )
            self._symbols.append( symbol )
            
            # For checking bonding
            if symbol not in symbol_types:
                symbol_types.append(symbol)
                
            # Masses
            self._masses.append( util.ATOMIC_MASS[ symbol ] )
            
            # Radii
            z = util.SYMBOL_TO_NUMBER[ symbol ]
            r = util.COVALENT_RADII[z] * util.BOHR2ANGSTROM
            self._atomRadii.append(r)
            # Remember the largest
            if r > self._maxAtomRadius:
                self._maxAtomRadius = r
            #print "ADDING R {} for label {}".format(r,label)
            

    def fromCarFile(self, carFile):
        """"Abbie did this.
        Gulp...
        """
        labels = []
        
        # numpy array
        coords = []
        
        reading = True
        with open( carFile, "r" ) as f:
            
            # First 4 lines just info - not needed
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            
            count=0
            while reading:
                line = f.readline()
                
                line = line.strip()
                fields = line.split()
                label = fields[0]
                
                # Check end of coordinates
                if label.lower() == "end":
                    reading=False
                    break
                     
                labels.append(label)
                
                coords.append( numpy.array(fields[1:4], dtype=numpy.float64) )
                
                count+=1
        
        self.fromLabelAndCoords( labels, coords )

    def fromXyzFile(self, xyzFile ):
        """"Jens did this.
        """
        
        labels = []
        
        # numpy array
        coords = []
        
        with open( xyzFile ) as f:
            
            # First line is number of atoms
            line = f.readline()
            natoms = int(line.strip())
            
            # Skip title
            line = f.readline()
            
            count = 0
            for _ in range(natoms):
                
                line = f.readline()
                line = line.strip()
                fields = line.split()
                label = fields[0]
                labels.append(label) 
                coords.append( numpy.array(fields[1:4], dtype=numpy.float64) )
                
                count += 1

        self.fromLabelAndCoords( labels, coords )
        
    def isEndGroup(self, idxAtom):
        """Return True if this atom is an endGroup - doesn't check if free"""
        return idxAtom in self._endGroups

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
            for m in self._masses:
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

    def rotate( self, axis, angle, center=None ):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
        if center==None:
            center = numpy.array([0,0,0])
        
        rmat = util.rotation_matrix( axis, angle )
        
        # loop through all _blocks and change all _coords
        # Need to check that the center bit is the best way of doing this-
        # am almost certainly doing more ops than needed
        for i,coord in enumerate( self._coords ):
            #self._coords[i] = numpy.dot( rmat, coord )
            coord = coord - center
            coord = numpy.dot( rmat, coord )
            self._coords[i] = coord + center
        
        self._changed = True

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # Make sure its a numpy array - is this needed?
        if isinstance(tvector,list):
            tvector = numpy.array( tvector )

        # Use len as we don't need to return the _coords
        for i in range( len (self._coords ) ):
            self._coords[i] += tvector
        
        self._changed = True
    
    def translateCentroid( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given coord
        """
        self.translate( position - self.centroid() )
        
        
    def update(self):
        # Recalculate the dynamic data
        self._calcCenters()
        self._calcRadius()

class TestFragment(unittest.TestCase):

    def makeCh4(self):
        """Create a CH4 molecule for testing"""
        
        ch4 = Fragment()
        ch4.fromXyzFile("../ch4.xyz")
        endGroups = [ 1, 2, 3, 4 ]
        self.assertEqual( endGroups, ch4._endGroups )
        self.assertEqual( ch4._endGroups, ch4._endGroups )
        angleAtoms = [ 0, 0, 0, 0, ]
        self.assertEqual( angleAtoms, ch4._angleAtoms )
        return ch4
    
    def makePaf(self):
        """Return the PAF molecule for testing"""
        
        paf = Fragment()
        paf.fromCarFile("../PAF_bb_typed.car")
        return paf
    
    def testAaaReadCar(self):
        """
        Test we can read a car file - needs to come first
        """
        
        paf = self.makePaf()
        self.assertTrue( paf._endGroups == [7, 12, 17, 22],
                         "Incorrect reading of endGroup contacts: {0}".format(paf._endGroups))

    def testAaaReadXyz(self):
        """
        Test we can read an xyz file - needs to come first
        """
        
        ch4 = self.makeCh4()
        self.assertTrue( ch4._endGroups == [ 1, 2, 3, 4 ],
                         "Incorrect reading of endGroup contacts: {0}".format(ch4._endGroups))

    def testCenterOfGeometry(self):
        """
        Test calculation of Center of Geometry
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        ch4 = self.makeCh4()
        cog = ch4.centroid()
        self.assertTrue( numpy.allclose( correct, cog, rtol=1e-9, atol=1e-6 ),
                         msg="testCenterOfGeometry incorrect COM.")

    def testCenterOfMass(self):
        """
        Test calculation of Center of Mass
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        ch4 = self.makeCh4()
        com = ch4.centerOfMass()
        self.assertTrue( numpy.allclose( correct, com, rtol=1e-9, atol=1e-7 ),
                         msg="testCenterOfMass incorrect COM.")

    def testMove(self):
        """Test we can move correctly"""
        
        paf = self.makePaf()
        m = paf.copy()
        m.translate( numpy.array( [5,5,5] ) )
        c = m.centroid()
        paf.translateCentroid( c )
        p = paf.centroid()
        
        self.assertTrue( numpy.allclose( p, c, rtol=1e-9, atol=1e-9 ), "simple move")
        
    def testRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = self.makeCh4()
        r = ch4.radius()
        #jmht - check...- old was: 1.78900031214
        # or maybe: 1.45942438719
        self.assertAlmostEqual(r, 1.79280605406, 7, "Incorrect radius: {}".format(str(r)) )
        
    def testMaxAtomRadius(self):
        """
        Test calculation of the radius
        """
        
        ch4 = self.makeCh4()
        r = ch4.maxAtomRadius()
        #jmht - check...- old was: 1.78900031214
        self.assertAlmostEqual(r, 0.70380574117, 7, "Incorrect radius: {}".format(str(r)) )
        
    def testRotate(self):
        """
        Test the rotation
        """
        
        ch4 = self.makeCh4()
        array1 = numpy.array([ -0.51336 ,  0.889165, -0.363 ])
        self.assertTrue( numpy.array_equal( ch4._coords[4], array1 ),
                         msg="testRotate arrays before rotation incorrect.")
        
        axis = numpy.array([1,2,3])
        angle = 2
        ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose( ch4._coords[4], array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")

    def makeH2C2(self):
        """
        Return H2C2 molecule for testing
        """
        coords = [ 
                  numpy.array([  0.000000,  0.000000,  0.000000 ] ),
                  numpy.array([  0.000000,  0.000000,  1.000000 ]),
                  numpy.array([  0.000000,  0.000000,  2.000000 ]),
                  numpy.array([  0.000000,  0.000000,  3.000000 ]),
                  ]

        labels = [ 'H', 'C', 'C', 'H' ]
        endGroups = [ 0, 3 ]
        angleAtoms = [ 1, 2 ]
        b = Fragment()
        b.createFromArgs( coords, labels, endGroups, angleAtoms )
        return b

    def XtestAlignBlocks(self):
        """Test we can align two _blocks correctly"""
    
        blockS = self.makePaf()
        block = blockS.copy()
        
        block.translateCentroid( [3,4,5] )
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

    def XtestPositionGrowBlock(self):
        
        blockS = self.makePaf()
        growBlock = blockS.copy()
        
        growBlock.translateCentroid( [3,4,5] )
        growBlock.randomRotate()
        
#        testb = blockS.copy()
#        testb += growBlock
#        testb.writeXyz("b4.xyz")
        
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
        
#        testb = blockS.copy()
#        testb += growBlock
#        testb.writeXyz("after.xyz")


if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
        