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
import util

class Block():
    '''
    classdocs
    '''

    def __init__( self, infile=None ):
        '''
        Constructor
        '''
        
        # List of all the _blocks that constitute this one ordered by endGroups
        # so a block bonded at the endGroup at index I is block I
        self._blocks = []
        
        # List of all _blocks contained or linked to this one - includes iterating over sub-_blocks
        # We add ourselves to this list
        self.allBlocks = [ self ]
        
        # The id of the parent block of this one (if we have one)
        self.parent = None
        
        # list of tuples of ( idxBlock, idxData ) mapping the block and index into the block array where the data lives
        self._dataMap = []
        
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
        
        # The index of the local endGroups
        self._myEndGroups = []
        
        # List of the indices of the atoms that define the angle for the _endGroups - need
        # to be in the same order as the _endGroups
        self._myAngleAtoms = []
        
        # Holds the center of mass for the combined block
        self._centerOfMass = numpy.zeros( 3 )
        
        # Holds the center of mass for this block
        self._myCenterOfMass = numpy.zeros( 3 )
        
        # Holds the centroid of the combined block
        self._centroid = numpy.zeros( 3 )
        
        # Holds the centroid of this block
        self._myCentroid = numpy.zeros( 3 )
        
        # The radius of the combined block assuming it is a circle centered on the centroid
        self._radius = 0
        
        # The radius of the combined block assuming it is a circle centered on the centroid
        self._myRadius = 0
        
        # maximum Atom Radius
        self._maxAtomRadius = 0
        self._myMaxAtomRadius = 0
        
        # Total mass of the combined block
        self._mass = None
        
        # Total mass of this block
        self._myMass = None
        
        # Flag to indicate if block has changed and parameters (such as centerOfMass)
        # need to be changed
        self._changed = True
        
        if infile:
            if infile.endswith(".car"):
                self.fromCarFile( infile )
            elif infile.endswith(".xyz"):
                self.fromXyzFile( infile )
            else:
                raise RuntimeError("Unrecognised file suffix: {}".format(infile) )

# new
    def atomCoord(self, idxAtom ):
        idxBlock, idxData = self._dataMap[ idxAtom ]
        return self.allBlocks[ idxBlock ]._coords[ idxData ]
    
    def atomLabel(self, idxAtom ):
        idxBlock, idxData = self._dataMap[ idxAtom ]
        return self.allBlocks[ idxBlock ]._labels[ idxData ]
    
    def atomMass(self, idxAtom ):
        idxBlock, idxData = self._dataMap[ idxAtom ]
        return self.allBlocks[ idxBlock ]._masses[ idxData ]
    
    def atomRadius(self, idxAtom ):
        idxBlock, idxData = self._dataMap[ idxAtom ]
        return self.allBlocks[ idxBlock ]._atomRadii[ idxData ]
    
    def atomSymbol(self, idxAtom):
        idxBlock, idxData = self._dataMap[ idxAtom ]
        return self.allBlocks[ idxBlock ]._symbols[ idxData ]
    
    def iterCoord(self):
        """Generator to return the coordinates"""
        for idx in range( len(self._dataMap) ):
            yield self.atomCoord( idx )

    def addBlock( self, idxAtom, addBlock, idxAddBlockAtom ):
        """ Add newBlock to this one
        """
        
        if addBlock.parent:
            raise RuntimeError, "addBlock has parent!"
        addBlock.parent = self
        
        # update endGroups for self 
        idxBlock, idxBlockAtom = self._dataMap[ idxAtom ]
        block = self.allBlocks[ idxBlock ]
        idxEndGroup = block._myEndGroups.index( idxBlockAtom )
        block._blocks[ idxEndGroup ] = addBlock
        
        # Add addBlock to the list of blocks contained in this one
        if not addBlock in self.allBlocks and not addBlock == self.parent:
            block.allBlocks += addBlock.allBlocks
        else:
            # Circular fragment
            raise RuntimeError, "addBlock  - addBlock was already in self"
        
        # update endGroups for addBlock
        idxBlock, idxBlockAtom = addBlock._dataMap[ idxAddBlockAtom ]
        block = addBlock.allBlocks[ idxBlock ]
        idxEndGroup = addBlock._myEndGroups.index( idxBlockAtom )
        block._blocks[ idxEndGroup ] = self
        
        self.update()

    def angleAtom( self, idxAtom ):
        """Return the coordinate for the angleAtom of the given atom"""
        idxBlock, idxBlockAtom = self._dataMap[ idxAtom ]
        block = self.allBlocks[ idxBlock ]
        idxEndGroup = block._myEndGroups.index( idxBlockAtom )
        idxAngleAtom = block._myAngleAtoms[ idxEndGroup ]
        return block._coords[ idxAngleAtom ]

    def alignBond(self, idxAtom, refVector ):
        """
        Align this block, so that the bond defined by idxBlockEG is aligned with
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

#        print "alignBlock BEFORE: {0} | {1}".format( endGroup, refVector )

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
            
#        endGroup =self._coords[ idxBlockEG ]
#        print "alignBlock AFTER: {0} | {1}".format( endGroup, refVector )
        return
        
    def createFromArgs(self, coords, labels, endGroups, angleAtoms ):
        """ Create from given arguments
        """
        # could check if numpy array here
        self._coords = coords
        self._labels = labels
        self._myEndGroups = endGroups
        self._myAngleAtoms = angleAtoms
        
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
  
    def _calcMyCenters(self):
        """Calculate the center of mass and geometry for this block
        """
        
        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        totalMass = 0.0
        for i, coord in enumerate( self._coords ):
            mass = self._masses[i]
            totalMass += mass
            sumG += coord
            sumM += mass * coord
        
        self._myCentroid = sumG / (i+1)
        self._myCenterOfMass = sumM / totalMass
        
    def _calcCenters(self): 
        """ Calc center for the combined block
        """
        
        self._centroid = 0
        self._centerOfMass = 0
        for block in self.allBlocks:
            block._calcMyCenters()
            self._centroid += block._myCentroid
            self._centerOfMass += block._myCenterOfMass
        
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
        for block in self.allBlocks:
            if self._myMaxAtomRadius > self._maxAtomRadius:
                self._maxAtomRadius = self._myMaxAtomRadius
            for coord in block._coords:
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
        
        # fill empty block list
        for _ in range( len( self._myEndGroups ) ):
            self._blocks.append( None )
        
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
            if r > self._myMaxAtomRadius:
                self._myMaxAtomRadius = r
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
            for block in self.allBlocks:
                for m in block._masses:
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

    def randomEndGroup(self):
        """Randomly select at atom that is an endGroup - return index in global array"""
        return random.choice( self._endGroups )

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

    def removeBlock( self, block ):
        """Remove the block from this one"""
        
        assert block in self.allBlocks
        
        block = self.allBlocks.pop( self.allBlocks.index( block ) )
        parent = block.parent
        assert block in parent._blocks
        idx = parent._blocks.index( block )
        parent._blocks[ idx ] = None
        self.update()
        return block

    def rotate( self, axis, angle, center=None ):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
        if center==None:
            center = numpy.array([0,0,0])
        
        rmat = util.rotation_matrix( axis, angle )
        
        # loop through all _blocks and change all _coords
        # Need to check that the center bit is the best way of doing this-
        # am almost certainly doing more ops than needed
        for block in self.allBlocks:
            for i,coord in enumerate( block._coords ):
                #self._coords[i] = numpy.dot( rmat, coord )
                coord = coord - center
                coord = numpy.dot( rmat, coord )
                block._coords[i] = coord + center
        
        self._changed = True

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # Make sure its a numpy array - is this needed?
        if isinstance(tvector,list):
            tvector = numpy.array( tvector )

        for block in self.allBlocks:
            # Use len as we don't need to return the _coords
            for i in range( len (block._coords ) ):
                block._coords[i] += tvector
        
        self._changed = True
    
    def translateCentroid( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given coord
        """
        self.translate( position - self.centroid() )
        
    def update( self ):
        if self.parent:
            self.parent.update()
        else:
            self._update()
    
    def _update( self  ):
        """Set the list of _endGroups & update data for new block"""
        
        # get list of all _blocks
        self.allBlocks = [ self ]
        for block in self._blocks:
            if block:
                assert block != self
                self.allBlocks += block.allBlocks
        
        self._dataMap = []
        self._endGroups = []
        count=0
        # Loop over every block contained in this and all subBlocks
        for i, block in enumerate( self.allBlocks ):
            for j in range( len(block._coords) ):
                self._dataMap.append( (i, j) )
                # Check if is a free endGroup
                if j in block._myEndGroups:
                    idx = block._myEndGroups.index( j )
                    if not block._blocks[ idx ]:
                        self._endGroups.append( count )
                count+=1
        
        # Fill atomCell list
        self.atomCell = [None]*len( self._dataMap )
                    
        # Recalculate the data for this new block
        self._calcCenters()
        self._calcRadius()
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
            
        
        coords = []
        labels = []
        
        for block in self.allBlocks:
            coords += block._coords
            labels += block._labels
            
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format( len(coords) ) )
            f.write( "id={}\n".format( str( id(self)  ) ) )
                             
            for i,c in enumerate(coords):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol(labels[i]), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
    def __str__(self):
        """
        Return a string representation of the molecule
        """
        
        mystr = ""
        mystr += "BlockID: {}\n".format(id(self))
        
        coords = []
        labels = []
        for block in self.allBlocks:
            coords += block._coords
            labels += block._labels
        mystr += "{}\n".format( len(coords) )
            
        for i,c in enumerate( coords ):
            #mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self._labels[i], c[0], c[1], c[2])
            mystr += "{0:5} {1:0< 15} {2:0< 15} {3:0< 15} \n".format( labels[i], c[0], c[1], c[2])
            
        #mystr += "radius: {}\n".format( self.radius() )
        mystr += "COM: {}\n".format( self._centerOfMass )
        mystr += "COG: {}\n".format( self._centroid )
        mystr += "_endGroups: {}\n".format( self._endGroups )
        mystr += "_myEndGroups: {}\n".format( self._myEndGroups )
        mystr += "_myAngleAtoms: {}\n".format( self._myAngleAtoms )
        mystr += "_blocks: {}\n".format( self._blocks )
        mystr += "allBlocks: {}\n".format( self.allBlocks )
        return mystr
    
    def __add__(self, other):
        """Add two _blocks - INCOMPLETE AND ONLY FOR TESTING"""
        self._coords += other._coords
        self._labels += other._labels
        self._symbols += other._symbols
        
        return self
        
        
class TestBuildingBlock(unittest.TestCase):

    def makeCh4(self):
        """Create a CH4 molecule for testing"""
        
#        _coords = [ numpy.array([  0.000000,  0.000000,  0.000000 ] ),
#        numpy.array([  0.000000,  0.000000,  1.089000 ]),
#        numpy.array([  1.026719,  0.000000, -0.363000 ]),
#        numpy.array([ -0.513360, -0.889165, -0.363000 ]),
#        numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]
#
#        #numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]
#        _labels = [ 'C', 'H', 'H', 'H', 'H' ]
#        
#        _endGroups = [ 1,2,3,4 ]
#        
#        endGroupContacts = { 1:0, 2:0, 3:0, 4:0 }
#        
#        ch4 = Block()
#        ch4.createFromArgs( _coords, _labels, _endGroups, endGroupContacts )

        ch4 = Block()
        ch4.fromXyzFile("../ch4.xyz")
        
        endGroups = [ 1, 2, 3, 4 ]
        self.assertEqual( endGroups, ch4._myEndGroups )
        self.assertEqual( ch4._endGroups, ch4._myEndGroups )
        angleAtoms = [ 0, 0, 0, 0, ]
        self.assertEqual( angleAtoms, ch4._myAngleAtoms )
        return ch4
    
    def testMakeCh4_3(self):
        
        ch4_1 = self.makeCh4()
        
        idxEG1 = 0
        idxEG2 = 1
        idxEG3 = 2
        idxEG4 = 3
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()
        
        #ch4_1.positionGrowBlock(idxEG1, ch4_2, idxEG2 )
        #ch4_2.positionGrowBlock(idxEG3, ch4_3, idxEG4 )
#        testb = ch4_1.copy()
#        testb += ch4_2
#        testb += ch4_3
#        testb.writeXyz("join.xyz")

        
        ch4_1.addBlock(1, ch4_2, 3)
        #print ch4_1
        
    def testAdd1(self):
        
        ch4_1 = self.makeCh4()
        ch4_2 = ch4_1.copy()
        
        idxAtom = 2
        idxGrowBlockAtom = 3
        ch4_1.addBlock( idxAtom, ch4_2, idxGrowBlockAtom )
        
        endGroups = [ 1, 3, 4, 6, 7, 9 ]
        self.assertEqual( endGroups, ch4_1._endGroups )
        
        blocks = [ None, ch4_2, None, None ]
        self.assertEqual( blocks, ch4_1._blocks )
        
        allBlocks = [ ch4_1, ch4_2 ]
        self.assertEqual( allBlocks, ch4_1.allBlocks )
        
    def testAdd2(self):
        
        ch4_1 = self.makeCh4()
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()
        
        idxEndAtom = 2
        idxGrowBlockAtom = 3
        ch4_1.positionGrowBlock( idxEndAtom, ch4_2, idxGrowBlockAtom )
        ch4_1.addBlock( idxEndAtom, ch4_2, idxGrowBlockAtom )
        
        
        idxEndAtom = 9
        idxGrowBlockAtom = 1
        ch4_1.positionGrowBlock(idxEndAtom, ch4_3, idxGrowBlockAtom )
        ch4_1.addBlock( idxEndAtom, ch4_3, idxGrowBlockAtom )
        
        endGroups = [ 1, 3, 4, 6, 7, 12, 13, 14 ]
        self.assertEqual( endGroups, ch4_1._endGroups )
        
        blocks = [ None, ch4_2, None, None ]
        self.assertEqual( blocks, ch4_1._blocks )
        
        allBlocks = [ ch4_1, ch4_2, ch4_3 ]
        self.assertEqual( allBlocks, ch4_1.allBlocks )
        
        #ch4_1.writeXyz("bonded.xyz")
        
        
    def testMultiAdd(self):
        
        ch4_1 = self.makeCh4()
        ch4_2 = ch4_1.copy()
        ch4_3 = ch4_1.copy()
        
        idxAtom = 2
        idxGrowBlockAtom = 4
        ch4_1.positionGrowBlock(idxAtom, ch4_2, idxGrowBlockAtom )
        ch4_1.addBlock( idxAtom, ch4_2, idxGrowBlockAtom )
        
        idxAtom = 1
        idxGrowBlockAtom = 4
        ch4_3.positionGrowBlock(idxAtom, ch4_1, idxGrowBlockAtom )
        ch4_3.addBlock( idxAtom, ch4_1, idxGrowBlockAtom )
        
        #ch4_3.writeXyz("bonded2.xyz")
        endGroups = [ 2, 3, 4, 6, 8, 11, 12, 13 ]
        self.assertEqual( endGroups, ch4_3._endGroups )
        
        blocks = [ ch4_1, None, None, None ]
        self.assertEqual( blocks, ch4_3._blocks )
        
        allBlocks = [ ch4_3, ch4_1, ch4_2 ]
        self.assertEqual( allBlocks, ch4_3.allBlocks )
        
    
    def makePaf(self):
        """Return the PAF molecule for testing"""
        
        paf = Block()
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
        
    def testAlignBlocks(self):
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

    def testPositionGrowBlock(self):
        
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


    def XtestBond(self):
        """Test we can correctly bond two _blocks at the given bond"""
        
        ch4 = self.makeCh4()
        m2 = ch4.copy()
        m2.translate( numpy.array( [2, 2, 2] ) )
        
        bond = (3,2)
        ch4.bond(m2, bond)
        
        self.assertTrue( ch4._endGroups == [1,2,4,6,8,9], "Incorrect calculation of _endGroups")
        self.assertTrue( ch4._angleAtoms == {1: 0, 2: 0, 4: 0, 6: 5, 8: 5, 9: 5}, "Incorrect calculation of endGroup contacts")
    

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
        b = Block()
        b.createFromArgs( coords, labels, endGroups, angleAtoms )
        return b
        
    def XtestNew(self):
        """
        test new coordinate class
        """
        
        b1 = self.makeH2C2()
        b2 = b1.copy()
        
        idxB1EG = 0
        idxB2EG = 1
        
        
        #print b1
        b1.positionGrowBlock( idxB1EG, b2, idxB2EG )
        
        #print b1
        #print b2
        
        #b1.join( idxB1EG, b2, idxB2EG )
        
        #b1.writeXyz()
        
        
        
        
        
if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
        