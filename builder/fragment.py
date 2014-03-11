'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import os
import types

import buildingBlock
import numpy
import util

class Fragment(object):
    '''
    classdocs
    '''

    def __init__( self, fragmentType ):
        '''
        Constructor
        '''

        # the coordinates for this block
        self._coords = []
        
        # ordered array of _labels
        self._labels = []
        
        # ordered array of atomTypes
        self._atomTypes = []
        
        # ordered array of _symbols (in upper case)
        self._symbols = []
        
        # ordered array of _masses
        self._masses = []
        
        # orderd array of atom radii
        self._atomRadii = []
        
        self._charges = []
        
        # These are calculated
        self._centroid = None
        self._centerOfMass = None
        self._totalMass = -1
        self._radius = -1
        self._maxAtomRadius = -1
        
        # Flag for when we've been moved and need to recalculate things
        self._changed = True
        
        # The type this fragment is (for bonding checks)
        self._fragmentType = fragmentType
        
        # The index in the list of block data where the data for this fragment starts
        self.blockIdx = None
        
        # A list of the endGroup objects
        self._endGroups = []
        #self._angleAtoms = []
        
        # List of internal fragment bonds
        self._bonds = []
        
        # a list of which body within this fragment each atom belongs to
        self._body = []
        
        # Additional rigid bodies attached to the fragment
        #self.processBodies( filepath )
        
        return

    def _calcCenters(self):
        """Calculate the center of mass and geometry for this fragment
        """
        
        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        self._totalMass = 0.0 # Not sure if it sensible to calculate each time - prob irrelevant
        for i, coord in enumerate( self._coords ):
            mass = self._masses[i]
            self._totalMass += mass
            sumG += coord
            sumM += mass * coord
        
        self._centroid = sumG / len( self._coords )
        self._centerOfMass = sumM / self._totalMass
        
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
        
        distances = []
        for i, coord in enumerate( self._coords ):
            #distances.append( numpy.linalg.norm(coord-cog) )
            distances.append( util.distance(self._centroid, coord) )
            
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

    def centroid(self):
        """
        Return or calculate the center of geometry for the building block.
        """
        if self._changed:
            self._calcProperties()
        return self._centroid

    def centerOfMass(self):
        """
        Return or calculate the center of mass for this building block.
        """
        if self._changed:
            self._calcProperties()
        return self._centerOfMass

    def body(self, atomIdx ):
        """Return the body in this fragment that the atom is in"""
        return self._body[ atomIdx ]
    
    def isEndGroup(self, atomIdx):
        """Return True if this atom is an endGroup - doesn't check if free"""
        #return atomIdx in self._endGroups
        return atomIdx in [ eg.fragmentEndGroupIdx for eg in self._endGroups ]
    
    def endGroupBonded( self, endGroup ):
        """Return a list of all atoms bonded to the given endGrop
        It excludes the capAtom.
        It does not however exclude atoms that could be capAtoms for other endGroups 
        """
        
        idxEndGroup = endGroup.fragmentEndGroupIdx
        idxCapAtom = endGroup.fragmentCapIdx
        
        # Could precalculate this if speed an issue
        bonded = []
        for b in self._bonds:
            if b[0] == idxEndGroup and b[1] != idxCapAtom:
                bonded.append( b[1] )
            elif b[1] == idxEndGroup and b[0] != idxCapAtom:
                bonded.append( b[0] )
        return bonded

    def fillData(self):
        """ Fill the data arrays from the label """
        
        self._totalMass = 0.0
        self._maxAtomRadius = 0.0
        for i in range( len(self._coords) ):
            
            symbol = self._symbols[ i ]
                
            # Masses
            mass = util.ATOMIC_MASS[ symbol ]
            self._masses.append( mass )
            self._totalMass += mass
            
            # Radii
            z = util.SYMBOL_TO_NUMBER[ symbol.upper() ]
            r = util.COVALENT_RADII[z] * util.BOHR2ANGSTROM
            self._atomRadii.append(r)
            
            # Remember the largest
            self._maxAtomRadius = max( r, self._maxAtomRadius )
            #print "ADDING R {} for label {}".format(r,label)
        
        return

    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius
    
    def mass(self):
        """Return the total mass for this block"""
        
        assert self._totalMass > 0
        return self._totalMass

    def processBodies(self, filepath ):
        """See if we split the fragment into bodies or not"""
        
        assert len(self._coords) > 0, "Coordinates must have been read before processing bodies"
        
        dirname, filename = os.path.split( filepath )
        basename, suffix = os.path.splitext( filename )
        
        bodyFile = os.path.join( dirname, basename+".ambody" )
        if os.path.isfile( bodyFile ):
            self._body = [ int( l.strip() ) for l in open( bodyFile ) ]
            assert len( self._body ) == len(self._coords), \
            "Must have as many bodies as coordinates: {0} - {1}!".format( len( self._body ), self._dataLen )
        else:
            # Just create an array with 0
            self._body = [ 0 ] * len( self._coords )
        return

    def radius(self):
        
        if self._changed:
            self._calcProperties()
        return self._radius

    def setData(self,
                coords    = None,
                labels    = None,
                symbols   = None,
                atomTypes = None,
                charges   = None,
                endGroups = None,
                capAtoms  = None,
                uwAtoms   = None
                ):
        
        self._coords = coords
        self._labels = labels
        self._symbols = symbols
        self._atomTypes = atomTypes
        self._charges = charges
        
        # Calculate anything we haven't been given
        self.fillData()
        
        # Specify internal bonds
        self._bonds = util.calcBonds( coords, symbols, maxAtomRadius=self.maxAtomRadius() )
        
        # Set up endGroups
        self.setEndGroups( endGroups, capAtoms, uwAtoms )
        
        # Recalculate all dynamic data
        #self._calcCenters()
        #self._calcRadius()
        #self._changed = False
        
        return

    def setEndGroups( self, endGroups, capAtoms, uwAtoms ):
        
        # Now set up the endGroup information
        self._endGroups =[]
        for i, e in enumerate( endGroups ):
            
            eg = buildingBlock.EndGroup()
            eg.fragment = self
            eg.fragmentEndGroupIdx = e
            eg.fragmentCapIdx = capAtoms[ i ]
            eg.fragmentUwIdx = uwAtoms[ i ]
            eg.fragmentBonded = self.endGroupBonded( eg )
            
            self._endGroups.append( eg )
        return

    #def rotate( self, axis, angle, center=None ):
    def rotate( self, rotationMatrix, center ):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
#         if center==None:
#             center = numpy.array([0,0,0])
#         rmat = util.rotation_matrix( axis, angle )
        
        # loop through all _blocks and change all _coords
        # Need to check that the center bit is the best way of doing this-
        # am almost certainly doing more ops than needed
        for i,coord in enumerate( self._coords ):
            coord = coord - center
            coord = numpy.dot( rotationMatrix, coord )
            self._coords[i] = coord + center
        
        self._changed = True
        
        return

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # FIX!!
        #assert isinstance( tvector, numpy.array )

        # Use len as we don't need to return the _coords
        for i in range( len (self._coords ) ):
            self._coords[i] += tvector
        
        self._changed = True
        
        return

    def type(self):
        return self._fragmentType

    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))

