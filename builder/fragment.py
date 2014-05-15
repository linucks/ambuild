'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import copy
import os
import types

import numpy
import util

class EndGroup(object):
    
    def __init__(self):
        
        self._isBonded            = False
        
        self._endGroupType        = None
        
        self.fragment             = None
        
        self.fragmentEndGroupIdx = None
        self.blockEndGroupIdx    = None
        
        self.fragmentCapIdx      = None
        self.blockCapIdx         = None
        
        self.fragmentDihedralIdx = -1
        self.blockDihedralIdx    = -1
        
        self.fragmentUwIdx       = -1
        self.blockUwIdx          = -1
        
        return 
    
    def block(self):
        return self.fragment.block

    def capIdx(self):
        """Return the index of the endGroup atom in external block indices"""
        return self.fragment.block._int2ext[self.blockCapIdx]
    
    def endGroupIdx(self):
        """Return the index of the endGroup atom in external block indices"""
        return self.fragment.block._int2ext[self.blockEndGroupIdx]

    def fragmentType(self):
        return self.fragment.fragmentType
    
    def isBonded(self):
        return self._isBonded
    
    def free(self):
        return not self.isBonded() and not self.saturated()
    
    def saturated(self):
        return self.fragment.endGroupSaturated( self.type() )
    
    def setBonded(self):
        self._isBonded = True
        self.fragment.addBond( self.type() )
        # Mask fragment cap and uw atoms now
        self.fragment._masked[ self.fragmentCapIdx ] = True
        if self.fragmentUwIdx != -1:
            self.fragment._masked[ self.fragmentUwIdx ] = True
        return
    
    def type(self):
        """The type of the endGroup"""
        return "{0}:{1}".format( self.fragment.fragmentType, self._endGroupType )
    
    def updateAncillaryIndices(self,cap2endGroup):
        """The block has been updated so we need to update our block indices based on where the
        fragment starts in the block"""
        if self.fragment.isMasked( self.fragmentCapIdx ):
            assert self.isBonded()
            # If this endGroup is involved in a bond, we want to get the block index of the
            # opposite endGroup as this has now become our cap atom
            self.blockCapIdx = cap2endGroup[ (self.fragment,self.fragmentCapIdx) ]
        else:
            self.blockCapIdx = self.fragmentCapIdx + self.fragment._blockIdx 
        
        # -1 means no dihedral or uw atom set
        if self.fragmentDihedralIdx != -1:
            if self.fragment.isMasked( self.fragmentDihedralIdx ):
                # Now work out which atom this is bonded to
                self.blockDihedralIdx = cap2endGroup[ (self.fragment,self.fragmentDihedralIdx) ]
            else:
                self.blockDihedralIdx = self.fragmentDihedralIdx + self.fragment._blockIdx
                
        if self.fragmentUwIdx != -1:
            #assert not self.fragment.isMasked( self.fragmentUwIdx )
            self.blockUwIdx = self.fragmentUwIdx + self.fragment._blockIdx 
        return
    
    def updateEndGroupIndex(self):
        """The block has been updated so we need to update our block indices based on where the
        fragment starts in the block"""
        self.blockEndGroupIdx = self.fragmentEndGroupIdx + self.fragment._blockIdx
        return
    
    def __str__(self):
        """List the data attributes of this object"""
#         me = {}
#         for slot in dir(self):
#             attr = getattr(self, slot)
#             if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
#               isinstance(attr, types.FunctionType) ):
#                 me[slot] = attr
#                 
#         t = []
#         for k in sorted(me):
#             t.append( str( ( k, me[k] ) ) )
#                                    
#         return "{0} : {1}".format(self.__repr__(), ",".join( t ) )

        s = "{0} {1}:{2}:{3} {4}->{5} {6}->{7}".format( self.type(), self.block().id,id(self.fragment),id(self),
                                                    self.fragmentEndGroupIdx, self.fragmentCapIdx,
                                                    self.blockEndGroupIdx,self.blockCapIdx)

        return s

class Fragment(object):
    '''
    classdocs
    '''

    def __init__( self, filePath=None, fragmentType=None ):
        '''
        Constructor
        '''

        # This first set of variables are shared by all fragments
        # When we copy a fragment the new fragment gets references to the variables
        # that were created for the first fragment (see copy)
        sharedAttrs = { 
            'radii'            : [],
            'types'            : [],
            'bodies'           : [], # a list of which body within this fragment each atom belongs to
            '_bonds'           : [], # List of internal fragment bonds
            '_bonded'          : [], # List of which atoms are bonded to which
            'charges'          : [],
            'fragmentType'     : fragmentType, 
            'labels'           : [],
            'masses'           : [],
            '_maxBonds'        : {},
            '_radius'          : -1,
            'symbols'          : [], # ordered array of symbols (in upper case)
            '_totalMass'       : -1,
            '_individualAttrs' : None,
            '_sharedAttrs'     : None,
            }
        
        #
        # The variables below here are specific to a fragment and change as the fragment moves
        # and is involved in bonds - each fragment gets its own copy of these
        individualAttrs = {
            'block'           : None,
            'coords'         : [],
            '_centroid'       : None,
            '_centerOfMass'   : None,
            '_maxAtomRadius'  : -1,
            '_changed'        : True, # Flag for when we've been moved and need to recalculate things
            '_blockIdx'       : None, # The index in the list of block data where the data for this fragment starts
            '_endGroups'      : [], # A list of the endGroup objects
            '_endGroupBonded' : [], # A list of the number of each endGroup that are used in bonds
            '_masked'         : [], # bool - whether the atoms is hidden (e.g. cap or uw atom)
            }
        
        # Set as attributes of self
        for a, v in sharedAttrs.iteritems():
            setattr( self, a, v )
            
        # Set as attributes of self
        for a, v in individualAttrs.iteritems():
            setattr( self, a, v )
            
        # Set these manually
        self._individualAttrs = individualAttrs
        self._sharedAttrs = sharedAttrs
        
        # Create from the file
        if filePath:
            self.fromFile(filePath)
        
        return
    
    def addBond(self, endGroupType):
        self._endGroupBonded[ endGroupType ] += 1
        return
    
    def bonds(self):
        """Return a list of bonds in block indices
        We exclude masked atoms"""
        bonds = []
        for b1, b2 in self._bonds:
            if not ( self._masked[b1] or self._masked[b2] ):
                # Map to internal and then add _blockIdx to get position in block
                b1 = b1 + self._blockIdx
                b2 = b2 + self._blockIdx
                bonds.append( (b1, b2) )
        return bonds

    def _calcBonded(self):
        assert len(self._bonds)
        # Create empty lists for all
        self._bonded =  [ [] for _ in xrange( len(self.coords) ) ]
        for b1,b2 in self._bonds:
            if b1 not in self._bonded[ b2 ]:
                self._bonded[ b2 ].append( b1 )
            if b2 not in self._bonded[ b1 ]:
                self._bonded[ b1 ].append( b2 )
        return

    def _calcCenters(self):
        """Calculate the center of mass and geometry for this fragment
        """
        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        self._totalMass = 0.0 # Not sure if it sensible to calculate each time - prob irrelevant
        for i, coord in enumerate( self.coords ):
            mass = self.masses[i]
            self._totalMass += mass
            sumG += coord
            sumM += mass * coord
        
        self._centroid = sumG / len( self.coords )
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
        for i, coord in enumerate( self.coords ):
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

#     def charge(self,atomIdx):
#         return self._charges[ atomIdx ]
# 
#     def coord(self,atomIdx):
#         return self._coords[ atomIdx ]
    
    def copy(self):
        """Create a copy of ourselves.
        Those attributes in shared are just copied as references as they do not change between fragments
        of the same fragmentType
        Those in single are deep-copied as each fragment has its own"""
        
        f = Fragment()
        for a in f.__dict__:
            if a in self._sharedAttrs.keys():
                setattr( f, a, getattr( self, a ) )
            elif a in self._individualAttrs.keys():
                setattr( f, a, copy.deepcopy( getattr( self, a ) ) )
            else:
                raise RuntimeError,"MISSING ATTRIBUTE {0}".format(a)
            
            # Update fragment references in the endGroups
            for e in f._endGroups:
                e.fragment = f
                
        return f

    def endGroups(self):
        return [ eg for eg in self._endGroups ]

    def endGroupSaturated(self, endGroupType):
        if self._maxBonds[ endGroupType ] is not None and \
        self._endGroupBonded[ endGroupType ] >= self._maxBonds[ endGroupType ]:
            return True
        return False
    
    def endGroupTypes(self):
        """Return a list of the endGroupTypes in this fragment"""
        return set(self._endGroupBonded.keys())

    def fillData(self):
        """ Fill the data arrays from the label """
        
        self._masked = [ False ] * len(self.coords)
        self._totalMass = 0.0
        self._maxAtomRadius = 0.0
        for i in range( len(self.coords) ):
            
            symbol = self.symbols[ i ]
                
            # Masses
            mass = util.ATOMIC_MASS[ symbol ]
            self.masses.append( mass )
            self._totalMass += mass
            
            # Radii
            z = util.SYMBOL_TO_NUMBER[ symbol.upper() ]
            r = util.COVALENT_RADII[z] * util.BOHR2ANGSTROM
            self.radii.append(r)
            
            # Remember the largest
            self._maxAtomRadius = max( r, self._maxAtomRadius )
            
        return
    
    def freeEndGroups(self):
        return [ eg for eg in self._endGroups if eg.free() ]

    def fromCarFile(self, carFile):
        """"Abbie did this.
        Gulp...
        """
        labels = []
        symbols = []
        atomTypes = []
        charges = []
        
        # numpy array
        coords = []
        
        reading = True
        with open( carFile, "r" ) as f:
            
            # skip first line
            f.readline()
            
            # 2nd states whether PBC: PBC=OFF
            pbc, state = f.readline().strip().split("=")
            assert pbc.strip() == "PBC"
            state=state.strip()
            nskip=3
            if state.upper() == "OFF":
                nskip=2 
            
            for i in range(nskip):
                f.readline()
            
            count=0
            while reading:
                line = f.readline()
                
                line = line.strip()
                if not line:
                    print "END OF CAR WITH NO END!!!"
                    break
                fields = line.split()
                label = fields[0]
                
                # Check end of coordinates
                if label.lower() == "end":
                    reading=False
                    break
                
                labels.append( label )
                coords.append( numpy.array(fields[1:4], dtype=numpy.float64 ) )
                atomTypes.append( fields[6] )
                symbols.append( fields[7] )
                charges.append( float( fields[8] ) )
                
                count+=1
        
        return  ( coords, labels, symbols, atomTypes, charges )

    def fromXyzFile(self, xyzFile ):
        """"Jens did this.
        """
        
        labels = []
        symbols = []
        atomTypes = [] # hack...
        charges = []
        
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
                symbol = util.label2symbol( label )
                symbols.append( symbol )
                atomTypes.append(  symbols )
                charges.append( 0.0 )
                
                count += 1
        
        return ( coords, labels, symbols, atomTypes, charges )
    
    def fromFile(self, filePath ):

        if filePath.endswith( ".car" ):
            ( coords, labels, symbols, atomTypes, charges ) = self.fromCarFile( filePath )
        elif filePath.endswith( ".xyz" ):
            ( coords, labels, symbols, atomTypes, charges ) = self.fromXyzFile( filePath )
        else:
            raise RuntimeError("Unrecognised file suffix: {}".format( filePath ) )
        
        # Get cap atoms and endgroups
        endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms = self.parseEndgroupFile( filePath )

        # Set the root fragment and its attributes
        self.setData( coords       = coords,
                     labels        = labels,
                     symbols       = symbols,
                     types         = atomTypes,
                     charges       = charges,
                     endGroupTypes = endGroupTypes,
                     endGroups     = endGroups,
                     capAtoms      = capAtoms,
                     dihedralAtoms = dihedralAtoms,
                     uwAtoms       = uwAtoms
                  )
        
        self.processBodies( filePath )
        
        self._calcProperties()
        
        return
    
    def isMasked(self,atomIdx):
        return self._masked[atomIdx]
    
#     def label(self,atomIdx):
#         return self._labels[ atomIdx ]
# 
#     def mass(self,atomIdx):
#         return self._masses[ atomIdx ]

    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius
    
    def totalMass(self):
        """Return the total mass for this block"""
        
        assert self._totalMass > 0
        return self._totalMass
    
    def numAtoms(self):
        return len(self.coords)

    def parseEndgroupFile(self, filePath ):
        
        dirname, filename = os.path.split( filePath )
        basename, suffix = os.path.splitext( filename )
        
        endGroupTypes = []
        endGroups     = []
        capAtoms      = []
        uwAtoms       = []
        dihedralAtoms = []
        egfile = os.path.join( dirname, basename+".ambi" )
        if os.path.isfile( egfile ):
            #raise RuntimeError,"Cannot find endgroup file {0} for car file {1}".format( egfile, filepath )
            egs = [ ( line.strip().split() ) for line in open( egfile, 'r') if line.strip() and not line.startswith("#") ]
            
            for eg in egs:
                # For now make sure first value is letter
                assert eg[0][0].isalpha(),"First column of ambi file needs to be a letter!"
                endGroupTypes.append( eg[0] )
                endGroups.append( int( eg[1] ) )
                capAtoms.append( int( eg[2] ) )
                if len(eg) > 3:
                    dihedralAtoms.append( int( eg[3] )  )
                else:
                    dihedralAtoms.append( -1 )
                if len(eg) > 4:
                    uwAtoms.append( int( eg[4] )  )
                else:
                    uwAtoms.append( -1 )
                    
                
            if self.fragmentType == 'cap' and len( endGroups ) != 1:
                raise RuntimeError, "Capfile had >1 endGroup specified!"
            
        return endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms

    def processBodies(self, filepath ):
        """See if we split the fragment into bodies or not"""
        
        assert len(self.coords) > 0, "Coordinates must have been read before processing bodies"
        
        dirname, filename = os.path.split( filepath )
        basename, suffix = os.path.splitext( filename )
        
        bodyFile = os.path.join( dirname, basename+".ambody" )
        if os.path.isfile( bodyFile ):
            self.bodies = [ int( l.strip() ) for l in open( bodyFile ) ]
            assert len( self.bodies ) == len(self.coords), \
            "Must have as many bodies as coordinates: {0} - {1}!".format( len( self.bodies ), self._dataLen )
        else:
            # Just create an array with 0
            self.bodies = [ 0 ] * len( self.coords )
        return

#     def radius(self,atomIdx):
#         return self._radii[ atomIdx ]

    def totalRadius(self):
        
        if self._changed:
            self._calcProperties()
        return self._radius

    def rotate( self, rotationMatrix, center ):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
#         if center==None:
#             center = numpy.array([0,0,0])
#         rmat = util.rotation_matrix( axis, angle )
        
        # loop through all _blocks and change all coords
        # Need to check that the center bit is the best way of doing this-
        # am almost certainly doing more ops than needed
        for i,coord in enumerate( self.coords ):
            coord = coord - center
            coord = numpy.dot( rotationMatrix, coord )
            self.coords[i] = coord + center
        
        self._changed = True
        
        return
    
    def setData(self,
                coords         = None,
                labels         = None,
                symbols        = None,
                types          = None,
                charges        = None,
                endGroupTypes  = None,
                endGroups      = None,
                capAtoms       = None,
                dihedralAtoms  = None,
                uwAtoms        = None
                ):
        
        self.coords  = coords
        self.labels  = labels
        self.symbols = symbols
        self.types   = types
        self.charges = charges
        
        # Calculate anything we haven't been given
        self.fillData()
        
        # Specify internal bonds - bond margin probably too big...
        self._bonds = util.calcBonds( coords,
                                      symbols,
                                      maxAtomRadius=self.maxAtomRadius(),
                                      bondMargin=0.25
                                       )
        
        # Create list of which atoms are bonded to each atom
        self._calcBonded()
        
        # Set up endGroups
        self.setEndGroups( endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms )
        
        return

    def setEndGroups( self, endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms ):
        
        # Now set up the endGroup information
        self._endGroups      = []
        self._maxBonds       = {}
        self._endGroupBonded = {}
        
        for i, e in enumerate( endGroups ):
            eg                     = EndGroup()
            eg.fragment            = self
            eg._endGroupType        = endGroupTypes[ i ]
            eg.fragmentEndGroupIdx = e
            eg.fragmentCapIdx      = capAtoms[ i ]
            eg.fragmentDihedralIdx = dihedralAtoms[ i ]
            eg.fragmentUwIdx       = uwAtoms[ i ]
            
            if eg.type() not in self._maxBonds:
                self._maxBonds[ eg.type() ]       = None
            if eg.type() not in self._endGroupBonded:
                self._endGroupBonded[ eg.type() ] = 0
            
            # sanity check
            assert eg.fragmentCapIdx in self._bonded[ eg.fragmentEndGroupIdx ],\
            "capAtom {0} is not bonded to endGroup {1}".format( eg.fragmentCapIdx,eg.fragmentEndGroupIdx )
            if eg.fragmentDihedralIdx != -1:
                assert eg.fragmentDihedralIdx in self._bonded[ eg.fragmentEndGroupIdx ],\
            "dihedral Atom {0} is not bonded to endGroup {1}".format( eg.fragmentDihedralIdx, eg.fragmentEndGroupIdx )
            if eg.fragmentUwIdx != -1:
                assert eg.fragmentUwIdx in self._bonded[ eg.fragmentEndGroupIdx ],\
            "uwAtom {0} is not bonded to endGroup {1}".format( eg.fragmentUwIdx, eg.fragmentEndGroupIdx )
            
            self._endGroups.append( eg )
        
        # sanity check - make sure no endGroup is the cap Atom for another endGroup
        eg = set( [ e.fragmentEndGroupIdx for e in self._endGroups ] )
        caps = set( [ e.fragmentCapIdx for e in self._endGroups ] )
        assert len(eg.intersection( caps )) == 0,"Cap atom for one endGroup is another endGroup!"
        
        return

    def setMaxBond(self, endGroupType, maxBond):
        assert endGroupType in self._maxBonds
        self._maxBonds[ endGroupType ] = maxBond
        return

#     def symbol(self,atomIdx):
#         return self._symbols[ atomIdx ]

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # FIX!!
        #assert isinstance( tvector, numpy.array )

        # Use len as we don't need to return the coords
        for i in range( len (self.coords ) ):
            self.coords[i] += tvector
        
        self._changed = True
        return

#     def type(self,atomIdx):
#         return self._types[ atomIdx ]

    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))

