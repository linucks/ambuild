'''
Created on Jan 15, 2013

@author: abbietrewin
'''

import copy
import os
import types

import buildingBlock
import numpy
import util

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
            '_atomRadii'       : [],
            '_atomTypes'       : [],
            '_body'            : [], # a list of which body within this fragment each atom belongs to
            '_bonds'           : [], # List of internal fragment bonds
            '_bonded'          : [], # List of which atoms are bonded to which
            '_charges'         : [],
            '_fragmentType'    : fragmentType, 
            '_labels'          : [],
            '_masses'          : [],
            '_radius'          : -1,
            '_symbols'         : [], # ordered array of _symbols (in upper case)
            '_totalMass'       : -1,
            '_individualAttrs' : None,
            '_sharedAttrs'     : None,
            }
        
        #
        # The variables below here are specific to a fragment and change as the framgent moves
        # and is involved in bonds - each fragment gets its own copy of these
        individualAttrs = {
            '_coords'        : [],
            '_centroid'      : None,
            '_centerOfMass'  : None,
            '_maxAtomRadius' : -1,
            '_changed'       : True, # Flag for when we've been moved and need to recalculate things
            '_blockIdx'      : None, # The index in the list of block data where the data for this fragment starts
            '_endGroups'     : [], # A list of the endGroup objects
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
    
    def bonded(self, idxAtom):
        return self._bonded[ idxAtom ]

    def _calcBonded(self):
        
        assert len(self._bonds)
        # Create empty lists for all
        self._bonded =  [ [] for _ in xrange( len(self._coords) ) ]
        for (b1,b2) in self._bonds:
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
    
    def copy(self):
        """Create a copy of ourselves.
        Those attributes in shared are just copied as references as they do not change between fragments
        of the same type
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

    def body(self, atomIdx ):
        """Return the body in this fragment that the atom is in"""
        return self._body[ atomIdx ]
    
    def isCapAtom(self, atomIdx):
        """Return True if this atom is a capAtom - doesn't check if bodned"""
        #return atomIdx in self._endGroups
        return atomIdx in [ eg.fragmentCapIdx for eg in self._endGroups ]
    
    def isEndGroup(self, atomIdx):
        """Return True if this atom is an endGroup - doesn't check if free"""
        #return atomIdx in self._endGroups
        return atomIdx in [ eg.fragmentEndGroupIdx for eg in self._endGroups ]

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
        endGroups, capAtoms, dihedralAtoms, uwAtoms = self.parseEndgroupFile( filePath )

        # Set the root fragment and its attributes
        self.setData( coords        = coords,
                     labels        = labels,
                     symbols       = symbols,
                     atomTypes     = atomTypes,
                     charges       = charges,
                     endGroups     = endGroups,
                     capAtoms      = capAtoms,
                     dihedralAtoms = dihedralAtoms,
                     uwAtoms       = uwAtoms
                  )
        
        self.processBodies( filePath )
        
        self._calcProperties()
        
        return

    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius
    
    def mass(self):
        """Return the total mass for this block"""
        
        assert self._totalMass > 0
        return self._totalMass

    def parseEndgroupFile(self, filePath ):
        
        dirname, filename = os.path.split( filePath )
        basename, suffix = os.path.splitext( filename )
        
        endGroups = []
        capAtoms = []
        uwAtoms = []
        dihedralAtoms = []
        egfile = os.path.join( dirname, basename+".ambi" )
        if os.path.isfile( egfile ):
            #raise RuntimeError,"Cannot find endgroup file {0} for car file {1}".format( egfile, filepath )
            egs = [ ( line.strip().split() ) for line in open( egfile, 'r') if line.strip() and not line.startswith("#") ]
            
            for eg in egs:
                endGroups.append( int( eg[0] ) )
                capAtoms.append( int( eg[1] ) )
                if len(eg) > 2:
                    dihedralAtoms.append( int( eg[2] )  )
                else:
                    dihedralAtoms.append( -1 )
                if len(eg) > 3:
                    uwAtoms.append( int( eg[3] )  )
                else:
                    uwAtoms.append( -1 )
                    
                
            if self._fragmentType == 'cap' and len( endGroups ) != 1:
                raise RuntimeError, "Capfile had >1 endGroup specified!"
            
        return endGroups, capAtoms, dihedralAtoms, uwAtoms

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
                coords         = None,
                labels         = None,
                symbols        = None,
                atomTypes      = None,
                charges        = None,
                endGroups      = None,
                capAtoms       = None,
                dihedralAtoms  = None,
                uwAtoms        = None
                ):
        
        self._coords = coords
        self._labels = labels
        self._symbols = symbols
        self._atomTypes = atomTypes
        self._charges = charges
        
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
        self.setEndGroups( endGroups, capAtoms, dihedralAtoms, uwAtoms )
        
        return

    def setEndGroups( self, endGroups, capAtoms, dihedralAtoms, uwAtoms ):
        
        # Now set up the endGroup information
        self._endGroups =[]
        for i, e in enumerate( endGroups ):
            
            eg                     = buildingBlock.EndGroup()
            eg.fragment            = self
            eg.fragmentEndGroupIdx = e
            eg.fragmentCapIdx      = capAtoms[ i ]
            eg.fragmentDihedralIdx = dihedralAtoms[ i ]
            eg.fragmentUwIdx       = uwAtoms[ i ]
            
            # sanity check
            assert eg.fragmentCapIdx in self.bonded( eg.fragmentEndGroupIdx ),\
            "capAtom {0} is not bonded to endGroup {1}".format( eg.fragmentCapIdx,eg.fragmentEndGroupIdx )
            assert eg.fragmentDihedralIdx in self.bonded( eg.fragmentEndGroupIdx ),\
            "dihedral Atom {0} is not bonded to endGroup {1}".format( eg.fragmentDihedralIdx, eg.fragmentEndGroupIdx )
            assert eg.fragmentUwIdx in self.bonded( eg.fragmentEndGroupIdx ),\
            "uwAtom {0} is not bonded to endGroup {1}".format( eg.fragmentUwIdx, eg.fragmentEndGroupIdx )
            
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

