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
import os
import copy
import math
import random
import unittest

# external imports
import numpy

# local imports
import util

class BuildingBlock():
    '''
    classdocs
    '''

    def __init__( self, infile=None ):
        '''
        Constructor
        Distance is the cells distance function used to calculate the distance between two vectors
        - passed in through constructor as we may be under periodic bounduary conditions.
        '''
        
        self.coords = []
        
        # ordered array of labels
        self.labels = []
        
        # ordered array of symbols (in upper case)
        self.symbols = []
        
        # ordered array of masses
        self.masses = []
        
        # orderd array of atom radii
        self.atom_radii = []
        
        # List of the cell (3-tuple) to which each atom belongs
        self.atomCell = []
        
        # Dict of dicts mapping atom types to their bond lengths for this block
        # used so we only check the types of atoms for this type of block
        # For now we assume that both blocks have the same atom types
        self._labels2bondlength = {}
        
        # A list of the indices of the atoms that are endGroups
        self.endGroups = []
        
        # Dictionary mapping the index of an endGroup to an atom bonded to it
        # Required for finding the angle when we check bonds
        self._endGroupContacts = {}
        
        # Holds the center of mass of the molecule
        self._centerOfMass = numpy.zeros( 3 )
        # Holds the center of geometry of the molecule
        self._centerOfGeometry = numpy.zeros( 3 )
        
        # The radius of the block assuming it is a circle centered on the COG
        self._radius = 0
        
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
                
    
    def bond (self, block, bond):
        """Bond the two blocks at the given bond - tuple is indices of self and other bond
        """
        
        print "Bonding Block"
        #print block
        #print "with"
        #print self
        
        # Needed to work out how much to add to the block indices
        lcoords = len(self.coords)
        
        self.coords.extend( block.coords )
        self.atom_radii.extend( block.atom_radii )
        self.labels.extend( block.labels )
        
        #jmht hack
        for i,l in enumerate(self.labels):
            if l[0] == 'H':
                self.labels[i] = 'X'+l[1:]
                print "CHANGING LABEL TO ",self.labels[i] 
        self.symbols.extend( block.symbols )
        self.masses.extend( block.masses )
        self.atomCell.extend( block.atomCell )
        
        # Need to remove the end groups used in the bond
        i = self.endGroups.index( bond[0] )
        self.endGroups.pop( i )
        i = block.endGroups.index( bond[1] )
        block.endGroups.pop( i )
        
        # Now add block to self, updating index
        for i in block.endGroups:
            self.endGroups.append( i+lcoords )
        
        del self._endGroupContacts[ bond[0] ]
        del block._endGroupContacts[ bond[1] ]
        
        for k,v in block._endGroupContacts.iteritems():
            self._endGroupContacts[ k+lcoords ] = v+lcoords
        
        self.update()

        
    def createFromArgs(self, coords, labels, endGroups, endGroupContacts ):
        """ Create from given arguments
        """
        # could check if numpy array here
        self.coords = coords
        self.labels = labels
        self.endGroups = endGroups
        self._endGroupContacts = endGroupContacts
        
        # Fill atomCell list
        self.atomCell = [None]*len(self.coords)
        
        self.fillData()
        
        self.update()
        
    def createFromLabelAndCoords(self, labels, coords):
        """ Given an array of labels and another of coords, create a block
        This requires determining the end groups and contacts from the label
        """
        
        # array of indexes
        endGroups = []
        
        endGroupContacts = {}
        
        # For tracking the mapping of the index of the endgroup to the label
        # used to mark the contat atom
        eGindex2label = {}
        
        # For tracking the mapping of the index used to mark the contact atom
        # to the true index of the atom
        eAlabel2index = {}
        
        for i, label in enumerate(labels):
            
            # End groups and the atoms which define their bond angles 
            # are of the form XX_EN for endgroups and XX_AN for the defining atoms
            # where N can be any number, but linking the two atoms together 
            # The indexing of the atoms starts from 1!!!!
            if "_" in label:
                _,ident = label.split("_")
                atype=ident[0]
                anum=int(ident[1:])
                
                if atype.upper() == "E":
                    # An endgroup
                    endGroups.append( i )
                    eGindex2label[ i ] = anum
                elif atype.upper() == "A":
                    eAlabel2index[ anum ] = i
                else:
                    raise RuntimeError,"Got a duff label! - {}".format(label)
        
        #Now match up the endgroups with their contact atoms
        for trueIndex, mapIndex in eGindex2label.iteritems():
            if trueIndex not in endGroups:
                raise RuntimeError,"Got a bad index for an endgroup!"
            endGroupContacts[ trueIndex ] = eAlabel2index[ mapIndex ]
        
        self.createFromArgs(coords, labels, endGroups, endGroupContacts)
  
    def calcCenterOfMassAndGeometry(self):
        """Calculate the center of mass and geometry
        """
        
        sumG = numpy.zeros( 3 )
        sumM = numpy.zeros( 3 )
        
        totalMass = 0.0
        for i, coord in enumerate( self.coords ):
            mass = self.masses[i]
            totalMass += mass
            sumG += coord
            sumM += mass * coord
        
        self._centerOfGeometry = sumG / (i+1)
        self._centerOfMass = sumM / totalMass
        
        
    def calcRadius(self):
        """
        Calculate a simple size metric of the block so that we can screen for whether
        two blocks are within touching distance
        
        First try a simple approach with a loop just to get a feel for things
        - Find the largest distance between any atom and the center of geometry
        - Get the covalent radius of that atom
        - return that distance + radius + buffer
        
        Should move to use scipy as detailed here:
        http://stackoverflow.com/questions/6430091/efficient-distance-calculation-between-n-points-and-a-reference-in-numpy-scipy
        """
        
        cog = self.centroid()
        
        distances = []
        for coord in self.coords:
            #distances.append( numpy.linalg.norm(coord-cog) )
            distances.append( util.distance(cog, coord) )
            
        imax = numpy.argmax( distances )
        dist = distances[ imax ]
        atomR = self.atom_radii[ imax ]
        
        # Set radius
        self._radius = dist + atomR
        
    def XcanBond( self, block, bondMargin=None, bondAngle=None ):
        """
        See if we can form a bond with the given block.
        Return the indicies of the two atoms (self and then other)
        or False if the molecules cannot bond but do not clash
        If the blocks clash (are close but cannot bond) we return
        the string "clash" to indicate a clash and that the move should be
        rejected
        """
        
        assert bondAngle
        
        # Might be able to bond so loop through all atoms and check
        # We only accept the case where just one pair of atoms is close enough to bond
        # more than one is considered a failure
        # NB ASSUMPION FOR BOND LENGTH CHECK IS BOTH BLOCKS HAVE SAME ATOM TYPES
        bond = None
        for i, c in enumerate( self.coords ):
            #i_radius = self.atom_radii[i]
            i_symbol = self.symbols[i]
            for j, b in enumerate( block.coords ):
                #j_radius = block.atom_radii[j]
                j_symbol = block.symbols[j]
                bond_length = self.bondLength( i_symbol, j_symbol )
                #bond_length = util.bondLength( i_symbol, j_symbol )
                #if ( numpy.linalg.norm( c - b ) < bond_length + bondMargin ):
                if ( self.distance( b,c ) < bond_length + bondMargin ):
                    # Check if both atoms are end groups
                    if not ( i in self.endGroups and j in block.endGroups ):
                        # 2 atoms close enough to bond but they are not end groups
                        return "clash"
                        
                    # Close enough to bond
                    if bond:
                        # Already got one so this is no good
                        return "clash"
                    bond = (i,j)
        
        if not bond:
            return False
        
        # Got a bond so check the angle
        icontact = self.endGroupContactIndex( bond[0] )
        contact = self.coords[ icontact ]
        sgatom = self.coords[ bond[0] ] 
        bgatom = block.coords[ bond[1] ] 
        angle = self.angle( contact, sgatom, bgatom )
        
#       angle = self.dihedral( ibcoord, igcoord, jgcoord, jbcoord )
        
        print "Got bond angle D: {}".format(angle)
        
        if ( bondAngle-15 < angle < bondAngle+15 ):
            return bond
        else:
            print "Cannot bond due to angle"
            return "clash"
        
    def centroid(self):
        """
        Return or calculate the center of geometry for the building block.
        """
        
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self._changed = False
        
        return self._centerOfGeometry
        
    def centerOfMass(self):
        """
        Return or calculate the center of mass for the building block.
        """
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self._changed = False
        
        return self._centerOfMass
    
    def Xclash( self, block, blockMargin=2.0, atomMargin=1.0 ):
        """ See if this molecule clashes (overlaps) with the given one
        by checking the atomic radii 
        """
        
        # First check if these two molecules are within range assuming they 
        # are circular, centered on the COG and with the given radius
        if not self.close(block, blockMargin=blockMargin):
            return False
        
        # Within range, so see if any atoms actually do overlap
        # Do this with scipy and no loops when we optimise
        for i, c in enumerate( self.coords ):
            i_radius = self.atom_radii[i]
            for j, b in enumerate( block.coords ):
                j_radius = block.atom_radii[j]
                #if ( numpy.linalg.norm( c - b ) < i_radius + j_radius + atomMargin ):
                if ( self.distance( b, c ) < i_radius + j_radius + atomMargin ):
                    return True
                
        return False
    
    def Xclose( self, block, blockMargin=2.0 ):
        """Return true of false depending on whether two blocks are close enough to bond/clash.
        Works from the overall radii of the two blocks
        Margin is allowed gap between their respective radii
        """
        #dist = numpy.linalg.norm( self.centroid() - block.centroid() )
        dist = self.distance( block.centroid(), self.centroid() )
        if ( dist < self.radius() + block.radius() + blockMargin ):
            return True
        else:
            return False
    
    def copy( self ):
        """
        Create a copy of ourselves and return it.
        OLD - Tried using deepcopy, but this took forever - I think this is because the distance functions'
        reference to the cell is causing the cell to be copied as well.
        """
        return copy.deepcopy(self)

    def fillData(self):
        """ Fill the data arrays from the label """
        
        if len(self.labels) < len(self.coords):
            raise RuntimeError("fillData needs labels filled!")
        
        symbol_types=[]
        for label in self.labels:
            
            # Symbols
            symbol = util.label2symbol( label )
            self.symbols.append( symbol )
            
            # For checking bonding
            if symbol not in symbol_types:
                symbol_types.append(symbol)
                
            # Masses
            self.masses.append( util.ATOMIC_MASS[ symbol ] )
            
            # Radii
            z = util.SYMBOL_TO_NUMBER[ symbol ]
            r = util.COVALENT_RADII[z] * util.BOHR2ANGSTROM
            self.atom_radii.append(r)
            #print "ADDING R {} for label {}".format(r,label)
            
            # Add bond lengths
            
        
        # make a dict mapping all atom types to their bonds to speed lookups
        # not sure if this actually saves anything...
        for i_symbol in symbol_types:
            self._labels2bondlength [ i_symbol ] = {}
            for j_symbol in symbol_types:
                bond_length = util.bondLength( i_symbol , j_symbol )
                self._labels2bondlength[ i_symbol ][ j_symbol ] = bond_length

    
    def XXXfindEndGroupContacts(self):
        """
        work out an atom connected to the endGroup so that we can define the angle 
        for bonding
        """
        
        MARGIN = 0.1 
        
        # Trundle through all the endGroups
        for e in self.endGroups:
            e_coord = self.coords[e]
            e_symbol = self.symbols[e]
            
            # Loop through all atoms
            for i, coord in enumerate( self.coords ):
                
                # Ignore self
                if i == e:
                    continue
                
                i_symbol = self.symbols[i]
                bond_length = util.bondLength( i_symbol , e_symbol )
                
                # See if these are bonded
                #if ( bond_length - MARGIN < numpy.linalg.norm( coord - e_coord ) < bond_length + MARGIN ):
                if ( bond_length - MARGIN < self.distance( e_coord, coord ) < bond_length + MARGIN ):
                    
                    # Found a bonded atom so add it to the list and move on
                    self._endGroupContacts[e] = i
                    break
    
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
        
        self.createFromLabelAndCoords( labels, coords )


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

        self.createFromLabelAndCoords( labels, coords )
    
    def angle(self,c1,c2,c3):
        """Return the angle in radians c1---c2---c3
        where c are the coordinates in a numpy array
        Taken from the CCP1GUI
        jmht - think about PBC
        """
        #p1 = a1.coord
        #p2 = a2.coord
        #p3 = a3.coord

        #r1 = cpv.distance(p1,p2)
        #r2 = cpv.distance(p2,p3)
        #r3 = cpv.distance(p1,p3)
        
        #r1 = numpy.linalg.norm( c1 - c2 )
        #r2 = numpy.linalg.norm( c2 - c3 )
        #r3 = numpy.linalg.norm( c1 - c3 )
        r1 = util.distance( c2, c1 )
        r2 = util.distance( c3, c2 )
        r3 = util.distance( c3, c1 )
        
        small = 1.0e-10
        #cnv   = 57.29577951
        if r1 + r2 - r3 < small:
            # printf("trig error %f\n",r3-r1-r2)
            # This seems to happen occasionally for 180 angles 
            theta = 180.0
        else:
            theta = math.acos( (r1*r1 + r2*r2  - r3*r3) / (2.0 * r1*r2) )
        return theta;

    def bondLength(self,symbola, symbolb):
        """Return the bond length between two atoms of the given type"""
        return self._labels2bondlength[symbola][symbolb]
    
    def dihedral(self, p1, p2, p3, p4):
        """ From the CCP1GUI"""

        #cnv=57.29577951

        vec_ij = p1 - p2
        vec_kj = p3 - p2
        vec_kl = p3 - p4

        # vec1 is the normal to the plane defined by atoms i, j, and k    
        vec1 = numpy.cross(vec_ij,vec_kj)
        magvec1 = numpy.dot(vec1,vec1)

        #  vec2 is the normal to the plane defined by atoms j, k, and l
        vec2 = numpy.cross(vec_kl,vec_kj)
        magvec2 = numpy.dot(vec2,vec2)

        # the definition of a dot product is used to find the angle between  
        # vec1 and vec2 and hence the angle between the planes defined by    
        # atoms i, j, k and j, k, l                                          
        #                                                                    
        # the factor of pi (180.0) is present since when we defined the      
        # vectors vec1 and vec2, one used the right hand rule while the      
        # other used the left hand rule                                      

        dotprod = numpy.dot(vec1,vec2)
        #print magvec1, magvec2
        #print type(magvec1), type(magvec2)
        fac = dotprod / math.sqrt(magvec1*magvec2)
        if(fac > 1.0):
            fac = 1.0
        if(fac < -1.0):
            fac = -1.0
        dihed = 180.0 - util.RADIANS2DEGREES * math.acos(fac )

        # the dot product between the bond between atoms i and j and the     
        # normal to the plane defined by atoms j, k, and l is used to        
        # determine whether or not the dihedral angle is clockwise or        
        # anti_clockwise                                                     
        #                                                                    
        # if the dot product is positive, the rotation is clockwise          

        sign_check = numpy.dot(vec_ij,vec2)
        if( sign_check > 0.0):
            dihed = dihed * -1.0

        return dihed

    def endGroupContactIndex(self, endGroupIndex ):
        """Return the index of the contact atom for this endGroup"""
        return self._endGroupContacts[ endGroupIndex ]
        

    def newBondPosition(self, targetEndGroupIndex, newBlock, newEndGroupIndex):
        """Return the coord where NewBlock would coord its bonding atom if joining
         to this one at targetEndGroupIndex by newEndGroupIndex
         I'm sure this algorithm is clunky in the extreme...
        """
        
        targetEndGroup = self.coords[ targetEndGroupIndex ]
        targetContact = self.endGroupContactIndex( targetEndGroupIndex )
        
        label = self.labels[ targetEndGroupIndex ]
        #jmht - use symbols!
        targetSymbol = util.label2symbol(label)
        
        label = newBlock.labels[ newEndGroupIndex ]
        newSymbol = util.label2symbol(label)
        
        # Get the bond length between these two atoms
        bondLength = util.bondLength(targetSymbol, newSymbol)
        
        # vector from target EndGroup contact to Endgroup
        bondVec =  targetEndGroup - targetContact 
        
        # Now extend it by Bondlength
        vLength = numpy.linalg.norm( bondVec )
        ratio = ( vLength + bondLength ) / vLength
        newVec = bondVec * ratio
        diff = newVec-bondVec
        newPosition = targetEndGroup + diff
        
        return newPosition
    
    def maxAtomRadius(self):
        """Return the maxium atom radius
        """
        
        rmax=0
        for r in self.atom_radii:
            if r>rmax:
                rmax=r
        return rmax
    
    def randomEndGroupIndex(self):
        """Return a randomly picked endgroup"""
        return random.choice( self.endGroups )
 
    def radius(self):
        
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self.calcRadius()
            self._changed = False
        
        return self._radius


    def rotate( self, axis, angle ):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
        rmat = self.rotation_matrix( axis, angle )
        
        # loop through and change all coords
        for i,coord in enumerate( self.coords ):
            self.coords[i] = numpy.dot( rmat, coord )
        
        self._changed = True

    def rotation_matrix(self, axis,angle):
        """
        http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
        """
        
        axis = axis/numpy.sqrt( numpy.dot(axis,axis) )
        a = numpy.cos(angle/2)
        b,c,d = -axis*numpy.sin(angle/2)
        return numpy.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
        
    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # Make sure its a numpy array
        if isinstance(tvector,list):
            tvector = numpy.array( tvector )

        for i in range( len(self.coords) ):
            self.coords[i] += tvector
        
        self._changed = True
    
    def translateCentroid( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given coord
        """
        self.translate( position - self.centroid() )
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
            
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.writeXyz( "{}\n".format(len(self.coords)) )
            f.writeXyz( "id={}\n".format(str(id(self))) )
                             
            for i,c in enumerate(self.coords):
                f.writeXyz("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol(self.labels[i]), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
    def update(self):
        """
        Run any methods that need to be done when the coordinates are changed
        """
        # Recalculate the data for this new block
        self.calcCenterOfMassAndGeometry()
        self.calcRadius()
        
    def __str__(self):
        """
        Return a string representation of the molecule
        """
        
        mystr = ""
        mystr += "BlockID: {}\n".format(id(self))
        mystr += "{}\n".format(len(self.coords))
        for i,c in enumerate(self.coords):
            #mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self.labels[i], c[0], c[1], c[2])
            mystr += "{0:5} {1:0< 15} {2:0< 15} {3:0< 15} \n".format( self.labels[i], c[0], c[1], c[2])
            
        mystr += "radius: {}\n".format( self.radius() )
        mystr += "COM: {}\n".format( self._centerOfMass )
        mystr += "COG: {}\n".format( self._centerOfGeometry )
        mystr += "endGroups: {}\n".format( self.endGroups )
        mystr += "endGroupContacts: {}\n".format( self._endGroupContacts )
        
        return mystr
        
        
class TestBuildingBlock(unittest.TestCase):

    def makeCh4(self):
        """Create a CH4 molecule for testing"""
        
#        coords = [ numpy.array([  0.000000,  0.000000,  0.000000 ] ),
#        numpy.array([  0.000000,  0.000000,  1.089000 ]),
#        numpy.array([  1.026719,  0.000000, -0.363000 ]),
#        numpy.array([ -0.513360, -0.889165, -0.363000 ]),
#        numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]
#
#        #numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]
#        labels = [ 'C', 'H', 'H', 'H', 'H' ]
#        
#        endGroups = [ 1,2,3,4 ]
#        
#        endGroupContacts = { 1:0, 2:0, 3:0, 4:0 }
#        
#        ch4 = BuildingBlock()
#        ch4.createFromArgs( coords, labels, endGroups, endGroupContacts )

        ch4 = BuildingBlock()
        ch4.fromXyzFile("../ch4.xyz")
        return ch4
    
    def makePaf(self):
        """Return the PAF molecule for testing"""
        
        paf = BuildingBlock()
        paf.fromCarFile("../PAF_bb_typed.car")
        return paf
    
    
    def testAaaReadCar(self):
        """
        Test we can read a car file - needs to come first
        """
        
        paf = self.makePaf()
        self.assertTrue( paf._endGroupContacts == {17: 3, 12: 2, 22: 4, 7: 1}, "Incorrect reading of endGroup contacts")

    def testAaaReadXyz(self):
        """
        Test we can read an xyz file - needs to come first
        """
        
        paf = self.makeCh4()
        self.assertTrue( paf._endGroupContacts == { 1:0, 2:0, 3:0, 4:0 }, "Incorrect reading of endGroup contacts")
        
        
        
    def testBond(self):
        """Test we can correctly bond two blocks at the given bond"""
        
        ch4 = self.makeCh4()
        m2 = ch4.copy()
        m2.translate( numpy.array( [2, 2, 2] ) )
        
        bond = (3,2)
        ch4.bond(m2, bond)
        
        self.assertTrue( ch4.endGroups == [1,2,4,6,8,9], "Incorrect calculation of endGroups")
        self.assertTrue( ch4._endGroupContacts == {1: 0, 2: 0, 4: 0, 6: 5, 8: 5, 9: 5}, "Incorrect calculation of endGroup contacts")
    

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
        
#        #m.writeXyz(name="after.xyz")
#        print m
#        print m.centerOfMass()
#        m.translate([1,1,1])
#        
#        #m.writeXyz(name="trans.xyz")
#        print m
#        print m.centerOfMass()

    def XtestClash(self):
        """
        Test we can spot a clash
        """
        
        
        paf = self.makePaf()
        
        block = paf.copy()
        self.assertTrue( paf.clash( block ) )
        
        radius = block.radius()
        block.translateCentroid( [0,0,radius*2+2.0] )
        self.assertFalse( paf.clash( block ) )
    
    def XtestClose(self):
        """Test we can check if molecules are close"""
        
        paf = self.makePaf()
        
        m = paf.copy()
        radius = paf.radius()
        m. translate( numpy.array( [0,0,radius*2+1.0] ) )
        
        self.assertFalse( paf.close( m, blockMargin=0.9 ), "Not close with 0.1 margin")
        
        self.assertTrue( paf.close( m, blockMargin=1.1 ), "Close with 0.1 margin")
        
        
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
        self.assertAlmostEqual(r, 1.78900031214, 7, "Incorrect radius: {}".format(str(r)) )
        
    def testRotate(self):
        """
        Test the rotation
        """
        
        ch4 = self.makeCh4()
        array1 = numpy.array([ -0.51336 ,  0.889165, -0.363 ])
        self.assertTrue( numpy.array_equal( ch4.coords[4], array1 ),
                         msg="testRotate arrays before rotation incorrect.")
        
        axis = numpy.array([1,2,3])
        angle = 2
        ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose( ch4.coords[4], array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")

        

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
        