'''
A block holding multiple subblocks
'''

# Python imports
import copy
import os
import random
import unittest

# external imports
import numpy

# local imports
import buildingBlock
import util

class CompositeBlock( buildingBlock.Block ):
    '''
    classdocs
    
    Parent must know about all child blocks
    
    call methods such as endGroup recursively - how to count how many _endGroups?
        
        
    
    def randomBlock(): - part of cell
    
    
    def newBlock(): - part of cell
        
    
    
    newBlock() - create a new block
    
    centroid
    copy
    flip
    maxAtomRadius
    mass
    positionGrowBlock
    
    
    '''

    def __init__( self, inputFile ):
        '''
        Constructor
        '''
        
        block = buildingBlock.Block( infile=inputFile )
        
        # List of all the blocks that this one contains - will be as many as there are _endGroups
        # For terminal blocks (those not containing other blocks we need to be careful we don't
        # iterate through these
        # Needs to filled with Nones for as many _endGroups as there are on init
        self.blocks = [ block ]
        
        # List of all blocks contained or linked to this one - includes iterating over sub-blocks
        self.allBlocks = [ block ]
        
        # The block containing this one or None if we are the top block/an isolated block
        self.parent = None
        
        # List of the atoms that are available to bond
        self._endGroups = []
        
        # List of the atoms that define the angle for the bonds (same order as _endGroups)
        self._angleAtoms = []
        
        
        
    def endGroup( self, idxEndGroup ):
        i,j = self._endGroups[ idxEndGroup ]
        block = self.blocks[i]
        return block.coords[ block._endGroups[j] ]
        
    
    def endGroupAngleAtom( self, idxEndGroup ):
        i,j = self._angleAtoms[ idxEndGroup ]
        block = self.blocks[i]
        return block.coords[ block._angleAtom[j] ]
        
    
    def bondBlock( self, idxEndGroup, newBlock, idxNewBlockEndGroup ):
        """ Add newBlock to this one
        """
        
        if newBlock not in self.blocks:
            self.blocks[ idxEndGroup ] = newBlock
        else:
            # Circular fragment
            print "Block was already in self"
            
        # Now add all subblocks to this one
        for block in newBlock.blocks:
            if block:
                self.allBlocks += block.allBlocks
            
        newBlock.blocks[ idxNewBlockEndGroup ] = self
        
        self.update( newBlock )
            
#    def removeBlock( self, block ):
#        """Delete the block from this one"""
#        self.blocks.pop( block )
#        self.bondedEndGroups.pop( ?? ) 
#        self.update()
        
    def update( self, newBlock ):
        if self.parent:
            self.parent.update( newBlock )
        else:
            self._update( newBlock )
    
    def _update( self, newBlock ):
        """Set the list of _endGroups"""
        
        if newBlock in self.blocks:
            raise RuntimeError, "Bonding block that is already in this - must fix..."
        
        # Loop over every block contained in this and all subBlocks
        for i, block in enumerate( self.allBlocks ):
            for j in range( len(block._endGroups) ):
                # Don't add _endGroups where there is a block attached
                if not block.blocks[j]:
                    self._endGroups.append(  (i, j) )
                    self._angleAtoms.append( (i, j) )

    def rotate(self, foo):
        """For any block"""
        for block in self.allBlocks:
            block._rotate( foo )
        
    def alignBond(self, idxBlockEG, refVector ):
        """
        Align this block, so that the bond defined by idxBlockEG is aligned with
        the refVector
        
        This assumes that the block has already been positioned with the contact atom at the origin
        """
        
        endGroup = self.coords[ idxBlockEG ]
        
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
            
#        endGroup =self.coords[ idxBlockEG ]
#        print "alignBlock AFTER: {0} | {1}".format( endGroup, refVector )
        return


    def join( self, idxEG, block, idxBlockEG ):
        """Bond block to this one at the _endGroups specified by the two indices
        endGroup = [ 7, 12, 17, 22 ]
        defAtom = [ 5, 6, 7, 8 ]
        
        idxEndGroup = randomEndGroupIndex()
        
        endGroup = block.coords[ block._endGroups[idxEndGroup] ]
        
        angleAtom = block.coords[ block._angleAtoms[idxEndGroup] ]
        
              H1
              |
              C0
          /  |  \
        H2   H3  H4
        
              H1             H6
              |              |
             C0             C5
          /  |  \        /  |  \
        H2   H3  H4 -- H7   H8  H9
        
        CH4 - 1
        _endGroups = [ 1,2,3,4 ]
        _angleAtoms = [ 0,0,0,0 ]
        bondedEndGroups = []
        
        CH4 - 2
        _endGroups = [ 1,2,3,4 ]
        _angleAtoms = [ 0,0,0,0 ]
        bondedEndGroups = []
        
        join 2,3 
        _endGroups = [ 1,2,4,6,7,8 ]
        _angleAtoms = [ 0,0,0,5,5,5 ]
        
        
        # Need to remember the _endGroups/_angleAtoms for each block
        # First is this block - although not used unless all blocks removed to leave us with ourselves
        blockEndGroups = [  [ 1,2,3,4 ], [ 1,2,3,4 ] ]
        blockAngleAtoms = [ [ 0,0,0,0 ], [ 0,0,0,0 ] ]
        
        So remove the bond we are making from both blocks and then add the new _endGroups/_angleAtoms
        # Remember the actual indices of the _endGroups for the parent block in the bondedEndGroups array of tuples
        bondedEndGroups.append(  ( _endGroups.pop( 2 ), _angleAtoms.pop( 2 ) ) )
        
        # Add the new free endgroups and their _angleAtoms to the array - we increment their indices accordingly
        for i, EG in enumerate( block._endGroups ):
            _endGroups.append( EG + start )
            _angleAtoms.append( block._angleAtoms[i] + start )
            
        
        # random
        0 1 2 3
        H-C=C-H
        _endGroups = [ 0,3 ]
        _angleAtoms = [ 1,2 ]
        
        0 1 2 3 4 5 6 7
        H-C=C-H_H-C=C-H
        _endGroups = [ 0, 7 ]
        _angleAtoms = [ 1, 6 ]
        
        0 1 2 3 4 5 6 7 8 9 1011
        H-C=C-H_H-C=C-H_H-C=C-H
        _endGroups = [ 0, 11 ]
        _angleAtoms = [ 1, 10 ]
        
        removeBlock( 0 )
        0 1 2 3 4 5 6 7
        H-C=C-H_H-C=C-H
        _endGroups = [ 0, 7 ]
        _angleAtoms = [ 1, 6 ]
        
        and 
        0 1 2 3
        H-C=C-H
        _endGroups = [ 0,3 ]
        _angleAtoms = [ 1,2 ]
        
        
        
        """
        
        # Needed to work out how much to add to the block indices
        start = len( self.coords )
        end = len( block.coords )
        
        # Keep track of where the new block starts and ends so we can remove it
        self.blocks.append[ (start, end)  ]
        
        # Add all the data structures together
        self.coords.extend( block.coords )
        self.atom_radii.extend( block.atom_radii )
        self.labels.extend( block.labels )
        self.symbols.extend( block.symbols )
        self.masses.extend( block.masses )
        self.atomCell.extend( block.atomCell )
        
        # Need to remove the end groups used in the bond
        self.bondedEndGroups.append( )
        
        i = self._endGroups.index( bond[0] )
        self._endGroups.pop( i )
        i = block._endGroups.index( bond[1] )
        block._endGroups.pop( i )
        
        # Now add block to self, updating index
        for i in block._endGroups:
            self._endGroups.append( i+start )
        
        del self._angleAtoms[ bond[0] ]
        del block._angleAtoms[ bond[1] ]
        
        for k,v in block._angleAtoms.iteritems():
            self._angleAtoms[ k+start ] = v+start
        
        self.update()

        
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
    
    def copy( self ):
        """
        Create a copy of ourselves and return it.
        Hit problem when passing in the distance function from the cell as the deepcopy
        also had a copy of the cell
        """
        return copy.deepcopy(self)


    def endGroupContactIndex(self, endGroupIndex ):
        """Return the index of the contact atom for this endGroup"""
        return self._angleAtoms[ endGroupIndex ]
        
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

    def maxAtomRadius(self):
        """Return the maxium atom radius
        """
        
        rmax=0
        for r in self.atom_radii:
            if r>rmax:
                rmax=r
        return rmax
    
    def mass(self):
        """Total mass of the molecule
        """
        if not self._mass:
            self._mass = 0.0
            for m in self.masses:
                self._mass += m
        return self._mass
                

    def positionGrowBlock( self, idxBlockEG, growBlock, idxGrowBlockEG ):
        """
        Position growBlock so it can bond to us, using the given _endGroups
        
        Arguments:
        idxBlockEG: the index of the endGroup to use for the bond
        growBlock: the block we are positioning
        idxGrowBlockEG: the index of the endGroup to use for the bond
        """

        # The vector we want to align along is the vector from the contact
        # to the endGroup
        endGroup =  self.coords[ idxBlockEG ]
        idxContact =  self.endGroupContactIndex( idxBlockEG )
        contact  =  self.coords[ idxContact ]
        refVector      =  endGroup - contact
        
        # Get the contact for the growBlock
        idxGrowBlockContact = growBlock.endGroupContactIndex( idxGrowBlockEG )
        growBlockContact = growBlock.coords[ idxGrowBlockContact ]
        
        # get the coord where the next block should bond
        # symbol of endGroup tells us the sort of bond we are making which determines
        # the bond length
        symbol = growBlock.symbols[ idxGrowBlockEG ]
        bondPos = self.newBondPosition( idxBlockEG, symbol )
        #print "got bondPos: {}".format( bondPos )
        
        # Shift block so contact at center, so the vector of the endGroup can be
        # aligned with the refVector
        growBlock.translate( -growBlockContact )
        
        # Align along the staticBlock bond
        growBlock.alignBond( idxGrowBlockEG, refVector )
        
        # Now turn the second block around
        # - need to get the new vectors as these will have changed due to the moves
        growBlockEG = growBlock.coords[ idxGrowBlockEG ]
        
        # Flip by 180 along bond axis
        growBlock.flip( growBlockEG )
        
        # Now need to place the endGroup at the bond coord - at the moment
        # the contact atom is on the origin so we need to subtract the vector to the endGroup
        # from the translation vector
        #jmht FIX!
        growBlock.translate( bondPos + growBlockEG )
        
        return

    def newBondPosition(self, idxTargetEG, symbol ):
        """Return the position where a bond to an atom of type 'symbol'
        would be placed if bonding to the target endgroup
         I'm sure this algorithm is clunky in the extreme...
        """
        
        targetEndGroup = self.coords[ idxTargetEG ]
        targetContactIndex = self.endGroupContactIndex( idxTargetEG )
        targetContact = self.coords[ targetContactIndex ]
        targetSymbol = self.symbols[ idxTargetEG ]
        
        # Get the bond length between these two atoms
        bondLength = util.bondLength(targetSymbol, symbol)
        
        # vector from target EndGroup contact to Endgroup
        bondVec =  targetEndGroup - targetContact 
        
        # Now extend it by Bondlength
        vLength = numpy.linalg.norm( bondVec )
        ratio = ( vLength + bondLength ) / vLength
        newVec = bondVec * ratio
        diff = newVec-bondVec
        newPosition = targetEndGroup + diff
        
        return newPosition
 
    def radius(self):
        
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self.calcRadius()
            self._changed = False
        
        return self._radius

    def randomEndGroupIndex(self):
        """Return a randomly picked endgroup"""
        return random.choice( self._endGroups )
    
    def rotate( self, axis, angle, center=None ):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        
        if center==None:
            center = numpy.array([0,0,0])
        
        rmat = util.rotation_matrix( axis, angle )
        
        # loop through all blocks and change all coords
        # Need to check that the center bit is the best way of doing this-
        # am almost certainly doing more ops than needed
        for block in self.blocks:
            for i,coord in enumerate( block.coords ):
                #self.coords[i] = numpy.dot( rmat, coord )
                coord = coord - center
                coord = numpy.dot( rmat, coord )
                block.coords[i] = coord + center
        
        self._changed = True

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        
        # Make sure its a numpy array - is this needed?
        if isinstance(tvector,list):
            tvector = numpy.array( tvector )

        for block in self.blocks:
            # Use len as we don't need to return hthe coords
            for i in range( len (block.coords ) ):
                block.coords[i] += tvector
        
        self._changed = True
    
    def translateCentroid( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given coord
        """
        self.translate( position - self.centroid() )
        
    def update(self):
        """
        Run any methods that need to be done when the coordinates are changed
        """
        # Recalculate the data for this new block
        self.calcCenterOfMassAndGeometry()
        self.calcRadius()
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
            
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format(len(self.coords)) )
            f.write( "id={}\n".format(str(id(self))) )
                             
            for i,c in enumerate(self.coords):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( util.label2symbol(self.labels[i]), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
        
    def __str__(self):
        """
        Return a string representation of the composite block
        """
        
        mystr = ""
        mystr += "BlockID: {} with {} blocks\n".format( id(self), len(self.blocks) )
        
        for block in self.blocks:
            mystr += str(block)
        
        return mystr
        
        
class TestCompositeBlock(unittest.TestCase):

    def makeH2C2(self):
        """
        Return H2C2 molecule for testing
        """
#        coords = [ 
#                  numpy.array([  0.000000,  0.000000,  0.000000 ] ),
#                  numpy.array([  0.000000,  0.000000,  1.000000 ]),
#                  numpy.array([  0.000000,  0.000000,  2.000000 ]),
#                  numpy.array([  0.000000,  0.000000,  3.000000 ]),
#                  ]
#
#        labels = [ 'H', 'C', 'C', 'H' ]
#        _endGroups = [ 0, 3 ]
#        endGroupContacts = { 0:1, 3:2 }
#        b = buildingBlock.Block()
#        b.createFromArgs( coords, labels, _endGroups, endGroupContacts )
        
        
        cb = CompositeBlock( inputFile="../h2c2.xyz" )
        return cb
    
    
    def testNew(self):
        """
        test new coordinate class
        """
        
        b1 = self.makeH2C2()
        print b1
        return
        b2 = b1.copy()
        
        idxB1EG = 0
        idxB2EG = 3
        
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
        