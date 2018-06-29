'''
Created on Jan 15, 2013

@author: abbietrewin
'''
import collections
import copy
import csv
import logging
import os

import numpy as np

# our imports
import xyz_core
import xyz_util

ENDGROUPBONDED = '*'
logger = logging.getLogger()

class RigidParticle(object):
    def __init__(self, body, bodyIdx, cell_dim=None, center=False):
        # Attributes of the central particle
        self.image = np.array([0, 0, 0])
        self.mass = None
        self.position = None
        self.principalMoments = None
        self.natoms = 0
        self.type = None
        # Attributes if the constituent particles
        self.b_charges = None # NB use index to prevent storing these twice
        self.b_diameters = None # NB use index to prevent storing these twice
        self.b_masses = None
        self.b_positions = None
        self.b_atomTypes = None
        
        self.fromBody(body, bodyIdx, cell_dim=cell_dim, center=center)
    
    def fromBody(self, body, bodyIdx, cell_dim=None, center=False):
        coords = body.coords
        self.natoms = coords.shape[0]
        com = xyz_core.centreOfMass(coords, body.masses)
        self.b_positions = coords - com # coordinate positions relative to com
        if cell_dim is not None:
            com, self.image = xyz_core.wrapCoord3(com, dim=cell_dim, center=center)
        self.position = com
        self.mass = np.sum(body.masses)
        self.principalMoments = xyz_core.principalMoments(coords, body.masses)
        self.type = "CP%d" % bodyIdx
        # Specify properties of consituent particles
        # To save memory can move to just saving the indices of the elements in the data array
        self.b_charges = body.charges
        self.b_diameters = body.diameters
        self.b_masses = body.masses
        self.b_atomTypes = body.atomTypes
        return self
    

class Body(object):
    """Proxy object for the body within a fragment"""
    def __init__(self, fragment, bodyIndex):
        self.fragment = fragment
        self.bodyIndex = bodyIndex
        self.indexes = self._setup_indexes()
        self.natoms = np.sum(self.indexes)
        self._centreOfMass = xyz_core.centreOfMass(self.coords, self.masses)
        return

    def _setup_indexes(self):
        return (np.array(self.fragment._bodies) == self.bodyIndex) &~ np.array(self.fragment.masked)

    @property
    def atomTypes(self):
        return list(np.compress(self.indexes, self.fragment._atomTypes, axis=0))

    @property
    def bodies(self):
        return [self.bodyIndex] * self.natoms
        
    def centreOfMass(self):
        return self._centreOfMass

    @property
    def charges(self):
        return list(np.compress(self.indexes, self.fragment._charges, axis=0))

    @property
    def coords(self):
        return np.compress(self.indexes, self.fragment._coords, axis=0)

    @property
    def diameters(self):
        return [ xyz_util.DUMMY_DIAMETER ] * self.natoms

    @property
    def masked(self):
        mask = []
        for i in [ i for i in self.fragment._ext2int.values() if self.fragment._bodies[i] == self.bodyIndex]:
            if (hasattr(self.fragment, 'unBonded') and self.fragment.unBonded[i]) or self.fragment._atomTypes[i].lower() == 'x':
                mask.append(True)
            else:
                mask.append(False)
        return mask

    @property
    def mass(self):
        return np.sum(self.masses())

    @property
    def masses(self):
        return np.compress(self.indexes, self.fragment._masses, axis=0)
    
    @property
    def principalMoments(self):
        return xyz_core.principalMoments(self.coords(), self.masses)
    
    @property
    def XrigidType(self):
        """return the type of this body based on the endGroup configuration"""
        return "{}{}{}".format(self.fragment.fragmentType, self.fragment.configStr, self.bodyIndex)
    
    def rigidParticle(self, bodyIdx, cell_dim=None, center=False):
        return RigidParticle(self, bodyIdx, cell_dim=cell_dim, center=center)

    @property
    def static(self):
        if hasattr(self.fragment, 'static') and self.fragment.static:
            v = True
        else:
            v = False
        return [ v ] * self.natoms

    @property
    def symbols(self):
        return list(np.compress(self.indexes, self.fragment._symbols, axis=0))

class EndGroup(object):

    def __init__(self):

        self.blocked = False # Used for indicating that this endGroup is unbonded but not free
        self.bonded = False

        self._endGroupType = None

        self.fragment = None

        self.fragmentEndGroupIdx = None
        self.blockEndGroupIdx = None

        self.fragmentCapIdx = None
        self.blockCapIdx = None
        self.capBondLength = None

        self.fragmentDihedralIdx = -1
        self.blockDihedralIdx = -1

        self.fragmentUwIdx = -1
        self.blockUwIdx = -1

        return

    def block(self):
        f = True
        b = True
        if self.fragment.block is None: b = False
        if self.fragment is None: f = False
        if not b and f:
            raise RuntimeError("None Block {0} Fragment {1}\n{2}".format(b, f, self))
        return self.fragment.block

    def capIdx(self):
        """Return the index of the endGroup atom in external block indices"""
        return self.blockCapIdx

    def bondedCatalyst(self):
        """Return True if this endGroup belongs to a catalyst that is bonded to another catalyst"""
        return self.fragment.catalyst and self._endGroupType.endswith(ENDGROUPBONDED)

    def coord(self,endGroup=True):
        """Need to think about an API for accessing coordinates for endGroups
        This just hacks in returning the endGroup.
        """
        return self.fragment._coords[self.fragmentEndGroupIdx]

    def dihedralIdx(self):
        """Return the index of the dihedral atom in external block indices"""
        return self.blockDihedralIdx

    def endGroupIdx(self):
        """Return the index of the endGroup atom in external block indices"""
        return self.blockEndGroupIdx

    def fragmentType(self):
        return self.fragment.fragmentType

    def free(self):
        return not self.bonded and not self.blocked

    def setBonded(self, bond):
        self.bonded = True
        self.fragment.addBond(self, bond)
        return

    def unBond(self, bondEndGroup):
        self.bonded = False
        # HACK WE REMOVE ALL SUFFIXES
        for eg in self.fragment.endGroups():
            if eg._endGroupType.endswith(ENDGROUPBONDED):
                logger.debug("unBond ENDGROUPBONDED")
                eg._endGroupType = eg._endGroupType.rstrip(ENDGROUPBONDED)
        self.fragment.delBond(self.type())
        # Unmask cap and set the coordinate to the coordinate of the last block atom
        # NOTE - NEED TO SCALE BY CORRECT LENGTH
        if hasattr(bondEndGroup, 'coord'):
            # Reposition the cap atom based on the bond vector
            #self.fragment._coords[self.fragmentCapIdx] = bondEndGroup.coord()
            # Get vector from this endGroup to the other endGroup
            egPos = self.fragment._coords[self.fragmentEndGroupIdx]
            v1 = bondEndGroup.coord() - egPos
            # Now get unit vector
            uv = v1 / np.linalg.norm(v1)
            # calculate noew position
            self.fragment._coords[self.fragmentCapIdx] = egPos + (uv * self.capBondLength)

        # Unhide the cap atom
        self.fragment.masked[self.fragmentCapIdx] = False
        # Mark capAtom as unBonded so that it won't be included in the optimisation
        self.fragment.unBonded[self.fragmentCapIdx] = True

        if self.fragmentUwIdx != -1:
            raise RuntimeError("Cannot unbond masked endGroups yet!")
            self.fragment.masked[ self.fragmentUwIdx ] = True
        self.fragment.update()
        return

    def type(self):
        """The type of the endGroup"""
        return "{0}:{1}".format(self.fragment.fragmentType, self._endGroupType)

    def updateAncillaryIndices(self, cap2endGroup):
        """The block has been updated so we need to update our block indices based on where the
        fragment starts in the block"""
        if self.fragment.masked[self.fragmentCapIdx]:
            assert self.bonded
            # If this endGroup is involved in a bond, we want to get the block index of the
            # opposite endGroup as this has now become our cap atom
            self.blockCapIdx = cap2endGroup[(self.fragment, self.fragmentCapIdx)]
        else:
            self.blockCapIdx = self.fragment._int2ext[self.fragmentCapIdx] + self.fragment.blockIdx

        # -1 means no dihedral or uw atom set
        if self.fragmentDihedralIdx != -1:
            if self.fragment.masked[self.fragmentDihedralIdx]:
                # Now work out which atom this is bonded to
                self.blockDihedralIdx = cap2endGroup[(self.fragment, self.fragmentDihedralIdx)]
            else:
                self.blockDihedralIdx = self.fragment._int2ext[self.fragmentDihedralIdx] + self.fragment.blockIdx

        # -1 means no uw atom and if it is masked then there will not be an external index
        if self.fragmentUwIdx != -1 and not self.fragment.masked[self.fragmentUwIdx]:
            # assert not self.fragment.isMasked( self.fragmentUwIdx )
            self.blockUwIdx = self.fragment._int2ext[self.fragmentUwIdx] + self.fragment.blockIdx
        return

    def updateEndGroupIndex(self):
        """The block has been updated so we need to update our block indices based on where the
        fragment starts in the block"""
        self.blockEndGroupIdx = self.fragment._int2ext[self.fragmentEndGroupIdx] + self.fragment.blockIdx
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
        #EndGroup:
        s = "{0} {1}:{2} {3}({4})->{5}({6}) {7}->{8}".format(self.type(), id(self.fragment), id(self),
                                                    self.fragmentEndGroupIdx, self.fragment._symbols[self.fragmentEndGroupIdx],
                                                    self.fragmentCapIdx,self.fragment._symbols[self.fragmentCapIdx],
                                                    self.blockEndGroupIdx, self.blockCapIdx)

        return s


class FragmentConfigManager(object):
    """Manages the configuration of each individual fragmentType
    
    As Fragments have their own class/instance managment machinary, we can't use standard class variables
    to track variables that are shared across all fragments of a given fragmentType. This class provides
    the functionality to track attributes that are shared by all fragments of a given type, but also need
    to be updated as the simulation progress.
    
    Need to think about how to reset this when running simulations and when unpickling
    """
    
    def __init__(self):
        self.configs = {}
    
    @staticmethod
    def calcConfigStr(i):
        import string
        alph = string.ascii_uppercase
        num_chars = 26 # letters in alphabet - unlikely to change...
        assert i < num_chars * num_chars, "Too many configurations!"
        return alph[int(i/num_chars)] + alph[i%num_chars]
    
    def reset(self):
        self.configs = {}

    def updateConfig(self, fragment):
        if fragment.fragmentType not in self.configs:
            self.configs[fragment.fragmentType] = []
        if fragment.config not in self.configs[fragment.fragmentType]:
            self.configs[fragment.fragmentType].append(fragment.config)
        idx = self.configs[fragment.fragmentType].index(fragment.config)
        return self.calcConfigStr(idx)
            
configManager = FragmentConfigManager()

class Fragment(object):
    '''
    classdocs
    '''
    def __init__(self,
                 filePath=None,
                 fragmentType=None,
                 solvent=False,
                 static=False,
                 markBonded=False,
                 catalyst=False
                 ):
        '''
        Constructor
        '''
        # This first set of variables are shared by all fragments
        # When we copy a fragment the new fragment gets references to the variables
        # that were created for the first fragment (see copy)
        sharedAttrs = {
            '_atomTypes'       : [],
            '_bodies'          : [],  # a list of which body within this fragment each atom belongs to
            '_bonds'           : [],  # List of internal fragment bonds
            '_bonded'          : [],  # List of which atoms are bonded to which
            'catalyst'         : catalyst, # Whether this block is a catalyst
            '_cellParameters'  : {},
            '_charges'         : [],
            'fragmentType'     : fragmentType,
            '_labels'          : [],
            '_masses'          : [],
            'markBonded'       : markBonded,
            'onbondFunction' : None,  # A function to be called when we bond an endGroup
            '_radii'           : [],
            '_maxBonds'        : {},
            '_radius'          :-1,
            'solvent'          : solvent,  # True if this fragment is solvent and should be excluded from clashChecks
            'static'           : static,
            '_symbols'         : [],  # ordered array of symbols (in upper case)
            '_totalMass'       :-1,
            '_individualAttrs' : None,
            '_sharedAttrs'     : None,
            }

        #
        # The variables below here are specific to a fragment and change as the fragment moves
        # and is involved in bonds - each fragment gets its own copy of these
        individualAttrs = {
            'block'           : None,
            'config'          : None,
            'configStr'     : None,
            '_coords'         : [],
            '_centroid'       : None,
            '_centerOfMass'   : None,
            '_ext2int'        : collections.OrderedDict(),
            '_int2ext'        : collections.OrderedDict(),
            '_centerOfMass'   : None,
            '_maxAtomRadius'  :-1,
            '_changed'        : True,  # Flag for when we've been moved and need to recalculate things
            'blockIdx'       : None,  # The index in the list of block data where the data for this fragment starts
            '_endGroups'      : [],  # A list of the endGroup objects
            '_endGroupBonded' : [],  # A list of the number of each endGroup that are used in bonds
            'masked'         : [],  # bool - whether the atoms is hidden (e.g. cap or uw atom)
            'unBonded'         : [],  # bool - whether an atom has just been unbonded
            }

        # Set as attributes of self
        for a, v in sharedAttrs.items():
            setattr(self, a, v)

        # Set as attributes of self
        for a, v in individualAttrs.items():
            setattr(self, a, v)

        # Set these manually
        self._individualAttrs = individualAttrs
        self._sharedAttrs = sharedAttrs

        # Create from the file
        if filePath:
            self.fromFile(filePath)
        return

    def _markBonded(self, endGroup, bond):
        """Append * to all endGroups in fragments that are involved in bonds to cat"""
        #if not endGroup.type() == 'cat:a': return
        if endGroup.type().endswith(ENDGROUPBONDED) or not (bond.endGroup1.fragment.catalyst or bond.endGroup2.fragment.catalyst):
            return
        logger.debug("_markBonded marking bonds for fragment {0}".format(self.fragmentType))
        for eg in self.endGroups():
            assert not eg._endGroupType.endswith(ENDGROUPBONDED),"Already got bonded endGroup"
            eg._endGroupType += ENDGROUPBONDED
        return

    def addBond(self, endGroup, bond):
        endGroupType = endGroup.type()

        # Mask fragment cap and uw atoms now
        self.masked[ endGroup.fragmentCapIdx ] = True
        if endGroup.fragmentUwIdx != -1:
            self.masked[ endGroup.fragmentUwIdx ] = True

        # Hack for starred endGroups
        if not endGroupType.endswith(ENDGROUPBONDED):
            self._endGroupBonded[ endGroupType] += 1
            # Handle maxBonds here
            if self._maxBonds[ endGroupType ] is not None:
                if self._endGroupBonded[ endGroupType ] >= self._maxBonds[ endGroupType ]:
                    # We have exceeded the bonding limit for these endGroupTypes, so we set any free ones
                    # of this type to blocked
                    for eg in self._endGroups:
                        if not eg.bonded and eg.type() == endGroupType:
                            eg.blocked = True

        # The user may have supplied a custom bonding function, so we call that here
        if self.onbondFunction:
            self.onbondFunction(endGroup)

        if  hasattr(self,'markBonded') and self.markBonded:
            self._markBonded(endGroup, bond)
        self.update()
        return

    def delBond(self, endGroupType):
        self._endGroupBonded[ endGroupType ] -= 1
        return

    def bodies(self):
        for bodyIdx in set(self._bodies):
            yield Body(self, bodyIdx)

    def body(self, idxAtom):
        return self._bodies[self._ext2int[idxAtom]]

    def bonds(self):
        """Return a list of bonds
        We exclude masked atoms"""
        bonds = []
        for b1, b2 in self._bonds:
            if not (self.masked[b1] or self.masked[b2]):
                # Map to internal and then add blockIdx to get position in block
                # b1 = b1 + self.blockIdx
                # b2 = b2 + self.blockIdx
                bonds.append((self._int2ext[b1], self._int2ext[b2]))
        return bonds

    def _calcBonded(self):
        assert len(self._bonds)
        # Create empty lists for all
        self._bonded = [ [] for _ in range(len(self._coords)) ]
        for b1, b2 in self._bonds:
            if b1 not in self._bonded[ b2 ]:
                self._bonded[ b2 ].append(b1)
            if b2 not in self._bonded[ b1 ]:
                self._bonded[ b1 ].append(b2)
        return

    def _calcCenters(self):
        """Calculate the center of mass and geometry for this fragment
        """
        self._centroid = np.sum(self._coords, axis=0) / np.size(self._coords, axis=0)
        self._totalMass = np.sum(self._masses)
        # Centre of mass is sum of the products of the individual coordinates multiplied by the mass, divded by the total mass
        self._centerOfMass = xyz_core.centreOfMass(self._coords, self._masses)
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
        distances = [xyz_core.distance(self._centroid, coord) for coord in self._coords]
        imax = np.argmax(distances)
        dist = distances[ imax ]
        # Add on the radius of the largest atom
        self._radius = dist + self.maxAtomRadius()
        return

    def _calcProperties(self):
        self._calcCenters()
        self._calcRadius()
        self._changed = False
        return

    def cellParameters(self):
        return self._cellParameters

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

    def charge(self, idxAtom):
        return self._charges[self._ext2int[idxAtom]]

    def coord(self, idxAtom, coord=None):
        if coord is not None:
            if isinstance(coord, list):
                coord = np.array(coord)
            self._coords[self._ext2int[idxAtom]] = coord
        return self._coords[self._ext2int[idxAtom]]

    def clearUnbonded(self):
        self.unBonded = [False] * len(self.unBonded)

    def copy(self):
        """Create a copy of ourselves.
        Those attributes in shared are just copied as references as they do not change between fragments
        of the same fragmentType
        Those in single are deep-copied as each fragment has its own"""

        f = Fragment()
        for a in f.__dict__:
            if a in self._sharedAttrs.keys():
                setattr(f, a, getattr(self, a))
            elif a in self._individualAttrs.keys():
                setattr(f, a, copy.deepcopy(getattr(self, a)))
            else:
                # HACKS FOR DEALING WITH OLD FILES
                msg = "Missing attribute in fragment copy: {0}".format(a)
                logger.critical(msg)
                raise RuntimeError(msg)
            # Update fragment references in the endGroups
            for e in f._endGroups:
                e.fragment = f
        return f

    def endGroups(self):
        return self._endGroups

    def endGroupTypes(self):
        """Return a list of the endGroupTypes in this fragment"""
        return set(self._endGroupBonded.keys())

    def fillData(self):
        """ Fill the data arrays from the label """

        self.masked = np.array([ False ] * len(self._coords))
        self.unBonded = [ False ] * len(self._coords)
        self._masses = np.array([ xyz_core.ATOMIC_MASS[ symbol ] for symbol in self._symbols ])
        self._totalMass = np.sum(self._masses)
        self._radii = np.array([ xyz_core.COVALENT_RADII[xyz_core.SYMBOL_TO_NUMBER[s.upper()]] * xyz_core.BOHR2ANGSTROM \
                                   for s in self._symbols])
        self._maxAtomRadius = np.max(self._radii)
        return

    def freeEndGroups(self):
        return [eg for eg in self._endGroups if eg.free()]

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
        with open(carFile, "r") as f:

            # skip first line
            f.readline()

            # 2nd states whether PBC: PBC=OFF
            pbc, state = f.readline().strip().split("=")
            assert pbc.strip() == "PBC"
            state = state.strip()

            # skip two lines
            f.readline()
            f.readline()
            if state.upper() == "ON":
                # Read in the PBC
                line = f.readline().strip()
                fields = line.split()
                assert fields[0].upper() == "PBC"
                self._cellParameters['A'] = float(fields[1])
                self._cellParameters['B'] = float(fields[2])
                self._cellParameters['C'] = float(fields[3])

            count = 0
            while reading:
                line = f.readline()

                line = line.strip()
                if not line:
                    logger.critical("END OF CAR WITH NO END!!!")
                    break
                fields = line.split()
                label = fields[0]

                # Check end of coordinates
                if label.lower() == "end":
                    reading = False
                    break

                labels.append(label)
                coords.append(np.array(fields[1:4], dtype=np.float64))
                atomTypes.append(fields[6])
                symbols.append(fields[7])
                charges.append(float(fields[8]))

                count += 1

        return  (coords, labels, symbols, atomTypes, charges)

    def fromXyzFile(self, xyzFile):
        """"Jens did this.
        """

        labels = []
        symbols = []
        atomTypes = []  # hack...
        charges = []

        # numpy array
        coords = []

        with open(xyzFile) as f:

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
                coords.append(np.array(fields[1:4], dtype=np.float64))
                symbol = xyz_util.label2symbol(label)
                symbols.append(symbol)
                atomTypes.append(symbols)
                charges.append(0.0)

                count += 1

        return (coords, labels, symbols, atomTypes, charges)

    def fromFile(self, filePath):

        if filePath.endswith(".car"):
            (coords, labels, symbols, atomTypes, charges) = self.fromCarFile(filePath)
        elif filePath.endswith(".xyz"):
            (coords, labels, symbols, atomTypes, charges) = self.fromXyzFile(filePath)
        else:
            raise RuntimeError("Unrecognised file suffix: {}".format(filePath))

        # Get cap atoms and endgroups
        endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms = self.parseEndgroupFile(filePath)

        # Set the root fragment and its attributes
        self.setData(coords=coords,
                     labels=labels,
                     symbols=symbols,
                     atomTypes=atomTypes,
                     charges=charges,
                     endGroupTypes=endGroupTypes,
                     endGroups=endGroups,
                     capAtoms=capAtoms,
                     dihedralAtoms=dihedralAtoms,
                     uwAtoms=uwAtoms
                  )

        self.processBodies(filePath)
        self.update()
        self._calcProperties()
        return

    def iterAtomTypes(self):
        """Generator to return the atomTypes"""
        for i in range(len(self._ext2int)):
            yield self.type(i)

    def iterCoord(self):
        """Generator to return the coordinates"""
        for i in range(len(self._ext2int)):
            yield self.coord(i)

    def label(self, idxAtom):
        return self._labels[self._ext2int[idxAtom]]

    def mass(self, idxAtom):
        return self._masses[self._ext2int[idxAtom]]

    def maxAtomRadius(self):
        assert self._maxAtomRadius > 0
        return self._maxAtomRadius

    def totalMass(self):
        """Return the total mass for this block"""
        assert self._totalMass > 0
        return self._totalMass

    def numAtoms(self):
        return len(self._ext2int)

    def numBondedEndGroups(self):
        """Return the total number of bonded endGroups in this fragment"""
        return sum([nbonded for nbonded in self._endGroupBonded.values() ])

    def parseEndgroupFile(self, filePath):

        dirname, filename = os.path.split(filePath)
        basename, suffix = os.path.splitext(filename)

        endGroupTypes = []
        endGroups = []
        capAtoms = []
        uwAtoms = []
        dihedralAtoms = []

        egfile = os.path.join(dirname, basename + ".csv")
        if not os.path.isfile(egfile):
            logger.critical("No endGroup definition file supplied for file: {0}".format(filePath))
            return endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms

        with open(egfile) as fh:
            csvreader = csv.reader(fh, delimiter=',', quotechar='"')
            for i, row in enumerate(csvreader):
                if i == 0:  # Header
                    if not len(row) >= 4 and row[0].lower() == "type" and \
                    row[1].lower() == "endgroup" and \
                    row[2].lower() == "capatom" and \
                    row[3].lower() == "delatom":
                        raise RuntimeError("First line of csv file must contain header line:\ntype,endgroup,capatom,dihedral,delatom")
                    continue

                # skip blank lines
                if not len(row):
                    continue

                # For now make sure first value is letter
                assert row[0][0].isalpha(), "First column of ambi file needs to be a letter!"
                endGroupTypes.append(row[0])
                endGroups.append(int(row[1]))
                capAtoms.append(int(row[2]))
                if len(row) >= 4 and row[3] and row[3] != -1:
                    dihedralAtoms.append(int(row[3]))
                else:
                    dihedralAtoms.append(-1)
                if len(row) >= 5 and row[4] and row[4] != -1:
                    uwAtoms.append(int(row[4]))
                else:
                    uwAtoms.append(-1)

            if self.fragmentType == 'cap' and len(endGroups) != 1:
                raise RuntimeError("Capfile had >1 endGroup specified!")

        return endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms

    def processBodies(self, filepath):
        """See if we split the fragment into bodies or not"""
        assert len(self._coords) > 0, "Coordinates must have been read before processing bodies"
        dirname, filename = os.path.split(filepath)
        basename, suffix = os.path.splitext(filename)
        bodyFile = os.path.join(dirname, basename + ".ambody")
        if os.path.isfile(bodyFile):
            self._bodies = np.array([ int(l.strip()) for l in open(bodyFile) ], dtype=np.int)
            assert len(self._bodies) == len(self._coords), \
            "Must have as many bodies as coordinates: {0} - {1}!".format(len(self.bodies), self._dataLen)
            assert self._bodies[0] == 0, "Bodies must start with zero!"
        else:
            # Just create an array with 0
            self._bodies = np.zeros(len(self._coords), dtype=np.int)
        return

    def radius(self, idxAtom):
        return self._radii[self._ext2int[idxAtom]]

    def totalRadius(self):

        if self._changed:
            self._calcProperties()
        return self._radius

    def rotate(self, rotationMatrix, center):
        """ Rotate the molecule about the given axis by the angle in radians
        """
        self._coords = self._coords - center
        #self._coords = np.array([ np.dot(rotationMatrix, c) for c in self._coords ])
        # I don't actually undestand why this works at all...
        # From: http://stackoverflow.com/questions/12148351/efficiently-rotate-a-set-of-points-with-a-rotation-matrix-in-numpy
        self._coords = np.dot(self._coords, rotationMatrix.T)
        self._coords = self._coords  + center
        self._changed = True
        return

    def setData(self,
                coords=None,
                labels=None,
                symbols=None,
                atomTypes=None,
                charges=None,
                endGroupTypes=None,
                endGroups=None,
                capAtoms=None,
                dihedralAtoms=None,
                uwAtoms=None
                ):

        self._charges = np.array([ c for c in charges])
        self._coords = np.array([ c for c in coords])
        self._labels = labels
        self._symbols = symbols
        self._atomTypes = atomTypes

        # Calculate anything we haven't been given
        self.fillData()

        # If under PBC we need to change how we calculate the bonds
        dim = None
        if self._cellParameters and self.static:
            dim = np.array([self._cellParameters['A'], self._cellParameters['B'], self._cellParameters['C']])

        # Specify internal bonds - bond margin probably too big...
        logger.debug("Calculating bonds for fragmentType: {0}".format(self.fragmentType))
        self._bonds = xyz_util.calcBonds(self._coords,
                                         atomTypes,
                                         dim=dim,
                                         maxAtomRadius=self.maxAtomRadius(),
                                         bondMargin=0.25)

        # Create list of which atoms are bonded to each atom
        self._calcBonded()

        # Set up endGroups
        self.setEndGroups(endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms)

        return

    def setEndGroups(self, endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms):

        # Now set up the endGroup information
        self._endGroups = []
        self._maxBonds = {}
        self._endGroupBonded = {}

        for i, e in enumerate(endGroups):
            eg = EndGroup()
            eg.fragment = self
            eg._endGroupType = endGroupTypes[ i ]
            eg.fragmentEndGroupIdx = e
            eg.fragmentCapIdx = capAtoms[ i ]
            eg.fragmentDihedralIdx = dihedralAtoms[ i ]
            eg.fragmentUwIdx = uwAtoms[ i ]

            eg.capBondLength =  xyz_core.distance(self._coords[eg.fragmentCapIdx], self._coords[eg.fragmentEndGroupIdx ])

            if eg.type() not in self._maxBonds:
                self._maxBonds[ eg.type() ] = None
            if eg.type() not in self._endGroupBonded:
                self._endGroupBonded[ eg.type() ] = 0

            # sanity check
            assert eg.fragmentCapIdx in self._bonded[ eg.fragmentEndGroupIdx ], \
            "capAtom {0} is not bonded to endGroup {1}".format(eg.fragmentCapIdx, eg.fragmentEndGroupIdx)
            if eg.fragmentDihedralIdx != -1:
                assert eg.fragmentDihedralIdx in self._bonded[ eg.fragmentEndGroupIdx ], \
            "dihedral Atom {0} is not bonded to endGroup {1}".format(eg.fragmentDihedralIdx, eg.fragmentEndGroupIdx)
            if eg.fragmentUwIdx != -1:
                assert eg.fragmentUwIdx in self._bonded[ eg.fragmentEndGroupIdx ], \
            "uwAtom {0} is not bonded to endGroup {1}".format(eg.fragmentUwIdx, eg.fragmentEndGroupIdx)

            self._endGroups.append(eg)

        # sanity check - make sure no endGroup is the cap Atom for another endGroup
        eg = set([ e.fragmentEndGroupIdx for e in self._endGroups ])
        caps = set([ e.fragmentCapIdx for e in self._endGroups ])
        assert len(eg.intersection(caps)) == 0, "Cap atom for one endGroup is another endGroup!"

        return

    def setMaxBond(self, endGroupType, maxBond):
        assert endGroupType in self._maxBonds
        self._maxBonds[ endGroupType ] = maxBond
        return

    def symbol(self, idxAtom):
        return self._symbols[self._ext2int[idxAtom]]

    def translate(self, tvector):
        """ translate the molecule by the given vector"""
        self._coords += tvector
        self._changed = True
        return

    def type(self, idxAtom):
        # Fix for old pkl files - _atomTypes used to be called _types
        if hasattr(self,'_atomTypes'):
            return self._atomTypes[self._ext2int[idxAtom]]
        else:
            return self._types[self._ext2int[idxAtom]]

    def update(self):
        """Update the mapping between the internal and external indices"""
        self._int2ext.clear()
        self._ext2int.clear()
        ecount = 0
        for i in range(len(self._coords)):
            if not self.masked[i]:
                self._ext2int[ecount] = i
                self._int2ext[i] = ecount
                ecount += 1
        self.config = [True if eg.free() else False for eg in self._endGroups]
        self.configStr = configManager.updateConfig(self)
        return

    def __str__(self):
        """List the data attributes of this object"""
#         me = {}
#         for slot in dir(self):
#             attr = getattr(self, slot)
#             if not slot.startswith("__") and not (isinstance(attr, types.MethodType) or
#               isinstance(attr, types.FunctionType)):
#                 me[slot] = attr
        me = "Fragment {0}\n".format(self.fragmentType)
        for i, c in enumerate(self._coords):
            me += "{0}  {1:5} {2:0< 15} {3:0< 15} {4:0< 15} \n".format(i, self._labels[i], c[0], c[1], c[2])
        return "{0} : {1}".format(self.__repr__(), me)

