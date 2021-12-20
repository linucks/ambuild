"""
Created on Jan 15, 2013

@author: abbietrewin
"""
import collections
import copy
import csv
import logging
import os
import warnings

import numpy as np

# our imports
from ambuild import ab_body
from ambuild import ab_endgroup
from ambuild import xyz_core
from ambuild import xyz_util

logger = logging.getLogger()


def fragmentFactory(
    fragmentType,
    filePath=None,
    catalyst=False,
    markBonded=False,
    solvent=False,
    static=False,
):
    """Dynamically create a fragment object.

    Dynamic classes allow us to have fragments where all fragments of the same type share
    class variables such as symbols that are the same between all fragments (to save memory),
    but have instance variables for holding attributes like coordinates.

    We need to define a __reduce__ method as otherwise the dynamic classes can't be pickled:
    https://stackoverflow.com/questions/11658511/pickling-dynamically-generated-classes
    """
    # Make sure fragments don't pollute the namespace - probably needs more thinking about
    assert not fragmentType in [
        "fragmentFactory",
        "Fragment",
    ], f"Unacceptable fragmentType: {fragmentType}"
    fobj = type(fragmentType, (Fragment,), dict())
    return fobj(
        filePath=filePath,
        fragmentType=fragmentType,
        solvent=solvent,
        markBonded=markBonded,
        catalyst=catalyst,
        static=static,
    )


class Fragment(object):
    """
    classdocs
    """

    #### These class variables are shared by all fragments of a particular type
    ## Public variables
    catalyst = None
    fragmentType = None
    markBonded = False
    onbondFunction = None  # A function to be called when we bond an endGroup
    solvent = (
        None  # True if this fragment is solvent and should be excluded from clashChecks
    )
    static = False

    ## Private variables
    _atomTypes = []
    _bodies = []  # a list of which body within this fragment each atom belongs to
    _bonds = []  # List of internal fragment bonds
    _bonded = []  # List of which atoms are bonded to which
    _cellParameters = {}
    _charges = []
    _coords = []
    _labels = []
    _masses = []
    _radii = []
    _maxBonds = {}
    _radius = -1
    _symbols = []  # ordered array of symbols (in upper case)
    _totalMass = -1

    def __init__(
        self,
        filePath=None,
        fragmentType=None,
        solvent=False,
        static=False,
        markBonded=False,
        catalyst=False,
    ):
        """
        Constructor
        """
        # The variables below here are specific to a fragment and change as the fragment moves
        # and is involved in bonds; each fragment gets its own copy of these
        self.block = None
        self._coords = []
        self._centroid = None
        self._centerOfMass = None
        self._ext2int = collections.OrderedDict()
        self._int2ext = collections.OrderedDict()
        self._centerOfMass = None
        self._maxAtomRadius = -1
        self._changed = (
            True  # Flag for when we've been moved and need to recalculate things
        )
        self.blockIdx = None  # The index in the list of block data where the data for this fragment starts
        self._endGroups = []  # A list of the endGroup objects
        self._endGroupBonded = (
            []
        )  # A list of the number of each endGroup that are used in bonds
        self.masked = []  # bool - whether the atoms is hidden (e.g. cap or uw atom)
        self.unBonded = []  # bool - whether an atom has just been unbonded

        # Init variables
        self.catalyst = catalyst
        self.fragmentType = fragmentType
        self.markBonded = markBonded
        self.solvent = solvent
        self.static = static

        # Create from the file
        if filePath:
            self.fromFile(filePath)
        return

    def _markBonded(self, endGroup, bond):
        """Append * to all endGroups in fragments that are involved in bonds to cat"""
        # if not endGroup.type() == 'cat:a': return
        if endGroup.type().endswith(ab_endgroup.ENDGROUPBONDED) or not (
            bond.endGroup1.fragment.catalyst or bond.endGroup2.fragment.catalyst
        ):
            return
        logger.debug(
            "_markBonded marking bonds for fragment {0}".format(self.fragmentType)
        )
        for eg in self.endGroups():
            assert not eg._endGroupType.endswith(
                ab_endgroup.ENDGROUPBONDED
            ), "Already got bonded endGroup"
            eg._endGroupType += ab_endgroup.ENDGROUPBONDED
        return

    def addBond(self, endGroup, bond):
        endGroupType = endGroup.type()

        # Mask fragment cap and uw atoms now
        self.masked[endGroup.fragmentCapIdx] = True
        if endGroup.fragmentUwIdx != -1:
            self.masked[endGroup.fragmentUwIdx] = True

        # Hack for starred endGroups
        if not endGroupType.endswith(ab_endgroup.ENDGROUPBONDED):
            self._endGroupBonded[endGroupType] += 1
            # Handle maxBonds here
            if self._maxBonds[endGroupType] is not None:
                if self._endGroupBonded[endGroupType] >= self._maxBonds[endGroupType]:
                    # We have exceeded the bonding limit for these endGroupTypes, so we set any free ones
                    # of this type to blocked
                    for eg in self._endGroups:
                        if not eg.bonded and eg.type() == endGroupType:
                            eg.blocked = True

        # The user may have supplied a custom bonding function, so we call that here
        if self.onbondFunction:
            self.onbondFunction(endGroup)

        if hasattr(self, "markBonded") and self.markBonded:
            self._markBonded(endGroup, bond)
        self.update()
        return

    def delBond(self, endGroupType):
        self._endGroupBonded[endGroupType] -= 1
        return

    def bodies(self):
        for bodyIdx in set(self._bodies):
            yield ab_body.Body(self, bodyIdx)

    def body(self, idxAtom):
        return self._bodies[self._ext2int[idxAtom]]

    def bonds(self):
        """Return a list of bonds
        We exclude masked atoms"""
        bonds = []
        for b1, b2 in self._bonds:
            if not (self.masked[b1] or self.masked[b2]):
                bonds.append((self._int2ext[b1], self._int2ext[b2]))
        return bonds

    def _calcBonded(self):
        if not len(self._bonds):
            warnings.warn("No bonds found in fragment %s" % self.fragmentType)
            return
        # Create empty lists for all
        self._bonded = [[] for _ in range(len(self._coords))]
        for b1, b2 in self._bonds:
            if b1 not in self._bonded[b2]:
                self._bonded[b2].append(b1)
            if b2 not in self._bonded[b1]:
                self._bonded[b1].append(b2)
        return

    def _calcConfigStr(self):
        """Specifies the type of this fragment based on its fragmentType and which endGroups are bonded"""
        return self.fragmentType + "".join(["1" if c else "0" for c in self.config])

    def _calcCenters(self):
        """Calculate the center of mass and geometry for this fragment"""
        self._centroid = np.sum(self._coords, axis=0) / np.size(self._coords, axis=0)
        self._totalMass = np.sum(self._masses)
        # Centre of mass is sum of the products of the individual coordinates multiplied by the mass, divded by the total mass
        self._centerOfMass = xyz_core.centreOfMass(self._coords, self._masses)
        return

    def _calcRadius(self):
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
        distances = [xyz_core.distance(self._centroid, coord) for coord in self._coords]
        imax = np.argmax(distances)
        dist = distances[imax]
        self._radius = (
            dist + self.maxAtomRadius()
        )  # Add on the radius of the largest atom
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
        f = copy.deepcopy(self)
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
        """Fill the data arrays from the label"""

        self.masked = np.array([False] * len(self._coords))
        self.unBonded = [False] * len(self._coords)
        self._masses = np.array(
            [xyz_core.ATOMIC_MASS[symbol] for symbol in self._symbols]
        )
        self._totalMass = np.sum(self._masses)
        self._radii = np.array(
            [
                xyz_core.COVALENT_RADII[xyz_core.SYMBOL_TO_NUMBER[s.upper()]]
                * xyz_core.BOHR2ANGSTROM
                for s in self._symbols
            ]
        )
        self._maxAtomRadius = np.max(self._radii)
        return

    # Can't use as in some cases the class name doesn't seem to get saved on pickling and
    # the fragment has the class 'Fragment'
    # @property
    # def fragmentType(self):
    #     return str(self.__class__).split(".")[-1].rstrip("'>")

    def freeEndGroups(self):
        return [eg for eg in self._endGroups if eg.free()]

    @staticmethod
    def fromCarFile(carFile):
        """ "Abbie did this."""
        labels = []
        symbols = []
        atomTypes = []
        charges = []
        coords = []
        box = None
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
                box = [float(fields[1]), float(fields[2]), float(fields[3])]
            count = 0
            while reading:
                line = f.readline().strip()
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
        return (coords, labels, symbols, atomTypes, charges, box)

    @staticmethod
    def fromXyzFile(xyzFile):
        """ "Jens did this."""
        labels = []
        symbols = []
        atomTypes = []
        charges = []
        coords = []
        box = None
        with open(xyzFile) as f:
            # First line is number of atoms
            line = f.readline()
            natoms = int(line.strip())
            # Skip title
            line = f.readline()
            count = 0
            for _ in range(natoms):
                line = f.readline().strip()
                fields = line.split()
                label = fields[0]
                labels.append(label)
                coords.append(np.array(fields[1:4], dtype=np.float64))
                symbol = xyz_util.label2symbol(label)
                symbols.append(symbol)
                atomTypes.append(symbols)
                charges.append(0.0)
                count += 1
        return (coords, labels, symbols, atomTypes, charges, box)

    def fromFile(self, filePath):

        if filePath.endswith(".car"):
            (coords, labels, symbols, atomTypes, charges, box) = self.fromCarFile(
                filePath
            )
        elif filePath.endswith(".xyz"):
            (coords, labels, symbols, atomTypes, charges, box) = self.fromXyzFile(
                filePath
            )
        else:
            raise RuntimeError("Unrecognised file suffix: {}".format(filePath))

        if box:
            self._cellParameters["A"] = box[0]
            self._cellParameters["B"] = box[1]
            self._cellParameters["C"] = box[2]

        # Get cap atoms and endgroups
        (
            endGroupTypes,
            endGroups,
            capAtoms,
            dihedralAtoms,
            uwAtoms,
        ) = self.parseEndgroupFile(filePath)

        # Set the root fragment and its attributes
        self.setData(
            coords=coords,
            labels=labels,
            symbols=symbols,
            atomTypes=atomTypes,
            charges=charges,
            endGroupTypes=endGroupTypes,
            endGroups=endGroups,
            capAtoms=capAtoms,
            dihedralAtoms=dihedralAtoms,
            uwAtoms=uwAtoms,
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
        return sum([nbonded for nbonded in self._endGroupBonded.values()])

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
            logger.critical(
                "No endGroup definition file supplied for file: {0}".format(filePath)
            )
            return endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms
        with open(egfile) as fh:
            csvreader = csv.reader(fh, delimiter=",", quotechar='"')
            for i, row in enumerate(csvreader):
                if i == 0:  # Header
                    if (
                        not len(row) >= 4
                        and row[0].lower() == "type"
                        and row[1].lower() == "endgroup"
                        and row[2].lower() == "capatom"
                        and row[3].lower() == "delatom"
                    ):
                        raise RuntimeError(
                            "First line of csv file must contain header line:\ntype,endgroup,capatom,dihedral,delatom"
                        )
                    continue
                # skip blank lines
                if not len(row):
                    continue
                # For now make sure first value is letter
                assert row[0][
                    0
                ].isalpha(), "First column of ambi file needs to be a letter!"
                endGroupTypes.append(row[0])
                endGroupIdx = int(row[1])
                endGroups.append(endGroupIdx)
                capAtomIdx = int(row[2])
                if capAtomIdx in capAtoms:
                    raise RuntimeError(
                        "multiple endGroups cannot share the same Cap Atom: {}".format(
                            capAtomIdx
                        )
                    )
                capAtoms.append(capAtomIdx)
                if len(row) >= 4 and row[3] and row[3] != -1:
                    dihedralAtoms.append(int(row[3]))
                else:
                    dihedralAtoms.append(-1)
                if len(row) >= 5 and row[4] and row[4] != -1:
                    uwAtoms.append(int(row[4]))
                else:
                    uwAtoms.append(-1)
            if self.fragmentType == "cap" and len(endGroups) != 1:
                raise RuntimeError("Capfile had >1 endGroup specified!")
        return endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms

    def processBodies(self, filepath):
        """See if we split the fragment into bodies or not"""
        assert (
            len(self._coords) > 0
        ), "Coordinates must have been read before processing bodies"
        dirname, filename = os.path.split(filepath)
        basename, suffix = os.path.splitext(filename)
        bodyFile = os.path.join(dirname, basename + ".ambody")
        if os.path.isfile(bodyFile):
            with open(bodyFile) as f:
                self._bodies = np.array([int(l.strip()) for l in f], dtype=int)
            assert len(self._bodies) == len(
                self._coords
            ), "Must have as many bodies as coordinates: {0} - {1}!".format(
                len(self.bodies), self._dataLen
            )
            assert self._bodies[0] == 0, "Bodies must start with zero!"
        else:
            # Just create an array with 0
            self._bodies = np.zeros(len(self._coords), dtype=int)
        return

    def radius(self, idxAtom):
        return self._radii[self._ext2int[idxAtom]]

    def totalRadius(self):

        if self._changed:
            self._calcProperties()
        return self._radius

    def rotate(self, rotationMatrix, center):
        """Rotate the molecule about the given axis by the angle in radians"""
        self._coords = self._coords - center
        # self._coords = np.array([ np.dot(rotationMatrix, c) for c in self._coords ])
        # I don't actually undestand why this works at all...
        # From: http://stackoverflow.com/questions/12148351/efficiently-rotate-a-set-of-points-with-a-rotation-matrix-in-numpy
        self._coords = np.dot(self._coords, rotationMatrix.T)
        self._coords = self._coords + center
        self._changed = True
        return

    def setData(
        self,
        coords=None,
        labels=None,
        symbols=None,
        atomTypes=None,
        charges=None,
        endGroupTypes=None,
        endGroups=None,
        capAtoms=None,
        dihedralAtoms=None,
        uwAtoms=None,
    ):

        self._charges = np.array(charges)
        self._coords = np.array(coords)
        self._labels = labels
        self._symbols = symbols
        self._atomTypes = atomTypes

        self.fillData()  # Calculate anything we haven't been explicitly given
        # If under PBC we need to change how we calculate the bonds
        dim = None
        if self._cellParameters and self.static:
            dim = np.array(
                [
                    self._cellParameters["A"],
                    self._cellParameters["B"],
                    self._cellParameters["C"],
                ]
            )
        # Specify internal bonds - bond margin probably too big...
        logger.debug(
            "Calculating bonds for fragmentType: {0}".format(self.fragmentType)
        )
        self._bonds = xyz_util.calcBonds(
            self._coords,
            atomTypes,
            dim=dim,
            maxAtomRadius=self.maxAtomRadius(),
            bondMargin=0.25,
        )
        self._calcBonded()
        self.setEndGroups(endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms)
        return

    def setEndGroups(self, endGroupTypes, endGroups, capAtoms, dihedralAtoms, uwAtoms):
        self._endGroups = []
        self._maxBonds = {}
        self._endGroupBonded = {}
        for i, e in enumerate(endGroups):
            eg = ab_endgroup.EndGroup()
            eg.fragment = self
            eg._endGroupType = endGroupTypes[i]
            eg.fragmentEndGroupIdx = e
            eg.fragmentCapIdx = capAtoms[i]
            eg.fragmentDihedralIdx = dihedralAtoms[i]
            eg.fragmentUwIdx = uwAtoms[i]
            eg.capBondLength = xyz_core.distance(
                self._coords[eg.fragmentCapIdx], self._coords[eg.fragmentEndGroupIdx]
            )
            if eg.type() not in self._maxBonds:
                self._maxBonds[eg.type()] = None
            if eg.type() not in self._endGroupBonded:
                self._endGroupBonded[eg.type()] = 0

            # sanity check
            assert (
                eg.fragmentCapIdx in self._bonded[eg.fragmentEndGroupIdx]
            ), "capAtom {0} is not bonded to endGroup {1}".format(
                eg.fragmentCapIdx, eg.fragmentEndGroupIdx
            )
            if eg.fragmentDihedralIdx != -1:
                assert (
                    eg.fragmentDihedralIdx in self._bonded[eg.fragmentEndGroupIdx]
                ), "dihedral Atom {0} is not bonded to endGroup {1}".format(
                    eg.fragmentDihedralIdx, eg.fragmentEndGroupIdx
                )
            if eg.fragmentUwIdx != -1:
                assert (
                    eg.fragmentUwIdx in self._bonded[eg.fragmentEndGroupIdx]
                ), "uwAtom {0} is not bonded to endGroup {1}".format(
                    eg.fragmentUwIdx, eg.fragmentEndGroupIdx
                )
            self._endGroups.append(eg)
        # sanity check - make sure no endGroup is the cap Atom for another endGroup
        eg = set([e.fragmentEndGroupIdx for e in self._endGroups])
        caps = set([e.fragmentCapIdx for e in self._endGroups])
        assert (
            len(eg.intersection(caps)) == 0
        ), "Cap atom for one endGroup is another endGroup!"

        self.setupCapTrilateration()
        return

    def setupCapTrilateration(self):
        NUM_ANCHORS = 4
        caps = set([e.fragmentCapIdx for e in self._endGroups])
        # Get 4 atoms evenly spaced across the list of all of them - this
        # should mean we get a good spacing for calculating distances
        indexes = [i for i in range(len(self._coords)) if i not in caps]
        if len(indexes) < NUM_ANCHORS:
            # Cannot do trilateration without 4 non-cap atoms
            return
        triAtoms = np.round(np.linspace(0, len(indexes) - 1, NUM_ANCHORS)).astype(int)
        for eg in self._endGroups:
            eg.triAtoms = triAtoms
            eg.triDistances = [
                xyz_core.distance(self._coords[eg.fragmentCapIdx], self._coords[a])
                for a in triAtoms
            ]
            # Check trilateration works for this set of atoms
            check_pos = xyz_util.trilaterate3D(
                eg.triDistances, [self._coords[p] for p in eg.triAtoms]
            )
            if not np.allclose(check_pos, self._coords[eg.fragmentCapIdx]):
                raise RuntimeError(
                    f"Incorrect trilateration check position: {check_pos} {self._coords[eg.fragmentCapIdx]}"
                )
        return

    def setMaxBond(self, endGroupType, maxBond):
        assert endGroupType in self._maxBonds
        self._maxBonds[endGroupType] = maxBond
        return

    def symbol(self, idxAtom):
        return self._symbols[self._ext2int[idxAtom]]

    def translate(self, tvector):
        """translate the molecule by the given vector"""
        self._coords += tvector
        self._changed = True
        return

    def type(self, idxAtom):
        # Fix for old pkl files - _atomTypes used to be called _types
        if hasattr(self, "_atomTypes"):
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
        self.config = [False if eg.free() else True for eg in self._endGroups]
        self.configStr = self._calcConfigStr()

    def updateCharges(self, charges):
        """Update the charges on all atoms"""
        if not len(charges) == len(self._charges):
            raise RuntimeError(
                "Charges arrays are of different lengths: %s vs %s"
                % (self._charges, charges)
            )
        self._charges = np.array(charges)

    def __reduce__(self):
        """Required to allow dynamically created classes to be pickled.
        See: https://docs.python.org/3/library/pickle.html#object.__reduce__
        """
        state = self.__dict__.copy()
        return (
            fragmentFactory,
            (self.fragmentType,),
            state,
        )

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
            me += "{0}  {1:5} {2:0< 15} {3:0< 15} {4:0< 15} \n".format(
                i, self._labels[i], c[0], c[1], c[2]
            )
        return "{0} : {1}".format(self.__repr__(), me)
