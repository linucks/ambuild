import logging
import math

from ambuild.ab_ffield import FFIELD
from ambuild.xyz_core import wrapCoord3

logger = logging.getLogger(__name__)
PARAM_SEP = ","


class DLPOLY(FFIELD):
    """
    write out CONFIG file
    write out FIELD file
    """

    def _writeCONFIG(
        self, cell, types, coords, fileName="CONFIG",
    ):
        """Write out DLPOLY CONFIG file

        DLPOLY assumes a centered cell.
        """

        with open(fileName, "w") as f:

            # header
            f.write("Ambuild CONFIG file\n")

            levcfg = 0  # just coordinates
            # imcon 1=cubic bounduary conditions, 2=orthorhombic boundary conditions
            imcon = 1 if cell[0] == cell[1] == cell[2] else 2
            f.write("{0:>10}{1:>10}\n".format(levcfg, imcon))

            # Now cell
            f.write("{0: > 20.6F}{1: > 20.6F}{2: > 20.6F}\n".format(cell[0], 0.0, 0.0))
            f.write("{0: > 20.6F}{1: > 20.6F}{2: > 20.6F}\n".format(0.0, cell[1], 0.0))
            f.write("{0: > 20.6F}{1: > 20.6F}{2: > 20.6F}\n".format(0.0, 0.0, cell[2]))

            # Loop through coordinates
            count = 1  # FORTRAN COUNTING
            for i, coord in enumerate(coords):
                # Remove atom index so we can be used by earlier DL-POLY versions
                # f.write("{0}    {1}\n".format(types[i], count))
                f.write("{0}\n".format(types[i]))
                x, y, z = coord
                f.write("{0: > 20.6F}{1: > 20.6F}{2: > 20.6F}\n".format(x, y, z))
                count += 1
        return

    def writeCONTROL(self):
        txt = """simulation NPT
temperature           zero
pressure              0.001
ensemble npt hoover   0.1 1.0
steps                 3
equilibration         2
print                 1
stats                 1
timestep              1.00E-4
cutoff                10.0
delr width            2.00E-01
rvdw cutoff           10.0
ewald precision       1.0E-06
job time              1.0E+04
close time            1.0E+02
finish
"""
        with open("CONTROL", "w") as w:
            w.write(txt)
        return

    def writeFIELDandCONFIG(
        self, cell, rigidBody=True, periodic=True, center=True, skipDihedrals=False
    ):
        angles = []
        angleTypes = []
        bonds = []
        bondTypes = []
        propers = []
        properTypes = []
        impropers = []
        improperTypes = []

        types = []
        coords = []
        images = []
        charges = []
        diameters = []
        masses = []
        symbols = []
        bodies = []
        bodies2 = []
        frozen = []

        atomCount = 0  # Global count in cell
        bodyCount = -1
        for block in cell.blocks.values():
            blockBonds = []
            blockBondTypes = []
            blockAngles = []
            blockAngleTypes = []
            blockPropers = []
            blockProperTypes = []
            blockImpropers = []
            blockImproperTypes = []
            # Do bonds first as counting starts from atomCount and goes up through the blocks
            # It turns out the counting is within a molecule, not global so need to fix.
            if not rigidBody:
                # add all bonds, angles and dihederals throughout the whole block
                # Add all bonds
                # blockBonds += [ (a1+atomCount, a2+atomCount) for a1, a2 in block.bonds() ]
                blockBonds += block.bonds()
                blockBondTypes += [
                    (block.type(a1), block.type(a2)) for a1, a2 in block.bonds()
                ]
                _angles, _propers, _impropers = block.anglesAndDihedrals()
                # Add all angles
                # blockAngles += [ (a1+atomCount, a2+atomCount, a3+atomCount) for a1, a2, a3 in _angles ]
                blockAngles += _angles
                blockAngleTypes += [
                    (block.type(a1), block.type(a2), block.type(a3))
                    for a1, a2, a3 in _angles
                ]
                # Add all propers
                # blockPropers += [ (a1+atomCount, a2+atomCount, a3+atomCount, a4+atomCount) \
                #                 for a1, a2, a3, a4 in _propers ]
                blockPropers += _propers
                blockProperTypes += [
                    (block.type(a1), block.type(a2), block.type(a3), block.type(a4))
                    for a1, a2, a3, a4 in _propers
                ]
                # Add all impropers
                # blockImpropers += [ (a1+atomCount, a2+atomCount, a3+atomCount, a4+atomCount) \
                #                   for a1, a2, a3, a4 in _impropers ]
                blockImpropers += _impropers
                blockImproperTypes += [
                    (block.type(a1), block.type(a2), block.type(a3), block.type(a4))
                    for a1, a2, a3, a4 in _impropers
                ]
            else:
                # Just add the bonds between blocks. Also add angles for all atoms connected to the bonds
                # we do this so that we can exclude them from VdW interactions in MD codes
                for b1, b2 in block.blockBonds():
                    # The bonds themselves
                    # blockBonds.append( (b1+atomCount,b2+atomCount) )
                    blockBonds.append((b1, b2))
                    blockBondTypes.append((block.type(b1), block.type(b2)))
                    _angles = set()
                    # Atoms connected to the endGroup that we need to specify as connected so we add as angles
                    for batom in block.atomBonded1(b1):
                        # The opposite endGroup is included in the list bonded to an endGroup so skip
                        if batom == b2:
                            continue
                        _angles.add((batom, b1, b2))
                    for batom in block.atomBonded1(b2):
                        # The opposite endGroup is included in the list bonded to an endGroup so skip
                        if batom == b1:
                            continue
                        _angles.add((b1, b2, batom))
                    # Add to overall lists
                    for a1, a2, a3 in _angles:
                        # blockAngles.append( (a1+atomCount, a2+atomCount, a3+atomCount))
                        blockAngles.append((a1, a2, a3))
                        blockAngleTypes.append(
                            (block.type(a1), block.type(a2), block.type(a3))
                        )
                    # Dihedrals
                    for dindices in block.dihedrals(b1, b2):
                        dihedral = (dindices[0], dindices[1], dindices[2], dindices[3])
                        blockPropers.append(dihedral)
                        blockProperTypes.append(
                            (
                                block.type(dindices[0]),
                                block.type(dindices[1]),
                                block.type(dindices[2]),
                                block.type(dindices[3]),
                            )
                        )
            blockTypes = []
            blockCoords = []
            blockImages = []
            blockCharges = []
            blockDiameters = []
            blockMasses = []
            blockSymbols = []
            blockBodies = []
            blockBodies2 = []
            blockFrozen = []
            # Now loop through fragments and coordinates
            if hasattr(block, "_fragments"):
                fragments = block._fragments
            else:
                fragments = block.fragments
            # for idxFrag,frag in enumerate(block.fragments): # need index of fragment in block
            molCount = 0
            for frag in fragments:
                # Body count always increments with fragment although it may go up within a fragment too
                bodyCount += 1
                lastBody = frag.body(0)
                bstart = molCount
                for i, coord in enumerate(frag.iterCoord()):
                    if periodic:
                        c, ic = wrapCoord3(coord, cell.dim, center=center)
                        blockCoords.append(c)
                        blockImages.append(ic)
                    else:
                        blockCoords.append(coord)
                        blockImages.append([0, 0, 0])
                    blockTypes.append(frag.type(i))
                    blockCharges.append(frag.charge(i))
                    blockDiameters.append(0.1)
                    blockMasses.append(frag.mass(i))
                    blockSymbols.append(frag.symbol(i))
                    # Work out which body this is in
                    b = frag.body(i)
                    if b != lastBody:
                        bodyCount += 1
                        lastBody = b
                        # blockBodies2.append((bstart,atomCount))
                        # bstart=atomCount+1
                        blockBodies2.append((bstart, molCount))
                        bstart = molCount + 1
                    blockBodies.append(bodyCount)
                    if (
                        frag.static and not rigidBody
                    ):  # Bit of a hack - can't have frozen atoms in rigid bodies
                        blockFrozen.append(1)
                    else:
                        blockFrozen.append(0)
                    # Increment global atom count
                    atomCount += 1  # NOT USED AS COUNTING IS INTER_MOLECULAR!
                    molCount += 1
                # Work out which fragment this is in
                # REM blockCount is NOT idxBlock in the dict - need to rationalise this.
                # d.tagIndices.append((idxBlock,idxFrag,atomCount-i,atomCount))
                # blockBodies2.append((bstart,atomCount))
                blockBodies2.append((bstart, molCount))
            # END loop through fragments
            # Now have all data pertaining to a particular block
            bonds.append(blockBonds)
            bondTypes.append(blockBondTypes)
            angles.append(blockAngles)
            angleTypes.append(blockAngleTypes)
            propers.append(blockPropers)
            properTypes.append(blockProperTypes)
            impropers.append(blockImpropers)
            improperTypes.append(blockImproperTypes)

            types.append(blockTypes)
            coords.append(blockCoords)
            images.append(blockImages)
            charges.append(blockCharges)
            diameters.append(blockDiameters)
            masses.append(blockMasses)
            symbols.append(blockSymbols)
            bodies.append(blockBodies)
            bodies2.append(blockBodies2)
            frozen.append(blockFrozen)

        # End block loop
        # First write out the CONFIG file
        # Unpack the coordinates and types from the blocks
        self._writeCONFIG(
            cell.dim, [j for i in types for j in i], [j for i in coords for j in i]
        )
        # Quick hack hijacking hoomdblue machinary
        # Check we have all the parameters we need
        self.bonds = set(["{0}-{1}".format(j[0], j[1]) for i in bondTypes for j in i])
        self.angles = set(
            ["{0}-{1}-{2}".format(j[0], j[1], j[2]) for i in angleTypes for j in i]
        )
        self.dihedrals = set(
            [
                "{0}-{1}-{2}-{3}".format(j[0], j[1], j[2], j[3])
                for i in properTypes
                for j in i
            ]
        )
        self.impropers = set(
            [
                "{0}-{1}-{2}-{3}".format(j[0], j[1], j[2], j[3])
                for i in improperTypes
                for j in i
            ]
        )
        self.atomTypes = set([j for i in types for j in i])
        # Pierre wants us to write things out even if there are missing dihedral parameters, but we need to know
        # if there are any valid parameters as that determines whether to add the relevant section
        self.checkParameters(skipDihedrals=skipDihedrals)
        # Now write out FIELD file
        # REM DLPOLY does FORTRAN counting so add 1 to everything

        # Each are organised in lists by molecule
        numMolecules = len(coords)
        with open("FIELD", "w") as f:
            # Header
            f.write("Ambuild FIELD file with {0} molecules\n".format(numMolecules))
            f.write("UNITS kcal\n")
            f.write("MOLECULES {0}\n".format(numMolecules))
            _types = set()
            for i in range(numMolecules):
                f.write("Molecule #{0}\n".format(i))
                f.write("NUMMOLS 1\n")
                f.write("ATOMS {0}\n".format(len(coords[i])))
                # Loop over each atom and add to set and write out at end - clumsy
                for j in range(len(coords[i])):
                    t = types[i][j]
                    _types.add(t)
                    f.write(
                        "{0:6}  {1:6}  {2:6}    1    {3}\n".format(
                            t, masses[i][j], charges[i][j], frozen[i][j]
                        )
                    )
                # Rigid bodies
                if rigidBody:
                    nb = len(bodies2[i])
                    f.write("RIGID {0}\n".format(nb))
                    for bstart, bend in bodies2[i]:
                        nsites = bend - bstart
                        s = "{0}    ".format(nsites)
                        # First line is length and up to 15 entries
                        for b in range(min(15, nsites)):
                            s += "{0}    ".format(bstart + b + 1)
                        s += "\n"
                        # Remaining lines can take 16 per line
                        if nsites > 15:
                            for b in range(nsites - 15):
                                if b % 16 == 0 and b > 0:
                                    s += "\n{0}    ".format(bstart + b + 15 + 1)
                                else:
                                    s += "{0}    ".format(bstart + b + 15 + 1)
                            s += "\n"
                        f.write(s)
                # Bonds
                if len(bonds[i]):
                    f.write("BONDS {0}\n".format(len(bonds[i])))
                    for j in range(len(bonds[i])):
                        b1, b2 = sorted(bondTypes[i][j])
                        b = "{0}-{1}".format(b1, b2)
                        param = self.ffield.bondParameter(b)
                        f.write(
                            "harm    {0}    {1}    {2}    {3}\n".format(
                                bonds[i][j][0] + 1,
                                bonds[i][j][1] + 1,
                                param["k"],
                                param["r0"],
                            )
                        )
                # Angles
                if len(angles[i]):
                    f.write("ANGLES {0}\n".format(len(angles[i])))
                    for j in range(len(angles[i])):
                        # a1,a2,a3 = sorted( angleTypes[i][j] )
                        a1, a2, a3 = angleTypes[i][j]
                        a = "{0}-{1}-{2}".format(a1, a2, a3)
                        param = self.ffield.angleParameter(a)
                        f.write(
                            "harm    {0}    {1}    {2}    {3}    {4}\n".format(
                                angles[i][j][0] + 1,
                                angles[i][j][1] + 1,
                                angles[i][j][2] + 1,
                                param["k"],
                                math.degrees(param["t0"]),
                            )
                        )
                if len(propers[i]):
                    # Check if any are valid
                    ok_propers = []
                    for j in range(len(propers[i])):
                        d1, d2, d3, d4 = properTypes[i][j]
                        d = "{0}-{1}-{2}-{3}".format(d1, d2, d3, d4)
                        if self.ffield.hasDihedral(d):
                            ok_propers.append(True)
                        else:
                            ok_propers.append(False)

                    if any(ok_propers):
                        f.write("DIHEDRALS {0}\n".format(sum(ok_propers)))
                        for j in range(len(propers[i])):
                            if not ok_propers[j]:
                                continue
                            # d1,d2,d3,d4 = sorted(properTypes[i][j] )
                            d1, d2, d3, d4 = properTypes[i][j]
                            d = "{0}-{1}-{2}-{3}".format(d1, d2, d3, d4)
                            param = self.ffield.dihedralParameter(d)
                            # A delta m - DLPOLY
                            # k d  n - hoomd
                            # d parameter should be 0 or 180
                            if param["d"] == -1:
                                dp = 180
                            elif param["d"] == 1:
                                dp = 0
                            f.write(
                                "cos  {0:6}  {1:6}  {2:6}  {3:6}  {4:6}  {5:6} {6:6}\n".format(
                                    propers[i][j][0] + 1,
                                    propers[i][j][1] + 1,
                                    propers[i][j][2] + 1,
                                    propers[i][j][3] + 1,
                                    param["k"] / 2,
                                    dp,
                                    param["n"],
                                )
                            )
                f.write("FINISH\n")
            # End of MOLECULE loop so write out non-bonded parameters
            _types = sorted(_types)
            p = []
            for i, atype in enumerate(_types):
                for j, btype in enumerate(_types):
                    if j >= i:
                        param = self.ffield.pairParameter(atype, btype)
                        p.append((atype, btype, param["epsilon"], param["sigma"]))
            f.write("VDW    {0}\n".format(len(p)))
            for a, b, e, s in p:
                f.write("{0:8} {1:6}  lj  {2:6}  {3:6}\n".format(a, b, e, s))
            f.write("CLOSE\n")
        return
