#!/usr/bin/env python
import collections
import csv
import logging
import math
import os

logger = logging.getLogger(__name__)
PARAM_SEP = ","
ANGLE_SEP = "-"

# Needs to be defined here to avoid circular import in util.py
def read_bond_params(bond_file):
    if not os.path.isfile(bond_file):
        raise RuntimeError("Cannot find bond parameter file: {0}".format(bond_file))
    BondParams = collections.namedtuple("BondParams", ["A", "B", "k", "r0", "comments"])
    params = []
    with open(bond_file) as fh:
        for i, bp in enumerate(
            map(BondParams._make, csv.reader(fh, delimiter=PARAM_SEP, quotechar='"'))
        ):
            if i == 0:  # check header ok
                if bp.A != "A" and bp.B != "B" and bp.k != "k" and bp.r0 != "r0":
                    raise RuntimeError(
                        "Incorrect header for bond param file: {0}".format(bond_file)
                    )
                continue
            params.append(bp)
    return params


class FfieldParameters(object):
    def __init__(self, paramsDir):
        # All parameters calculated to fit PCFF adapted forcefield Holden et al j.phys.chem.c 2012
        # from combining rules calculated using dlpoly-prep (Willock)
        # bonds adapted from quartic PCFF
        self.bonds = self._processBonds(
            os.path.abspath(os.path.join(paramsDir, "bond_params.csv"))
        )
        self.angles = self._processAngles(
            os.path.abspath(os.path.join(paramsDir, "angle_params.csv"))
        )
        self.dihedrals = self._processDihedrals(
            os.path.abspath(os.path.join(paramsDir, "dihedral_params.csv"))
        )
        self.impropers = self._processImpropers(
            os.path.abspath(os.path.join(paramsDir, "improper_params.csv"))
        )
        self.pairs = self._readPairs(
            os.path.abspath(os.path.join(paramsDir, "pair_params.csv"))
        )
        self.paramsDir = os.path.abspath(paramsDir)
        return

    def _processBonds(self, bond_file):
        bond_params = read_bond_params(bond_file)
        if not len(bond_params):
            raise RuntimeError(
                "Could not read any parameters from bond file: {0}".format(bond_file)
            )
        bonds = {}
        for p in bond_params:
            try:
                bonds[(p.A, p.B)] = {"k": float(p.k), "r0": float(p.r0)}
            except Exception as e:
                logger.critical(
                    "Error reading angle parameters from file: {0}".format(bond_file)
                )
                logger.critical("Error occured with parameter {0}".format(p))
                raise e
        return bonds

    def _processAngles(self, angle_file):
        # REM: angles are stored in degrees in the parameter file, but we convert to radians for our uses
        if not os.path.isfile(angle_file):
            raise RuntimeError(
                "Cannot find angle parameter file: {0}".format(angle_file)
            )
        header = ["angle", "k", "t0", "comments"]
        angles = {}
        with open(angle_file) as f:
            reader = csv.reader(f, delimiter=PARAM_SEP, quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header:
                        raise RuntimeError(
                            "Header for {0} file, should be: {1}".format(
                                angle_file, ",".join(header)
                            )
                        )
                    else:
                        continue
                if row[0].startswith("#"):
                    continue
                try:
                    angle = row[0]
                    k = float(row[1])
                    t0 = math.radians(float(row[2]))
                    angles[angle] = {"k": k, "t0": t0}
                    # Angles are symmetrical so we also need to add the opposite notation
                    angle_reversed = ANGLE_SEP.join(angle.split(ANGLE_SEP)[::-1])
                    angles[angle_reversed] = {"k": k, "t0": t0}
                except Exception as e:
                    logger.critical(
                        "Error reading angle parameters from file: {0}".format(
                            angle_file
                        )
                    )
                    logger.critical(
                        "Error occured on line {0}: {1}".format(
                            i + 1, PARAM_SEP.join(row)
                        )
                    )
                    raise e
        return angles

    def _processDihedrals(self, dihedral_file):
        if not os.path.isfile(dihedral_file):
            raise RuntimeError(
                "Cannot find dihedral parameter file: {0}".format(dihedral_file)
            )
        header = ["dihedral", "k", "d", "n", "comments"]
        dihedrals = {}
        with open(dihedral_file) as f:
            reader = csv.reader(f, delimiter=PARAM_SEP, quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header:
                        raise RuntimeError(
                            "Header for {0} file, should be: {1}".format(
                                dihedral_file, ",".join(header)
                            )
                        )
                    else:
                        continue
                if row[0].startswith("#"):
                    continue
                try:
                    dihedral = row[0]
                    k = float(row[1])
                    d = int(row[2])
                    n = int(row[3])
                    dihedrals[dihedral] = {"k": k, "d": d, "n": n}
                except Exception as e:
                    logger.critical(
                        "Error reading dihedral parameters from file: {0}".format(
                            dihedral_file
                        )
                    )
                    logger.critical(
                        "Error occured on line {0}: {1}".format(
                            i + 1, PARAM_SEP.join(row)
                        )
                    )
                    raise e
        return dihedrals

    def _processImpropers(self, improper_file):
        if not os.path.isfile(improper_file):
            raise RuntimeError(
                "Cannot find dihedral parameter file: {0}".format(improper_file)
            )
        header = ["improper", "k", "chi", "comments"]
        impropers = {}
        with open(improper_file) as f:
            reader = csv.reader(f, delimiter=PARAM_SEP, quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header:
                        raise RuntimeError(
                            "Header for {0} file, should be: {1}".format(
                                improper_file, ",".join(header)
                            )
                        )
                    else:
                        continue
                if row[0].startswith("#"):
                    continue
                try:
                    improper = row[0]
                    k = float(row[1])
                    chi = float(row[2])
                    impropers[improper] = {"k": k, "chi": chi}
                except Exception as e:
                    logger.critical(
                        "Error reading improper parameters from file: {0}".format(
                            improper_file
                        )
                    )
                    logger.critical(
                        "Error occured on line {0}: {1}".format(
                            i + 1, PARAM_SEP.join(row)
                        )
                    )
                    raise e
        return impropers

    def _readPairs(self, pairs_file):
        if not os.path.isfile(pairs_file):
            raise RuntimeError(
                "Cannot find pair parameter file: {0}".format(pairs_file)
            )
        header = ["atom1", "atom2", "epsilon", "sigma", "comments"]
        pairs = {}
        with open(pairs_file) as f:
            reader = csv.reader(f, delimiter=PARAM_SEP, quotechar='"')
            for i, row in enumerate(reader):
                if i == 0:
                    if row != header:
                        raise RuntimeError(
                            "Header for {0} file, should be: {1}".format(
                                pairs_file, ",".join(header)
                            )
                        )
                    else:
                        continue
                if row[0].startswith("#"):
                    continue
                try:
                    atom1 = row[0]
                    atom2 = row[1]
                    epsilon = float(row[2])
                    sigma = float(row[3])
                    pairs[(atom1, atom2)] = {"epsilon": epsilon, "sigma": sigma}
                except Exception as e:
                    logger.critical(
                        "Error reading pair parameters from file: {0}".format(
                            pairs_file
                        )
                    )
                    logger.critical(
                        "Error occured on line {0}: {1}".format(
                            i + 1, PARAM_SEP.join(row)
                        )
                    )
                    raise e
        return pairs

    def angleParameter(self, angle):
        return self.angles[angle]

    def hasAngle(self, angle):
        return angle in self.angles.keys()

    def bondParameter(self, bond):
        a, b = bond.split("-")
        if (a, b) in self.bonds.keys():
            return self.bonds[(a, b)]
        elif (b, a) in self.bonds.keys():
            return self.bonds[(b, a)]
        else:
            raise RuntimeError("Missing bond parameter for: {0} {1}".format(a, b))

    def hasBond(self, bond):
        a, b = bond.split("-")
        return (a, b) in self.bonds.keys() or (b, a) in self.bonds.keys()

    def dihedralParameter(self, dihedral):
        return self.dihedrals[dihedral]

    def hasDihedral(self, dihedral):
        return dihedral in self.dihedrals.keys()

    def improperParameter(self, improper):
        return self.impropers[improper]

    def hasImproper(self, improper):
        return improper in self.impropers.keys()

    def pairParameter(self, p1, p2):
        """Return whichever pair is defined"""
        if (p1, p2) in self.pairs:
            return self.pairs[(p1, p2)]
        if (p2, p1) in self.pairs:
            return self.pairs[(p2, p1)]

        assert False, "Could not find {0} {1}".format(p1, p2)
        return

    def hasPair(self, p1, p2):
        """dummy atoms treated specially"""
        if p1.lower() == "x" or p2.lower() == "x":
            return True
        if (p1, p2) in self.pairs or (p2, p1) in self.pairs:
            return True
        return False


class FFIELD(object):
    def __init__(self, paramsDir):
        self.ffield = FfieldParameters(paramsDir)
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        self.impropers = None
        self.atomTypes = None
        return

    def checkParameters(self, skipDihedrals=False):

        assert self.ffield
        assert self.atomTypes

        ok = True
        missingBonds = []
        for bond in self.bonds:
            if not self.ffield.hasBond(bond):
                ok = False
                missingBonds.append(bond)
        missingAngles = []
        for angle in self.angles:
            if not self.ffield.hasAngle(angle):
                ok = False
                missingAngles.append(angle)
        missingDihedrals = []
        missingImpropers = []
        if not skipDihedrals:
            for dihedral in self.dihedrals:
                if not self.ffield.hasDihedral(dihedral):
                    ok = False
                    missingDihedrals.append(dihedral)
            for improper in self.impropers:
                if not self.ffield.hasImproper(improper):
                    ok = False
                    missingImpropers.append(improper)

        missingPairs = []
        for i, atype in enumerate(self.atomTypes):
            for j, btype in enumerate(self.atomTypes):
                if j >= i:
                    if not self.ffield.hasPair(atype, btype):
                        ok = False
                        missingPairs.append((atype, btype))

        if not ok:
            msg = "The following parameters could not be found:\n"
            if missingBonds:
                msg += "Bonds: {0}\n".format(missingBonds)
            if missingAngles:
                msg += "Angles: {0}\n".format(missingAngles)
            if missingDihedrals:
                msg += "Dihedrals: {0}\n".format(missingDihedrals)
            if missingImpropers:
                msg += "Impropers: {0}\n".format(missingImpropers)
            if missingPairs:
                msg += "Pairs: {0}\n".format(missingPairs)

            msg += "Please add these to the files in the directory: {0}\n".format(
                self.ffield.paramsDir
            )
            raise RuntimeError(msg)
        return
