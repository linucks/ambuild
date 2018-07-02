#!/usr/bin/env python

import logging
import numpy as np
import math

logger = logging.getLogger()


# degrees to radians
BOHR2ANGSTROM = 0.529177249

# double check these values...
# hvd values obtained from http://www.webelements.com/ and recorded to their
#    known accuracy.
ATOMIC_MASS = {
   'H'  :   1.00794,
   'He' :   4.002602,
   'HE' :   4.002602,
   'Li' :   6.941,
   'LI' :   6.941,
   'Be' :   9.012182,
   'BE' :   9.012182,
   'B'  :  10.811,
   'C'  :  12.0107,
   'N'  :  14.0067,
   'O'  :  15.9994,
   'F'  :  18.9984032,
   'Ne' :  20.1797,
   'NE' :  20.1797,
   'Na' :  22.989770,
   'NA' :  22.989770,
   'Mg' :  24.3050,
   'MG' :  24.3050,
   'Al' :  26.981538,
   'AL' :  26.981538,
   'Si' :  28.0855,
   'SI' :  28.0855,
   'P'  :  30.973761,
   'S'  :  32.065,
   'Cl' :  35.453,
   'CL' :  35.453,
   'Ar' :  39.948,
   'AR' :  39.948,
   'K'  :  39.0983,
   'Ca' :  40.078,
   'CA' :  40.078,
   'Sc' :  44.955910,
   'SC' :  44.955910,
   'Ti' :  47.867,
   'TI' :  47.867,
   'V'  :  50.9415,
   'Cr' :  51.9961,
   'CR' :  51.9961,
   'Mn' :  54.938049,
   'MN' :  54.938049,
   'Fe' :  55.845,
   'FE' :  55.845,
   'Co' :  58.933200,
   'CO' :  58.933200,
   'Ni' :  58.6934,
   'NI' :  58.6934,
   'Cu' :  63.546,
   'CU' :  63.546,
   'Zn' :  65.39,
   'ZN' :  65.39,
   'Ga' :  69.723,
   'GA' :  69.723,
   'Ge' :  72.64,
   'GE' :  72.64,
   'As' :  74.92160,
   'AS' :  74.92160,
   'Se' :  78.96,
   'SE' :  78.96,
   'Br' :  79.904,
   'BR' :  79.904,
   'Kr' :  83.80,
   'KR' :  83.80,
   'Rb' :  85.4678,
   'RB' :  85.4678,
   'Sr' :  87.62,
   'SR' :  87.62,
   'Y'  :  88.90585,
   'Zr' :  91.224,
   'ZR' :  91.224,
   'Nb' :  92.90638,
   'NB' :  92.90638,
   'Mo' :  95.94,
   'MO' :  95.94,
   'Tc' :  98,
   'TC' :  98,
   'Ru' : 101.07,
   'RU' : 101.07,
   'Rh' : 102.90550,
   'RH' : 102.90550,
   'Pd' : 106.42,
   'PD' : 106.42,
   'Ag' : 107.8682,
   'AG' : 107.8682,
   'Cd' : 112.411,
   'CD' : 112.411,
   'In' : 114.818,
   'IN' : 114.818,
   'Sn' : 118.710,
   'SN' : 118.710,
   'Sb' : 121.760,
   'SB' : 121.760,
   'Te' : 127.60,
   'TE' : 127.60,
   'I'  : 126.90447,
   'Xe' : 131.293,
   'XE' : 131.293,
   'Cs' : 132.90545,
   'CS' : 132.90545,
   'Ba' : 137.327,
   'BA' : 137.327,
   'La' : 138.9055,
   'LA' : 138.9055,
   'Ce' : 140.116,
   'CE' : 140.116,
   'Pr' : 140.90765,
   'PR' : 140.90765,
   'Nd' : 144.24,
   'ND' : 144.24,
   'Pm' : 145,
   'PM' : 145,
   'Sm' : 150.36,
   'SM' : 150.36,
   'Eu' : 151.964,
   'EU' : 151.964,
   'Gd' : 157.25,
   'GD' : 157.25,
   'Tb' : 158.92534,
   'TB' : 158.92534,
   'Dy' : 162.50,
   'DY' : 162.50,
   'Ho' : 164.93032,
   'HO' : 164.93032,
   'Er' : 167.259,
   'ER' : 167.259,
   'Tm' : 168.93421,
   'TM' : 168.93421,
   'Yb' : 173.04,
   'YB' : 173.04,
   'Lu' : 174.967,
   'LU' : 174.967,
   'Hf' : 178.49,
   'HF' : 178.49,
   'Ta' : 180.9479,
   'TA' : 180.9479,
   'W'  : 183.84,
   'Re' : 186.207,
   'RE' : 186.207,
   'Os' : 190.23,
   'OS' : 190.23,
   'Ir' : 192.217,
   'IR' : 192.217,
   'Pt' : 195.078,
   'PT' : 195.078,
   'Au' : 196.96655,
   'AU' : 196.96655,
   'Hg' : 200.59,
   'HG' : 200.59,
   'Tl' : 204.3833,
   'TL' : 204.3833,
   'Pb' : 207.2,
   'PB' : 207.2,
   'Bi' : 208.98038,
   'BI' : 208.98038,
   'Po' : 208.98,
   'PO' : 208.98,
   'At' : 209.99,
   'AT' : 209.99,
   'Rn' : 222.02,
   'RN' : 222.02,
   'Fr' : 223.02,
   'FR' : 223.02,
   'Ra' : 226.03,
   'RA' : 226.03,
   'Ac' : 227.03,
   'AC' : 227.03,
   'Th' : 232.0381,
   'TH' : 232.0381,
   'Pa' : 231.03588,
   'PA' : 231.03588,
   'U'  : 238.02891,
   'Np' : 237.05,
   'NP' : 237.05,
   'Pu' : 244.06,
   'PU' : 244.06,
   'Am' : 243.06,
   'AM' : 243.06,
   'Cm' : 247.07,
   'CM' : 247.07,
   'Bk' : 247.07,
   'BK' : 247.07,
   'Cf' : 251.08,
   'CF' : 251.08,
   'Es' : 252.08,
   'ES' : 252.08,
   'Fm' : 257.10,
   'FM' : 257.10,
   'Md' : 258.10,
   'MD' : 258.10,
   'No' : 259.10,
   'NO' : 259.10,
   'Lr' : 262.11,
   'LR' : 262.11,
   'Rf' : 261.11,
   'RF' : 261.11,
   'Db' : 262.11,
   'DB' : 262.11,
   'Sg' : 266.12,
   'SG' : 266.12,
   'Bh' : 264.12,
   'BH' : 264.12,
   'Hs' : 269.13,
   'HS' : 269.13,
   'Mt' : 268.14,
   'MT' : 268.14,
   # jmht - dummy atom
   'X' : 0,
}

# mapping symbols to Z
SYMBOL_TO_NUMBER = {
'H'  : 1,
'HE' : 2,
'LI' : 3,
'BE' : 4,
'B'  : 5,
'C'  : 6,
'N'  : 7,
'O'  : 8,
'F'  : 9,
'NE' : 10,
'NA' : 11,
'MG' : 12,
'AL' : 13,
'SI' : 14,
'P'  : 15,
'S'  : 16,
'CL' : 17,
'AR' : 18,
'K'  : 19,
'CA' : 20,
'SC' : 21,
'TI' : 22,
'V'  : 23,
'CR' : 24,
'MN' : 25,
'FE' : 26,
'CO' : 27,
'NI' : 28,
'CU' : 29,
'ZN' : 30,
'GA' : 31,
'GE' : 32,
'AS' : 33,
'SE' : 34,
'BR' : 35,
'KR' : 36,
'RB' : 37,
'SR' : 38,
'Y'  : 39,
'ZR' : 40,
'NB' : 41,
'MO' : 42,
'TC' : 43,
'RU' : 44,
'RH' : 45,
'PD' : 46,
'AG' : 47,
'CD' : 48,
'IN' : 49,
'SN' : 50,
'SB' : 51,
'TE' : 52,
'I'  : 53,
'XE' : 54,
'CS' : 55,
'BA' : 56,
'LA' : 57,
'CE' : 58,
'PR' : 59,
'ND' : 60,
'PM' : 61,
'SM' : 62,
'EU' : 63,
'AD' : 64,
'TB' : 65,
'DY' : 66,
'HO' : 67,
'ER' : 68,
'TM' : 69,
'YB' : 70,
'LU' : 71,
'HF' : 72,
'TA' : 73,
'W'  : 74,
'RE' : 75,
'OS' : 76,
'IR' : 77,
'PT' : 78,
'AU' : 79,
'HG' : 80,
'TL' : 81,
'PB' : 82,
'BI' : 83,
'PO' : 84,
'AT' : 85,
'RN' : 86,
'FR' : 87,
'RA' : 88,
'AC' : 89,
'TH' : 90,
'PA' : 91,
'U'  : 92,
'NP' : 93,
'PU' : 94,
'AM' : 95,
'CM' : 96,
'BK' : 97,
'CF' : 98,
'ES' : 99,
'FM' : 100,
'MD' : 101,
'NO' : 102,
'LR' : 103,
'UNA' : 104,
'UNP' : 105,
'X'  : 0,
}


# the display table (aus)
# H=.7 is just for appearance
#
#
# note that indexing with -1 accesses the last element (0 copied at start and edn)
# UNITS ARE IN BOHR!!!
COVALENT_RADII = [
0,
0.7, 3.80,
2.76, 1.99, 1.62, 1.33, 1.23, 1.14, 0.95, 3.80,
3.42, 2.85, 2.38, 2.09, 1.90, 1.90, 1.90, 3.80,
4.18, 3.42,
     3.04, 2.66, 2.57, 2.66, 2.66, 2.66, 2.57, 2.57, 2.57, 2.57,
                          2.47, 2.38, 2.18, 2.18, 2.18, 3.80,
4.46, 3.80,
     3.42, 2.94, 2.76, 2.76, 2.57, 2.47, 2.57, 2.66, 3.04, 2.94,
                          2.94, 2.76, 2.76, 2.66, 2.66, 3.80,
4.94, 4.09,
     3.71,
     3.52, 3.52, 3.52, 3.52, 3.52, 3.52, 3.42, 3.33, 3.33, 3.33, 3.33, 3.33, 3.33, 3.33,
          2.94, 2.76, 2.57, 2.57, 2.47, 2.57, 2.57, 2.57, 2.85,
                         3.61, 3.42, 3.04, 3.61, 3.61, 3.80,
4.94, 4.09,
     3.71,
     3.42, 3.42, 3.33, 3.33, 3.33, 3.33, 3.23, 3.13, 3.13, 3.13, 3.13, 3.13, 3.13, 3.13, 1., 1., 1., 1., 1., 1., 0, ]


# joe lennards table angstroms
# 1.0 added as index [0] and [-1] see above
# UNITS ARE IN BOHR!!!
VDW_RADII = [  1.0,
  1.20, 1.40,
  1.82, 1.78, 1.74, 1.70, 1.55, 1.52, 1.47, 1.54,
  2.27, 2.22, 2.16, 2.10, 1.80, 1.80, 1.75, 1.88,
  2.75, 2.57,
       2.56, 2.54, 2.52, 2.50, 2.48, 2.46, 2.44, 2.42, 2.41, 2.40,
            2.40, 2.10, 1.85, 1.90, 1.85, 2.02,
  3.10, 2.80,
       2.77, 2.74, 2.71, 2.68, 2.65, 2.62, 2.59, 2.56, 2.53, 2.51,
            2.50, 2.20, 2.10, 2.06, 1.98, 2.16, 1., 1., 1.]


"""
# Characteristic lengths of single bonds.
# Reference: CRC Handbook of Chemistry and Physics, 87th edition, (2006), Sec. 9 p. 46
   As   Br   C    Cl   F    Ge   H    I    N    O    P    S    Sb   Se   Si
As 2.10
Br 2.32 2.28
C  1.96 1.94 1.53
Cl 2.17 2.14 1.79 1.99
F  1.71 1.76 1.39 1.63 1.41
Ge    2.30 1.95 2.15 1.73 2.40
H  1.51 1.41 1.09 1.28 0.92 1.53 0.74
I    2.47 2.13 2.32 1.91 2.51 1.61 2.67
N         1.46 1.90 1.37     1.02       1.45
O         1.42 1.70 1.42     0.96       1.43 1.48
P       2.22 1.85 2.04 1.57     1.42       1.65      2.25
S       2.24 1.82 2.05 1.56     1.34            1.69      2.00
Sb              2.33           1.70
Se           1.95      1.71      1.47                               2.33
Si      2.21 1.87 2.05 1.58      1.48 2.44      1.63      2.14           2.33
Sn           2.14 2.28           1.71 2.67
Te                     1.82      1.66
"""
# REM - symbols should be in lower case!
# UNITS ARE IN ANGSTROM!!!
ELEMENT_TYPE_BOND_LENGTHS = {}
ELEMENT_TYPE_BOND_LENGTHS['AS'] = { 'AS' : 2.10,
                                   'BR' : 2.32,
                                   'C'  : 1.96,
                                   'CL' : 2.17,
                                   'F'  : 1.71,
                                   'H'  : 1.51 }

ELEMENT_TYPE_BOND_LENGTHS['BR'] = { 'BR': 2.28,
                                   'C' : 1.94,
                                   'CL': 2.14,
                                   'F' : 1.76,
                                   'GE': 2.30,
                                   'H' : 1.41,
                                   'I' : 2.47,
                                   'P' : 2.22,
                                   'S' : 2.24,
                                   'SI' : 2.21 }

ELEMENT_TYPE_BOND_LENGTHS['C'] = { 'C' : 1.53,
                                  'CL': 1.79,
                                  'F' : 1.39,
                                  'GE': 1.95,
                                  'H' : 1.09,
                                  'I' : 2.13,
                                  'O' : 1.42,
                                  # 'N' : 1.46,
                                  'N' : 1.4,  # jmht changed for abbie
                                  'P' : 1.85,
                                  'S' : 1.82,
                                  'SE': 1.95,
                                  'SI': 1.87,
				  'NI': 2.75,
                                  'SN': 2.14,
                                  'Cu': 1.98 }

ELEMENT_TYPE_BOND_LENGTHS['CL'] = { 'CL' : 1.99,
                                   'F'  : 1.63,
                                   'GE' : 2.15,
                                   'H'  : 1.28,
                                   'I'  : 2.32,
                                   'N'  : 1.90,
                                   'O'  : 1.70,
                                   'P'  : 2.04,
                                   'S'  : 2.05,
                                   'SB' : 2.33,
                                   'SI' : 2.05,
                                   'SN' : 2.28 }

ELEMENT_TYPE_BOND_LENGTHS['F'] = { 'F'  : 1.41,
                                  'GE' : 1.73,
                                  'H'  : 0.92,
                                  'I'  : 1.91,
                                  'N'  : 1.37,
                                  'O'  : 1.42,
#                                  'P'  : 1.57, # changed for abbie
                                  'P'  : 1.63,
                                  'S'  : 1.56,
                                  'SE' : 1.71,
                                  'SI' : 1.58,
                                  'TE' : 1.82 }

ELEMENT_TYPE_BOND_LENGTHS['GE'] = { 'GE' : 2.40,
                                   'H'  : 1.53,
                                   'I'  : 2.51 }

ELEMENT_TYPE_BOND_LENGTHS['H'] = { 'H' : 0.74,
                                  'I' : 1.61,
                                  'N' : 1.02,
                                  'O' : 0.96,
                                  'P' : 1.42,
                                  'S' : 1.34,
                                  'SB': 1.70,
                                  'SE': 1.47,
                                  'SI': 1.48,
                                  'SN': 1.71,
				  'NI': 2.75,
                                  'TE': 1.66 }

ELEMENT_TYPE_BOND_LENGTHS['I'] = { 'I'  : 2.67,
                                   'SI' : 2.44,
                                   'SN' : 2.67 }

ELEMENT_TYPE_BOND_LENGTHS['N'] = { 'N' : 1.45,
                                   'O' : 1.43,
                                   'P' : 1.65,
                                   'ZN' : 2.166,  # Added for abbie
                                   }

ELEMENT_TYPE_BOND_LENGTHS['O'] = { 'O'  : 1.48,
                                   'SI' : 1.63 }

ELEMENT_TYPE_BOND_LENGTHS['P'] = { 'P' : 2.25 }

ELEMENT_TYPE_BOND_LENGTHS['S'] = { 'S' : 2.00,
                                   'O' : 1.69 }

ELEMENT_TYPE_BOND_LENGTHS['SE'] = { 'SE' : 2.33 }

ELEMENT_TYPE_BOND_LENGTHS['SI'] = { 'SI' : 2.33 }
ELEMENT_TYPE_BOND_LENGTHS['LI'] = { 'LI' : 2.67 }

ELEMENT_TYPE_BOND_LENGTHS['Cu'] = { 'C' : 1.98 }


def angle(c1, c2, c3, dim=None, pbc=None):
    """Return the angle in radians c1---c2---c3
    where c are the coordinates in a numpy array
    """
    r1 = distance(c2, c1, dim=dim, pbc=pbc)
    r2 = distance(c3, c2, dim=dim, pbc=pbc)
    r3 = distance(c3, c1, dim=dim, pbc=pbc)
    x = (r1 * r1 + r2 * r2 - r3 * r3) / (2.0 * r1 * r2)
    assert not np.isnan(x)
    # print "r1: {0}, r2: {1}, r3: {2}, x: {3}".format( r1, r2, r3, x )
    if np.allclose(x, 1.0):
        theta = 0.0
    elif np.allclose(x, -1.0):
        theta = math.pi
    else:
        theta = np.arccos(x)
    return theta


def centroid(coords):
    return np.sum(coords, axis=0) / np.size(coords, axis=0)


def centreOfMass(coords, masses):
    return np.sum(coords * masses[:,np.newaxis], axis=0) / np.sum(masses)


def dihedral(p1, p2, p3, p4, dim=None, pbc=None):
    """Dihedral angle between 4 points"""
    vec_ij = vecDiff(p1, p2, dim=dim, pbc=pbc)
    vec_kj = vecDiff(p3, p2, dim=dim, pbc=pbc)
    vec_kl = vecDiff(p3, p4, dim=dim, pbc=pbc)
    # vec1 is the normal to the plane defined by atoms i, j, and k
    vec1 = np.cross(vec_ij, vec_kj)
    magvec1 = np.dot(vec1, vec1)
    #  vec2 is the normal to the plane defined by atoms j, k, and l
    vec2 = np.cross(vec_kl, vec_kj)
    magvec2 = np.dot(vec2, vec2)
    # the definition of a dot product is used to find the angle between
    # vec1 and vec2 and hence the angle between the planes defined by
    # atoms i, j, k and j, k, l
    #
    # the factor of pi (180.0) is present since when we defined the
    # vectors vec1 and vec2, one used the right hand rule while the
    # other used the left hand rule
    dotprod = np.dot(vec1, vec2)
    # print type(magvec1), type(magvec2)
    fac = dotprod / math.sqrt(magvec1 * magvec2)
    # We get a nan when all the atoms are in line - for these we return zero
    if np.isnan(fac):
        #raise RuntimeError("Nan in dihedral!")
        return 0.0
    fac = min(fac, 1.0)
    fac = max(fac, -1.0)
    # dihed = 180.0 - math.degrees( math.acos(fac ) )
    dihed = math.pi - math.acos(fac)
    # the dot product between the bond between atoms i and j and the
    # normal to the plane defined by atoms j, k, and l is used to
    # determine whether or not the dihedral angle is clockwise or
    # anti_clockwise
    #
    # if the dot product is positive, the rotation is clockwise
    sign_check = np.dot(vec_ij, vec2)
    if(sign_check > 0.0):
        dihed = dihed * -1.0
    return dihed


def distance(v1, v2, dim=None, pbc=[True,True,True]):
    """Distance with numpy taking PBC into account
    This works either with 2 points or a vector of any number of points
    Adapted from: http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co
    Changed so that it can cope with distances across more than one cell
    """
    return np.sqrt((vecDiff(v1, v2, dim=dim, pbc=pbc) ** 2).sum(axis=-1))


def momentOfInertia(coords, masses):
    """Moment of Inertia Tensor"""
    totalMass = np.sum(masses)
    centreOfMass = np.sum(coords * masses[:,np.newaxis], axis=0) / totalMass
    coords = coords - centreOfMass # Coords relative to centre of mass
    if False:
        # Below from Pierre but doesn't seem to replicate results of doing it 'manually'
        I = np.dot(coords.transpose(), coords)
    else:
        x = 0
        y = 1
        z = 2
        I = np.zeros(shape=(3, 3))
        I[x, x] = np.sum((np.square(coords[:, y]) + np.square(coords[:, z])) * masses)
        I[y, y] = np.sum((np.square(coords[:, x]) + np.square(coords[:, z])) * masses)
        I[z, z] = np.sum((np.square(coords[:, x]) + np.square(coords[:, y])) * masses)
        I[x, y] = np.sum(coords[:, x] * coords[:, y] * masses)
        I[y, x] = I[x, y]
        I[y, z] = np.sum(coords[:, y] * coords[:, z] * masses)
        I[z, y] = I[y, z]
        I[x, z] = np.sum(coords[:, x] * coords[:, z] * masses)
        I[z, x] = I[x, z]
    return I


def principalMoments(coords, masses):
    """http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node67.html
    """
    I = momentOfInertia(coords, masses)
    eigval, eigvec = np.linalg.eig(I)
    return np.sort(eigval)

def rigid_rotate(A, B):
    """Return rotation matrix to rotate A to B. Assumes both sets of points are already centred.
    
    Uses SVD to calculate the rotation.
    
    Taken from: http://nghiaho.com/?page_id=671"""
    A = np.matrix(A)
    B = np.matrix(B)
    assert len(A) == len(B)
    # dot is matrix multiplication for array
    H = A.T * B
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T * U.T
    # special reflection case
    if np.linalg.det(R) < 0:
        # Reflection detected
        Vt[2,:] *= -1
        R = Vt.T * U.T
    return R


def rotation_matrix(axis, angle):
    """
    Return the rotation matrix to rotate a vector by the given angle about the
    axis.

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """

    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2)
    b, c, d = -axis * np.sin(angle / 2)
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])


def rotation_matrix_to_quaternion(M):
    """Return the quaternion of the rotation matrix
    
    https://math.stackexchange.com/questions/893984/conversion-of-rotation-matrix-to-quaternion
    http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    """
    tr = np.trace(M)
    if (tr > 0):
        S = np.sqrt(tr + 1.0) * 2 # S=4*qw 
        qw = 0.25 * S
        qx = (M[2, 1] - M[1, 2]) / S
        qy = (M[0, 2] - M[2, 0]) / S
        qz = (M[1, 0] - M[0, 1]) / S
    elif ((M[0, 0] > M[1, 1])&(M[0, 0] > M[2, 2])):
        S = np.sqrt(1.0 + M[0, 0] - M[1, 1] - M[2, 2]) * 2 # S=4*qx 
        qw = (M[2, 1] - M[1, 2]) / S
        qx = 0.25 * S
        qy = (M[0, 1] + M[1, 0]) / S 
        qz = (M[0, 2] + M[2, 0]) / S 
    elif (M[1, 1] > M[2, 2]):
        S = np.sqrt(1.0 + M[1, 1] - M[0, 0] - M[2, 2]) * 2 # S=4*qy
        qw = (M[0, 2] - M[2, 0]) / S
        qx = (M[0, 1] + M[1, 0]) / S 
        qy = 0.25 * S
        qz = (M[1, 2] + M[2, 1]) / S 
    else:
        S = np.sqrt(1.0 + M[2, 2] - M[0, 0] - M[1, 1]) * 2 # S=4*qz
        qw = (M[1, 0] - M[0, 1]) / S
        qx = (M[0, 2] + M[2, 0]) / S
        qy = (M[1, 2] + M[2, 1]) / S
        qz = 0.25 * S
    return np.array([qw, qx, qy, qz])
    

def vecDiff(v1, v2, dim=None, pbc=[True,True,True]):
    """Difference between vectors with numpy taking PBC into account
    This works either with 2 points or a vector of any number of points
    Adapted from: http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co
    Changed so that it can cope with distances across more than one cell
   Args:
   dim - 3 element array with dimensions of unit cell
   pbc - 3 element boolean array indicating if this dimension has periodic boundaries
    """
    # Currently (e.g. cell.dataDict return a list of coordinates rather than a numpy array, so we need to check if we have been given a list
    # or numpy array and convert accordingly
    if type(v1) is list:
        delta = np.array(v1) - np.array(v2)
    else:
        delta = v1 - v2
    if dim is not None:
        assert type(dim) is np.ndarray
        
        # Account for multiple cells
        if delta.ndim == 1:
            if pbc[0]: delta[0] = np.remainder(delta[0], dim[0])
            if pbc[1]: delta[1] = np.remainder(delta[1], dim[1])
            if pbc[2]: delta[2] = np.remainder(delta[2], dim[2])
        else:
            if pbc[0]: delta[:,0] = np.remainder(delta[:,0], dim[0])
            if pbc[1]: delta[:,1] = np.remainder(delta[:,1], dim[1])
            if pbc[2]: delta[:,2] = np.remainder(delta[:,2], dim[2])
        # Could use below as we don't really need a true cell in that direction - change if slow 
        #delta = np.remainder(delta, dim)
            
        # Wherever the distance is > half the cell length we subtract the cell length
        # we multiply by the pbc array to only make change where we have pbc
        delta = np.where(np.abs(delta) > 0.5 * dim,
                            delta - (np.copysign(dim, delta) * pbc),
                            delta)
    return delta


def vectorAngle(v1, v2):
    """ Calculate the angle between two vectors
    Return value in Radians
    A . B = |A|*|B|*cos(theta)
    so: theta = arccos( X.Y / |X||Y| )
    Stolen from: http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """

    # print "v1: {}".format(v1)
    # print "v1 norm: {}".format(np.linalg.norm(v1))
    v1_u = v1 / np.linalg.norm(v1)
    # print "v1_u: {}".format(v1_u)
    v2_u = v2 / np.linalg.norm(v2)

    angle = np.arccos(np.dot(v1_u, v2_u))
    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle


def unWrapCoord3(coord, image, ldim, centered=False, inplace=False):
    """Unwrap a coordinate back into a cell
    """
    if not inplace:
        coord = np.copy(coord)
    if centered:
        # Put it back with origin at corner
        coord += ldim / 2
        assert np.all((coord >= 0.0) & (coord <= ldim)), "Coords out of cell: {} : {}".format(coord, ldim)
    coord = coord + ldim * image # Move according to image
    return coord


def wrapCoord3(coord, dim, center=False, inplace=False):
    """Wrap coodinate triple into a cell of length dim
    return the wrapped coordinates and the images
    """
    if not inplace:
        coord = np.copy(coord)
    image = np.floor(coord/dim).astype(int)
    if np.any(image < 0):
        coord += -image * dim # Make positive so the modulo operator works
    wcoord = np.fmod(coord, dim)
    assert np.all(wcoord >= 0.0), "-ve Coord!: {0} -> {1} : {2}".format(coord, wcoord, image)
    # Change the coord so the origin is at the center of the box (we start from the corner)
    if center:
        wcoord -= dim / 2
    return wcoord, image
