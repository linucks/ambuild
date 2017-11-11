#!/usr/bin/env python
'''
Created on Feb 3, 2013

@author: jmht

Utility functions
'''
import cPickle
import logging
import os
import numpy
import math
import sys
import warnings
import xml.etree.ElementTree as ET

from opt import read_bond_params
from paths import PARAMS_DIR

logger = logging.getLogger()

# Bits stolen from the CCP1GUI: http://sourceforge.net/projects/ccp1gui/
# However, I wrote bits of that so I assume its ok

DUMMY_DIAMETER = 0.1

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


class BondLength(object):
    """Class to hold calculate bond lengths.
    
    This is required because we use data stored in the ATOM_TYPE_BOND_LENGTHS dict
    which is calculated from the parameter files, which is read at run time. Therefore
    we initialise this object using a parameter file and set it as a module object for 
    use by the bondLength function
    """
    def __init__(self, bond_param_file):
        self.ATOM_TYPE_BOND_LENGTHS = {}
        for p in read_bond_params(bond_param_file):
            if p.A not in self.ATOM_TYPE_BOND_LENGTHS:
                if p.B in self.ATOM_TYPE_BOND_LENGTHS:
                    self.ATOM_TYPE_BOND_LENGTHS[p.B][p.A] = float(p.r0)
                else:
                    self.ATOM_TYPE_BOND_LENGTHS[p.A] = { p.B : float(p.r0) }
            else:
                self.ATOM_TYPE_BOND_LENGTHS[p.A][p.B] = float(p.r0)
    
    def bondLength(self, atomType1, atomType2):
        """ Get the characteristic lengths of single bonds as defined in:
            Reference: CRC Handbook of Chemistry and Physics, 87th edition, (2006), Sec. 9 p. 46
        """
        # We first see if we can find the bond length in the ATOM_TYPE_BOND_LENGTHS table
        # If not we fall back to using the bonds calculated from element types
        if self.ATOM_TYPE_BOND_LENGTHS.has_key(atomType1) and self.ATOM_TYPE_BOND_LENGTHS[ atomType1 ].has_key(atomType2):
            #print "ATOM TYPE"
            return self.ATOM_TYPE_BOND_LENGTHS[ atomType1 ][ atomType2 ]
        elif self.ATOM_TYPE_BOND_LENGTHS.has_key(atomType2) and self.ATOM_TYPE_BOND_LENGTHS[ atomType2 ].has_key(atomType1):
            #print "ATOM TYPE"
            return self.ATOM_TYPE_BOND_LENGTHS[ atomType2 ][ atomType1 ]
            
        symbol1 = label2symbol(atomType1).upper()
        symbol2 = label2symbol(atomType2).upper()
        if ELEMENT_TYPE_BOND_LENGTHS.has_key(symbol1) and ELEMENT_TYPE_BOND_LENGTHS[ symbol1 ].has_key(symbol2):
            #print "ELEMENT TYPE"
            return ELEMENT_TYPE_BOND_LENGTHS[ symbol1 ][ symbol2 ]
        elif ELEMENT_TYPE_BOND_LENGTHS.has_key(symbol2) and ELEMENT_TYPE_BOND_LENGTHS[ symbol2 ].has_key(symbol1):
            #print "ELEMENT TYPE"
            return ELEMENT_TYPE_BOND_LENGTHS[ symbol2 ][ symbol1 ]
        warnings.warn('No data for bond length for %s-%s' % (atomType1, atomType2))
        return 1.0

# This needs to be set to the bondLength function of the BondLength class 
# See cell.Cell_setUtilBondLength
# We set this to raise an error if unset
def __STOP(x,y): raise NotImplementedError("Need to set the bondLength function - see setModuleBondLength")
bondLength = __STOP
    
def setModuleBondLength(paramFile):
    """Set the bondLength function of the util module"""
    global bondLength
    BL = BondLength(paramFile)
    bondLength = BL.bondLength
    return

def angle(c1, c2, c3, dim=None, pbc=None):
    """Return the angle in radians c1---c2---c3
    where c are the coordinates in a numpy array
    """
    r1 = distance(c2, c1, dim=dim, pbc=pbc)
    r2 = distance(c3, c2, dim=dim, pbc=pbc)
    r3 = distance(c3, c1, dim=dim, pbc=pbc)
    x = (r1 * r1 + r2 * r2 - r3 * r3) / (2.0 * r1 * r2)
    assert not numpy.isnan(x)
    # print "r1: {0}, r2: {1}, r3: {2}, x: {3}".format( r1, r2, r3, x )
    if numpy.allclose(x, 1.0):
        theta = 0.0
    elif numpy.allclose(x, -1.0):
        theta = math.pi
    else:
        theta = numpy.arccos(x)
    return theta

def calcBondsHACK(coords, symbols, maxAtomRadius=None, bondMargin=0.2, boxMargin=1.0):
    """HACK FOR NETWORK"""

    # print "s ",symbols
    bonds = []
    for i, coord1 in enumerate(coords):
        symbol1 = symbols[ i ]
        for j, coord2 in enumerate(coords):
            if j < i:
                continue
            symbol2 = symbols[ j ]
            dist = distance(coord1, coord2)
            if symbol1 == 'j' and symbol2 == 'C':
                bond_length = 4.31
            elif symbol2 == 'j' and symbol1 == 'C':
                bond_length = 4.31
            else:
                bond_length = bondLength(symbol1, symbol2)
            # print "GOT ",symbol1,symbol2,bond_length, dist
            if  bond_length - bondMargin < dist < bond_length + bondMargin:
                # print "BONDING"
                bonds.append((i, j))
    return bonds

def calcBonds(coords, symbols, dim=None, maxAtomRadius=None, bondMargin=0.2, boxMargin=1.0):
    """Calculate the bonds for the fragments. This is done at the start when the only coordinates
    are those in the fragment.
    symbols can be chemical elements or atomTypes
    If supplied cell is a list/numpy array with the dimensions of the simulation cell, in which case
    PBC will be applied
    """
    close = closeAtoms(coords, symbols, dim, maxAtomRadius, boxMargin)
    v1 = []
    v2 = []
    for idxAtom1, idxAtom2 in close:
        v1.append(coords[idxAtom1])
        v2.append(coords[idxAtom2])
    
    distances = distance(numpy.array(v1), numpy.array(v2), dim=dim)
    bonds = []
    
    for i, (idxAtom1, idxAtom2) in enumerate(close):
        bond_length = bondLength(symbols[idxAtom1], symbols[idxAtom2])
        if bond_length < 0: continue
        logger.debug("Dist:length {1}({0})-{3}({2}) {4} {5}".format( symbols[idxAtom1], idxAtom1, symbols[idxAtom2], idxAtom2, distances[i], bond_length ))
        if  bond_length - bondMargin < distances[i] < bond_length + bondMargin:
            bonds.append((idxAtom1, idxAtom2))
    return bonds

def closeAtoms(coords, symbols, dim=None, maxAtomRadius=None, boxMargin=1.0):
    """Return a list of which atoms in the cell are close to each other.
    Close is defined by the maxAtomRadius and boxMargin"""
    if maxAtomRadius is None:
        for s in set(symbols):
            maxAtomRadius = max(COVALENT_RADII[SYMBOL_TO_NUMBER[s.upper()]] * BOHR2ANGSTROM, maxAtomRadius)

    boxSize = (maxAtomRadius * 2) + boxMargin
    # If we are under PBC calculate the number of boxes in each dimenson
    boxNum = None
    if dim is not None:
        if type(dim) is list:
            dim = numpy.array(dim)
        boxNum = [ int(math.ceil(dim[0] / boxSize)),
                 int(math.ceil(dim[1] / boxSize)),
                 int(math.ceil(dim[2] / boxSize))]

    atomCells = []  # List of which cell each atom is in - matches coords array
    # For cells and hCells, the key is a triple of the indices of the cell position (a,b,c)
    cells = {}  # Dictionary of the cells, each containing a list of atoms in that cell
    hCells = {}  # Dictionary keyed by cell with a list of the cells that surround a particular cell

    # Work out which box each atom is in and the surrounding boxes
    for idxAtom1, coord in enumerate(coords):
        key = getCell(coord, boxSize, dim=dim)
        atomCells.append(key)
        if cells.has_key(key):
            cells[ key ].append(idxAtom1)
        else:
            # Add to main list
            cells[ key ] = [ (idxAtom1) ]
            # Map surrounding boxes
            hCells[ key ] = haloCells(key, boxNum=boxNum)

    # Now calculate the bonding
    _closeAtoms = []
    for idxAtom1 in range(len(coords)):
        key = atomCells[ idxAtom1 ]
        # Loop through all cells surrounding this one
        for scell in hCells[ key ]:
            # Check if we have a cell with anything in it
            try: alist = cells[ scell ]
            except KeyError: continue  # Trigger exception to save searching through the keys
            for idxAtom2 in alist:
                if idxAtom2 > idxAtom1:  # Skip atoms we've already processed
                    _closeAtoms.append((idxAtom1, idxAtom2))
    return _closeAtoms

def getCell(coord, boxSize, dim=None, pbc=[True, True, True]):
    """Return the cell that the coord is in under periodic boundaries"""
    # Periodic Boundaries
    if dim is not None:
        x = coord[0] % dim[0] if pbc[0] else coord[0]
        y = coord[1] % dim[1] if pbc[1] else coord[1]
        z = coord[2] % dim[2] if pbc[2] else coord[2]
    else:
        x, y, z = coord[0], coord[1], coord[2]

    # Calculate which cell the atom is in
    a = int(math.floor(x / boxSize))
    b = int(math.floor(y / boxSize))
    c = int(math.floor(z / boxSize))
    return (a, b, c)

def cellFromPickle(pickleFile, paramsDir=None):
    """Recreate a cell from a pickled file and apply any hacks so that we can work with older versions"""
    def fixFragment(fragment):
        # Need to make sure coords and masses are numpy arrays
        if type(fragment._coords) is list:
            fragment._coords = numpy.array(fragment._coords)
            fragment._masses = numpy.array(fragment._masses)
    
        if not hasattr(fragment,'_atomTypes'):
            fragment._atomTypes = fragment._types
            fragment._sharedAttrs['_atomTypes'] = None
        if not hasattr(fragment,'solvent'):
            # Solvent is a new attribute so we set to false
            fragment.solvent = False
            fragment._sharedAttrs['solvent'] = None
        if not hasattr(fragment,'onbondFunction'):
            fragment.onbondFunction = None
            fragment._sharedAttrs['onbondFunction'] = None
        if not hasattr(fragment,'markBonded'):
            fragment.markBonded = None
            fragment._sharedAttrs['markBonded'] = None
        if not hasattr(fragment,'unBonded'):
            fragment.unBonded = [ False ] * len(fragment._coords)
            fragment._individualAttrs['unBonded'] = None
        if not hasattr(fragment,'catalyst'):
            fragment.catalyst = False
            fragment._sharedAttrs['catalyst'] = False
            
        for e in fragment._endGroups:
            # More horrible hacks for old versions
            if hasattr(e,'_isBonded'):
                e.bonded = e._isBonded
            if not hasattr(e,'blocked'): e.blocked = False
        return
        
    with open(pickleFile) as f:
        myCell = cPickle.load(f)
        
    # Need to hack to work with older versions
    if not hasattr(myCell, 'dim'):
        myCell.dim = numpy.array([myCell.A, myCell.B, myCell.C])
        myCell.numBoxes = [myCell.numBoxA, myCell.numBoxB, myCell.numBoxC]
    if not type(myCell.dim) is numpy.ndarray:
        myCell.dim = numpy.array(myCell.dim)
    if not hasattr(myCell, 'pbc'):
        myCell.pbc = [True, True, True]
        myCell.walls = [False, False, False]
    
    # Set parameter Directory
    if paramsDir is None:
        if hasattr(myCell, 'paramsDir'):
            paramsDir = myCell.paramsDir
        else:
            paramsDir = PARAMS_DIR
    if not os.path.isdir(paramsDir):
        raise RuntimeError("Cannot find cell paramsDir: {0}".format(paramsDir))
    logger.info("Getting parameter files from directory: {0}".format(paramsDir))
    setModuleBondLength(os.path.join(paramsDir,'bond_params.csv'))
    
    # Fix all the fragments
    for fragment in myCell._fragmentLibrary.values():
        fixFragment(fragment)
        
    for block in myCell.blocks.values():
        for fragment in block.fragments:
            fixFragment(fragment)
        
    return myCell

def dihedral(p1, p2, p3, p4, dim=None, pbc=None):
    """ From the CCP1GUI
    """
#     if cell is not None:
#         # We need to fix all distances between vectors for PBC
#         dimensions = numpy.array(cell)
#         vec_ij = numpy.remainder(p1 - p2, dimensions)
#         vec_kl = numpy.remainder(p3 - p4, dimensions)
#         vec_kj = numpy.remainder(p3 - p2, dimensions)
#         vec_ij = numpy.where(numpy.abs(vec_ij) > 0.5 * dimensions, vec_ij - numpy.copysign(dimensions, vec_ij), vec_ij)
#         vec_kj = numpy.where(numpy.abs(vec_kj) > 0.5 * dimensions, vec_kj - numpy.copysign(dimensions, vec_kj), vec_kj)
#         vec_kl = numpy.where(numpy.abs(vec_kl) > 0.5 * dimensions, vec_kl - numpy.copysign(dimensions, vec_kl), vec_kl)
#     else:
#         vec_ij = p1 - p2
#         vec_kj = p3 - p2
#         vec_kl = p3 - p4
    
    vec_ij = vecDiff(p1, p2, dim=dim, pbc=pbc)
    vec_kj = vecDiff(p3, p2, dim=dim, pbc=pbc)
    vec_kl = vecDiff(p3, p4, dim=dim, pbc=pbc)

    # vec1 is the normal to the plane defined by atoms i, j, and k
    vec1 = numpy.cross(vec_ij, vec_kj)
    magvec1 = numpy.dot(vec1, vec1)

    #  vec2 is the normal to the plane defined by atoms j, k, and l
    vec2 = numpy.cross(vec_kl, vec_kj)
    magvec2 = numpy.dot(vec2, vec2)

    # the definition of a dot product is used to find the angle between
    # vec1 and vec2 and hence the angle between the planes defined by
    # atoms i, j, k and j, k, l
    #
    # the factor of pi (180.0) is present since when we defined the
    # vectors vec1 and vec2, one used the right hand rule while the
    # other used the left hand rule
    dotprod = numpy.dot(vec1, vec2)
    # print type(magvec1), type(magvec2)
    fac = dotprod / math.sqrt(magvec1 * magvec2)
    # We get a nan when all the atoms are in line - for these we return zero
    if numpy.isnan(fac):
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
    sign_check = numpy.dot(vec_ij, vec2)
    if(sign_check > 0.0):
        dihed = dihed * -1.0
    return dihed

def distance(v1, v2, dim=None, pbc=[True,True,True]):
    """Distance with numpy taking PBC into account
    This works either with 2 points or a vector of any number of points
    Adapted from: http://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co
    Changed so that it can cope with distances across more than one cell
    """
    return numpy.sqrt((vecDiff(v1, v2, dim=dim, pbc=pbc) ** 2).sum(axis=-1))

def XdistanceP(self, v1, v2):
    """
    Calculate the distance between two vectors in the cell
    under periodic boundary conditions - from wikipedia entry
    """

    # return numpy.linalg.norm(v1-v2)
    dx = v2[0] - v1[0]
    if math.fabs(dx) > self.A[0] * 0.5:
        dx = dx - math.copysign(self.A[0], dx)
    dy = v2[1] - v1[1]
    if math.fabs(dy) > self.B[1] * 0.5:
        dy = dy - math.copysign(self.B[1], dy)
    dz = v2[2] - v1[2]
    if math.fabs(dz) > self.C[2] * 0.5:
        dz = dz - math.copysign(self.C[2], dz)

    return math.sqrt(dx * dx + dy * dy + dz * dz)

def dumpPkl(pickleFile, split=None, nonPeriodic=False):

    fpath = os.path.abspath(pickleFile)
    logger.info("Dumping pkl file: {0}".format(fpath))
    dname, fname = os.path.split(fpath)
    prefix = os.path.splitext(fname)[0]

    mycell = cellFromPickle(pickleFile)
    if split == "fragments":
        for t in mycell.fragmentTypes().keys():
            data = mycell.dataDict(fragmentType=t)
            mycell.writeXyz("{0}_{1}_P.xyz".format(prefix, t),
                            data=data,
                            periodic=True)
            mycell.writeCml("{0}_{1}_PV.cml".format(prefix, t),
                            data=data,
                            periodic=True,
                            pruneBonds=True)
    elif split == "blocks":
        periodic = True
        for i, b in enumerate(mycell.blocks.values()):
            # Write out each block to a separate file
            if periodic:
                b.writeXyz("{0}_block{1}.xyz".format(prefix, i),
                           cell=mycell.dim)
                b.writeCml("{0}_block{1}.cml".format(prefix, i),
                           cell=mycell.dim)
            else:
                b.writeXyz("{0}_block{1}.xyz".format(prefix, i))
                b.writeCml("{0}_block{1}.cml".format(prefix, i))
    else:
        if nonPeriodic:
            data = mycell.dataDict(rigidBody=False, periodic=False)
            mycell.writeCml(prefix + "_NP.cml", data,
                            periodic=False, pruneBonds=False)
            mycell.writeXyz(prefix + "_NP.xyz", data=data, periodic=False)
            mycell.writeXyz(prefix + "_NP_types.xyz", data=data, periodic=False, atomTypes=True)
        else:
            data = mycell.dataDict(rigidBody=False)
            mycell.writeXyz(prefix + "_P.xyz", data=data, periodic=True)
            mycell.writeXyz(prefix + "_P_types.xyz", data=data, periodic=True, atomTypes=True)
            # self.writeCar(prefix+"_P.car",data=data,periodic=True)
            mycell.writeCml(prefix + "_PV.cml", data=data, periodic=True, pruneBonds=True)
    return

def dumpDLPOLY(pickleFile, rigidBody=False, skipDihedrals=False):
    fpath = os.path.abspath(pickleFile)
    logger.info("Dumping DLPOLY files from pkl file: {0}".format(fpath))
    mycell = cellFromPickle(pickleFile)

    # Need to do this here or else hoomdblue gets the command line arguments on import of the module
    import opt

    d = opt.DLPOLY()
    d.writeCONTROL()
    d.writeFIELDandCONFIG(mycell, rigidBody=rigidBody, skipDihedrals=skipDihedrals)
    return

def haloCells(key, boxNum=None, pbc=[True, True, True]):
    """Returns the list of cells surrounding a cell
    boxNum is the number of boxes in each dimension of the containing cell if we
    are under PBC
    wall is an array with 0 where there is a wall along that axis
    """
    a, b, c = key
    # cells = set()
    cells = set()
    for  i in [ 0, -1, +1 ]:
        for j in [ 0, -1, +1 ]:
            for k in [ 0, -1, +1 ]:
                ai = a + i
                bj = b + j
                ck = c + k
                if boxNum:
                    # Impose periodic boundaries 
                    if ai < 0 and pbc[0]:
                        ai = boxNum[0] - 1
                    elif ai > boxNum[0] - 1 and pbc[0]:
                        ai = 0
                    if bj < 0 and pbc[1]:
                        bj = boxNum[1] - 1
                    elif bj > boxNum[1] - 1 and pbc[1]:
                        bj = 0
                    if ck < 0 and pbc[2]:
                        ck = boxNum[2] - 1
                    elif ck > boxNum[2] - 1 and pbc[2]:
                        ck = 0
                skey = (ai, bj, ck)
                # print "sKey ({},{},{})->({})".format(a,b,c,skey)
                # cells.add(skey)
                cells.add(skey)

    return list(cells)

def frange(start, stop, step):
    """
    Range function that works with floating points
    http://stackoverflow.com/questions/4189766/python-range-with-step-of-type-float
    """
    while start < stop:
        yield start
        start += step

def label2symbol(name):
    """ Determine the element type of an atom from its name, e.g. Co_2b -> Co
        Returns a capitalised element name
        Originally written by Jens Thomas in the CCP1GUI
    """

    origName = name
    name = name.strip().upper()

    # Determine the element from the first 2 chars of the name
    if len(name) > 2:
        name = name[0:2]

    if len(name) == 2 and name[0].isalpha() and name[1].isalpha():
        # 2 Character name, so see if it matches any 2-character elements
        sym2c = filter(lambda x: len(x) == 2, SYMBOL_TO_NUMBER.keys())
        # HACK: NEED TO REMOVE NP and NB
        sym2c.remove('NP')
        sym2c.remove('NB')
        if name in sym2c:
            return name.capitalize()

    # If it was a valid 2 character symbol we should have picked it up so now only 1 symbol
    name = name[0]
    if not name.isalpha():
        raise RuntimeError, "label2symbol first character of name is not a character: {0}".format(origName)

    # Hack - for x return x
    if name.lower() == 'x': return 'x'

    # Get 1 character element names
    sym1c = filter(lambda x: len(x) == 1 and x != 'X', SYMBOL_TO_NUMBER.keys())

    if name in sym1c: return name.capitalize()

    raise RuntimeError, "label2symbol cannot convert name {0} to symbol!".format(origName)
    return

def newFilename(filename, separator="_"):

    dname, name = os.path.split(filename)

    # Create a new filename using _1 etc
    name, suffix = os.path.splitext(name)

    try:
        basename, num = name.split(separator)
    except ValueError:
        # No separator so assume is an un-numbered file
        return os.path.join(dname, name + separator + "1" + suffix)

    num = int(num) + 1
    return os.path.join(dname, basename + separator + str(num) + suffix)

def pickleObj(obj, fileName):
    """Pickle an object - required as we can't pickle in the cell as otherwise the open filehandle
    is within the cell which is the object we are trying to pickle..."""

    with open(fileName, 'w') as pfile:
        cPickle.dump(obj , pfile)

    return

def readMol2(filename):

    coords = []
    symbols = []
    # bonds=[]
    with open(filename) as f:
        captureAtom = False
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                captureAtom = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                captureAtom = False
                break
                captureBond = False
                continue
            if captureAtom:
                f = line.split()
                symbols.append(label2symbol(f[1]))
                coords.append([float(f[2]), float(f[3]), float(f[4])])
#             if captureBond:
#                 f=line.split()

    return coords, symbols

def rotation_matrix(axis, angle):
    """
    Return the rotation matrix to rotate a vector by the given angle about the
    axis.

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """

    axis = axis / numpy.sqrt(numpy.dot(axis, axis))
    a = numpy.cos(angle / 2)
    b, c, d = -axis * numpy.sin(angle / 2)
    return numpy.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])

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
        delta = numpy.array(v1) - numpy.array(v2)
    else:
        delta = v1 - v2
    if dim is not None:
        assert type(dim) is numpy.ndarray
        
        # Account for multiple cells
        if delta.ndim == 1:
            if pbc[0]: delta[0] = numpy.remainder(delta[0], dim[0])
            if pbc[1]: delta[1] = numpy.remainder(delta[1], dim[1])
            if pbc[2]: delta[2] = numpy.remainder(delta[2], dim[2])
        else:
            if pbc[0]: delta[:,0] = numpy.remainder(delta[:,0], dim[0])
            if pbc[1]: delta[:,1] = numpy.remainder(delta[:,1], dim[1])
            if pbc[2]: delta[:,2] = numpy.remainder(delta[:,2], dim[2])
        # Could use below as we don't really need a true cell in that direction - change if slow 
        #delta = numpy.remainder(delta, dim)
            
        # Wherever the distance is > half the cell length we subtract the cell length
        # we multiply by the pbc array to only make change where we have pbc
        delta = numpy.where(numpy.abs(delta) > 0.5 * dim,
                            delta - (numpy.copysign(dim, delta) * pbc),
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
    # print "v1 norm: {}".format(numpy.linalg.norm(v1))
    v1_u = v1 / numpy.linalg.norm(v1)
    # print "v1_u: {}".format(v1_u)
    v2_u = v2 / numpy.linalg.norm(v2)

    angle = numpy.arccos(numpy.dot(v1_u, v2_u))
    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return numpy.pi
            # return numpy.pi/RADIANS2DEGREES
    # return angle/RADIANS2DEGREES
    return angle

def unWrapCoord(coord, image, ldim, centered=False):
    """Unwrap a coordinate back into a cell
    """
    if centered:
        # Put it back with origin at corner
        coord += ldim / 2
        assert 0.0 <= coord <= ldim, "Coord outside box: {0}".format(coord)
    # Now move it according to its image
    coord += ldim * float(image)
    return coord

def wrapCoord(coord, ldim, center=False):
    """Wrap coodinate into a cell of length ldim
    return the wrapped coordinate and the image index
    """
    image = int(math.floor(coord / ldim))
    # Make the coordinate positive so the math modulo operator works
    if image < 0: coord += -image * ldim

    # Use fmod to avoid overflow problems with python modulo operater - see stackexchange
    wcoord = math.fmod(coord, ldim)

    # Should never be negative
    assert wcoord >= 0.0, "Coord {0} -> {1} : {2}".format(coord, wcoord, image)

    # Change the coord so the origin is at the center of the box (we start from the corner)
    if center: wcoord -= ldim / 2
    
    return wcoord, image

def wrapCoord3(coord, dim, center=False):
    """Wrap coodinate triple into a cell of length dim
    return the wrapped coordinates and the images
    """
    image = numpy.floor(coord/dim).astype(int)
    if numpy.any(image < 0): coord += -image * dim # Make positive so the modulo operator works
    wcoord = numpy.fmod(coord, dim)
    assert numpy.all(wcoord >= 0.0), "-ve Coord!: {0} -> {1} : {2}".format(coord, wcoord, image)
    # Change the coord so the origin is at the center of the box (we start from the corner)
    if center: wcoord -= dim / 2
    return wcoord, image

def writeCml(cmlFilename,
             coords,
             symbols,
             bonds=None,
             atomTypes=None,
             cell=None,
             pruneBonds=False):

    assert len(coords) == len(symbols)
    if pruneBonds:
        assert cell is not None
        pcoords = []  # need to save periodic coordinates

    root = ET.Element('molecule')
    root.attrib['xmlns'] = "http://www.xml-cml.org/schema"
    root.attrib['xmlns:cml'] = "http://www.xml-cml.org/dict/cml"
    root.attrib['xmlns:units'] = "http://www.xml-cml.org/units/units"
    root.attrib['xmlns:xsd'] = "http://www.w3c.org/2001/XMLSchema"
    root.attrib['xmlns:iupac'] = "http://www.iupac.org"
    root.attrib['id'] = "mymolecule"

    if cell is not None:
        # First set up the cell
        crystal = ET.SubElement(root, "crystal")

        crystalANode = ET.SubElement(crystal, "scalar")
        crystalBNode = ET.SubElement(crystal, "scalar")
        crystalCNode = ET.SubElement(crystal, "scalar")
        crystalAlphaNode = ET.SubElement(crystal, "scalar")
        crystalBetaNode = ET.SubElement(crystal, "scalar")
        crystalGammaNode = ET.SubElement(crystal, "scalar")

        crystalANode.attrib["title"] = "a"
        crystalBNode.attrib["title"] = "b"
        crystalCNode.attrib["title"] = "c"
        crystalAlphaNode.attrib["title"] = "alpha"
        crystalBetaNode.attrib["title"] = "beta"
        crystalGammaNode.attrib["title"] = "gamma"

        crystalANode.attrib["units"] = "units:angstrom"
        crystalBNode.attrib["units"] = "units:angstrom"
        crystalCNode.attrib["units"] = "units:angstrom"
        crystalAlphaNode.attrib["units"] = "units:degree"
        crystalBetaNode.attrib["units"] = "units:degree"
        crystalGammaNode.attrib["units"] = "units:degree"

        crystalANode.text = str(cell[0])
        crystalBNode.text = str(cell[1])
        crystalCNode.text = str(cell[2])

        # Only support orthorhombic? cells
        crystalAlphaNode.text = "90"
        crystalBetaNode.text = "90"
        crystalGammaNode.text = "90"

    if atomTypes:
        assert len(atomTypes) == len(coords)
        # Need to collate atomTypes
        for atype in set(atomTypes):
            atomTypeNode = ET.SubElement(root, "atomType")
            atomTypeNode.attrib['name'] = atype
            atomTypeNode.attrib['title'] = atype

    # Now atom data
    atomArrayNode = ET.SubElement(root, "atomArray")
    for i, coord in enumerate(coords):
        atomNode = ET.SubElement(atomArrayNode, "atom")
        atomNode.attrib['id'] = "a{0}".format(i)
        atomNode.attrib['elementType'] = symbols[i]
        if cell is None:
            x = coord[0]
            y = coord[1]
            z = coord[2]
        else:
            x, ix = wrapCoord(coord[0], cell[0], center=False)
            y, iy = wrapCoord(coord[1], cell[1], center=False)
            z, iz = wrapCoord(coord[2], cell[2], center=False)
            if pruneBonds:
                pcoords.append(numpy.array([x, y, z]))

        atomNode.attrib['x3'] = str(x)
        atomNode.attrib['y3'] = str(y)
        atomNode.attrib['z3'] = str(z)

        if atomTypes:
            # Now add atomType as child node referring to the atomType
            atomTypeNode = ET.SubElement(atomNode, "atomType")
            atomTypeNode.attrib['ref'] = atomTypes[ i ]

    if len(bonds):
        # Now do bonds
        if pruneBonds:
            # Hack to get vis working
            # Calculate all bond distances
            distances = distance([ pcoords[b1] for b1, b2 in bonds],
                                  [ pcoords[b2] for b1, b2 in bonds])

            # Complete hack - just see if it's longer then 0.5 the cell A distance - assumes cubic cell
            bonds = [ b for i, b in enumerate(bonds) if distances[i] <= cell[0] * 0.5 ]

        bondArrayNode = ET.SubElement(root, "bondArray")
        for b in bonds:
            bondNode = ET.SubElement(bondArrayNode, "bond")
            bondNode.attrib['atomRefs2'] = "a{0} a{1}".format(b[0], b[1])
            bondNode.attrib['order'] = "1"

    tree = ET.ElementTree(root)
    # ET.dump(tree)

    cmlFilename = os.path.abspath(cmlFilename)
    # tree.write(file_or_filename, encoding, xml_declaration, default_namespace, method)
    tree.write(cmlFilename, encoding="utf-8", xml_declaration=True)

    return cmlFilename

def writeXyz(fileName, coords, symbols, cell=None):
    """Write out an xyz file to fileName
    if cell is a list of 3 floats, write out a periodic cell using the 3 floats as the cell parameters
    """

    if cell is not None: assert len(cell) == 3

    xyz = "{}\n".format(len(coords))
    if cell is None:
        xyz += "ambuld xyz file\n"
    else:
        xyz += "ambuld xyz file. Axes:{}:{}:{}\n".format(cell[0], cell[1], cell[2])

    for i, coord in enumerate(coords):
        if cell is None:
            x, y, z = coord
        else:
            x, ix = wrapCoord(coord[0], cell[0], center=False)
            y, iy = wrapCoord(coord[1], cell[1], center=False)
            z, iz = wrapCoord(coord[2], cell[2], center=False)

        xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format(symbols[i], x, y, z)

    with open(fileName, 'w') as f:
        fpath = os.path.abspath(f.name)
        f.writelines(xyz)

    return fpath

def hoomdCml(xmlFilename):

    tree = ET.parse(xmlFilename)
    root = tree.getroot()

    coords = []
    x = root.findall(".//position")
    ptext = x[0].text
    for line in ptext.split(os.linesep):
        line = line.strip()
        if line:
            x, y, z = line.split()
            coords.append(numpy.array([ float(x), float(y), float(z) ]))

    symbols = []
    atext = root.findall(".//type")[0].text
    for line in atext.split(os.linesep):
        atomType = line.strip()
        if atomType:
            symbols.append(label2symbol(atomType))

    bonds = []
    x = root.findall(".//bond")
    ptext = x[0].text
    for line in ptext.split(os.linesep):
        line = line.strip()
        if line:
            label, b1, b2 = line.split()
            bonds.append((b1, b2))

    writeCml(xmlFilename + ".cml",
             coords,
             symbols,
             bonds=bonds,
             atomTypes=None,
             cell=None,
             pruneBonds=False)

    return

def hoomdContacts(xmlFilename):

    tree = ET.parse(xmlFilename)
    root = tree.getroot()

    coords = []
    x = root.findall(".//position")
    ptext = x[0].text
    for line in ptext.split(os.linesep):
        line = line.strip()
        if line:
            x, y, z = line.split()
            coords.append(numpy.array([ float(x), float(y), float(z) ]))

    symbols = []
    atext = root.findall(".//type")[0].text
    for line in atext.split(os.linesep):
        atomType = line.strip()
        if atomType:
            symbols.append(label2symbol(atomType))


    # Strip x-atoms
    toGo = []
    for i, s in enumerate(symbols):
        if s.lower() == 'x':
            toGo.append(i)

    gone = 0
    for i in toGo:
        coords.pop(i - gone)
        symbols.pop(i - gone)
        gone += 1

    assert len(coords) == len(symbols)
    bonds, md = _calcBonds(coords, symbols)

    return

def xyzContacts(xyzFile):

    symbols = []
    coords = []

    with open(xyzFile, 'r') as f:
        natoms = int(f.readline().strip())
        f.readline()
        line = f.readline()
        while line:
            s, x, y, z = line.strip().split()
            symbols.append(s)
            coords.append(numpy.array([ float(x), float(y), float(z) ]))
            line = f.readline()

    assert len(coords) == natoms

    bonds, md = _calcBonds(coords, symbols)

    return

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('pkl_file', type=str, metavar='pickle_file.pkl', help='Ambuild pickle file')
    p.add_argument('-f', '--split_fragments', action='store_true', default=False, help="Split the cell into fragments")
    p.add_argument('-b', '--split_blocks', action='store_true', default=False, help="Split the cell into blocks")
    p.add_argument('-da', '--dlpoly_allatom', action='store_true', default=False, help="Create all-atom DLPOLY CONFIG and FIELD files")
    p.add_argument('-dr', '--dlpoly_rigid', action='store_true', default=False, help="Create rigid-body DLPOLY CONFIG and FIELD files")
    p.add_argument('-np', '--non_periodic', action='store_true', default=False, help="Dump a non-periodic system")
    p.add_argument('-id', '--ignore_dihedrals', action='store_true', default=False, help="Ignore missing dihedrals when creating DL-POLY files.")
    args=p.parse_args()
    
    split = None
    dlpoly = False
    rigid = True
    nonPeriodic = False
    skipDihedrals = False
    
    if args.split_fragments:
        split='fragments'
    elif args.split_blocks:
        split='blocks'
        
    if args.dlpoly_allatom:
        dlpoly = True
        rigid = False
    elif args.dlpoly_rigid:
        dlpoly = True
        rigid = True
        
    if args.non_periodic: nonPeriodic = True
    if args.ignore_dihedrals: skipDihedrals = True     

    # Need to reset sys.argv as otherwise hoomdblue eats it and complains
    sys.argv = [sys.argv[0]]

    if dlpoly:
        dumpDLPOLY(args.pkl_file, rigidBody=rigid, skipDihedrals=skipDihedrals)
    else:
        dumpPkl(args.pkl_file, split=split, nonPeriodic=nonPeriodic)
