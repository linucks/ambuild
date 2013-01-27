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
import os
import re
import copy
import unittest

import numpy


# This stolen from the CCP1GUI: http://sourceforge.net/projects/ccp1gui/


# double check these values...
#hvd values obtained from http://www.webelements.com/ and recorded to their
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
# note that indexing with -1 accesses the last element (0.4 copied at start and edn)
#
COVALENT_RADII = [
0.4,
0.7,                                   3.80,    
2.76,1.99,                1.62,1.33,1.23,1.14,0.95,3.80,
3.42,2.85,                2.38,2.09,1.90,1.90,1.90,3.80,
4.18,3.42,
     3.04,2.66,2.57,2.66,2.66,2.66,2.57,2.57,2.57,2.57,
                          2.47,2.38,2.18,2.18,2.18,3.80,
4.46,3.80,               
     3.42,2.94,2.76,2.76,2.57,2.47,2.57,2.66,3.04,2.94,
                          2.94,2.76,2.76,2.66,2.66,3.80,
4.94,4.09,
     3.71,
     3.52,3.52,3.52,3.52,3.52,3.52,3.42,3.33,3.33,3.33,3.33,3.33,3.33,3.33,
          2.94,2.76,2.57,2.57,2.47,2.57,2.57,2.57,2.85,
                         3.61,3.42,3.04,3.61,3.61,3.80,
4.94,4.09,
     3.71,
     3.42,3.42,3.33,3.33,3.33,3.33,3.23,3.13,3.13,3.13,3.13,3.13,3.13,3.13, 1., 1., 1.,1.,1.,1., 0.4,]
     

# joe lennards table angstroms
# 1.0 added as index [0] and [-1] see above
VDW_RADII = [  1.0,
  1.20,                              1.40,
  1.82,1.78,1.74,1.70,1.55,1.52,1.47,1.54,
  2.27,2.22,2.16,2.10,1.80,1.80,1.75,1.88,
  2.75,2.57,
       2.56,2.54,2.52,2.50,2.48,2.46,2.44,2.42,2.41,2.40,
            2.40,2.10,1.85,1.90,1.85,2.02,
  3.10,2.80,
       2.77,2.74,2.71,2.68,2.65,2.62,2.59,2.56,2.53,2.51,
            2.50,2.20,2.10,2.06,1.98,2.16, 1., 1., 1.]


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
S       2.24 1.82 2.05 1.56     1.34                     2.00
Sb              2.33           1.70
Se           1.95      1.71      1.47                               2.33
Si      2.21 1.87 2.05 1.58      1.48 2.44      1.63      2.14           2.33
Sn           2.14 2.28           1.71 2.67
Te                     1.82      1.66
"""

# REM - symbols should be in lower case!
BOND_LENGTHS = {}
BOND_LENGTHS['AS'] = { 'AS' : 2.10,
                       'BR' : 2.32,
                       'C'  : 1.96,
                       'CL' : 2.17,
                       'F'  : 1.71,
                       'H'  : 1.51 }
                       
BOND_LENGTHS['BR'] = { 'BR': 2.28,
                       'C' : 1.94,
                       'CL': 2.14,
                       'F' : 1.76,
                       'GE': 2.30,
                       'H' : 1.41,
                       'I' : 2.47,
                       'P' : 2.22,
                       'S' : 2.24,
                       'SI' : 2.21 }

BOND_LENGTHS['C'] = { 'C' : 1.53,
                      'CL': 1.79,
                      'F' : 1.39,
                      'GE': 1.95,
                      'H' : 1.09,
                      'I' : 2.13,
                      'N' : 1.46,
                      'O' : 1.42,
                      'P' : 1.85,
                      'S' : 1.82,
                      'SE': 1.95,
                      'SI': 1.87,
                      'SN': 2.14 }

BOND_LENGTHS['CL'] = { 'CL' : 1.99,
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

BOND_LENGTHS['F'] = { 'F'  : 1.41,
                      'GE' : 1.73,
                      'H'  : 0.92,
                      'I'  : 1.91,
                      'N'  : 1.37,
                      'O'  : 1.42,
                      'P'  : 1.57,
                      'S'  : 1.56,
                      'SE' : 1.71,
                      'SI' : 1.58,
                      'TE' : 1.82 }

BOND_LENGTHS['GE'] = { 'GE' : 2.40,
                       'H'  : 1.53,
                       'I'  : 2.51 }
                       
BOND_LENGTHS['H'] = { 'H' : 0.74,
                      'I' : 1.61,
                      'N' : 1.02,
                      'O' : 0.96,
                      'P' : 1.42,
                      'S' : 1.34,
                      'SB': 1.70,
                      'SE': 1.47,
                      'SI': 1.48,
                      'SN': 1.71,
                      'TE': 1.66 }

BOND_LENGTHS['I'] = { 'I'  : 2.67,
                      'SI' : 2.44,
                      'SN' : 2.67 }

BOND_LENGTHS['N'] = { 'N' : 1.45,
                      'O' : 1.43,
                      'P' : 1.65 }

BOND_LENGTHS['O'] = { 'O'  : 1.48,
                      'SI' : 1.63 }

BOND_LENGTHS['P'] = { 'P' : 2.25 }

BOND_LENGTHS['S'] = { 'S' : 2.00 }

BOND_LENGTHS['SE'] = { 'SE' : 2.33 }

BOND_LENGTHS['SI'] = { 'SI' : 2.33 }

def get_bond_length( symbol1,symbol2 ):
    """ Get the characteristic lengths of single bonds as defined in:
        Reference: CRC Handbook of Chemistry and Physics, 87th edition, (2006), Sec. 9 p. 46
        If we can't find one return 1.0 as a default 
    """
    global BOND_LENGTHS

    symbol1 = symbol1.upper()
    symbol2 = symbol2.upper()

    #print "Getting bond length for %s-%s" % ( symbol1, symbol2 )

    if BOND_LENGTHS.has_key( symbol1 ):
        if BOND_LENGTHS[ symbol1 ].has_key( symbol2 ):
            return BOND_LENGTHS[ symbol1 ][ symbol2 ]
        
    if BOND_LENGTHS.has_key( symbol2 ):
        if BOND_LENGTHS[ symbol2 ].has_key( symbol1 ):
            return BOND_LENGTHS[ symbol2 ][ symbol1 ]
    
    print 'No data for bond length for %s-%s' % (symbol1,symbol2)
    return 1.0


class BuildingBlock():
    '''
    classdocs
    '''

    def __init__( self, infile = None ):
        '''
        Constructor
        '''
        
        # An python array of numpy.array
        self.coords = []
        
        # ordered array of labels
        self.labels = []
        
        # ordered array of symbols (in upper case)
        self.symbols = []
        
        # ordered array of masses
        self.masses = []
        
        # orderd array of atom radii
        self.atom_radii = []
        
        # Groups at end that can particpate in bonds
        self.endGroups = []
        
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
        for i, coord in enumerate(block.coords):
            self.coords.append( coord )
            self.atom_radii.append( block.atom_radii[i] )
            self.labels.append( block.labels[i] )
            self.symbols.append( block.symbols[i] )
            self.masses.append( block.masses[i] )
            self.endGroups.append( block.endGroups[i] )
            
        # Now add bond
        
        
    def createFromArgs(self, coords, labels, endGroups ):
        """ Create from given arguments
        """
        # could check if numpy array here
        self.coords = coords
        self.labels = labels
        self.endGroups = endGroups
        
        self.fillData()
        
        self.calcCenterOfMassAndGeometry()
        self.calcRadius()
    
  
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
        
        cog = self.centerOfGeometry()
        
        distances = []
        for coord in self.coords:
            distances.append(  numpy.linalg.norm(coord-cog) )
            
        imax = numpy.argmax( distances )
        dist = distances[ imax ]
        atomR = self.atom_radii[ imax ]
        
        # Set radius
        self._radius = dist + atomR
        
    def canBond( self, block ):
        """See if we can form a bond with the given block.
        Return the indicies of the two atoms (self and then other)
        or False if the molecules cannot bond but do not clash
        If the blocks clash (are close but cannot bond) we return
        the string "clash" to indicate a clash and that the move should be
        rejected
        """
        # First see if we are close enough to consider bonding
        #jmht- check longest possible bond - possibly precalculate for each molecule?
        if not self.close(block, margin=2.0):
            return False
        
        # Might be able to bond so loop through all atoms and check
        # We only accept the case where just one pair of atoms is close enough to bond
        # more than one is considered a failure
        global get_bond_lengths
        MARGIN = 0.2
        bond = None
        for i, c in enumerate( self.coords ):
            i_radius = self.atom_radii[i]
            i_symbol = self.symbols[i]
            for j, b in enumerate( block.coords ):
                j_radius = block.atom_radii[j]
                j_symbol = block.symbols[j]
                bond_length = get_bond_length( i_symbol , j_symbol )
                if ( numpy.linalg.norm( c - b ) > bond_length - MARGIN and
                     numpy.linalg.norm( c - b ) < bond_length + MARGIN):
                    # Check if both atoms are end groups
                    if not ( self.endGroups[i] and block.endGroups[j] ):
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
        return bond
        
    def centerOfMass(self):
        """
        Return or calculate the center of mass for the building block.
        """
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self._changed = False
        
        return self._centerOfMass
    
    def clash( self, block ):
        """ See if this molecule clashes (overlaps) with the given one """
        
        # First check if these two molecules are within range assuming they 
        # are circular, centered on the COG and with the given radius
        if not self.close(block):
            return False
        
        # Within range, so see if any atoms actually do overlap
        # Do this with scipy and no loops when we optimise
        MARGIN = 1
        for i, c in enumerate( self.coords ):
            i_radius = self.atom_radii[i]
            for j, b in enumerate( block.coords ):
                j_radius = block.atom_radii[j]
                if ( numpy.linalg.norm( c - b ) < i_radius + j_radius + MARGIN ):
                    return True
                
        return False
    
    def close( self, block, margin=1.0 ):
        """Return true of false depending on whether two blocks are close enough to bond/clash.
        Works from the overall radii of the two blocks
        Margin is allowed gap between their respective radii
        """
        dist = numpy.linalg.norm( self.centerOfGeometry() - block.centerOfGeometry() )
        if ( dist < self.radius() + block.radius() + margin ):
            return True
        else:
            return False
    
    def copy( self ):
        """Create a copy of ourselves and return it"""
        return copy.deepcopy(self)

    def centerOfGeometry(self):
        """
        Return or calculate the center of geometry for the building block.
        """
        
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self._changed = False
        
        return self._centerOfGeometry

    def fillData(self):
        """ Fill the data arrays from the label """
        
        global SYMBOL_TO_NUMBER, COVALENT_RADII

        if len(self.labels) < len(self.coords):
            raise RuntimeError("fillData needs labels filled!")
        
        for label in self.labels:
            # Symbols
            symbol = self.labelToSymbol( label )
            self.symbols.append( symbol )
            
            # Masses
            self.masses.append( ATOMIC_MASS[ symbol ] )
            
            # Radii
            z = SYMBOL_TO_NUMBER[ symbol ]
            r = COVALENT_RADII[z]
            self.atom_radii.append(r)
            #print "ADDING R {} for label {}".format(r,label)
            
    
    def fromCarFile(self, carFile):
        """"Abbie did this.
        Gulp...
        """
        labels = []
        
        # numpy array
        coords = []
        
        # array of bools - could change to an array of the indexes
        endGroups = []
        
        reading = True
        with open( carFile ) as f:
            
            count = 0
            # First 4 lines just info - not needed
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            
            while reading:
                line = f.readline()
                
                line = line.strip()
                fields = line.split()
                label = fields[0]
                if label.lower() == "end":
                    reading=False
                    break
                     
                labels.append(label) 
                
                # End groups are denoted by an underscore at the end
                # of the label name
                if fields[0].endswith('_'):
                    endGroups.append(True)
                else:
                    endGroups.append(False)
                coords.append( numpy.array(fields[1:4], dtype=numpy.float64) )
                
        self.createFromArgs(coords, labels, endGroups)


    def fromXyzFile(self, xyzFile):
        """"Jens did this.
        """
        
        labels = []
        
        # numpy array
        coords = []
        
        # array of bools - could change to an array of the indexes
        endGroups = []
        
        with open( xyzFile ) as f:
            
            # First line is number of atoms
            line = f.readline()
            natoms = int(line.strip())
            
            # Skip title
            line = f.readline()
            
            for i in range(natoms):
                
                line = f.readline()
                line = line.strip()
                fields = line.split()
                labels.append(fields[0]) 
                
                # End groups are denoted by an underscore at the end
                # of the label name
                if fields[0].endswith('_'):
                    endGroups.append(True)
                else:
                    endGroups.append(False)
                    
                coords.append( numpy.array(fields[1:4], dtype=numpy.float64) )
                
        self.createFromArgs(coords, labels, endGroups)
                
    def labelToSymbol( self, name ):
        """ Determine the element type of an atom from its name, e.g. Co_2b -> Co
            Returns a capitalised element name
            Originally written by Jens Thomas in the CCP1GUI
        """
    
        # Determine the element from the first 2 chars of the name
        if ( len( name ) == 1 ):
            if not re.match( '[a-zA-Z]', name ):
                print "Error converting name to symbol for atom %s!" % name
                element = 'XX'
            else:
                element = name
        else:
            # See if 2nd char is a character - if so use 1st 2 chars as symbol
            if re.match( '[a-zA-Z]', name[1] ):
                element = name[0:2]
            else:
                element = name[0]
    
        return element.capitalize()
    
    def radius(self):
        
        if self._changed:
            self.calcCenterOfMassAndGeometry()
            self.calcRadius()
            self._changed = False
        
        return self._radius


    def rotate( self, axis, angle ):
        """ Rotate the molecule about the given axis by the angle in radians - I think...
        """
        
        rmat = self.rotation_matrix( axis, angle )
        
        # loop through and change all coords
        for i,coord in enumerate( self.coords ):
            self.coords[i] = numpy.dot( rmat,coord )
        
        self._changed = True

    def rotation_matrix(self, axis,angle):
        """
        http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
        """
        axis = axis/numpy.sqrt(numpy.dot(axis,axis))
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
    
    def translateCenterOfGeometry( self, position ):
        """Translate the molecule so the center of geometry moves
        to the given position
        """
        self.translate( position - self.centerOfGeometry() )
        
    

        
        
    def writeXyz(self,name=None):
        
        if not name:
            name = str(id(self))+".xyz"
            
        with open(name,"w") as f:
            fpath = os.path.abspath(f.name)
            f.write( "{}\n".format(len(self.coords)) )
            f.write( "id={}\n".format(str(id(self))) )
                             
            for i,c in enumerate(self.coords):
                f.write("{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( self.labelToSymbol(self.labels[i]), c[0], c[1], c[2]))
        
        print "Wrote file: {0}".format(fpath)
    
    def __add__(self):
        """Add two blocks - concatenate the coords and recalculate the variables"""
        
        
        

        
    def __str__(self):
        """
        Return a string representation of the molecule
        """
        
        mystr = ""
        for i,c in enumerate(self.coords):
            #mystr += "{0:4}:{1:5} [ {2:0< 15},{3:0< 15},{4:0< 15} ]\n".format( i+1, self.labels[i], c[0], c[1], c[2])
            mystr += "{0:5} {1:0< 15} {2:0< 15} {3:0< 15} \n".format( self.labels[i], c[0], c[1], c[2])
            
        mystr += "radius: {}\n".format( self.radius() )
        mystr += "COM: {}\n".format( self._centerOfMass )
        mystr += "COG: {}\n".format( self._centerOfGeometry )
        
        return mystr
        
        
class TestBuildingBlock(unittest.TestCase):

    def setUp(self):
        """Create a methane molecule for testing"""
        
        coords = [ numpy.array([  0.000000,  0.000000,  0.000000 ] ),
        numpy.array([  0.000000,  0.000000,  1.089000 ]),
        numpy.array([  1.026719,  0.000000, -0.363000 ]),
        numpy.array([ -0.513360, -0.889165, -0.363000 ]),
        numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]

        #numpy.array([ -0.513360,  0.889165, -0.363000 ]) ]
        labels = [ 'C', 'H', 'H', 'H', 'H' ]
        
        endGroups = [1,2,3,4]
        
        self.ch4 = BuildingBlock()
        self.ch4.createFromArgs( coords, labels, endGroups )
        
    def testRotate(self):
        """
        Test the rotation
        """
        
        array1 = numpy.array([ -0.51336 ,  0.889165, -0.363 ])        
        self.assertTrue( numpy.array_equal(self.ch4.coords[4], array1 ),
                         msg="testRotate arrays before rotation incorrect.")
        
        axis = numpy.array([1,2,3])
        angle = 2
        self.ch4.rotate(axis, angle)
        
        array2 = numpy.array([  1.05612011, -0.04836936, -0.26113713 ])

        # Need to use assertTrue as we get a numpy.bool returned and need to test this will
        # bool - assertIs fails
        self.assertTrue( numpy.allclose(self.ch4.coords[4], array2, rtol=1e-9, atol=1e-8 ),
                         msg="testRotate arrays after rotation incorrect.")

    def testCenterOfMass(self):
        """
        Test calculation of Center of Mass
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        com = self.ch4.centerOfMass()
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

    def testCenterOfGeometry(self):
        """
        Test calculation of Center of Geometry
        """
        
        correct = numpy.array([  0.000000,  0.000000,  0.000000 ])
        cog = self.ch4.centerOfGeometry()
        self.assertTrue( numpy.allclose( correct, cog, rtol=1e-9, atol=1e-6 ),
                         msg="testCenterOfGeometry incorrect COM.")
        
    def testCalcRadius(self):
        """
        Test calculation of the radius
        """
        
        r = self.ch4.radius()
        self.assertAlmostEqual(r, 1.78900031214, 7, "Incorrect radius: {}".format(str(r)) )
        
    def testClash(self):
        """
        Test we can spot a clash
        """
        
        block = self.ch4.copy()
        self.assertTrue( self.ch4.clash( block ) )
        
        block.translateCenterOfGeometry( [3,3,3] )
        self.assertFalse( self.ch4.clash( block ) )
        

if __name__ == '__main__':
    """
    Run the unit tests
    """
    unittest.main()
        