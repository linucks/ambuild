#!/usr/bin/env python
'''
Created on Feb 3, 2013

@author: jmht

Utility functions
'''
import logging
import os
import numpy as np
import math
import warnings
import xml.etree.ElementTree as ET
import xml.dom.minidom

from ambuild import xyz_core
from ambuild.ab_ffield import read_bond_params

logger = logging.getLogger()

DUMMY_DIAMETER = 0.1

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
        if atomType1 in self.ATOM_TYPE_BOND_LENGTHS and atomType2 in self.ATOM_TYPE_BOND_LENGTHS[ atomType1 ]:
            #print "ATOM TYPE"
            return self.ATOM_TYPE_BOND_LENGTHS[ atomType1 ][ atomType2 ]
        elif atomType2 in self.ATOM_TYPE_BOND_LENGTHS and atomType1 in self.ATOM_TYPE_BOND_LENGTHS[ atomType2 ]:
            #print "ATOM TYPE"
            return self.ATOM_TYPE_BOND_LENGTHS[ atomType2 ][ atomType1 ]
        
        symbol1 = label2symbol(atomType1).upper()
        symbol2 = label2symbol(atomType2).upper()
        if symbol1 in xyz_core.ELEMENT_TYPE_BOND_LENGTHS and symbol2 in xyz_core.ELEMENT_TYPE_BOND_LENGTHS[ symbol1 ]:
            #print "ELEMENT TYPE"
            return xyz_core.ELEMENT_TYPE_BOND_LENGTHS[ symbol1 ][ symbol2 ]
        elif symbol2 in xyz_core.ELEMENT_TYPE_BOND_LENGTHS and symbol1 in xyz_core.ELEMENT_TYPE_BOND_LENGTHS[ symbol2 ]:
            #print "ELEMENT TYPE"
            return xyz_core.ELEMENT_TYPE_BOND_LENGTHS[ symbol2 ][ symbol1 ]
        warnings.warn('No data for bond length for %s-%s' % (atomType1, atomType2))
        return 1.0

# This needs to be set to the bondLength function of the BondLength class 
# See cell.Cell_setUtilBondLength
# We set this to raise an error if unset
def __STOP(x,y):
    raise NotImplementedError("Need to set the bondLength function - see setModuleBondLength")
bondLength = __STOP
    

def setModuleBondLength(paramFile):
    """Set the bondLength function of the util module"""
    global bondLength
    BL = BondLength(paramFile)
    bondLength = BL.bondLength
    return


def calcBondsHACK(coords, symbols, bondMargin=0.2):
    """HACK FOR NETWORK"""
    bonds = []
    for i, coord1 in enumerate(coords):
        symbol1 = symbols[ i ]
        for j, coord2 in enumerate(coords):
            if j < i:
                continue
            symbol2 = symbols[ j ]
            dist = xyz_core.distance(coord1, coord2)
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
    
    distances = xyz_core.distance(np.array(v1), np.array(v2), dim=dim)
    bonds = []
    
    for i, (idxAtom1, idxAtom2) in enumerate(close):
        bond_length = bondLength(symbols[idxAtom1], symbols[idxAtom2])
        if bond_length < 0: continue
        logger.debug("Dist:length {1}({0})-{3}({2}) {4} {5}".format( symbols[idxAtom1], idxAtom1, symbols[idxAtom2], idxAtom2, distances[i], bond_length ))
        if  bond_length - bondMargin < distances[i] < bond_length + bondMargin:
            bonds.append((idxAtom1, idxAtom2))
    
    # We sort the bonds on return so that results are deterministic and can be tested
    return sorted(bonds)


def closeAtoms(coords, symbols, dim=None, maxAtomRadius=None, boxMargin=1.0):
    """Return a list of which atoms in the cell are close to each other.
    Close is defined by the maxAtomRadius and boxMargin"""
    if maxAtomRadius is None:
        maxAtomRadius = max([xyz_core.COVALENT_RADII[xyz_core.SYMBOL_TO_NUMBER[s.upper()]] * xyz_core.BOHR2ANGSTROM for s in set(symbols)])

    boxSize = (maxAtomRadius * 2) + boxMargin
    # If we are under PBC calculate the number of boxes in each dimenson
    boxNum = None
    if dim is not None:
        if type(dim) is list:
            dim = np.array(dim)
        boxNum = [int(math.ceil(dim[0] / boxSize)),
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
        if key in cells:
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
            try:
                alist = cells[ scell ]
            except KeyError:
                continue  # Trigger exception to save searching through the keys
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


def label2symbol(name):
    """ Determine the element type of an atom from its name, e.g. Co_2b -> Co
        Returns a capitalised element name
    """
    origName = name
    name = name.strip().upper()
    # Determine the element from the first 2 chars of the name
    if len(name) > 2:
        name = name[0:2]
    if len(name) == 2 and name[0].isalpha() and name[1].isalpha():
        # 2 Character name, so see if it matches any 2-character elements
        sym2c = list(filter(lambda x: len(x) == 2, xyz_core.SYMBOL_TO_NUMBER.keys()))
        # HACK: NEED TO REMOVE NP and NB
        sym2c.remove('NP')
        sym2c.remove('NB')
        if name in sym2c:
            return name.capitalize()
    # If it was a valid 2 character symbol we should have picked it up so now only 1 symbol
    name = name[0]
    if not name.isalpha():
        raise RuntimeError("label2symbol first character of name is not a character: {0}".format(origName))
    # Hack - for x return x
    if name.lower() == 'x':
        return 'x'
    # Get 1 character element names
    sym1c = list(filter(lambda x: len(x) == 1 and x != 'X', xyz_core.SYMBOL_TO_NUMBER.keys()))
    if name in sym1c:
        return name.capitalize()
    raise RuntimeError("label2symbol cannot convert name {0} to symbol!".format(origName))
    return


def writeCml(cmlFilename,
             coords,
             symbols,
             bonds=None,
             atomTypes=None,
             cell=None,
             pruneBonds=False,
             prettyPrint=False):
    
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
        # Need to collate atomTypes - sorted to ensure consistency for testing
        for atype in sorted(set(atomTypes)):
            atomTypeNode = ET.SubElement(root, "atomType")
            atomTypeNode.attrib['name'] = atype
            atomTypeNode.attrib['title'] = atype
    # Now atom data
    atomArrayNode = ET.SubElement(root, "atomArray")
    for i, coord in enumerate(coords):
        atomNode = ET.SubElement(atomArrayNode, "atom")
        atomNode.attrib['id'] = "a{0}".format(i)
        atomNode.attrib['elementType'] = symbols[i]
        if cell is not None:
            coord, _ = xyz_core.wrapCoord3(coord, cell, center=False)
            if pruneBonds:
                pcoords.append(coord)
        x, y, z = coord
        atomNode.attrib['x3'] = "{:.12}".format(x)
        atomNode.attrib['y3'] = "{:.12}".format(y)
        atomNode.attrib['z3'] = "{:.12}".format(z)
        if atomTypes:
            # Now add atomType as child node referring to the atomType
            atomTypeNode = ET.SubElement(atomNode, "atomType")
            atomTypeNode.attrib['ref'] = atomTypes[ i ]
    if len(bonds):
        # Now do bonds
        if pruneBonds:
            # Hack to get vis working
            # Calculate all bond distances
            distances = xyz_core.distance([ pcoords[b1] for b1, b2 in bonds],
                                  [ pcoords[b2] for b1, b2 in bonds])

            # Complete hack - just see if it's longer then 0.5 the cell A distance - assumes cubic cell
            bonds = [ b for i, b in enumerate(bonds) if distances[i] <= cell[0] * 0.5 ]

        bondArrayNode = ET.SubElement(root, "bondArray")
        for b in bonds:
            bondNode = ET.SubElement(bondArrayNode, "bond")
            bondNode.attrib['atomRefs2'] = "a{0} a{1}".format(b[0], b[1])
            bondNode.attrib['order'] = "1"

    cmlFilename = os.path.abspath(cmlFilename)
    if prettyPrint:
        estring = xml.dom.minidom.parseString(ET.tostring(root))
        with open(cmlFilename, 'w') as w:
            w.write(estring.toprettyxml())
    else:
        tree = ET.ElementTree(root)
        tree.write(cmlFilename, encoding="utf-8", xml_declaration=True)
    return cmlFilename

def writeXyz(fileName, coords, symbols, cell=None):
    """Write out an xyz file to fileName
    if cell is a list of 3 floats, write out a periodic cell using the 3 floats as the cell parameters
    """
    if cell is not None:
        assert len(cell) == 3

    xyz = "{}\n".format(len(coords))
    if cell is None:
        xyz += "ambuld xyz file\n"
    else:
        xyz += "ambuld xyz file. Axes:{}:{}:{}\n".format(cell[0], cell[1], cell[2])

    for i, coord in enumerate(coords):
        if cell is not None:
            coord, _ = xyz_core.wrapCoord3(coord, cell, center=False)
        x, y, z = coord
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
            coords.append(np.array([ float(x), float(y), float(z) ]))

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
            _, b1, b2 = line.split()
            bonds.append((b1, b2))

    writeCml(xmlFilename + ".cml",
             coords,
             symbols,
             bonds=bonds,
             atomTypes=None,
             cell=None,
             pruneBonds=False)

    return
