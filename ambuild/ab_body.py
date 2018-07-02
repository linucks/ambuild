'''
Created on Jan 15, 2013

@author: abbietrewin
'''
import logging

import numpy as np
# our imports
import xyz_core
import xyz_util


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
