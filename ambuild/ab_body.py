'''
Created on Jan 15, 2013

@author: abbietrewin
'''
import logging

import numpy as np
# our imports
import ab_rigidparticle
import xyz_core
import xyz_util


logger = logging.getLogger()


class Body(object):
    """Proxy object for the body within a fragment"""
    def __init__(self, fragment, bodyIdx):
        self.fragment = fragment
        self.bodyIdx = bodyIdx
        self.indexes = self._setup_indexes()
        self.natoms = np.sum(self.indexes)
        self._centreOfMass = xyz_core.centreOfMass(self.coords, self.masses)
        return

    def _setup_indexes(self):
        return (np.array(self.fragment._bodies) == self.bodyIdx) &~ np.array(self.fragment.masked)

    @property
    def atomTypes(self):
        return list(np.compress(self.indexes, self.fragment._atomTypes, axis=0))

    @property
    def bodies(self):
        return [self.bodyIdx] * self.natoms
        
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
        for i in [ i for i in self.fragment._ext2int.values() if self.fragment._bodies[i] == self.bodyIdx]:
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
    def orientation(self):
        return self.fragment.cell.rigidParticleMgr.orientation(self)
    
    @property
    def principalMoments(self):
        return xyz_core.principalMoments(self.coords(), self.masses)
    
    @property
    def rigidType(self):
        """return the type of this body based on the endGroup configuration"""
        return "{}{}".format(self.bodyIdx, self.fragment.configStr)
    
    def rigidParticle(self):
        return ab_rigidparticle.RigidParticle(self)

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
