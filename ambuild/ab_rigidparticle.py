'''
Created on Jan 15, 2013

@author: abbietrewin
'''
import copy
import logging

import numpy as np
# our imports
from ambuild import xyz_core


logger = logging.getLogger()

class RigidParticle(object):
    def __init__(self, body):
        # Attributes of the central particle
        self.image = np.array([0, 0, 0])
        self.natoms = body.natoms
        self.position = body.centreOfMass
        self.b_positions = body.coordsRelativeToCom
        self.mass = body.mass
        self.principalMoments = body.principalMoments
        # These are set by the particle manager
        self.type = None
        self.orientation = None
        # Specify properties of consituent particles
        self.b_rigidConfigStr = body.rigidConfigStr # for debugging
        self.b_charges = body.charges
        self.b_diameters = body.diameters
        self.b_masses = body.masses
        self.b_atomTypes = body.atomTypes   
        return

class RigidParticleManager(object):
    """Manages the configuration of the rigid particles
    
    As Fragments have their own class/instance managment machinary, we can't use standard class variables
    to track variables that are shared across all fragments of a given fragmentType. This class provides
    the functionality to track attributes that are shared by all fragments of a given type, but also need
    to be updated as the simulation progress.
    
    Need to think about how to reset this when running simulations and when unpickling
    
    need to track
    * bodyConfigStr
    * rigidParticleType
    * originalOrientation
    
    lookup is by bodyConfigStr, so
    bodyConfigStr -> rigidParticleType
    """
    PREFIX_STR = ')' # Sometthing that hopefully won't turn up starting an atom type
    def __init__(self):
        # Needs to maintain list of configs and the original rigid particle so we
        # can calculate a relative orientation
        self._positions = {}
        self._types = {}
        self._configStr = {}
    
    def calcConfigStr(self, idx):
        """Return a 2-letter id str for the index idx"""
        import string
        alph = string.ascii_uppercase
        num_chars = 26 # letters in alphabet - unlikely to change...
        assert idx < num_chars * num_chars, "Too many configurations!"
        return self.PREFIX_STR + alph[int(idx/num_chars)] + alph[idx % num_chars]
    
    def configStr(self, body):
        return self._configStr[body.rigidConfigStr]
    
    def createParticle(self, body):
        self.updateConfig(body)
        rigidParticle = RigidParticle(body)
        refCoords = self._positions[body.rigidConfigStr]     
        rigidParticle.orientation = xyz_core.orientationQuaternion(refCoords, body.coordsRelativeToCom)
        rigidParticle.type = self._configStr[body.rigidConfigStr]
        return rigidParticle
    
    def checkConfigStrClashes(self, atomTypes):
        overlap = set(atomTypes).intersection(set(self._configStr.values()))
        if overlap:
            raise RuntimeError("Clashing atomType and rigidPaticle type identifiers: {}".format(overlap))
    
    @property
    def referenceParticles(self):
        for k in self._configStr.keys():
            yield self._configStr[k], self._positions[k], self._types[k]
    
    def reset(self):
        self._positions.clear()
        self._types.clear()
        self._configStr.clear()

    def updateConfig(self, body):
        cstr = body.rigidConfigStr
        if cstr not in self._configStr:
            self._positions[cstr] = body.coordsRelativeToCom.copy()
            self._types[cstr] = copy.copy(body.atomTypes)
            self._configStr[cstr] = self.calcConfigStr(len(self._positions) - 1)

