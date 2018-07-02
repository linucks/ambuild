'''
Created on Jan 15, 2013

@author: abbietrewin
'''
import logging

import numpy as np
# our imports
import xyz_core


logger = logging.getLogger()

class RigidParticle(object):
    def __init__(self, body):
        # Attributes of the central particle
        self.image = np.array([0, 0, 0])
        self.mass = None
        self.orientation = None
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
        
        self.fromBody(body)
    
    def fromBody(self, body):
        coords = body.coords
        self.natoms = coords.shape[0]
        self.position = xyz_core.centreOfMass(coords, body.masses)
        self.b_positions = coords - self.position # coordinate positions relative to com
        self.mass = np.sum(body.masses)
        self.orientation = body.orientation
        self.principalMoments = xyz_core.principalMoments(coords, body.masses)
        self.type = body.rigidType
        # Specify properties of consituent particles
        self.b_charges = body.charges
        self.b_diameters = body.diameters
        self.b_masses = body.masses
        self.b_atomTypes = body.atomTypes
        return self

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
    
    def __init__(self):
        # Needs to maintain list of configs and the original rigid particle so we
        # can calculate a relative orientation
        self._configs = {}
        self._configStr = {}
    
    @staticmethod
    def calcConfigStr(idx):
        """Return a 2-letter id str for the index idx"""
        import string
        alph = string.ascii_uppercase
        num_chars = 26 # letters in alphabet - unlikely to change...
        assert idx < num_chars * num_chars, "Too many configurations!"
        return alph[int(idx/num_chars)] + alph[idx % num_chars]
    
    def configStr(self, body):
        return self._configStr[body.rigidType]
    
    def orientation(self, body):
        refCoords = self._configs[body.rigidType]
        M = xyz_core.rigid_rotate(body.coords, refCoords)
        return xyz_core.rotation_matrix_to_quaternion(M)
    
    def reset(self):
        self._configs.clear()

    def updateConfig(self, fragment):
        for body in fragment.bodies():
            if body.rigidType not in self._configs:
                self._configs[body.rigidType] = body.coords.copy()
                self._configStr[body.rigidType] = self.calcConfigStr(len(self._configs) - 1)
        return