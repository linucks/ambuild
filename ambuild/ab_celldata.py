class CellData(object):
    def __init__(self):
        self.cell = []
        self.natoms = 0
        self.atomTypes = []
        self.bodies = []
        self.coords = []
        self.charges = []
        self.diameters = []
        self.images = []
        self.masses = []
        self.masked = [] # Atoms that are to be ignored in MD/optimisation
        self.symbols = []
        self.static = []  # If this atom is part of a group that isn't to be moved
        # multi-particle properties
        self.bonds = []
        self.bondLabels = []
        self.angles = []
        self.angleLabels = []
        self.propers = []
        self.properLabels = []
        self.impropers = []
        self.improperLabels = []
        # for computing block/fragment enegies
        self.tagIndices = []
        # Central particles for hoomd-blue rigid bodies
        self.rigidParticles = []
        return
