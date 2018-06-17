class CellData(object):
    def __init__(self):
        self.cell = []
        
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
        self.rigid_body = []
        self.rigid_image = []
        self.rigid_mass = []
        self.rigid_centre = []
        self.rigid_type = []
        self.rigid_moment_inertia = []
        self.rigid_orientation = []
        self.rigid_fragments = {} # key fragmentType -> {'atomTypes': list, 'coords' : list}
        return
