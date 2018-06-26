#!/usr/bin/env python
import sys
sys.path.insert(0, '../ambuild')
from xyz_util import writeXyz
import itertools
import numpy as np
import hoomd
import hoomd.md

BOX_WIDTH = 20.0

def setup():
    from ab_cell import Cell
    from ab_block import Block
    
    boxDim = [BOX_WIDTH, BOX_WIDTH, BOX_WIDTH]
    mycell = Cell(boxDim)
    mycell.libraryAddFragment(filename='../blocks/ch4.car', fragmentType='A')
    mycell.addBondType('A:a-A:a')
    # Create block manually
    b1 = Block(filePath='../blocks/ch4.car', fragmentType='A')
    b2 = Block(filePath='../blocks/ch4.car', fragmentType='A')
    # Position block so that it's aligned along x-axis
    b1.alignAtoms(0, 1, [ 1, 0, 0 ])
    b1.translateCentroid([ mycell.dim[0] / 2, mycell.dim[1] / 2, mycell.dim[2] / 2 ])
    endGroup1 = b1.freeEndGroups()[ 0 ]
    endGroup2 = b2.freeEndGroups()[ 0 ]
    b1.positionGrowBlock(endGroup1, endGroup2)
    b1.translate([-2.0, 0.0, 0.0])
    
    #symbols = ['C', 'H', 'H', 'H', 'H']
    b = list(b1.fragments[0].bodies())[0]
    print b.symbols()
    print b.coords()
    
    b = list(b2.fragments[0].bodies())[0]
    print b.symbols()
    print b.coords()
#         print "1 ",b.coords()
#         print "2 ",b.coords(bodyFrame=True)
#         writeXyz('f1.xyz', b.coords(), symbols)
#         writeXyz('f2.xyz', b.coords(bodyFrame=True), symbols)
    sys.exit()
    
    mycell.addBlock(b1)
    mycell.addBlock(b2)
    mycell.dataDict(self, rigidBody=True, periodic=False, center=True)
    #mycell.dump()
    #mycell.runMD(rigidBody=False, dump=True)


#boxdim = np.array([BOX_WIDTH, BOX_WIDTH, BOX_WIDTH])
#setup()
#sys.exit(1)

class CentralParticle(object):
    def __init__(self, molecule):
        self.mass = None
        self.position = None
        self.principalMoments = None
        self.type = None
        
        self.m_positions = None
        self.m_types = None
        
        self.fromMolecule(molecule)
    
    def fromMolecule(self, molecule):
        self.mass = np.sum(molecule.masses)
        self.position = centreOfMass(molecule.positions, molecule.masses)
        self.principalMoments = principalMoments(molecule.positions, molecule.masses)
        
        self.m_positions = molecule.positions
        self.m_types = molecule.types
        return self

class Molecule(object):
    def __init__(self):
        self.positions = None
        self.masses = None
        self.types = None
    
    def centralParticle(self):
        return CentralParticle(self)


def wrapBox(positions, boxdim):
    # Move into centre of cell
    positions = np.fmod(positions, boxdim)
    # Change the coord so the origin is at the center of the box (we start from the corner)
    positions -= boxdim / 2
    return positions

def centreOfMass(coords, masses):
    totalMass = np.sum(masses)
    return np.sum(coords * masses[:,np.newaxis], axis=0) / totalMass

def momentOfInertia(coords, masses):
    """Moment of Inertia Tensor"""
    # positions relative to the centre of mass
    coords = coords - centreOfMass(coords, masses)
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
    I = momentOfInertia(coords, masses)
    eigval, eigvec = np.linalg.eig(I)
    return np.sort(eigval)


# Create 2 molecules
mol1 = Molecule()
mol1.positions = np.array([[ -2.00000000e+00,   0.00000000e+00,  -2.00000001e-07],
                           [ -9.11000000e-01,   0.00000000e+00,  -2.00000001e-07],
                           [ -2.36300000e+00,   0.00000000e+00,  -1.02671920e+00],
                           [ -2.36300000e+00,  -8.89165000e-01,   5.13359800e-01],
                           [ -2.36300000e+00,   8.89165000e-01,   5.13359800e-01]])
mol1.types = ['C', 'H', 'H', 'H', 'H']
mol1.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])

mol2 = Molecule()
mol2.positions = np.array([[  3.53000000e+00,   0.00000000e+00,  -2.00000001e-07],
                           [  2.44100000e+00,   0.00000000e+00,  -2.00000001e-07],
                           [  3.89300000e+00,   0.00000000e+00,   1.02671880e+00],
                           [  3.89300000e+00,  -8.89165000e-01,  -5.13360200e-01],
                           [  3.89300000e+00,   8.89165000e-01,  -5.13360200e-01]])
mol2.types = ['C', 'H', 'H', 'H', 'H']
mol2.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])

molecules = [mol1, mol2]

# Now create the central particles
centralParticles = []
for i, mol in enumerate(molecules):
    cp = mol.centralParticle()
    cp.type = "CP%d" % i
    centralParticles.append(cp)
    
# Need to set the positions of the molcules particles to be relative to those of the central particle
for i, cp in enumerate(centralParticles):
    mol = molecules[i]
    mol.positions = mol.positions - cp.position
    print "BETTER WAY TO DO THIS"
    cp.m_positions = mol.positions

# Get total number of particles and the set of types
particle_types = set()
nparticles = 0
exclusions = []
for cp in centralParticles:
    nparticles += 1
    particle_types.add(cp.type)
    exclusions.append(cp.type)
for mol in molecules:
    nparticles += mol.positions.shape[0]
    particle_types.update(set(mol.types))
particle_types = list(particle_types)
#
# Start of hoomd code
#
hoomd.context.initialize()

lx = BOX_WIDTH
ly = BOX_WIDTH
lz = BOX_WIDTH
snapshot = hoomd.data.make_snapshot(N=nparticles,
                                    box=hoomd.data.boxdim(Lx=lx, Ly=ly, Lz=lz),
                                    particle_types=particle_types)

for i, cp in enumerate(centralParticles):
    snapshot.particles.body[i] = i
    snapshot.particles.mass[i] = cp.mass
    snapshot.particles.position[i] = cp.position
    snapshot.particles.typeid[i] = snapshot.particles.types.index(cp.type)
    snapshot.particles.moment_inertia[i] = cp.principalMoments

idx = i + 1
for i, mol in enumerate(molecules):
    for j in range(mol.positions.shape[0]):
        snapshot.particles.body[idx] = i
        snapshot.particles.mass[idx] = mol.masses[j]
        snapshot.particles.position[idx] = mol.positions[j]
        snapshot.particles.typeid[idx] = snapshot.particles.types.index(mol.types[j])
        idx += 1
        
# print "BODIES ",snapshot.particles.body
# print "1 ",snapshot.particles.typeid
# print "2 ",snapshot.particles.types
# print "3 ",snapshot.particles.position
# for i, cp in enumerate(centralParticles):
#     print cp.m_positions
# sys.exit()
        
#writeXyz('foo.xyz', snapshot.particles.position, [snapshot.particles.types[t] for t in snapshot.particles.typeid])
    
system = hoomd.init.read_snapshot(snapshot)
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=5.0, nlist=nl)
lj.pair_coeff.set('C', 'C', epsilon=0.0968, sigma=3.4)
lj.pair_coeff.set('C', 'H', epsilon=0.1106, sigma=3.7736)
lj.pair_coeff.set('H', 'H', epsilon=0.1106, sigma=1.7736)
# Ignore all interactions with the central particles
for atype, btype in itertools.combinations_with_replacement(particle_types, 2):
    if atype in exclusions or btype in exclusions:
        epsilon = 0.0
        sigma = 1.0
        lj.pair_coeff.set(atype, btype, epsilon=epsilon, sigma=sigma)

nl.reset_exclusions(exclusions=['bond', '1-3', '1-4', 'angle', 'dihedral', 'body'])

# Set up the rigid bodies
rigid = hoomd.md.constrain.rigid()
for cp in centralParticles:
    rigid.set_param(cp.type,
                    positions=cp.m_positions,
                    types=cp.m_types)
rigid.validate_bodies()

#group_all = hoomd.group.all()
group_all = hoomd.group.rigid_center()

dt=0.0001
Nmin=5
alpha_start=0.1
ftol=1e-2
Etol=1e-5
finc=1.1
fdec=0.5
fire = hoomd.md.integrate.mode_minimize_fire(dt=dt,
                                             Nmin=Nmin,
                                             alpha_start=alpha_start,
                                             ftol=ftol,
                                             Etol=Etol,
                                             finc=finc,
                                             fdec=fdec)

integrate_nve = hoomd.md.integrate.nve(group=group_all)

dgsd = hoomd.dump.gsd(filename="opt.gsd",
                      period=1,
                      group=group_all,
                      overwrite=True)

opt_cycles=1000
hoomd.run(opt_cycles,
          callback=lambda x:-1 if fire.has_converged() else 0,
          callback_period=1)
 
