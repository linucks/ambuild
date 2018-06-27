#!/usr/bin/env python
import sys
sys.path.insert(0, '../ambuild')
from xyz_util import writeXyz
import ab_bond
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
    
    if False:
        b = list(b1.fragments[0].bodies())[0]
        print b.symbols()
        print b.coords()
        b = list(b2.fragments[0].bodies())[0]
        print b.symbols()
        print b.coords()
    else:
        b3 = b2.copy() # Copy before bonding so can move away from bonded pair 
        bond = ab_bond.Bond(endGroup1, endGroup2)
        bond.engage()
        b3.translate([3.0, 0.0, 0.0])
    
    for f in b1.fragments:
        for b in f.bodies():
            print b.symbols()
            print b.coords()
    for f in b3.fragments:
        for b in f.bodies():
            print b.symbols()
            print b.coords()


def wrapBox(positions, boxdim):
    # Move into centre of cell
    positions = np.fmod(positions, boxdim)
    # Change the coord so the origin is at the center of the box (we start from the corner)
    positions -= boxdim / 2
    return positions

#boxdim = np.array([BOX_WIDTH, BOX_WIDTH, BOX_WIDTH])
# setup()
# sys.exit(1)

class RigidParticle(object):
    def __init__(self, molecule):
        # Attributes of the central particle
        self.mass = None
        self.position = None
        self.principalMoments = None
        self.type = None
        # Attributes if the constituent particles
        self.m_positions = None
        self.m_types = None
        
        self.fromMolecule(molecule)
    
    def fromMolecule(self, molecule):
        self.mass = np.sum(molecule.masses)
        self.position = centreOfMass(molecule.positions, molecule.masses)
        self.principalMoments = principalMoments(molecule.positions, molecule.masses)
        # Specify types and position of consituent particles
        # The positions of the constituent particles need to be relative to the central particle
        self.m_positions = molecule.positions - self.position
        self.m_types = molecule.types
        self.m_masses = molecule.masses
        return self

class Molecule(object):
    def __init__(self):
        self.positions = None
        self.masses = None
        self.types = None
    
    def rigidParticle(self):
        return RigidParticle(self)

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


bonded = True
if not bonded:
    # Create 2 methane molecules
    mol1 = Molecule()
    mol1.positions = np.array([[ -2.00000000e+00,   0.00000000e+00,  -2.00000001e-07],
                               [ -9.11000000e-01,   0.00000000e+00,  -2.00000001e-07],
                               [ -2.36300000e+00,   0.00000000e+00,  -1.02671920e+00],
                               [ -2.36300000e+00,  -8.89165000e-01,   5.13359800e-01],
                               [ -2.36300000e+00,   8.89165000e-01,   5.13359800e-01]])
    mol1.types = ['C', 'H', 'H', 'H', 'H']
    mol1.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])
    
    mol2 = Molecule()
    mol2.positions = np.array([[  1.53000000e+00,   0.00000000e+00,  -2.00000001e-07],
                               [  4.41000000e-01,   0.00000000e+00,  -2.00000001e-07],
                               [  1.89300000e+00,   0.00000000e+00,   1.02671880e+00],
                               [  1.89300000e+00,  -8.89165000e-01,  -5.13360200e-01],
                               [  1.89300000e+00,   8.89165000e-01,  -5.13360200e-01]])
    mol2.types = ['C', 'H', 'H', 'H', 'H']
    mol2.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])
    
    molecules = [mol1, mol2]
else:
    # 2 bonded, 1 not
    mol1 = Molecule()
    mol1.positions = np.array([[  0.00000000e+00,   0.00000000e+00,  -2.00000001e-07],
                               [ -3.63000000e-01,   0.00000000e+00,  -1.02671920e+00],
                               [ -3.63000000e-01,  -8.89165000e-01,   5.13359800e-01],
                               [ -3.63000000e-01,   8.89165000e-01,   5.13359800e-01]])
    mol1.types = ['C', 'H', 'H', 'H']
    mol1.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794])
    
    mol2 = Molecule()
    mol2.positions = np.array([[  1.53000000e+00,   0.00000000e+00,  -2.00000001e-07],
                               [  1.89300000e+00,   0.00000000e+00,   1.02671880e+00],
                               [  1.89300000e+00,  -8.89165000e-01,  -5.13360200e-01],
                               [  1.89300000e+00,   8.89165000e-01,  -5.13360200e-01]])
    mol2.types = ['C', 'H', 'H', 'H']
    mol2.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794])
    
    
    mol3 = Molecule()
    mol3.positions = np.array([[  4.53000000e+00,   0.00000000e+00,  -2.00000001e-07],
                               [  3.44100000e+00,   0.00000000e+00,  -2.00000001e-07],
                               [  4.89300000e+00,   0.00000000e+00,   1.02671880e+00],
                               [  4.89300000e+00,  -8.89165000e-01,  -5.13360200e-01],
                               [  4.89300000e+00,   8.89165000e-01,  -5.13360200e-01]])
    mol3.types = ['C', 'H', 'H', 'H', 'H']
    mol3.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])    
    
    molecules = [mol1, mol2, mol3]


# boxdim = np.array([BOX_WIDTH, BOX_WIDTH, BOX_WIDTH])
# for m in molecules:
#     m.positions = wrapBox(m.positions, boxdim)
#     print repr(m.positions)


# Create the central particles
rigidParticles = []
for i, mol in enumerate(molecules):
    rp = mol.rigidParticle()
    rp.type = "CP%d" % i
    rigidParticles.append(rp)

# Get total number of particles and the set of types
particle_types = set()
nparticles = 0
exclusions = [] # to stop the central particles ineracting directly
for rp in rigidParticles:
    nparticles += 1
    particle_types.add(rp.type)
    exclusions.append(rp.type)
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
bond_types = None
if bonded:
    bond_types = ['C-C']
snapshot = hoomd.data.make_snapshot(N=nparticles,
                                    box=hoomd.data.boxdim(Lx=lx, Ly=ly, Lz=lz),
                                    particle_types=particle_types,
                                    bond_types=bond_types)
# Add centreal particles at start of snapshot
for i, rp in enumerate(rigidParticles):
    snapshot.particles.body[i] = i
    snapshot.particles.mass[i] = rp.mass
    snapshot.particles.position[i] = rp.position
    snapshot.particles.typeid[i] = snapshot.particles.types.index(rp.type)
    snapshot.particles.moment_inertia[i] = rp.principalMoments
# Then add in the constituent molecule particles
idx = i + 1
for i, rp in enumerate(rigidParticles):
    for j in range(rp.m_positions.shape[0]):
        snapshot.particles.body[idx] = i
        snapshot.particles.mass[i] = rp.m_masses[i] # Only need to display consituent particles in gsd file
        snapshot.particles.position[idx] = rp.m_positions[j]
        snapshot.particles.typeid[idx] = snapshot.particles.types.index(rp.m_types[j])
        idx += 1
        
if bonded:
    # Create bond between two C-atoms of first 2 molecules
    snapshot.bonds.resize(1)
    b0 = 3
    b1 = 7
    snapshot.bonds.group[0] = [b0, b1]
    snapshot.bonds.typeid[0] = 0
    
# print "BODIES ",snapshot.particles.body
# print "1 ",snapshot.particles.typeid
# print "2 ",snapshot.particles.types
# print "3 ",snapshot.particles.position
# for i, rp in enumerate(rigidParticles):
#     print rp.m_positions
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

if bonded:
    bond_harmonic = hoomd.md.bond.harmonic(name="bond_harmonic")
    bond_harmonic.bond_coeff.set(bond_types[0], k=606.2, r0=1.535)

nl.reset_exclusions(exclusions=['bond', '1-3', '1-4', 'angle', 'dihedral', 'body'])

# Set up the rigid bodies
rigid = hoomd.md.constrain.rigid()
for rp in rigidParticles:
    rigid.set_param(rp.type,
                    positions=rp.m_positions,
                    types=rp.m_types)
rigid.validate_bodies()

groupAll = hoomd.group.all()
groupRigid = hoomd.group.rigid_center()

dgsd = hoomd.dump.gsd(filename="trajectory.gsd",
                      period=1,
                      group=groupAll,
                      overwrite=True)

optimise = False
ncycles = 1000
if optimise:
    dt = 0.0001
    Nmin = 5
    alpha_start = 0.1
    ftol = 1e-2
    Etol = 1e-5
    finc = 1.1
    fdec  =0.5
    fire = hoomd.md.integrate.mode_minimize_fire(dt=dt,
                                                 Nmin=Nmin,
                                                 alpha_start=alpha_start,
                                                 ftol=ftol,
                                                 Etol=Etol,
                                                 finc=finc,
                                                 fdec=fdec)
    
    integrate_nve = hoomd.md.integrate.nve(group=groupRigid)
    hoomd.run(ncycles,
              callback=lambda x:-1 if fire.has_converged() else 0,
              callback_period=1)
else:
    T = 1.0
    tau = 0.5
    dt = 0.0005
    integrator_mode = hoomd.md.integrate.mode_standard(dt=dt)
    integ = hoomd.md.integrate.nvt(group=groupRigid, kT=T, tau=tau)
    hoomd.run(ncycles)

