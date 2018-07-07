#!/usr/bin/env python
import sys
import numpy as np
import hoomd
import hoomd.md

class RigidParticle(object):
    def __init__(self, molecule):
        # Attributes of the central particle
        self.mass = None
        self.orientation = np.array([1.0, 0.0, 0.0, 0.0])
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


def centreOfMass(positions, masses):
    totalMass = np.sum(masses)
    return np.sum(positions * masses[:,np.newaxis], axis=0) / totalMass


def momentOfInertia(positions, masses):
    """Moment of Inertia Tensor"""
    positions = positions - centreOfMass(positions, masses) # positions relative to the centre of mass
    x = 0
    y = 1
    z = 2
    I = np.zeros(shape=(3, 3))
    I[x, x] = np.sum((np.square(positions[:, y]) + np.square(positions[:, z])) * masses)
    I[y, y] = np.sum((np.square(positions[:, x]) + np.square(positions[:, z])) * masses)
    I[z, z] = np.sum((np.square(positions[:, x]) + np.square(positions[:, y])) * masses)
    I[x, y] = np.sum(positions[:, x] * positions[:, y] * masses)
    I[y, x] = I[x, y]
    I[y, z] = np.sum(positions[:, y] * positions[:, z] * masses)
    I[z, y] = I[y, z]
    I[x, z] = np.sum(positions[:, x] * positions[:, z] * masses)
    I[z, x] = I[x, z]
    return I


def quaternion_from_matrix(M):
    """Return the quaternion of the rotation matrix
    """
    tr = np.trace(M)
    if (tr > 0):
        S = np.sqrt(tr + 1.0) * 2 # S=4*qw 
        qw = 0.25 * S
        qx = (M[2, 1] - M[1, 2]) / S
        qy = (M[0, 2] - M[2, 0]) / S
        qz = (M[1, 0] - M[0, 1]) / S
    elif ((M[0, 0] > M[1, 1])&(M[0, 0] > M[2, 2])):
        S = np.sqrt(1.0 + M[0, 0] - M[1, 1] - M[2, 2]) * 2 # S=4*qx 
        qw = (M[2, 1] - M[1, 2]) / S
        qx = 0.25 * S
        qy = (M[0, 1] + M[1, 0]) / S 
        qz = (M[0, 2] + M[2, 0]) / S 
    elif (M[1, 1] > M[2, 2]):
        S = np.sqrt(1.0 + M[1, 1] - M[0, 0] - M[2, 2]) * 2 # S=4*qy
        qw = (M[0, 2] - M[2, 0]) / S
        qx = (M[0, 1] + M[1, 0]) / S 
        qy = 0.25 * S
        qz = (M[1, 2] + M[2, 1]) / S 
    else:
        S = np.sqrt(1.0 + M[2, 2] - M[0, 0] - M[1, 1]) * 2 # S=4*qz
        qw = (M[1, 0] - M[0, 1]) / S
        qx = (M[0, 2] + M[2, 0]) / S
        qy = (M[1, 2] + M[2, 1]) / S
        qz = 0.25 * S
    return np.array([qw, qx, qy, qz])


def orientationQuaternion(A, B):
    """Return the quaternion to rotate A to B"""
    M = rigid_rotate(A, B)
    return quaternion_from_matrix(M)


def principalMoments(positions, masses):
    I = momentOfInertia(positions, masses)
    eigval, eigvec = np.linalg.eig(I)
    return np.sort(eigval)


def rigid_rotate(A, B):
    """Return rotation matrix to rotate A to B. Assumes both sets of points are already centred.
    Uses SVD to calculate the rotation.
    """
    A = np.matrix(A)
    B = np.matrix(B)
    assert len(A) == len(B)
    H = A.T * B
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T * U.T
    # special reflection case
    if np.linalg.det(R) < 0:
        # Reflection detected
        Vt[2,:] *= -1
        R = Vt.T * U.T
    B_r = np.dot(A, R.T) # Check that rotating A by the matrix gives us B again
    assert np.allclose(B_r, B, atol=0.0001), "{}\n{}".format(B_r, B)
    return R


BOX_WIDTH = 20.0
ORIENTED_PARTICLES = False

# Create the two molecules that will be separated by a single bond
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

molecules = [mol1, mol2]

# Create the central particles
RIGID_TYPE = 'A'
rigidParticles = []
ref_orientation = None
for i, mol in enumerate(molecules):
    rp = mol.rigidParticle()
    if ORIENTED_PARTICLES:
        if i == 0:
            # Use first molecule as reference orientation
            ref_orientation = rp.m_positions
        rp.type = RIGID_TYPE
        rp.orientation = orientationQuaternion(ref_orientation, rp.m_positions)
    else:
        rp.type = "CP%d" % i
    rigidParticles.append(rp)

# Get total number of particles and the set of types
particle_types = set()
nparticles = 0
for rp in rigidParticles:
    nparticles += 1
    particle_types.add(rp.type)
for mol in molecules:
    nparticles += mol.positions.shape[0]
    particle_types.update(set(mol.types))
particle_types = list(particle_types)
#
# Start of HOOMD-blue code
#
hoomd.context.initialize()
lx = BOX_WIDTH
ly = BOX_WIDTH
lz = BOX_WIDTH
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
    snapshot.particles.orientation[i] = rp.orientation
    
# Then add in the constituent molecule particles
idx = i + 1
for i, rp in enumerate(rigidParticles):
    for j in range(rp.m_positions.shape[0]):
        snapshot.particles.body[idx] = i
        snapshot.particles.mass[i] = rp.m_masses[i] # Only need to display consituent particles in gsd file
        snapshot.particles.position[idx] = rp.m_positions[j]
        snapshot.particles.typeid[idx] = snapshot.particles.types.index(rp.m_types[j])
        idx += 1
        
#writeXyz('foo2.xyz', snapshot.particles.position, [snapshot.particles.types[t] for t in snapshot.particles.typeid])

# Create bond between two C-atoms of first 2 molecules
snapshot.bonds.resize(1)
b0 = 2
b1 = 6
snapshot.bonds.group[0] = [b0, b1]
snapshot.bonds.typeid[0] = 0

system = hoomd.init.read_snapshot(snapshot)

# bond potential
bond_harmonic = hoomd.md.bond.harmonic(name="bond_harmonic")
bond_harmonic.bond_coeff.set(bond_types[0], k=606.2, r0=1.535)

# pair potentials
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=5.0, nlist=nl)
lj.pair_coeff.set('C', 'C', epsilon=0.0968, sigma=3.4)
lj.pair_coeff.set('C', 'H', epsilon=0.1106, sigma=3.7736)
lj.pair_coeff.set('H', 'H', epsilon=0.1106, sigma=1.7736)
# Ignore all interactions with the central particles
if ORIENTED_PARTICLES:
    lj.pair_coeff.set(RIGID_TYPE, [RIGID_TYPE, 'C', 'H'], epsilon=0.0, sigma=0.0, r_cut=False)
else:
    lj.pair_coeff.set('CP0', ['CP0', 'CP1', 'C', 'H'], epsilon=0.0, sigma=0.0, r_cut=False)
    lj.pair_coeff.set('CP1', ['CP1', 'CP0', 'C', 'H'], epsilon=0.0, sigma=0.0, r_cut=False)
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
# Dump out the trajectory
dgsd = hoomd.dump.gsd(filename="trajectory.gsd",
                      period=1,
                      group=groupAll,
                      overwrite=True)

# Run the MD
ncycles = 1000
T = 1.0
tau = 0.5
dt = 0.0005
integrator_mode = hoomd.md.integrate.mode_standard(dt=dt)
integ = hoomd.md.integrate.nvt(group=groupRigid, kT=T, tau=tau)
hoomd.run(ncycles)
