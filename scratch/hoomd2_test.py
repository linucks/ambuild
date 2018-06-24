#!/usr/bin/env python
import sys
sys.path.insert(0, '../ambuild')
from xyz_util import writeXyz
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

#setup()
#sys.exit(1)

class CentralParticle(object):
    def __init__(self, molecule):
        self.position = None
        self.principalMoments = None
        self.type = None
        self.fromMolecule(molecule)
    
    def fromMolecule(self, molecule):
        self.position = centreOfMass(molecule.positions, molecule.masses)
        self.principalMoments = principalMoments(molecule.positions, molecule.masses)
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
mol1.types = np.array(['C', 'H', 'H', 'H', 'H'])
mol1.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])

mol2 = Molecule()
mol2.positions = np.array([[  3.53000000e+00,   0.00000000e+00,  -2.00000001e-07],
                           [  2.44100000e+00,   0.00000000e+00,  -2.00000001e-07],
                           [  3.89300000e+00,   0.00000000e+00,   1.02671880e+00],
                           [  3.89300000e+00,  -8.89165000e-01,  -5.13360200e-01],
                           [  3.89300000e+00,   8.89165000e-01,  -5.13360200e-01]])
mol2.types = np.array(['C', 'H', 'H', 'H', 'H'])
mol2.masses = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])

molecules = [mol1, mol2]

# Now create the central particles
centralParticles = []
for i, mol in enumerate(molecules):
    cp = mol.centralParticle()
    cp.type = "CP%d" % i
    centralParticles.append(cp)
    
for cp in centralParticles:
    print cp.type
    
sys.exit()
    


centralParticle1 = centreOfMass(positions1, masses1)
# Get positions relative to central particle
positions1 = positions1 - centralParticle1
prinMom1 = principalMoments(coords1, masses1)


centralParticle2 = centreOfMass(positions2, masses2)
# Get positions relative to central particle
positions2 = positions2 - centralParticle2
prinMom2 = principalMoments(coords2, masses2)


types = np.concatenate([types1, types2])
positions = np.concatenate([positions1, positions2])
masses = np.concatenate([masses1, masses2])
bodies = np.concatenate([bodies1, bodies2])

boxdim = np.array([BOX_WIDTH, BOX_WIDTH, BOX_WIDTH])


writeXyz('foo.xyz', positions, types)
sys.exit()



# Start of hoomd code
hoomd.context.initialize()

nparticles = positions.shape[0]
box_width = 20.0
lx = box_width
ly = box_width
lz = box_width
particle_types = list(set((types)))
snapshot = hoomd.data.make_snapshot(N=nparticles,
                                    box=hoomd.data.boxdim(Lx=lx, Ly=ly, Lz=lz),
                                    particle_types=particle_types)

#for i, p in enumerate(snapshot.particles):
#snap.particles.moment_inertia[i] = data.rigid_moment_inertia[i]
for i in range(nparticles):
    snapshot.particles.mass[i] = masses[i]
    snapshot.particles.position[i] = positions[i]
    snapshot.particles.body[i] = bodies[i]
    snapshot.particles.typeid[i] = snapshot.particles.types.index(types[i])
    
    
    
rigid = hoomd.md.constrain.rigid()
for ftype, fdata in data.rigid_fragments.iteritems():
    rigid.set_param(ftype,
                    types='TYPES',
                    positions='FOO')
    rigid.validate_bodies()

    
system = hoomd.init.read_snapshot(snapshot)
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=5.0, nlist=nl)
lj.pair_coeff.set('C', 'C', epsilon=0.0968, sigma=3.4)
lj.pair_coeff.set('C', 'H', epsilon=0.1106, sigma=3.7736)
lj.pair_coeff.set('H', 'H', epsilon=0.1106, sigma=1.7736)

nl.reset_exclusions(exclusions=['bond', '1-3', '1-4', 'angle', 'dihedral', 'body'])

group_all = hoomd.group.all()

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
 
