#!/usr/bin/env python
import sys
sys.path.insert(0, '../ambuild')
from xyz_util import writeXyz
import numpy as np
import hoomd
import hoomd.md

def setup():
    from ab_cell import Cell
    from ab_block import Block
    
    boxDim = [20.0, 20.0, 20.0]
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
box_width = 20.0
types1 = np.array(['C', 'H', 'H', 'H', 'H'])
masses1 = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])
positions1 = np.array([[  8.,         10.,          9.9999998],
                       [  9.089,      10.,          9.9999998],
                       [  7.637,      10.,          8.9732808],
                       [  7.637,       9.110835,   10.5133598],
                       [  7.637,      10.889165,   10.5133598]])


types2 = np.array(['C', 'H', 'H', 'H', 'H'])
masses2 = np.array([12.0107, 1.00794, 1.00794, 1.00794, 1.00794])
positions2 = np.array([[ 11.53,       10.,          9.9999998],
                       [ 10.441,      10.,          9.9999998],
                       [ 11.893,      10.,         11.0267188],
                       [ 11.893,       9.110835,    9.4866398],
                       [ 11.893,      10.889165,    9.4866398]])

types = np.concatenate([types1, types2])
positions = np.concatenate([positions1, positions2])
masses = np.concatenate([masses1, masses2])


boxdim = np.array([20.0, 20.0, 20.0])

# Move into centre of cell
positions = np.fmod(positions, boxdim)
# Change the coord so the origin is at the center of the box (we start from the corner)
positions -= boxdim / 2
#writeXyz('foo.xyz', positions, types)

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
for i in range(nparticles):
    snapshot.particles.mass[i] = masses[i]
    snapshot.particles.position[i] = positions[i]
    snapshot.particles.typeid[i] = snapshot.particles.types.index(types[i])
    
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
 
