#!/Users/jmht/miniconda2/envs/hoomd2/bin/python
# source activate hoomd2

import itertools
import logging
import sys

# 3rd-party imports
import hoomd
import hoomd.md
import numpy as np

# Our imports
from ab_ffield import FfieldParameters
import xyz_core

logger = logging.getLogger(__name__)


class Hoomd2(object):
    """
    TODO in 2
    * fix masked atoms
    * handle reading data back into cell - need to exclude central particles

    TODO in 1
    * change checkParameters to use self.particleTypes etc
    * change setBonds etc to match hoomd2 so creation of oboject in routine


    Each rigid body in the cell is a type of rigid body - define fragmentType:bodyCount
    - calc centre for the body -> need coords and atomTypes for each body in a fragment

    """
    def __init__(self, paramsDir):
        self.ffield = FfieldParameters(paramsDir)
        self.debug = False
        self.rCut = 5.0
        self.system = None
        self.rigidBody = False
        self.exclusions = set()  # particle tags to be ignored in pair-pair interactions

    def checkParameters(self, skipDihedrals=False):
        assert self.ffield
        assert self.particleTypes
        ok = True
        missingBonds = []
        for bond in self.bond_types:
            if not self.ffield.hasBond(bond):
                ok = False
                missingBonds.append(bond)
        missingAngles = []
        for angle in self.angle_types:
            if not self.ffield.hasAngle(angle):
                ok = False
                missingAngles.append(angle)
        missingDihedrals = []
        missingImpropers = []
        if not skipDihedrals:
            for dihedral in self.dihedral_types:
                if not self.ffield.hasDihedral(dihedral):
                    ok = False
                    missingDihedrals.append(dihedral)
#             for improper in self.impropers:
#                 if not self.ffield.hasImproper(improper):
#                     ok = False
#                     missingImpropers.append(improper)
        missingPairs = []
        for atype, btype in itertools.combinations_with_replacement(self.particleTypes, 2):
            if atype in self.exclusions or btype in self.exclusions:
                continue
            if not self.ffield.hasPair(atype, btype):
                ok = False
                missingPairs.append((atype, btype))
        if not ok:
            msg = "The following parameters could not be found:\n"
            if missingBonds:
                msg += "Bonds: {0}\n".format(missingBonds)
            if missingAngles:
                msg += "Angles: {0}\n".format(missingAngles)
            if missingDihedrals:
                msg += "Dihedrals: {0}\n".format(missingDihedrals)
            if missingImpropers:
                msg += "Impropers: {0}\n".format(missingImpropers)
            if missingPairs:
                msg += "Pairs: {0}\n".format(missingPairs)
            msg += "Please add these to the files in the directory: {0}\n".format(self.ffield.paramsDir)
            raise RuntimeError(msg)
        return

    def createSnapshot(self, data, doCharges=True, doDihedral=True):
        """Create a populated snapshot with all the particles.

        Rigid body will require creating the central particles so we'll do this later
        """
        # Reset exclusions here for time being
        self.exclusions = set()
        # Create snapshot
        # snap attributes: angles, bonds, box, constraints, dihedrals, impropers, pairs, particles
        if self.rigidBody:
            # Combined size includes rigid centre and constituent particles
            nRigidParticles = len(data.rigidParticles)
            nparticles = data.natoms + nRigidParticles
            rigidCenters = set([r.type for r in data.rigidParticles])
            atomTypes = set()
            for r in data.rigidParticles:
                atomTypes.update(r.b_atomTypes)
            self.particleTypes = list(atomTypes.union(rigidCenters))
            self.exclusions = set(rigidCenters)
        else:
            nparticles = len(data.coords)
            self.particleTypes = list(set(data.atomTypes))

        assert nparticles > 0, "Simulation needs some particles!"
        # NEED TO THINK ABOUT WHAT TO DO ABOUT MASKED ATOMS - set particleTypes?
        # self.masked = data.masked
        self.bond_types = list(set(data.bondLabels)) if len(data.bonds) else []
        self.angle_types = list(set(data.angleLabels)) if len(data.angles) else []
        self.dihedral_types = list(set(data.properLabels)) if len(data.propers) and doDihedral else []
        snapshot = hoomd.data.make_snapshot(N=nparticles,
                                            box=hoomd.data.boxdim(Lx=data.cell[0], Ly=data.cell[1], Lz=data.cell[2]),
                                            particle_types=self.particleTypes,
                                            bond_types=self.bond_types,
                                            angle_types=self.angle_types,
                                            dihedral_types=self.dihedral_types)
        # Add Bonds
        if len(self.bond_types):
            snapshot.bonds.resize(len(data.bonds))
            for i, b in enumerate(data.bonds):
                if self.rigidBody:  # center particles are at the front of the arrays, so everything gets shifted up
                    b0 = b[0] + nRigidParticles
                    b1 = b[1] + nRigidParticles
                else:
                    b0, b1 = b
                snapshot.bonds.group[i] = [b0, b1]
                snapshot.bonds.typeid[i] = self.bond_types.index(data.bondLabels[i])
        # Add Angles
        if len(self.angle_types):
            snapshot.angles.resize(len(data.angles))
            for i, a in enumerate(data.angles):
                if self.rigidBody:  # center particles are at the front of the arrays, so everything gets shifted up
                    a0 = a[0] + nRigidParticles
                    a1 = a[1] + nRigidParticles
                    a2 = a[2] + nRigidParticles
                else:
                    a0, a1, a2 = a
                snapshot.angles.group[i] = [a0, a1, a2]
                snapshot.angles.typeid[i] = self.angle_types.index(data.angleLabels[i])
        # Add Dihedrals
        if doDihedral and len(self.dihedral_types):
            snapshot.dihedrals.resize(len(data.propers))
            for i, d in enumerate(data.propers):
                if self.rigidBody:  # center particles are at the front of the arrays, so everything gets shifted up
                    d0 = d[0] + nRigidParticles
                    d1 = d[1] + nRigidParticles
                    d2 = d[2] + nRigidParticles
                    d3 = d[3] + nRigidParticles
                else:
                    d0, d1, d2, d3 = d
                snapshot.dihedrals.group[i] = [d0, d1, d2, d3]
                snapshot.dihedrals.typeid[i] = self.dihedral_types.index(data.properLabels[i])
        # Populate  particle data
        if self.rigidBody:
            # Central particles need to be first in the list before any constituent particles.
            for i, rp in enumerate(data.rigidParticles):
                snapshot.particles.body[i] = i
                # Wrap central paticle into the cell. We use the image of this for all the constituent particles
                position, rp_image = xyz_core.wrapCoord3(rp.position, dim=data.cell, center=True)
                snapshot.particles.position[i] = position
                snapshot.particles.image[i] = rp_image
                snapshot.particles.mass[i] = rp.mass
                snapshot.particles.orientation[i] = rp.orientation
                snapshot.particles.typeid[i] = snapshot.particles.types.index(rp.type)
                snapshot.particles.moment_inertia[i] = rp.principalMoments
            # Then add in the constituent molecule particles
            idx = i + 1
            for i, rp in enumerate(data.rigidParticles):
                for j in range(rp.natoms):
                    snapshot.particles.body[idx] = i # to match central particle
                    if doCharges:
                        snapshot.particles.charge[idx] = rp.b_charges[j]
                    snapshot.particles.diameter[idx] = rp.b_diameters[j]
                    snapshot.particles.image[idx] = rp_image
                    snapshot.particles.mass[idx] = rp.b_masses[j]
                    snapshot.particles.position[idx] = rp.b_positions[j]
                    snapshot.particles.typeid[idx] = snapshot.particles.types.index(rp.b_atomTypes[j])
                    idx += 1
        else:
            for i in range(nparticles):
                if doCharges:
                    snapshot.particles.charge[i] = data.charges[i]
                snapshot.particles.diameter[i] = data.diameters[i]
                snapshot.particles.image[i] = data.images[i]
                snapshot.particles.mass[i] = data.masses[i]
                snapshot.particles.position[i] = data.coords[i]
                snapshot.particles.typeid[i] = self.particleTypes.index(data.atomTypes[i])
        return snapshot

    def optimiseGeometry(self,
                          data,
                          rigidBody=True,
                          doDihedral=False,
                          doImproper=False,
                          doCharges=True,
                          rCut=None,
                          quiet=None,
                          walls=None,
                          wallAtomType=None,
                          **kw):
        if rCut is not None:
            self.rCut = rCut
        if doDihedral and doImproper:
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.rigidBody = rigidBody
        self.setupContext(quiet=quiet)
        snapshot = self.createSnapshot(data, doCharges=doCharges, doDihedral=doDihedral)
        self.setupSimulation(snapshot, data, walls=walls, wallAtomType=wallAtomType)
        hlog = self._createLog('geomopt.tsv')
        optimised = self._optimiseGeometry(**kw)
        # Extract the energy
        if 'd' in kw and kw['d'] is not None:
            for i in ['potential_energy' ]:
                kw['d'][ i ] = hlog.query(i)
        return True
        return optimised

    def _optimiseGeometry(self,
                          optCycles=1000000,
                          dump=False,
                          dumpPeriod=1,
                          dt=0.005,
                          Nmin=5,
                          alpha_start=0.1,
                          ftol=1e-2,
                          Etol=1e-5,
                          finc=1.1,
                          fdec=0.5,
                          max_tries=3,
                          retries_on_error=3,
                          **kw):
        """Optimise the geometry with hoomdblue"""

        assert max_tries > 0
        if retries_on_error < 1:
            retries_on_error = 0
        retries_on_error += 1

        for i in range(retries_on_error):
            # Try optimising and lower the timestep and increasing the number of cycles each time
            try:
                # Create the integrator with the values specified
                fire = hoomd.md.integrate.mode_minimize_fire(dt=dt,
                                                             Nmin=Nmin,
                                                             alpha_start=alpha_start,
                                                             ftol=ftol,
                                                             Etol=Etol,
                                                             finc=finc,
                                                             fdec=fdec)
                integrate_nve = hoomd.md.integrate.nve(group=self.groupActive)
                if dump:
                    dgsd = hoomd.dump.gsd(filename="opt.gsd",
                                          period=dumpPeriod,
                                          group=self.groupAll,
                                          overwrite=True)
                optimised = False
                for j in range(max_tries):
                    logger.info("Running {0} optimisation cycles in macrocycle {1}".format(optCycles, j))
                    hoomd.run(optCycles,
                              callback=lambda x:-1 if fire.has_converged() else 0,
                              callback_period=1)
                    if fire.has_converged():
                        logger.info("Optimisation converged on macrocycle {0}".format(j))
                        optimised = True
                        break
                    else:
                        logger.info("Optimisation failed to converge on macrocycle {0}".format(j))
                        if j + 1 < max_tries:
                            logger.info("Attempting another optimisation macrocycle")
                break  # Break out of try/except loop
            except RuntimeError as e:
                logger.info("Optimisation step {0} failed!\n{1}".format(i, e))
                fire.reset()
                integrate_nve.disable()
                if i + 1 < retries_on_error:
                    dt_old = dt
                    dt = dt_old * 0.1
                    logger.info("Rerunning optimisation changing dt {0} -> {1}".format(dt_old, dt))
                else:
                    # If we're not going to try again we re-raise the last exception that we caught
                    raise
        # Delete variables before we return to stop memory leaks
        # del fire
        if dump and False:
            dgsd.disable()
            del dgsd
        return optimised

    def runMD(self,
              data,
              doDihedral=False,
              doImproper=False,
              doCharges=True,
              rigidBody=False,
              rCut=None,
              quiet=None,
              walls=None,
              wallAtomType=None,
              **kw):
        if rCut is not None:
            self.rCut = rCut
        if doDihedral and doImproper:
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.rigidBody = rigidBody
        self.setupContext(quiet=quiet)
        snapshot = self.createSnapshot(data, doCharges=doCharges, doDihedral=doDihedral)
        self.setupSimulation(snapshot, data, walls=walls, wallAtomType=wallAtomType)
        hlog = self._createLog('runmd.log')
        self._runMD(**kw)
        # Extract the energy
        if 'd' in kw and kw['d'] is not None:
            for i in ['potential_energy' ]:
                kw['d'][ i ] = hlog.query(i)
        return True

    def _runMD(self,
               mdCycles=100000,
               integrator='nvt',
               T=1.0,
               tau=0.5,
               P=1,
               tauP=0.5,
               dt=0.0005,
               dump=False,
               dumpPeriod=100,
               **kw):
        # Added **kw arguments so that we don't get confused by arguments intended for the optimise
        # when MD and optimiser run together
        integrator_mode = hoomd.md.integrate.mode_standard(dt=dt)
        if integrator == 'nvt':
            integ = hoomd.md.integrate.nvt(group=self.groupActive, kT=T, tau=tau)
        elif integrator == 'npt':
            integ = hoomd.md.integrate.npt(group=self.groupActive, kT=T, tau=tau, P=P, tauP=tauP)
        else:
            raise RuntimeError("Unrecognised integrator: {0}".format(integrator))

        if dump:
            dgsd = hoomd.dump.gsd(filename="runmd.gsd",
                                  period=dumpPeriod,
                                  group=self.groupAll,
                                  overwrite=True)
        # run mdCycles time steps
        hoomd.run(mdCycles)
        integ.disable()
        if dump:
            dgsd.disable()
            del dgsd
        del integ
        del integrator_mode
        return

    def setBonds(self):
        if not self.bond_types:
            return
        bond_harmonic = hoomd.md.bond.harmonic(name="bond_harmonic")
        for bond in self.bond_types:
            param = self.ffield.bondParameter(bond)
            if self.debug:
                logger.info("DEBUG: bond_harmonic.bond_coeff.set( '{0}',  k={1}, r0={2} )".format(bond,
                                                                                                  param['k'],
                                                                                                  param['r0']))
            bond_harmonic.bond_coeff.set(bond, k=param['k'], r0=param['r0'])
        return

    def setAngles(self):
        if not self.angle_types: return
        angle_harmonic = hoomd.md.angle.harmonic()
        for angle in self.angle_types:
            param = self.ffield.angleParameter(angle)
            if self.debug:
                logger.info("DEBUG: angle_harmonic.angle_coeff.set( '{0}',  k={1}, t0={2} )".format(angle,
                                                                                                    param['k'],
                                                                                                    param['t0']))
            angle_harmonic.angle_coeff.set(angle, k=param['k'], t0=param['t0'])
        return

    def setDihedrals(self):
        if not self.dihedral_types: return
        dihedral_harmonic = hoomd.md.dihedral.harmonic()
        for dihedral in self.dihedral_types:
            param = self.ffield.dihedralParameter(dihedral)
            if self.debug:
                logger.info("DEBUG: dihedral_harmonic.dihedral_coeff.set('{0}',  k={1}, d={2}, n={3})".format(dihedral,
                                                                                                              param['k'],
                                                                                                              param['d'], param['n']))
            dihedral_harmonic.dihedral_coeff.set(dihedral, k=param['k'], d=param['d'], n=param['n'])
        return

    def setPairs(self):
        nl = hoomd.md.nlist.cell()
        lj = hoomd.md.pair.lj(r_cut=self.rCut, nlist=nl)
        for atype, btype in itertools.combinations_with_replacement(self.particleTypes, 2):
#             # dummy atoms treated specially
#             #if atype.lower() in ['x','hn'] or btype.lower() in ['x','hn']:
#             if self.masked[i] or self.masked[j]:
#                 epsilon = 0.0
#                 sigma = 0.0
#             else:
            if atype in self.exclusions or btype in self.exclusions:
                epsilon = 0.0
                sigma = 1.0
            else:
                param = self.ffield.pairParameter(atype, btype)
                epsilon = param['epsilon']
                sigma = param['sigma']
            lj.pair_coeff.set(atype, btype, epsilon=epsilon, sigma=sigma)
            if self.debug:
                logger.info("DEBUG: lj.pair_coeff.set('{0}', '{1}', epsilon={2}, sigma={3} )".format(atype,
                                                                                                     btype,
                                                                                                     epsilon,
                                                                                                     sigma))        
        # Don't think we need to include body any more for rigid bodies, as these are already excluded by default?
        # nl.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle', 'body'])
        # nl.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle'])
        nl.reset_exclusions(exclusions=['bond', '1-3', '1-4', 'angle', 'dihedral', 'body'])
        return

    def setupContext(self, quiet=False):
        hoomd.context.initialize()
        if quiet:
            hoomd.util.quiet_status()
            hoomd.option.set_notice_level(0)
        return

    def setupGroups(self, data):
        self.groupAll = hoomd.group.all()
        self.groupRigid = hoomd.group.rigid_center()
        # All atoms that are part of static fragments
        if any(data.static):
            self.groupStatic = hoomd.group.tag_list(name="static",
                                                      tags=[i for i, s in enumerate(data.static) if s])
            if self.rigidBody:
                self.groupActive = hoomd.group.difference(name="active",
                                                            a=self.groupRigid,
                                                            b=self.groupStatic)
            else:
                self.groupActive = hoomd.group.difference(name="active",
                                                            a=self.groupAll,
                                                            b=self.groupStatic)
        else:
            if self.rigidBody:
                self.groupActive = self.groupRigid
            else:
                self.groupActive = self.groupAll
        return

    def setupSimulation(self, snapshot, data, walls=None, wallAtomType=None):
        self.system = hoomd.init.read_snapshot(snapshot)
        self.setupRigidBody(data)
        self.checkParameters()
        self.setBonds()
        self.setAngles()
        self.setDihedrals()
        self.setPairs()
        self.setupGroups(data)
        self.setupWalls(walls, wallAtomType)
        return

    def setupRigidBody(self, data):
        if not self.rigidBody:
            return
        rigid = hoomd.md.constrain.rigid()
#         for ftype, fdata in data.rigid_fragments.iteritems():
#             rigid.set_param(ftype,
#                             types=fdata['atomTypes'],
#                             positions=data.coords[fdata['coord_idxs'][0] : fdata['coord_idxs'][1]])
        for rp in data.rigidParticles:
            rigid.set_param(rp.type,
                            positions=rp.b_positions,
                            types=rp.b_atomTypes)
        rigid.validate_bodies()
        return

    def setupWalls(self, walls, wallAtomType):
        """Set up walls for the simulation

        I think that we require two walls. One on the front side with the potential facing in,
        the other on the back wall with the potential facing back towards the other potential.
        The origin of the wall is the centre of plane but then back half a cell along the axis
        that isn't part of the wall.
        """
        if walls is None:
            return
        assert len(walls) is 3  # array of three booleans - one per wall
        if not any(walls):
            return
        assert wallAtomType is not None, "Need to set a wallAtomType!"
        wtypes = ['XOY', 'XOZ', 'YOZ']  # Order of walls
        wallstructure = None
        for do, wtype in zip(walls, wtypes):
            if do:
                logger.info('Setting wall in HOOMDBLUE for {0} of atomType {1}'.format(wtype, wallAtomType))
                if wtype is 'XOY':
                    originFront = (0, 0, -self.system.box.Lz / 2)
                    originBack = (0, 0, self.system.box.Lz / 2)
                    normal = (0, 0, 1)
                elif wtype is 'XOZ':
                    originFront = (0, -self.system.box.Ly / 2, 0)
                    originBack = (0, self.system.box.Ly / 2, 0)
                    normal = (0, 1, 0)
                elif wtype is 'YOZ':
                    originFront = (-self.system.box.Lx / 2, 0, 0)
                    originBack = (self.system.box.Lx / 2, 0, 0)
                    normal = (1, 0, 0)
                else:
                    raise RuntimeError("Unrecognised Wall Type! {0}".format(wtype))
                if not wallstructure:
                    # We only create the wall and the LJ potentials once as they are used
                    # by all subsequent walls in the group
                    try:
                        wallstructure = hoomd.md.wall.group()
                    except AttributeError:
                        raise RuntimeError('HOOMD-blue wall does not have a group attribute. You may need to update your version of HOOMD-Blue in order to use walls')
                    lj = hoomd.md.wall.lj(wallstructure, r_cut=self.rCut)
                    for atype in self.particleTypes:
                        param = self.ffield.pairParameter(atype, wallAtomType)
                        lj.force_coeff.set(atype, epsilon=param['epsilon'], sigma=param['sigma'])

                # Add the two walls
                # Front
                wallstructure.add_plane(origin=originFront, normal=normal, inside=True)
                # Back
                wallstructure.add_plane(origin=originBack, normal=normal, inside=False)
        return

    def updateCell(self, cell):
        """Reset the particle positions from hoomdblue system"""
        box = np.array([self.system.box.Lx, self.system.box.Ly, self.system.box.Lz])
        snapshot = self.system.take_snapshot()
        if self.rigidBody:
            atomIdx = len(hoomd.group.rigid_center())
        else:
            atomIdx = 0
        for block in cell.blocks.itervalues():
            for i in range(block.numAtoms()):
                coord = snapshot.particles.position[atomIdx]
                coord = xyz_core.unWrapCoord3(coord,
                                              snapshot.particles.image[atomIdx],
                                              box,
                                              centered=True)
                block.coord(i, coord)
                atomIdx += 1
        if atomIdx != snapshot.particles.N:
            raise RuntimeError("Read {0} positions but there were {1} particles!".format(atomIdx, len(self.system.particles)))

        # If we are running (e.g.) an NPT simulation, the cell size may have changed. In this case we need to update
        # our cell parameters. Repopulate cells will then update the halo cells and add the new blocks
        if not np.allclose(box, cell.dim):
            logger.info("Changing cell dimensions after HOOMD-blue simulation from: {0} to: {1}".format(cell.dim, box))
            cell.dim = box
        # Now have the new coordinates, so we need to put the atoms in their new cells
        cell.repopulateCells()
        return

    def _createLog(self, filename):
        return hoomd.analyze.log(filename=filename,
                                     quantities=[ 'num_particles',
                                                 'pair_lj_energy',
                                                 'potential_energy',
                                                 'kinetic_energy',
                                                 ],
                                     period=100,
                                     header_prefix='#',
                                     overwrite=True)


def snap2xyz(snapshot, fpath='foo.xyz'):
    with open(fpath, 'w') as w:
        w.write("%d\n" % snapshot.particles.N)
        w.write("JENS\n")
        types = snapshot.particles.types
        for i in range(snapshot.particles.N):
            w.write("%s    %f    %f    %f\n" % (types[snapshot.particles.typeid[i]],
                                                snapshot.particles.position[i][0],
                                                snapshot.particles.position[i][1],
                                                snapshot.particles.position[i][2]))
    print("JMHT WROTE SNAPSHOT: %s" % fpath)
    return

if __name__ == "__main__":
    from ab_paths import PARAMS_DIR
    import ab_util
    mycell = ab_util.cellFromPickle(sys.argv[1])
    rigidBody = True
    data = mycell.dataDict(periodic=True, center=True, rigidBody=rigidBody)
    hmd = Hoomd2(PARAMS_DIR)
    if False:
        ok = hmd.optimiseGeometry(data,
                                  rigidBody=rigidBody,
                                  doDihedral=True,
                                  doImproper=False,
                                  doCharges=True,
                                  dump=True,
                                  optCycles=1000,
                                  max_tries=1,
                                  retries_on_error=0)
    else:
        ok = hmd.runMD(data,
                       rigidBody=rigidBody,
                       doDihedral=True,
                       doImproper=False,
                       doCharges=True,
                       dump=True)
