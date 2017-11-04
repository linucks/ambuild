#!/Users/jmht/miniconda2/envs/hoomd2/bin/python
# source activate hoomd2

import itertools
import logging
import os
import sys


# 3rd-party imports
import hoomd
import hoomd.md

# Our imports
from opt import FFIELD, FfieldParameters
import util

logger = logging.getLogger(__name__)

class Hoomd2(object):
    """
    TODO in 2
    * fix masked atoms
    
    TODO in 1
    * change checkParameters to use self.particle_types etc
    * change setBonds etc to match hoomd2 so creation of oboject in routine
    
    """
    def __init__(self, paramsDir):
        self.ffield = FfieldParameters(paramsDir)
        self.debug = True

    def checkParameters(self, skipDihedrals=False):

        assert self.ffield
        assert self.particle_types

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
        for atype, btype in itertools.combinations(self.particle_types, 2):
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

    def createSnapshot(self, data, rigidBody=False):
        """Create a populated snapshot with all the particles.
        
        Rigid body will require creating the central particles so we'll do this later
        """
        
        if not rigidBody:
            # Create snapshot
            nparticles = len(data.coords)
            assert nparticles > 0,"Simulation needs some particles!"
            # NEED TO THINK ABOUT WHAT TO DO ABOUT MASKED ATOMS - set particle_types?
            #self.masked = data.masked

            self.particle_types = list(set(data.atomTypes))
    
            # now bond etc attributes
            # snap attributes: angles, bonds, box, constraints, dihedrals, impropers, pairs, particles
            self.bond_types = list(set(data.bondLabels)) if len(data.bonds) else None
            self.angle_types = list(set(data.angleLabels)) if len(data.angles) else None
            self.dihedral_types = list(set(data.properLabels)) if len(data.propers) else None
                
            snap = hoomd.data.make_snapshot(N=nparticles,
                                            box=hoomd.data.boxdim(Lx=data.cell[0], Ly=data.cell[1], Lz=data.cell[2]),
                                            particle_types = self.particle_types,
                                            bond_types = self.bond_types,
                                            angle_types = self.angle_types,
                                            dihedral_types = self.dihedral_types,
                                            #pair_types = pair_types
                                            )
            
            # Add Bonds
            if len(self.bond_types):
                snap.bonds.resize(len(data.bonds))
                for i, b in enumerate(data.bonds):
                    snap.bonds.group[i] = [b[0], b[1]]
                    snap.bonds.typeid[i] = self.bond_types.index(data.bondLabels[i])
                    
            # Add Angles
            if len(self.angle_types):
                snap.angles.resize(len(data.angles))
                for i, a in enumerate(data.angles):
                    snap.angles.group[i] = [a[0], a[1], a[2]]
                    snap.angles.typeid[i] = self.angle_types.index(data.angleLabels[i])
    
            # Add Dihedrals
            if len(self.dihedral_types):
                snap.dihedrals.resize(len(data.propers))
                for i, d in enumerate(data.propers):
                    snap.dihedrals.group[i] = [d[0], d[1], d[2], d[3]]
                    snap.dihedrals.typeid[i] = self.dihedral_types.index(data.properLabels[i])
            
            # Populate  particle data
            # particle attributes:  acceleration, angmom, body, charge, diameter, image, is_accel_set, mass, 
            # moment_inertia, orientation, position, typeid, types, velocity
            for i in range(nparticles):
                snap.particles.body[i] = data.bodies[i]
                snap.particles.charge[i] = data.charges[i]
                snap.particles.diameter[i] = data.diameters[i]
                snap.particles.image[i] = data.images[i]
                snap.particles.mass[i] = data.masses[i]
                snap.particles.position[i] = data.coords[i]
                snap.particles.typeid[i] = self.particle_types.index(data.atomTypes[i])
        return snap
    
    def setBonds(self):
        if not self.bond_types: return
        bond_harmonic = hoomd.md.bond.harmonic(name="bond_harmonic")
        for bond in self.bond_types:
            param = self.ffield.bondParameter(bond)
            if self.debug: logger.info("DEBUG: bond_harmonic.bond_coeff.set( '{0}',  k={1}, r0={2} )".format(bond, param['k'], param['r0']))
            bond_harmonic.bond_coeff.set(bond, k=param['k'], r0=param['r0'])
        return

    def setAngles(self):
        if not self.angle_types: return
        angle_harmonic = hoomd.md.angle.harmonic()
        for angle in self.angle_types:
            param = self.ffield.angleParameter(angle)
            if self.debug: logger.info("DEBUG: angle_harmonic.angle_coeff.set( '{0}',  k={1}, t0={2} )".format(angle, param['k'], param['t0']))
            angle_harmonic.angle_coeff.set(angle, k=param['k'], t0=param['t0'])
        return
    
    def setDihedrals(self):
        if not self.dihedral_types: return
        dihedral_harmonic = hoomd.md.dihedral.harmonic()
        for dihedral in self.dihedral_types:
            param = self.ffield.dihedralParameter(dihedral)
            if self.debug:logger.info("DEBUG: dihedral_harmonic.dihedral_coeff.set('{0}',  k={1}, d={2}, n={3})".format(dihedral, param['k'], param['d'], param['n']))
            dihedral_harmonic.dihedral_coeff.set(dihedral, k=param['k'], d=param['d'], n=param['n'])
        return

    def setPairs(self, rCut, rigidBody):
        nl = hoomd.md.nlist.cell()
        lj = hoomd.md.pair.lj(r_cut=rCut, nlist=nl)
        for atype, btype in itertools.combinations(self.particle_types, 2):
#             # dummy atoms treated specially
#             #if atype.lower() in ['x','hn'] or btype.lower() in ['x','hn']:
#             if self.masked[i] or self.masked[j]:
#                 epsilon = 0.0
#                 sigma = 0.0
#             else:
            param = self.ffield.pairParameter(atype, btype)
            epsilon = param['epsilon']
            sigma = param['sigma']
            lj.pair_coeff.set(atype, btype, epsilon=epsilon, sigma=sigma)
            if self.debug: logger.info("DEBUG: lj.pair_coeff.set( '{0}', '{1}', epsilon={2}, sigma={3} )".format(atype, btype, epsilon, sigma))

        if rigidBody:
            nl.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle', 'body'])
        else:
            nl.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle'])

        return

    def setupGroups(self, data):
        self.groupAll = hoomd.group.all()
        # All atoms that are part of static fragments
        if any(data.static):
            self.groupStatic = hoomd.group.tag_list(name="static",
                                                      tags=[i for i, s in enumerate(data.static) if s])
            self.groupActive = hoomd.group.difference(name="active",
                                                        a=self.groupAll,
                                                        b=self.groupStatic)
        else:
            self.groupActive = self.groupAll
        return

    def setupSimulation(self, data, snap, rCut, rigidBody):

        # Check we have all the parameters for this snapshot
        self.checkParameters()
        
        # Init the sytem from the snapshot
        system = hoomd.init.read_snapshot(snap)

#         # Below disables pretty much all output
#         if quiet:
#             logger.info("Disabling HOOMD-Blue output!")
#             hoomdblue.globals.msg.setNoticeLevel(0)

        # Set the parameters
        self.setBonds()
        self.setAngles()
        self.setDihedrals()
        self.setPairs(rCut, rigidBody)
        
        # Specify the groups
        self.setupGroups(data)
        
        # Add any walls
        #self.setupWalls(walls, system, wallAtomType, rCut)
        return

class HoomdOptimiser(FFIELD):

    def __init__(self, paramsDir):

        self.ffield = FfieldParameters(paramsDir)
        self.masked = None
        self.rCut = 5.0
        
        self.groupActive = None
        self.groupAll = None
        self.groupStatic = None
        
        self.debug = False

        # kB = 8.310 * 10**-23 Angstroms**2 g mole**-1 s**-2 K**-1
        # Incoming T is in Kelvin so we multiply by kB
        self.CONVERSIONFACTOR = 8.310E-23
        return

    def optimiseGeometry(self,
                          data,
                          xmlFilename="hoomdOpt.xml",
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
            raise RuntimeError, "Cannot have impropers and dihedrals at the same time"

        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       rigidBody=rigidBody,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       walls=walls,
                                       wallAtomType=wallAtomType,
                                       )

        # Logger - we don't write anything but query the final value - hence period 0 and overwrite
        self.hlog = hoomdblue.analyze.log(filename='geomopt.tsv',
                                     quantities=[
                                                 'num_particles',
                                                 'pair_lj_energy',
                                                 'potential_energy',
                                                 'kinetic_energy',
                                                 ],
                                     period=100,
                                     header_prefix='#',
                                     overwrite=True
                                     )

        optimised = self._optimiseGeometry(rigidBody=rigidBody,
                                            **kw)

        # Extract the energy
        if 'd' in kw and kw['d'] is not None:
            for i in ['potential_energy' ]:
                kw['d'][ i ] = self.toStandardUnits(self.hlog.query(i))

        return optimised

    def _optimiseGeometry(self,
                          carOut="hoomdOpt.car",
                          rigidBody=True,
                          optCycles=1000000,
                          dump=False,
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
        
        assert max_tries > 0 and retries_on_error > 0

        for i in range(retries_on_error):
            # Try optimising and lower the timestep and increasing the number of cycles each time
            try:
                # Create the integrator with the values specified
                if rigidBody:
                    fire = hoomdblue.integrate.mode_minimize_rigid_fire(group=self.groupActive,
                                                                         dt=dt,
                                                                         Nmin=Nmin,
                                                                         alpha_start=alpha_start,
                                                                         ftol=ftol,
                                                                         Etol=Etol,
                                                                         finc=finc,
                                                                         fdec=fdec)
                else:
                    fire = hoomdblue.integrate.mode_minimize_fire(group=self.groupActive,
                                                                dt=dt,
                                                                Nmin=Nmin,
                                                                alpha_start=alpha_start,
                                                                ftol=ftol,
                                                                Etol=Etol,
                                                                finc=finc,
                                                                fdec=fdec)
                if dump:
                    # For tracking the optimsation
                    xmld = hoomdblue.dump.xml(filename="runopt.xml", vis=True)
                    dcdd = hoomdblue.dump.dcd(filename="runopt.dcd",
                                              period=1,
                                              unwrap_full=True,
                                              overwrite=True)
                optimised = False
                for j in range(max_tries):
                    logger.info("Running {0} optimisation cycles in macrocycle {1}".format(optCycles,j))
                    hoomdblue.run(optCycles,
                                   callback=lambda x:-1 if fire.has_converged() else 0,
                                   callback_period=1)
                    if fire.has_converged():
                        logger.info("Optimisation converged on macrocycle {0}".format(j))
                        optimised = True
                        break
                    else:
                        logger.info("Optimisation failed to converge on macrocycle {0}".format(j))
                        if j+1 < max_tries:  logger.info("Attempting another optimisation macrocycle")
                break # Break out of try/except loop
            except RuntimeError as e:
                logger.info("Optimisation step {0} failed!\n{1}".format(i,e))
                if i+1 < retries_on_error:
                    dt_old = dt
                    dt = dt_old * 0.1
                    logger.info("Rerunning optimisation changing dt {0} -> {1}".format(dt_old, dt))
                else:
                    # If we're not going to try again we re-raise the last exception that we caught
                    raise
            # if float( self.hlog.query( 'potential_energy' ) ) < 1E-2:
            #    print "!!!!!!!!!!!HACK CONVERGENCE CRITERIA!!!!!!"
            #    optimised=True
            #    break
        # Delete variables before we return to stop memory leaks
        # del fire
        if dump and False:
            xmld.disable()
            dcdd.disable()
            del xmld
            del dcdd
#         del harmonic
#         del aharmonic
#         del improper
#         del lj
        # Write out a car file so we can see what happened
        # self.writeCar( system=self.system, filename=carOut, unwrap=True )
        return optimised

    def runMD(self,
              data,
              xmlFilename="hoomdMD.xml",
              doDihedral=False,
              doImproper=False,
              doCharges=True,
              rigidBody=True,
              rCut=None,
              quiet=None,
              walls=None,
              wallAtomType=None,
              **kw):

        if rCut is not None:
            self.rCut = rCut

        if doDihedral and doImproper:
            raise RuntimeError, "Cannot have impropers and dihedrals at the same time"

        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       rigidBody=rigidBody,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       walls=walls,
                                       wallAtomType=wallAtomType,
                                       )

        # Logger - we don't write anything but query the final value - hence period 0 and overwrite
        hlog = hoomdblue.analyze.log(filename='runmd.tsv',
                                     quantities=[
                                                 'num_particles',
                                                 'pair_lj_energy',
                                                 'potential_energy',
                                                 'kinetic_energy',
                                                 ],
                                     period=100,
                                     header_prefix='#',
                                     overwrite=True
                                     )

        self._runMD(rigidBody=rigidBody, **kw)

        # Extract the energy
        if 'd' in kw and kw['d'] is not None:
            for i in ['potential_energy' ]:
                kw['d'][ i ] = self.toStandardUnits(hlog.query(i))

        return True

    def runMDAndOptimise(self,
                         data,
                         xmlFilename="hoomdMdOpt.xml",
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
            raise RuntimeError, "Cannot have impropers and dihedrals at the same time"

        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       rigidBody=rigidBody,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       walls=walls,
                                       wallAtomType=wallAtomType,
                                       )

        # Logger - we don't write anything but query the final value - hence period 0 and overwrite
        self.hlog = hoomdblue.analyze.log(filename='md_geomopt.tsv',
                                     quantities=[
                                                 'num_particles',
                                                 'pair_lj_energy',
                                                 'potential_energy',
                                                 'kinetic_energy',
                                                 ],
                                     period=0,
                                     header_prefix='#',
                                     overwrite=True
                                     )

        # pre-optimise for preopt steps to make sure the sytem is sane - otherwise the MD
        # blows up
        preOptCycles = 5000
        self._optimiseGeometry(optCycles=preOptCycles,
                               dt=0.0001)
        # Now run the MD steps
        self._runMD(**kw)

        # Finally do a full optimisation
        optimised = self._optimiseGeometry(**kw)

        # Extract the energy
        if 'd' in kw and kw['d'] is not None:
            for i in ['potential_energy' ]:
                kw['d'][ i ] = self.toStandardUnits(self.hlog.query(i))

        return optimised

    def _runMD(self, mdCycles=100000,
               rigidBody=True,
               integrator='nvt',
               T=1.0,
               tau=0.5,
               P=1,
               tauP=0.5,
               dt=0.0005,
               dump=False,
               dumpPeriod=100, **kw):

        # Added **kw arguments so that we don't get confused by arguments intended for the optimise
        # when MD and optimiser run together
        # T= 0.1

        # Convert T
        # kB = 8.310 * 10**-23 Angstroms**2 g mole**-1 s**-2 K**-1
        # Incoming T is in Kelvin so we multiply by kB
        T = self.fromStandardUnits(T)
        integrator_mode = hoomdblue.integrate.mode_standard(dt=dt)

        if integrator == 'nvt':
            if rigidBody:
                # nvt = hoomdblue.integrate.nvt_rigid(group=hoomdblue.group.rigid(), T=T, tau=tau )
                integ = hoomdblue.integrate.nvt_rigid(group=self.groupActive, T=T, tau=tau)
            else:
                integ = hoomdblue.integrate.nvt(group=self.groupActive, T=T, tau=tau)
        elif integrator == 'npt':
            if rigidBody:
                integ = hoomdblue.integrate.npt_rigid(group=self.groupActive, T=T, tau=tau, P=P, tauP=tauP)
            else:
                integ = hoomdblue.integrate.npt(group=self.groupActive, T=T, tau=tau, P=P, tauP=tauP)
        else:
            raise RuntimeError("Unrecognised integrator: {0}".format(integrator))

        if dump:
            xmld = hoomdblue.dump.xml(filename="runmd.xml", vis=True)
            dcdd = hoomdblue.dump.dcd(filename="runmd.dcd",
                                      period=dumpPeriod,
                                      unwrap_full=True,
                                      overwrite=True)

        # run mdCycles time steps
        hoomdblue.run(mdCycles)

        integ.disable()
        if dump:
            #xmld.disable()
            dcdd.disable()
            del xmld
            del dcdd
        del integ
        del integrator_mode
        return

    def setupSystem(self,
                    data,
                    xmlFilename,
                    doDihedral=False,
                    doImproper=False,
                    doCharges=True,
                    rigidBody=True,
                    rCut=None,
                    quiet=False,
                    maxE=False,
                    walls=None,
                    wallAtomType=None,
                     ):

        self.writeXml(data,
                      xmlFilename=xmlFilename,
                      rigidBody=rigidBody,
                      doCharges=doCharges,
                      doDihedral=doDihedral,
                      doImproper=doImproper)

        # Set the masked array
        self.masked = data.masked
        # Read parameters from file, set them as  attributes & check them
        self.setAttributesFromFile(xmlFilename)
        self.checkParameters()

        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()

        # Init the sytem from the file
        system = hoomdblue.init.read_xml(filename=xmlFilename)

        # Below disables pretty much all output
        if quiet:
            logger.info("Disabling HOOMD-Blue output!")
            hoomdblue.globals.msg.setNoticeLevel(0)

        # Set the parameters
        harmonic = None
        if len(self.bonds):
            harmonic = hoomdblue.bond.harmonic()
            self.setBond(harmonic)

        aharmonic = None
        if len(self.angles):
            aharmonic = hoomdblue.angle.harmonic()
            self.setAngle(aharmonic)

        dharmonic = improper = None
        if doDihedral and len(self.dihedrals):
                dharmonic = hoomdblue.dihedral.harmonic()
                self.setDihedral(dharmonic)
        elif doImproper and len(self.dihedrals):
            improper = hoomdblue.improper.harmonic()
            self.setImproper(improper)

        lj = hoomdblue.pair.lj(r_cut=rCut)
        self.setPair(lj)
        
        # Specify the groups
        self.setupGroups(data)
        
        # Add any walls
        self.setupWalls(walls, system, wallAtomType, rCut)

        if rigidBody:
            hoomdblue.globals.neighbor_list.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle', 'body'])
        else:
            hoomdblue.globals.neighbor_list.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle'])

        # For calculating groups
        # self.labels=[]
        if maxE:
            for idxBlock, idxFragment, start, end in data.tagIndices:
                # create the group
                l = "{0}:{1}".format(start, end)
                g = hoomdblue.group.tag_list(name=l, tags=range(start, end))
                print "SET GROUP ", start, end
                # create the compute for this group
                c = hoomdblue.compute.thermo(group=g)
                # self.labels.append(l)


# Do we need to think about saving the references and deleting them?
#         del harmonic
#         del aharmonic
#         del improper
#         del lj

        return system
    
    def setupWalls(self, walls, system, wallAtomType, rCut):
        """Set up walls for the simulation
        
        I think that we require two walls. One on the front side with the potential facing in,
        the other on the back wall with the potential facing back towards the other potential.
        The origin of the wall is the centre of plane but then back half a cell along the axis 
        that isn't part of the wall.
        """
        if walls is None: return
        assert len(walls) is 3 # array of three booleans - one per wall
        if not any(walls): return
        assert wallAtomType is not None,"Need to set a wallAtomType!"
        wtypes = ['XOY', 'XOZ', 'YOZ'] # Oroder of walls
        wallstructure = None
        for do, wtype in zip(walls, wtypes):
            if do:
                logger.info('Setting wall in HOOMDBLUE for {0} of atomType {1}'.format(wtype, wallAtomType))
                if wtype is 'XOY':
                    originFront = (0, 0, -system.box.Lz/2)
                    originBack  = (0, 0,  system.box.Lz/2)
                    normal = (0, 0, 1)
                elif wtype is 'XOZ':
                    originFront = (0, -system.box.Ly/2, 0)
                    originBack  = (0,  system.box.Ly/2, 0)
                    normal = (0, 1, 0)
                elif wtype is 'YOZ':
                    originFront = (-system.box.Lx/2, 0, 0)
                    originBack  = ( system.box.Lx/2, 0, 0)
                    normal = (1, 0, 0)
                else:
                    raise RuntimeError("Unrecognised Wall Type! {0}".format(wtype))
            
                if not wallstructure:
                    # We only create the wall and the LJ potentials once as they are used
                    # by all subsequent walls in the group
                    try: 
                        wallstructure = hoomdblue.wall.group()
                    except AttributeError:
                        raise RuntimeError('HOOMD-blue wall does not have a group attribute. You may need to update your version of HOOMD-Blue in order to use walls')
                    lj = hoomdblue.wall.lj(wallstructure, r_cut=rCut)
                    for atype in self.atomTypes:
                        param = self.ffield.pairParameter(atype, wallAtomType)
                        lj.force_coeff.set(atype, epsilon=param['epsilon'], sigma=param['sigma'])
                
                # Add the two walls
                # Front
                wallstructure.add_plane(origin=originFront, normal=normal, inside=True)
                # Back
                wallstructure.add_plane(origin=originBack, normal=normal, inside=False)
        return

    def toStandardUnits(self, value):
        # return float(value) / self.CONVERSIONFACTOR
        return float(value)

    def writeCar(self, system, filename, unwrap=True, pbc=True):
        """Car File
        """

        car = "!BIOSYM archive 3\n"
        car += "PBC=ON\n"

        car += "ambuild generated car file\n"
        tstr = time.strftime("%a %b %d %H:%M:%S %Y", time.gmtime())
        car += "!DATE {0}\n".format(tstr)

        xdim = system.box[0]
        ydim = system.box[1]
        zdim = system.box[2]
        if pbc:
            car += "PBC  {0: < 9.4F} {1: < 9.4F} {2: < 9.4F}  90.0000   90.0000   90.0000 (P1)\n".format(xdim,
                                                                                                      ydim,
                                                                                                      zdim)

        for p in system.particles:

                label = atype = p.type.strip()

                # Treat x-atoms differently
                if label[0].lower() == 'x':
                    symbol = 'x'
                else:
                    symbol = util.label2symbol(label)

                x, y, z = p.position
                ix, iy, iz = p.image
                charge = float(p.charge)

                if unwrap:
                    x = util.unWrapCoord(x, ix, xdim, centered=False)
                    y = util.unWrapCoord(y, iy, ydim, centered=False)
                    z = util.unWrapCoord(z, iz, zdim, centered=False)
                else:
                    # Put back with origin at corner
                    x = x + (xdim / 2)
                    y = y + (ydim / 2)
                    z = z + (zdim / 2)

                car += "{0: <5} {1: >15.10} {2: >15.10} {3: >15.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format(label, x, y, z, atype, symbol, charge)

        car += "end\nend\n\n"

        with open(filename, 'w') as f:
            f.writelines(car)

        return

    def writeXml(self,
                 data,
                 xmlFilename="hoomd.xml",
                 rigidBody=True,
                 doCharges=True,
                 doDihedral=False,
                 doImproper=False
                 ):
        """Write out a HOOMD Blue XML file.
        """
        d = data
        #
        # Got data so now put into xml
        #
        body = "\n" + "\n".join(map(str, d.bodies)) + "\n"
        charge = "\n" + "\n".join(map(str, d.charges)) + "\n"
        diameter = "\n" + "\n".join(map(str, d.diameters)) + "\n"
        mass = "\n" + "\n".join(map(str, d.masses)) + "\n"
        ptype = "\n" + "\n".join(map(str, d.atomTypes)) + "\n"

        image = "\n"
        position = "\n"
        for i in range(len(d.coords)):
            image += "{0} {1} {2}\n".format(d.images[i][0],
                                             d.images[i][1],
                                             d.images[i][2])
            position += "{0} {1} {2}\n".format(d.coords[i][0],
                                                d.coords[i][1],
                                                d.coords[i][2])

        # Now do all angles and bonds
        bond = False
        if len(d.bonds):
            bond = "\n"
            for i, b in enumerate(d.bonds):
                bond += "{0} {1} {2}\n".format(d.bondLabels[i], b[0], b[1])

        angle = False
        if len(d.angles):
            angle = "\n"
            for i, a in enumerate(d.angles):
                angle += "{0} {1} {2} {3}\n".format(d.angleLabels[i], a[0], a[1], a[2])

        dihedral = False
        if len(d.propers):
            dihedral = "\n"
            for i, dh in enumerate(d.propers):
                dihedral += "{0} {1} {2} {3} {4}\n".format(d.properLabels[i], dh[0], dh[1], dh[2], dh[3])

        root = ET.Element('hoomd_xml', version="1.4")
        config = ET.SubElement(root, "configuration", timestep="0")

        # e = ET.SubElement( config, "box",
        ET.SubElement(config, "box",
                        Lx=str(d.cell[0]),
                        Ly=str(d.cell[1]),
                        Lz=str(d.cell[2]))

        e = ET.SubElement(config, "position")
        e.text = position
        e = ET.SubElement(config, "image")
        e.text = image

        if rigidBody:
            e = ET.SubElement(config, "body")
            e.text = body
        if doCharges:
            e = ET.SubElement(config, "charge")
            e.text = charge
        e = ET.SubElement(config, "diameter")
        e.text = diameter
        e = ET.SubElement(config, "type")
        e.text = ptype
        e = ET.SubElement(config, "mass")
        e.text = mass

        if bond:
            e = ET.SubElement(config, "bond")
            e.text = bond
        if angle:
            e = ET.SubElement(config, "angle")
            e.text = angle
        if dihedral:
            if doDihedral:
                e = ET.SubElement(config, "dihedral")
                e.text = dihedral
            elif doImproper:
                e = ET.SubElement(config, "improper")
                e.text = dihedral

        tree = ET.ElementTree(root)

        # ET.dump(tree)

        # tree.write(file_or_filename, encoding, xml_declaration, default_namespace, method)
        tree.write(xmlFilename)

        return True

    def writeXyz(self, system, filename=None, unwrap=True):
        """Car File
        """

        xyz = "{0}\n".format(len(system.particles))
        xyz += "opt.xyz file\n"

        xdim = system.box[0]
        ydim = system.box[1]
        zdim = system.box[2]

        for p in system.particles:

                label = util.label2symbol(p.type.strip())
                x, y, z = p.position
                ix, iy, iz = p.image

                if unwrap:
                    x = util.unWrapCoord(x, ix, xdim, centered=False)
                    y = util.unWrapCoord(y, iy, ydim, centered=False)
                    z = util.unWrapCoord(z, iz, zdim, centered=False)
                else:
                    # Put back with origin at corner
                    x = x + (xdim / 2)
                    y = y + (ydim / 2)
                    z = z + (zdim / 2)
                xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format(label, x, y, z)

        with open(filename, 'w') as f:
            f.writelines(xyz)

        return


if __name__ == "__main__":
    from paths import PARAMS_DIR
    mycell = util.cellFromPickle(sys.argv[1])
    rigidBody = False
    data = mycell.dataDict(periodic=True, center=True, rigidBody=rigidBody)
    
    opt = Hoomd2(PARAMS_DIR)
    hoomd.context.initialize()
    snap = opt.createSnapshot(data, rigidBody)
    rCut = 5.0
    opt.setupSimulation(data, snap, rCut, rigidBody)
#     ok = opt.optimiseGeometry(data,
#                             xmlFilename="opt.xml",
#                             rigidBody=rigidBody,
#                             doDihedral=True,
#                             doImproper=False,
#                             doCharges=True
#                             )

