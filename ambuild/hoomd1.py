#!/usr/bin/env python
import sys
from hoomd_script import *
# Need to avoid name clash with dump now that we've removed the hooomdblue namespace we had earlier
hdump = dump
del dump
if __name__ != '__main__' and 'context' in locals().keys():
    # Later versions of hoomd-blue have a context that needs to be updated
    # Context parses the command-line, so we need to make sure there
    # are no arguments or else we lose the ability to display our own help
    _argv = sys.argv
    sys.argv = [sys.argv[0]]
    context.initialize()
    sys.argv = _argv
import logging
import numpy as np
import os
import time
import xml.etree.ElementTree as ET
# our imports
from ab_ffield import FFIELD, FfieldParameters
import xyz_core
import xyz_util

logger = logging.getLogger(__name__)

class Hoomd1(FFIELD):

    def __init__(self, paramsDir):
        self.ffield = FfieldParameters(paramsDir)
        self.system = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        self.impropers = None
        self.atomTypes = None
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

    def fragMaxEnergy(self,
                      data,
                      xmlFilename,
                      doDihedral=False,
                      doImproper=False,
                      rigidBody=True,
                      rCut=None,
                      quiet=None,
                       **kw):
        if rCut is not None:
            self.rCut = rCut
        if doDihedral and doImproper:
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       maxE=True)
        quantities = ['potential_energy', 'kinetic_energy']
        # see setupSystem for where the groups are created
        for idxBlock, idxFragment, start, end in data.tagIndices:
            quantities.append("potential_energy_{0}:{1}".format(start, end))
        self.hlog = analyze.log(filename='mylog1.csv',
                                quantities=quantities,
                                period=1,
                                header_prefix='#',
                                overwrite=True )
        optimised = self._optimiseGeometry(rigidBody=rigidBody, optCycles=10, **kw)

        maxe = -10000
        maxi = -1
        i=0 # hack
        for i, (idxBlock, idxFragment, start, end) in enumerate(data.tagIndices):
            label = "potential_energy_{0}:{1}".format(start, end)
            e = self.toStandardUnits(self.hlog.query(label))
            if e > maxe:
                maxe = e
                maxi = i
        assert i != -1
        idxBlock, idxFragment, start, end = data.tagIndices[maxi]
        return maxe, idxBlock, idxFragment
        return float(value)

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
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       rigidBody=rigidBody,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       walls=walls,
                                       wallAtomType=wallAtomType)
        # Logger - we don't write anything but query the final value - hence period 0 and overwrite
        self.hlog = analyze.log(filename='geomopt.tsv',
                                quantities=[
                                    'num_particles',
                                    'pair_lj_energy',
                                    'potential_energy',
                                    'kinetic_energy',
                                    ],
                                period=100,
                                header_prefix='#',
                                overwrite=True)
        optimised = self._optimiseGeometry(rigidBody=rigidBody, **kw)
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
        assert max_tries > 0
        if retries_on_error < 1:
            retries_on_error = 0
        retries_on_error += 1
        for i in range(retries_on_error):
            # Try optimising and lower the timestep and increasing the number of cycles each time
            try:
                # Create the integrator with the values specified
                if rigidBody:
                    fire = integrate.mode_minimize_rigid_fire(group=self.groupActive,
                                                              dt=dt,
                                                              Nmin=Nmin,
                                                               alpha_start=alpha_start,
                                                               ftol=ftol,
                                                               Etol=Etol,
                                                               finc=finc,
                                                               fdec=fdec)
                else:
                    fire = integrate.mode_minimize_fire(group=self.groupActive,
                                                        dt=dt,
                                                        Nmin=Nmin,
                                                        alpha_start=alpha_start,
                                                        ftol=ftol,
                                                        Etol=Etol,
                                                        finc=finc,
                                                        fdec=fdec)
                if dump:
                    # For tracking the optimsation
                    xmld = hdump.xml(filename="runopt.xml", vis=True)
                    dcdd = hdump.dcd(filename="runopt.dcd",
                                    period=1,
                                    unwrap_full=True,
                                    overwrite=True)
                optimised = False
                for j in range(max_tries):
                    logger.info("Running {0} optimisation cycles in macrocycle {1}".format(optCycles,j))
                    run(optCycles,
                        callback=lambda x:-1 if fire.has_converged() else 0,
                        callback_period=1)
                    if fire.has_converged():
                        logger.info("Optimisation converged on macrocycle {0}".format(j))
                        optimised = True
                        break
                    else:
                        logger.info("Optimisation failed to converge on macrocycle {0}".format(j))
                        if j+1 < max_tries:
                            logger.info("Attempting another optimisation macrocycle")
                break # Break out of try/except loop
            except RuntimeError as e:
                logger.info("Optimisation step {0} failed!\n{1}".format(i,e))
                del xmld
                del dcdd
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
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       rigidBody=rigidBody,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       walls=walls,
                                       wallAtomType=wallAtomType)

        # Logger - we don't write anything but query the final value - hence period 0 and overwrite
        hlog = analyze.log(filename='runmd.tsv',
                           quantities=[
                               'num_particles',
                               'pair_lj_energy',
                               'potential_energy',
                               'kinetic_energy',
                               ],
                           period=100,
                           header_prefix='#',
                           overwrite=True)
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
            raise RuntimeError("Cannot have impropers and dihedrals at the same time")
        self.system = self.setupSystem(data,
                                       xmlFilename=xmlFilename,
                                       rigidBody=rigidBody,
                                       doDihedral=doDihedral,
                                       doImproper=doImproper,
                                       rCut=self.rCut,
                                       quiet=quiet,
                                       walls=walls,
                                       wallAtomType=wallAtomType)

        # Logger - we don't write anything but query the final value - hence period 0 and overwrite
        self.hlog = analyze.log(filename='md_geomopt.tsv',
                                quantities=[
                                    'num_particles',
                                    'pair_lj_energy',
                                    'potential_energy',
                                    'kinetic_energy',
                                    ],
                                period=0,
                                header_prefix='#',
                                overwrite=True)
        # pre-optimise for preopt steps to make sure the sytem is sane - otherwise the MD
        # blows up
        preOptCycles = 5000
        self._optimiseGeometry(optCycles=preOptCycles, dt=0.0001)
        self._runMD(**kw)
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
        integrator_mode = integrate.mode_standard(dt=dt)
        if integrator == 'nvt':
            if rigidBody:
                # nvt = integrate.nvt_rigid(group=group.rigid(), T=T, tau=tau )
                integ = integrate.nvt_rigid(group=self.groupActive, T=T, tau=tau)
            else:
                integ = integrate.nvt(group=self.groupActive, T=T, tau=tau)
        elif integrator == 'npt':
            if rigidBody:
                integ = integrate.npt_rigid(group=self.groupActive, T=T, tau=tau, P=P, tauP=tauP)
            else:
                integ = integrate.npt(group=self.groupActive, T=T, tau=tau, P=P, tauP=tauP)
        else:
            raise RuntimeError("Unrecognised integrator: {0}".format(integrator))

        if dump:
            xmld = hdump.xml(filename="runmd.xml", vis=True)
            dcdd = hdump.dcd(filename="runmd.dcd",
                            period=dumpPeriod,
                            unwrap_full=True,
                            overwrite=True)
        run(mdCycles)
        integ.disable()
        if dump:
            #xmld.disable()
            dcdd.disable()
            del xmld
            del dcdd
        del integ
        del integrator_mode
        return

    def setAngle(self, anglePotential):
        for angle in self.angles:
            param = self.ffield.angleParameter(angle)
            if self.debug:
                logger.info("DEBUG: angle set_coeff( '{0}',  k={1}, t0={2} )".format(angle, param['k'], param['t0']))
            anglePotential.set_coeff(angle, k=param['k'], t0=param['t0'])
        return

    def setBond(self, bondPotential):
        for bond in self.bonds:
            param = self.ffield.bondParameter(bond)
            if self.debug:
                logger.info("DEBUG: bond_coeff.set( '{0}',  k={1}, r0={2} )".format(bond, param['k'], param['r0']))
            bondPotential.bond_coeff.set(bond, k=param['k'], r0=param['r0'])
        return

    def setDihedral(self, dihedralPotential):
        for dihedral in self.dihedrals:
            param = self.ffield.dihedralParameter(dihedral)
            if self.debug:
                logger.info("DEBUG: dihedral set_coeff.('{0}',  k={1}, d={2}, n={3})".format(dihedral, param['k'], param['d'], param['n']))
            dihedralPotential.set_coeff(dihedral, k=param['k'], d=param['d'], n=param['n'])
        return

    def setImproper(self, improperPotential):
        for improper in self.impropers:
            param = self.ffield.improperParameter(improper)
            improperPotential.set_coeff(improper, k=param['k'], chi=param['chi'])
        return

    def setPair(self, pairPotential):
        for i, atype in enumerate(self.atomTypes):
            for j, btype in enumerate(self.atomTypes):
                if j >= i:
                    # dummy atoms treated specially
                    #if atype.lower() in ['x','hn'] or btype.lower() in ['x','hn']:
                    if self.masked[i] or self.masked[j]:
                        epsilon = 0.0
                        sigma = 0.0
                    else:
                        param = self.ffield.pairParameter(atype, btype)
                        epsilon = param['epsilon']
                        sigma = param['sigma']

                    pairPotential.pair_coeff.set(atype, btype, epsilon=epsilon, sigma=sigma)
                    if self.debug:
                        logger.info("DEBUG: pair_coeff.set( '{0}', '{1}', epsilon={2}, sigma={3} )".format(atype,
                                                                                                           btype,
                                                                                                           epsilon,
                                                                                                           sigma))
        return

    def setAttributesFromFile(self, xmlFilename):
        """Parse the xml file to extract the bonds, angles etc."""
        tree = ET.parse(xmlFilename)
        root = tree.getroot()
        bonds = []
        x = root.findall(".//bond")
        if len(x):
            btext = x[0].text
            for line in btext.split(os.linesep):
                line = line.strip()
                if line:
                    bond = line.split()[0]
                    if bond not in bonds:
                        bonds.append(bond)
        self.bonds = bonds

        angles = []
        x = root.findall(".//angle")
        if len(x):
            atext = x[0].text
            for line in atext.split(os.linesep):
                line = line.strip()
                if line:
                    angle = line.split()[0]
                    if angle not in angles:
                        angles.append(angle)
        self.angles = angles

        dihedrals = []
        dn = root.findall(".//dihedral")
        if len(dn):
            dtext = dn[0].text
            for line in dtext.split(os.linesep):
                line = line.strip()
                if line:
                    dihedral = line.split()[0]
                    if dihedral not in dihedrals:
                        dihedrals.append(dihedral)
        self.dihedrals = dihedrals

        impropers = []
        dn = root.findall(".//improper")
        if len(dn):
            dtext = dn[0].text
            for line in dtext.split(os.linesep):
                line = line.strip()
                if line:
                    improper = line.split()[0]
                    if improper not in impropers:
                        impropers.append(improper)
        self.impropers = impropers

        atomTypes = []
        atext = root.findall(".//type")[0].text
        for line in atext.split(os.linesep):
            atomType = line.strip()
            if atomType and atomType not in atomTypes:
                atomTypes.append(atomType)
        self.atomTypes = atomTypes
        return

    def setupGroups(self, data):
        self.groupAll = group.all()
        # All atoms that are part of static fragments
        if any(data.static):
            self.groupStatic = group.tag_list(name="static",
                                                      tags=[i for i, s in enumerate(data.static) if s])
            self.groupActive = group.difference(name="active",
                                                        a=self.groupAll,
                                                        b=self.groupStatic)
        else:
            self.groupActive = self.groupAll
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

        if init.is_initialized(): init.reset()

        # Init the sytem from the file
        system = init.read_xml(filename=xmlFilename)

        # Below disables pretty much all output
        if quiet:
            logger.info("Disabling HOOMD-Blue output!")
            globals.msg.setNoticeLevel(0)

        # Set the parameters
        harmonic = None
        if len(self.bonds):
            harmonic = bond.harmonic()
            self.setBond(harmonic)

        aharmonic = None
        if len(self.angles):
            aharmonic = angle.harmonic()
            self.setAngle(aharmonic)

        dharmonic = improper = None
        if doDihedral and len(self.dihedrals):
                dharmonic = dihedral.harmonic()
                self.setDihedral(dharmonic)
        elif doImproper and len(self.dihedrals):
            improper = improper.harmonic()
            self.setImproper(improper)

        lj = pair.lj(r_cut=rCut)
        self.setPair(lj)

        # Specify the groups
        self.setupGroups(data)

        # Add any walls
        self.setupWalls(walls, system, wallAtomType, rCut)

        if rigidBody:
            globals.neighbor_list.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle', 'body'])
        else:
            globals.neighbor_list.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'angle'])

        # For calculating groups
        # self.labels=[]
        if maxE:
            for idxBlock, idxFragment, start, end in data.tagIndices:
                # create the group
                l = "{0}:{1}".format(start, end)
                g = group.tag_list(name=l, tags=range(start, end))
                # create the compute for this group
                c = compute.thermo(group=g)
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
                        wallstructure = wall.group()
                    except AttributeError:
                        raise RuntimeError('HOOMD-blue wall does not have a group attribute. You may need to update your version of HOOMD-Blue in order to use walls')
                    lj = wall.lj(wallstructure, r_cut=rCut)
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

    def updateCell(self, cell):
        """Reset the particle positions from hoomdblue system"""  # Should really check HOOMD version but...
        if hasattr(self.system.box, "Lx"):
            box = np.array([self.system.box.Lx, self.system.box.Ly, self.system.box.Lz])
        else:
            box = np.array(self.system.box)
        # Read back in the particle positions
        atomCount = 0
        for block in cell.blocks.itervalues():
            for k in range(block.numAtoms()):
                p = self.system.particles[ atomCount ]
                coord = xyz_core.unWrapCoord3(p.position, p.image, box, centered=True)
                block.coord(k, coord)
                atomCount += 1

        if atomCount != len(self.system.particles):
            raise RuntimeError("Read {0} positions but there were {1} particles!".format(atomCount, len(self.system.particles)))

        # If we are running (e.g.) an NPT simulation, the cell size may have changed. In this case we need to update
        # our cell parameters. Repopulate cells will then update the halo cells and add the new blocks
        if not np.allclose(box, cell.dim):
            logger.info("Changing cell dimensions after HOOMD-blue simulation from: {0} to: {1}".format(cell.dim,box))
            cell.dim = box

        # Now have the new coordinates, so we need to put the atoms in their new cells
        cell.repopulateCells()
        return

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
                    symbol = xyz_util.label2symbol(label)

                charge = float(p.charge)

                if unwrap:
                    coord = xyz_core.unWrapCoord3(p.position, p.image, system.box, centered=False)
                else:
                    pass
                    # Put back with origin at corner
#                     x = x + (xdim / 2)
#                     y = y + (ydim / 2)
#                     z = z + (zdim / 2)

                car += "{0: <5} {1: >15.10} {2: >15.10} {3: >15.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format(label, coord[0], coord[1], coord[2], atype, symbol, charge)
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
                label = xyz_util.label2symbol(p.type.strip())
                x, y, z = p.position
                ix, iy, iz = p.image
                if unwrap:
                    x, y, z = xyz_core.unWrapCoord3(p.position, p.image, system.box, centered=False)
                else:
                    # Put back with origin at corner
                    x = x + (xdim / 2)
                    y = y + (ydim / 2)
                    z = z + (zdim / 2)
                xyz += "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format(label, x, y, z)
        with open(filename, 'w') as f:
            f.writelines(xyz)
        return

def xml2xyz(xmlFilename, xyzFilename):
    """Convert a hoomdblue xml file to xyz"""
    tree = ET.parse(xmlFilename)
    root = tree.getroot()
    atext = root.findall(".//type")[0].text
    atomTypes = [ line.strip() for line in atext.split(os.linesep) if line.strip() ]
    ptext = root.findall(".//position")[0].text
    positions = [ line.strip().split() for line in ptext.split(os.linesep) if line.strip() ]
    assert len(atomTypes) == len(positions)
    # Convert atom types to symbols
    symbols = []
    for at in atomTypes:
        if at.lower()[0] == 'x':
            symbols.append('x')
        else:
            symbols.append(xyz_util.label2symbol(at))
    # Now write out xyz
    xyz_util.writeXyz(xyzFilename, positions, symbols)
#     with open( xyzFilename, 'w') as o:
#         o.write("{0}\n".format( len( positions ) ) )
#         o.write("XYZ file created from: {0}\n".format( xmlFilename ) )
#         for i, symbol in enumerate( symbols ):
#             o.write( "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( symbol,
#                                                                            float( positions[i][0] ),
#                                                                            float( positions[i][1] ),
#                                                                            float( positions[i][2] )
#                                                                            ) )
#
#         o.write("\n")
    logger.info("Wrote file: {0}".format(xyzFilename))
    return

if __name__ == "__main__":
    from ab_paths import PARAMS_DIR
    import util
    mycell = util.cellFromPickle(sys.argv[1])
    rigidBody = True
    data = mycell.dataDict(periodic=True, center=True, rigidBody=rigidBody)

    hmd = Hoomd1(PARAMS_DIR)
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
