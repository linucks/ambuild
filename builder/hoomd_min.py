#!/Applications/HOOMD-blue.app/Contents/MacOS/hoomd

# ---- init_xml.py ----
"""

"""
import math
import sys

from hoomd_script import *

filename="hoomd.xml"

# read in the file
init.read_xml(filename=filename)

def calculateR0( ri, rj, chiI, chiJ, bondorder ):
    # precompute the equilibrium geometry
    # From equation 3
    rbo = -0.1332 * (ri + rj) * math.log(bondorder)
    # From equation 4
    
    dchi = math.sqrt(chiI) - math.sqrt(chiJ)
    ren = ri * rj * dchi * dchi / (chiI * ri + chiJ * rj)
    
    # From equation 2
    # NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
    # There is a typo in the published paper
    return (ri + rj + rbo - ren)

# here we fold the 1/2 into the kij from equation 1a
# Otherwise, this is equation 6 from the UFF paper.
# kb = KCAL332 * parA.dVal[PAR_Z] * parB.dVal[PAR_Z] / (r0 * r0 * r0);

# example4.xml defines a single polymer: use the same force field as in example 3
# force field setup
harmonic = bond.harmonic()
harmonic.bond_coeff.set('C-C', k=330.0, r0=5.84)

if False:
    # dump every few steps
    dump.mol2(filename="jens_h2o.mol2", period=10)

    # integrate NVT for a bunch of time steps
    all = group.all()
    integrate.mode_standard(dt=0.005)
    integrate.nvt(group=all, T=1.2, tau=0.5)
    run(100)

#fire=integrate.mode_minimize_fire( group=group.all(), dt=0.05, ftol=1e-2, Etol=1e-7)
fire=integrate.mode_minimize_rigid_fire( group=group.all(), dt=0.05, ftol=1e-2, Etol=1e-7)

while not(fire.has_converged()):
    #xml = dump.mol2(filename="dump",period=10)
    xml = dump.dcd(filename="trajectory.dcd",period=10)
    run(100)
