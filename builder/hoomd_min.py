#!/Applications/HOOMD-blue.app/Contents/MacOS/hoomd

# ---- init_xml.py ----
"""

"""
import cPickle
import math
import sys

from hoomd_script import *


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

def write_xyz( system, filename ):
    
    bonds = [ 63, 89, 48, 135, 176, 227, 268, 299, 350, 391, 340, 427, 627, 678, 622, 714, 709, 745, 637, 791 ]
    angles =  [ 45, 42, 166, 250, 332, 330, 617, 616, 699, 619 ]
    
    f = open( filename, 'w')
    f.write( "{0}\n".format( len( system.particles ) ) )
    f.write("Hoomdblue XYZ file\n")
    for i, p in enumerate( system.particles ):
        
        if i in bonds:
            ptype = 'Cl'
        elif i in angles:
            ptype = 'F'
        else:
            ptype = p.type
        f.write( "{0:5} {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( ptype, p.position[0],p.position[1], p.position[2] ) )
    f.close()
    

def pickle_positions( system, filename="positions.pkl" ):
    positions = []
    for p in system.particles:
        positions.append( p.position )
    f = open( filename, 'w' )
    cPickle.dump(positions, f )
    f.close()
    return
        

filename="hoomd.xml"

# read in the file
system = init.read_xml(filename=filename)

write_xyz( system, "before.xyz")

bharmonic = bond.harmonic()
bharmonic.bond_coeff.set('C-C', k=330.0, r0=5.84)

aharmonic = angle.harmonic()
aharmonic.set_coeff('C-C-C', k=330.0, t0=math.pi)

# simple lennard jones potential
lj = pair.lj(r_cut=10.0)
lj.pair_coeff.set('C', 'C', epsilon=0.15, sigma=4.00)
lj.pair_coeff.set('C', 'H', epsilon=0.0055, sigma=3.00)
lj.pair_coeff.set('H', 'H', epsilon=0.02, sigma=2.00)

#fire=integrate.mode_minimize_fire( group=group.all(), dt=0.05, ftol=1e-2, Etol=1e-7)
fire = integrate.mode_minimize_rigid_fire( group=group.all(), dt=0.05, ftol=1e-2, Etol=1e-7)

# dcd = dump.dcd(filename="trajectory.dcd",period=100)
# mol2 = dump.mol2(filename="trajectory",period=10000)
# run(100000)

count = 0
while not(fire.has_converged()):
    #dcd = dump.dcd(filename="trajectory.dcd",period=100,unwrap_full=True,unwrap_rigid=True)
    #mol2 = dump.mol2(filename="trajectory",period=1000)
    run(1000)
    count += 1
    if count > 20:
        print "TOO MANY ITERATIONS!"
        break

write_xyz( system, "after.xyz")

pickle_positions( system )
