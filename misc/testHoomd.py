from hoomd_script import *

xmlFilename="hoomdOpt.xml"
system = init.read_xml( filename=xmlFilename )

harmonic = bond.harmonic()
harmonic.bond_coeff.set( 'cp-cp',  k=1550.0, r0=1.384 )

dharmonic = dihedral.harmonic()
dharmonic.set_coeff( 'c-cp-cp-c', k=30, d=-1, n=2  )

lj = pair.lj( r_cut=5.0)
lj.pair_coeff.set( 'c', 'c', epsilon=0.0933, sigma=3.7736 )
lj.pair_coeff.set( 'c', 'cp', epsilon=0.1016, sigma=3.7736 )
lj.pair_coeff.set( 'c', 'br', epsilon=0.1016, sigma=3.7736 )
lj.pair_coeff.set( 'cp', 'cp', epsilon=0.1106, sigma=1.7736 )
lj.pair_coeff.set( 'cp', 'br', epsilon=0.1106, sigma=3.7736 )
lj.pair_coeff.set( 'br', 'br', epsilon=0.0376, sigma=3.4893 )

globals.neighbor_list.reset_exclusions(exclusions = ['1-2', '1-3', '1-4', 'angle', 'body'] )

hlog = analyze.log(filename='mylog.csv',
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

fire = integrate.mode_minimize_rigid_fire( group=group.all(),
                                           dt=0.005,
                                           Nmin=5,
                                           alpha_start=0.1,
                                           ftol=1e-2,
                                           Etol=1e-4,
                                           finc=1.1,
                                           fdec=0.5 )
                                           
xmld = dump.xml(filename="runopt.xml", vis=True)
dcdd = dump.dcd(filename="runopt.dcd",
               period=1, 
               unwrap_full=True,
               overwrite=True,
                )

for i in range(5):
     run( 100000 )
     if fire.has_converged():
         break
 
