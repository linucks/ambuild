#!/usr/bin/env python
import os
import sys

ambuild_home = "/opt/ambuild.git"
sys.path.insert(0, ambuild_home)

# This imports the builder cell module - this is the only module that should be required
from ambuild import ab_cell

#
# Initialising the system
#
boxDim=[20,20,20] # Create a variable to hold the cell dimensions - this is the first argument to create the cell

# Create a cell object called mycell and specify the parameters that will dictate how close the atoms are permitted
# and what bonds are acceptable
mycell = ab_cell.Cell(boxDim, atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15)

#
# Create two variables that hold the path to the .car files with the coordinates. A corresponding .csv file with
# the endGroup information will need to be found in the same location, so for the amine, there will need to be a
# file called: /Users/jmht/Documents/abbie/AMBI/ambuild/blocks/amine_typed.csv"
fragA=os.path.join(ambuild_home,"blocks/amine_typed.car")
fragB=os.path.join(ambuild_home,"blocks/triquin_typed.car")

# Add the two fragments to the cell and name the fragment type - can be anything you want that doesn't contain a "-"
# or ":" character
mycell.libraryAddFragment(filename=fragA, fragmentType='amine')
mycell.libraryAddFragment(filename=fragB, fragmentType='triquin')

# We now specify which endGroups in which fragments can bond to which. The amine csv file only contains one endGroup
# type (the first column in the csv file) and this is 'a', the triquin csv file also only contains one endGroup type
# the the first column specifies this as a 'b', so we say that amine a can bond to triquin b.
mycell.addBondType( 'amine:a-triquin:b' )

# Even though the amine fragments have 4 endGroups, we will only allow two bonds to each fragment
mycell.setMaxBond('amine:a',2)

# We now start to fill the cell. We seed with 5 amine's, with the first block placed at the center.
mycell.seed(5, fragmentType='amine', center=True )

# Next we add 5 triquin
# The functions seed, growBlocks and joinBlocks all return the number of blocks they have added, so in the below
# case we check that the seed did as we expected and if not we stop the script, saving the last state as a pickle file
toAdd=5
added = mycell.seed( toAdd, fragmentType='triquin')
if added != toAdd:
    print("Failed to seed enough triquin fragments!")
    mycell.dump() # Save the pickle file
    sys.exit() # Exit the script

# We would like to save the state at this point, so we write out a pickle file
# The pickle file will be called step_1.pkl (as this is the first time we have saved the state). To extract the
# coordinates, we use the util.py script in the builder directory and call it with the path to the pkl file as
# the first argument. As everything in this run is in the builder directory, we can run the command:
# ./util.py step_1.pkl
mycell.dump()

# We now grow 5 new triquin fragments onto the amine blocks in the cell. We specify that the dihedral angle of the
# dihedral atoms (as specified in the two csv files) around the bond will be 0 degrees.
# if we had just called:
# mycell.growBlocks(5)
# then 5 fragments of random types would be added to random blocks in the cell at random angles
toAdd=5
added = mycell.growBlocks(toAdd, cellEndGroups=['amine:a'], libraryEndGroups=['triquin:b'], dihedral=0)
if added != toAdd:
    print("Failed to grow enough blocks!")
    mycell.dump() # Save the pickle file
    sys.exit() # Exit the script

# Save the state after growing the blocks - this will be called step_2.pkl
mycell.dump()

# We now do a rigid-body optimisation of the geometry of all the blocks. We will include dihedral terms in the forcefield
# and will attempt 100000 cycles of optimisation - although the optimisation will stop as soon as convergence is reached.
# All the other keyword arguments are just passed through to the mode_minimize_rigid_fire optimiser of hoomdblue. To see
# what arguments are acceptable, see:
# http://codeblue.umich.edu/hoomd-blue/doc/classhoomd__script_1_1integrate_1_1mode__minimize__rigid__fire.html
#
# The optimiser returns True or False depending on whether it converged or not
ok = mycell.optimiseGeometry(rigidBody=True,
                        doDihedral=True,
                        optCycles = 1000000,
                        rCut=5.0,
                        dt=0.005,
                        Etol=1e-5)

# Save the result of the optimisation regardless
mycell.dump()

# check if it worked
if not ok:
    print("Failed to optimise!")
    sys.exit() # Exit the script

# We now create a loop where we run some MD to jiggle the system about, see if we can zip any blocks together
# and optimise the result if we make any bonds
for i in range(3):
    print("Running zip loop: {0}".format(i))

    # Run a few cycles of MD
    mycell.runMD( doDihedral=True,
                  quiet=True,
                  rCut=5.0,
                  mdCycles=100,
                  T=1.0,
                  tau=0.5,
                  dt=0.0005)

    # We now see if any endGroups have come close enough to bond
    # Set clashCheck to True to disallow the bond if any atoms come close to the bond vector
    added = mycell.zipBlocks(bondMargin=5, bondAngleMargin=45, clashCheck=True)
    if added > 0:
        print("Added {0} bonds in zip step".format(added))
        # We've made some bonds! We now optimise the geometry of the new structure
        ok = mycell.optimiseGeometry(rigidBody=True,
                        doDihedral=True,
                        optCycles = 1000000,
                        rCut=5.0,
                        dt=0.005,
                        Etol=1e-5)
        # Save the result
        mycell.dump()

        # check if it worked
        if not ok:
            print("Failed to optimise!")
            sys.exit() # Exit the script

        # End of loop, so go around again

# The loop has finished so dump out the last state
mycell.dump()
