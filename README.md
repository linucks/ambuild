# Amorphous Builder (Ambuild)
Ambuild is a python program for creating polymeric molecular structures.

## Concepts
Ambuild works with fragments, blocks, endGroups and cells. The cell is a cubic volume of space that will be filled
with molecular blocks. The blocks themselves are containers for fragments. Fragments are molecules that can
contain one or more endGroups that may be used to form bonds. Fragments never exist by themselves in the cell,
they are always contained within blocks. A block may contain one or more fragments of the same or different types.
The fragments within a block are bonded to each other via their endGroups, and rules can be specified for which
endGroups can bond to which, the dihedral angle around the bond that the fragment should adopt, and any atoms that
may need to be removed on bonding.

### Specifying fragments and endGroups
The first stage is to specify the molecular fragments. This is done with two files:

* a ".car" file that contains the coordinates, atom-types etc.
* a ".csv" file that defines which atoms are endGroups, capAtoms etc.

The .car file is a standard Insight II file, which can be prepared with programs such as Avogadro or Materials Studio.
More information can be found about the format [here](http://chem5.nchc.org.tw/software/document2007/insightII/doc/formats980/File_Formats_1998.html#781840).

The .csv file should be named identically to the .car file, but with a .csv suffix instead of .car. The .csv
file is a standard comma-separated file that can be created and edited with any spreadsheet program. The first
line is a header that names the columns (case isn't important, but the field names must be correct):

```
Type,EndGroup,CapAtom,Dihedral,DelAtom
```

Each subsequent line defines an endGroup, with the columns separated by a comma. The first three columns are
required, the other two are optional.

NB: all indexes in the file start from 0 - this means that the first atom has an index of 0, the second 1, the third
2 etc.

The meanings of the columns are:

* **Type:** this is just a string (starting with a letter) that names the endGroup and is used to specify bonding rules
  and to control which endGroups in which building blocks are bonded in Grow or Join steps (see later).

* **EndGroup:** this is the index of the atom in the .car file that is to be the endGroup. The endGroups are the atoms
  that bond the fragments together. The indexing starts from 0, so if the first atom in the car file is to be the
  endGroup, the endGroup column should be zero. If it was the 6th atom, the column should have 5 in it.

* **CapAtom:** this is the atom that caps the endGroup and defines the vector that the bond will be be made along. The
  capAtom is removed when the bond is formed and the corresponding endGroup in the other block will replace the capAtom,
  although the other endGroup will probably not be in the same position, as length of the bond will be dependent on
  the types of the two endGroups.

* **Dihedral:** this is an optional column that specifies the atom that defines the dihedral angle around the bond between
  two endGroups. This can be used to specify the orientation of molecules when they are attached with a Grow step.
* **DelAtom:** this is an optional column that specifies an atom that will be removed when the endGroup forms a bond this
  can be used to unsaturate an atom on bonding. So for example, if the endGroup is the N of an NH2, then one of the H
  atoms could be the capAtom and the other the delatom, so that the on bonding to another NH2 endGroup a N=N bond would
  be formed.

The following sections describe the various operations that can be performed by ambuild. For an annotated example script,
see the builder.py file in the builder directory.

### Setting up the system
Ambuild creates its structures in a cubic cell. The first step is to create an empty cell object. This is done with
code like the following:

```
boxDim = [20,20,20]
mycell = cell.Cell(boxDim, atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15, doLog=False)
```

This creates an empty cell object called mycell (we call it mycell so it doesn't clash with cell, which is the name
we import the Cell code module with). The boxDim argument sets the size of the A,B and C cell axes.

The arguments to cell are:
Name | Description
boxDim | a list with three numbers specifying the size of the cell A,B and C dimensions (angstroms). Dimensions are from 0 - A, B or C
filePath | a (.car) file with the cell dimensions and the coordinates of molecules that will be kept static throught the simulation.
atomMargin | the additional distance that will be added on to the VdW radii of two atoms to determine if they are close enough to clash.
bondMargin | two atoms are considered close enough to bond if they are within the bond length defined for the two atoms +/- the bondMargin.
bondAngleMargin | the tolerance (in degrees) from the ideal of 180 that defines an acceptable bond
doLog | True/False - specifies if a log will be created - not recommended as it generates lots of data and slows the program.
                
Alternatively, the cell can be read from a .car file using the filePath argument instead of the boxDim argument. In
this case, the cell dimensions will be taken from the .car file PBC line, e.g.:

```
PBC   40.0000   30.0000   30.0000   90.0000   90.0000   90.0000 (P1)
```

And the molecule will be read in as a static structure - i.e. it will not be moved in any of the MD or geometry optimisation
steps. The name of this static will be the named after the .car file.

With the cell defined, the next stage is to add the fragments to the cell library so that they will become available
for use. This is done with the following command:

```
mycell.libraryAddFragment( filename=fragA, fragmentType='A' )
```

filename - the path to the .car file. There will need to be a corresponding .csv file that defines
         - the endGroups, capAtoms etc.
fragmentType - a name that will be used to identify the fragment - cannot contain the ":" or "-" characters
solvent - specify that this fragmentType is solvent and so won't be clash-checked in Zip steps

As many fragments can be added as desired, the only requirement is that they are named differently.

With the fragments added, we specify how they may bond. This is done with the following command:

```
mycell.addBondType( 'A:a-B:a' )
```

endGroups are defined by the fragmentType they belong to (which is set by the fragmentType argument
to libraryAddFragment), together with the identifier for that endGroup (which is specified by the first column
in the .csv file). These are separated by a colon, so an endGroup identifier is of the form:

```
FRAGMENT1:ENDGROUP1
```

A bond is defined by two endGroups, separated by a hyphen, so a bond identifier has the form:

```
FRAGMENT1:ENDGROUP1-FRAGMENT2:ENDGROUP2
```

In the above example, the block was added with a fragmentType of 'A' and the first column of the endGroup would have
been 'a'.

Sometimes it is necessary to limit the number of bonds of a particular type to an individual fragment. For example if
a fragment has three nitrogen endGroups, but once one is used, the others become unavailable. This is achieved with
the setMaxBonds argument as shown below:

```
mycell.addBondType( bondType='A:a-B:a', count=1 )
```

The arguments are:

bondType - the bondType (FRAGMENT1:ENDGROUP1-FRAGMENT2:ENDGROUP2) as was specified with the call
           to addBondType
count - the maximum number of permissible bonds for a single fragment.

This defines the system. We can then go about filling it with blocks in the following section.

### Creating a polymer
This is the simplest and most automatic way to create a structure:

```
growPolymer(monomers=['A','B;], ratio=[1,1], length=10, random=False, center=False)
```

Ths will create a linear polymer.

Arguments:
monomers - list of ambuild fragmentTypes (as specified when adding the fragment to the cell with libraryAddFragment)
           that will be joined to create a subunit
ratio - list of integers specifying the number of monomers of each type (needs to be the same length as monomers) in
        the subunit
length - the number of subumits that will be created
random - True/False - whether to build up the polymer deterministically (as per monomers/ratio) or stochastically, but
       - where the final ratio of fragments will approach that specified in ratio.
center - place the first monomer polymer in the center of the cell (if possible)

The fragment and ratio lists define the construction of the subunit. E.g.:
```
monomers = ['A', 'B']
ratio = [1,2]
gives: [ABB]
```

We could even have:
```
monomers = ['A', 'B', 'A', 'B']
ratio = [1,2,3,1]
gives: [ABBAAAB]
```

### Filling the cell with blocks
The first stage is to seed the system with molecular building blocks. This is done with the seed command as shown
below:

```
seed(100, fragmentType='A', maxTries=500, center=False )
```

The arguments to seed are:

nblocks - the number of blocks to add.
fragmentType - the type of fragments to add. If fragment is None, or omitted, then fragments will be randomly
               chosen from the library.
maxTries - the number of attempts to make when adding a block before the seed step is fails and returns.
center - True/False - if True, place the first block in the center of the cell.
point - a list of three floats defiting a point around which the centroids of the blocks will be seeded
        (requires the radius argument).
radius - a single float specifying the radius of a sphere around point, within which the centroids of 
         the blocks will be seeded (requires point argument).
zone - a list of 6 floats specifying a box within the cell within which the centroids of the blocks will
       be seeded.

Once the system is seeded, Grow or Join steps are used either add new blocks to the system, or join existing blocks
together.

Growing a block involves adding a fragment to an existing block in the cell. This is done with the command:

```
mycell.growBlocks( 10, cellEndGroups=['A:a'], libraryEndGroups=['B:a'], dihedral=90, maxTries=500 )
```

The arguments to growBlocks are:

toGrow: number of blocks to add
cellEndGroups: a list of the endGroup types in the cell that the new blocks will be bonded to.
              If more than one endGroup type is supplied, the endGroup will be randomly chosen from
              that list.
libraryEndGroups: a list of the endGroup types from the library that will be used form the bonds.
              If more than one endGroup type is supplied, the endGroup will be randomly chosen from
              that list.
dihedral: the dihedral angle about the bond (3rd column in csv file)
maxTries: number of attempts to make before giving up


### Joining blocks in the cell together
#### joinBlocks

To join existing blocks in the cell together use the joinBlocks command:

```
mycell.joinBlocks( 10, cellEndGroups=['A:a','B:a'], dihedral=90, maxTries=500 )
```

This will randomly select two blocks (with endGroups as specified in the list supplied to cellEndGroups, or just
two random blocks if cellEndGroups is None) and attempt to bond them by removing the second block from the cell
and attaching it to the first.

If the bond fails (due to clashes etc), the second block will be returned to its original position.

The arguments to joinBlocks are:

toJoin - number of blocks to join
cellEndGroups - a list of the different endGroupTypes that should be bonded. If this is None
                randomly chosen endGroups will be used.
dihedral: the dihedral angle about the bond (3rd column in csv file)
maxTries - the maximum number of moves to try when joining

#### zipBlocks
Join existing blocks in the cell by changing the bondMargin and bondAngleMargin parameters that were
specified when the cell was created, and then looping over all the free endGroups to see if any can bond
with the new parameters. The blocks are not moved in this step.

The arguments to zipBlock are:
bondMargin - the new bondMargin [degrees]
bondAngleMargin - the new bondAngleMargin [degrees] 
clashCheck - True/False check for clashes between the bond and any atoms that fall within a cylinder
             of radius clashDist (default=1.6A) centered on the bond axis.
clashDist  - a float specifying the perpendicular distance from the bond axis that determines if an atom 
             is clashing with the bond. clashDist needs to be < the the cell box size as otherwise we won't
             see the atom.
selfBond  - boolean to specify if zip will allow a block to bond to itself (True) or not (False) [default: True]

zipBlocks is called thus:

```
zipBlocks(bondMargin=5, bondAngleMargin=30, clashTest=True, clashDist=1.6)
```

### Deleting Blocks
Blocks can be removed from the cell with two commands:

deleteBlocksIndices
Arguments:
indices: list of indices of the blocks to be deleted
save: save all deleted blocks so they can be readded with the restoreBlocks command

deleteBlocksTypes
delete numBlocks of blocks of the specified fragment type

Arguments:
fragmentType: the fragmentType, or list of fragmentTypes of the blocks to remove
numBlocks: optional - the number of blocks to remove, otherwise all blocks of the specified types will be removed
multiple: if True remove blocks that contain > 1 fragment, else only single-fragment blocks
save: save all deleted blocks so they can be readded with the restoreBlocks command

### Optimising the geometry of the blocks in the cell and running molecular dynamics simulations
Once the cell has been filled with blocks, it is possible to either optimise the geometry or run a molecular dynamics
simulation.

#### Geometry Optimisation
Geometry optimisations and MD runs are carried out using [HOOMD-blue](http://codeblue.umich.edu/hoomd-blue)

Rigid body optimisation is carried out using the mode_minimize_rigid_fire optimiser. The optimiser is run with a command
similar to the following:

```
mycell.optimiseGeometry( rigidBody=True,
                         doDihedral=True,
                         quiet=True,
                         rCut=5.0,
                         optCycles = 1000000,
                         dt=0.005,
                         Etol=1e-5,
                         )
```

The arguments accepted by optimiseGeometry are:

rigidBody - True/False - do rigid body or all-atom optimisation
doDihderal - True/False - include dihedral terms
doImproper - True/False - include improper terms
rCut - the VdW cutoff to use [angstroms]
optCycles - the number of hoomdblue optimisation cycles to run.
quiet - True/False - don't print out the normal hoomdblue runtime output to the screen.

All other arguments accepted by the mode_minimize_rigid_fire are passed through, so for example, you can set fdec, alpha_start
etc. in the call to optimiseGeometry. The arguments accepted by mode_minimize_rigid_fire are listed in the hoomdblue [documentation](http://codeblue.umich.edu/hoomd-blue/doc/classhoomd__script_1_1integrate_1_1mode__minimize__rigid__fire.html)

#### Molecular Dynamics
MD runs are run using the nvt, nvt_rigid, npt, or npt_rigid integrators within hoomdblue.

The documentation for the integrators is [here](http://codeblue.umich.edu/hoomd-blue/doc/page_command_list.html#sec_index_integration)

To run MD, a call like the following can be used:

```
mycell.runMD(doDihedral=True,
             quiet=False,
             rCut=5.0,
             mdCycles=100,
             T=1.0,
             tau=0.5,
             P=1.0,
             tauP=0.5,
             dt=0.0005,
             dump=True,
             dumpPeriod=20,
             integrator='npt',
            )
```

The integrator option can either by 'npt' or 'nvt'.

#### MD and Geometry optimisation combined
In order to save copying data between ambuild and hoomdblue, a combined MD and Optimisation command has been created. This runs
a very short optimisation loop to stabilise the structure, then runs an MD simulation, before running a full geometry optimisation
to convergence (if possible).

An example command to do this is:

```
mycell.runMDAndOptimise( doDihedral=True )
```

All the arguments accepted by optimiseGeometry and runMD are accepted by runMDAndOptimise.

### Saving Results and extracting coordinates
By default, ambuild doesn't save any data or coordinates. The dump command is used to save the current state of the cell.

This is done as follows:

```
mycell.dump()
```

This will create a python pickle file (a serialised representation of the entire cell object) that will be named step_X.pkl
where X is a number that will increment each time that dump is called.

In order to extract coordinates the util.py script in the builder directory can be used. The script should be called with
the path to the pkl file as the first argument, so if the build was taking place in the same directory, the command would be:

```
./util.py step_4.pkl
```

The util script will then create an XYZ file of the coordinates and a CML file suitable for viewing with Avogadro. The CML
file is periodic, but the bonds that extend across the periodic boundary have been omitted so this file SHOULD NOT be used
for extracting structural data - it is just for visualisation.

The util script also accepts a "--fragments" or "-f" argument, which tells the script to write out files for each fragment separately.
This is done as follows:

```
./util.py -f step_4.pkl
```

The util script also accepts a "--blocks" or "-b" argument, which tells the script to write out separate files for each block in
the cell. This is done as follows:

```
./util.py -b step_4.pkl
```



