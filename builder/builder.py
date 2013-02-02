'''
Created on Jan 15, 2013

@author: abbietrewin



'''

import buildingBlock
import cell
import os

def newFile(filename):
    # Create a new filename using _1 etc
    name,suffix = os.path.splitext( filename )
    
    count=0
    while True:
        try:
            int(name[count-1])
        except ValueError:
            break
        count-=1
    
    nstr=name[count:]
            
    if not name[count-1] == "_" or count==0:
        raise RuntimeError,"Filename needs to be of the form: NAME_1.xyz"
    
    n=int(nstr)
    n=n+1
    name = name[:count]+str(n)
    
    return name+suffix


INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/pyrene_typed.car" 
INFILE = "/Users/jmht/Dropbox/Amorphousbuilder/pyrene_typed.car" 
#INFILE = "/Users/jmht/1.xyz" 
#INFILE = "/Users/abbietrewin/Dropbox/Amorphousbuilder/builder/ch4.xyz"
#OUTFILE1 = "/Users/abbietrewin/Dropbox/Amorphousbuilder/abbie_output/cell_1.xyz"
OUTFILE = "/Users/jmht/cell_1.xyz"
#OUTFILE2 = "/Users/abbietrewin/Dropbox/Amorphousbuilder/abbie_output/cell_2.xyz"
#OUTFILE2 = "/Users/jmht/cell_2.xyz"
nblocks = 45 
CELLA = [ 30,  0,  0 ]
CELLB = [ 0, 30,  0 ]
CELLC = [ 0,  0, 30 ]
STEPS = 100
NMOVES=50

BONDANGLE=180


# Create building block and read in car file
firstBlock = buildingBlock.BuildingBlock( infile = INFILE )

# Create Cell and seed it with the blocks
cell = cell.Cell( CELLA, CELLB, CELLC )
cell.seed (nblocks, firstBlock )
cell.write( OUTFILE )

# Loop through as many shimmy stages as required
while True:
    OUTFILE=newFile(OUTFILE)
    cell.shimmy( nsteps=STEPS, nmoves=NMOVES, bondAngle=BONDANGLE  )
    cell.write( OUTFILE )
    response = raw_input('Do we have to do this _again_? (y/n)')
    if response.lower() == 'n':
        break
    
            
#if __name__ == '__main__':
#    pass