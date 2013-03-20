'''
Created on Mar 7, 2013

@author: abbietrewin
'''
# This is the energy input file taken directly as outputted from dlmonte 
INENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3a/ENERGY_150_80"
#This outputs the total energy for the framework and the nitrogen as a function of iteration
OUTENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3a/TOT_ENERGY_OUT_150_80.csv"
#This is the first output file with the core framework atoms removed so that only the N2 for each step is included.
#ENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/ENERGY.out"
ATOMS = ['n=','c1','c2', 'cp_1', 'cp_2', 'hc_1', 'hc_2', 'c=1', 'NQ','COM', 'COM']
STEP = ['ITERATION']


# This reads in ENERGY.000 file and writes out the iteration and the total energy (frame and N2)

m = open( INENERGY, "r")
k = open( OUTENERGY, "w")
energies = []

iteration=None
inserts=None
for line in m:
    
    if line.startswith(" ITERATION"):
        if len(energies):
            tote = sum(energies)
            k.write("{0},{1},{2},{3}\n".format(iteration,inserts,iteration,tote))
            energies=[]
        
        # reset iteration
        fields = line.strip().split()
        iteration = fields[1]
        inserts = fields[2]
        continue
        
    if line.startswith(" MOLECULE"):
        continue
    
    fields = line.strip().split()
    energies.append( float(fields[5]) )

m.close()
k.close()

