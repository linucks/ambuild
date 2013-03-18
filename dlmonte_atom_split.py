'''
Created on Feb 20, 2013

@author: abbietrewin
'''
# This is the energy input file taken directly as outputted from dlmonte 
INENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/ENERGY.000"
#This is the first output file with the core framework atoms removed so that only the N2 for each step is included.
ENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/ENERGY.out"
#This is the output file for the iteration that is specified in STEP_2 outputted as an xyz file 
STRUCTURE = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/LOADING.xyz"
# The following files are the N and COM positions split into four files dependent upon their energies 
STRUCTURE_E_1 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/EBAND-1.xyz"
STRUCTURE_E_2 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/EBAND_2.xyz"
STRUCTURE_E_3 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/EBAND_3.xyz"
STRUCTURE_E_4 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/EBAND_4.xyz"
ATOMS_CORE = ['n=','c1','c2', 'cp_1', 'cp_2', 'hc_1', 'hc_2', 'c=1']
STEP = ['ITERATION']
#Change step_2 for which iteration you want out
STEP_2 = ['15000000']
ATOMS_N2 = [ 'NQ','COM' ]
ATOMS_N2_COM = ['COM' ]

# This section reads in the full ENERGY file output as is from dlmonte and returns the positions of the 
# N atoms and COM as N's (COM returned as N's so that it can be viewed). writes to ENERGY.out
f = open( INENERGY, "r" )
o = open( ENERGY, "w" )
          
for line in f:
    fields = line.split()
    label = fields[0]
    if label in STEP:
        newline = label + "   " + fields[1] + "   " + fields[2] + "   " + fields[3] + "\n" 
        o.write( newline )
    elif label in ATOMS_N2:
        label = 'N' 
        newline = label + "   " + fields[2] + "   " + fields[3] + "   " + fields[4] + "   " + fields[5] + "\n"  
        o.write( newline )

f.close()
o.close() 

# This section takes the ENERGY.out file and checks for the specified iteration and writes this to a LOADING file

m = open( ENERGY, "r")
k = open( STRUCTURE, "w")
  

# Set flags and counts here
wrote=0
writing=False

for line in m:
    fields = line.split()
    label = fields[0]
    label_2 = fields[1]
    
    if label in STEP and label_2 in STEP_2:
        # Found the iteration we are after - work out how many lines to read
        count = int( fields[2] ) * 3 - 3
 # http://docs.python.org/2/library/string.html#format-examples
        newline = "{}\n{}  {}  {}  {}\n".format(count, fields[0], fields[1], fields[2], fields[3] )
        k.write( newline )
        # set writing to true so we know that we are in writing mode
        writing=True
        # continue skips to the next line as we don't want to write out the header line - I think
        continue
    
    # This bit will only be executed when we are writing the bits we want
    if writing:
        # So, need to convert string numbers to float
        # space after the colon: "indicates that a leading space should be used on positive numbers, and a minus sign on negative numbers."
        # 10 to force a minimum width of 10 digits, the 6 to use 6 decimal places and F to write out in decimal notation (E for scientific)
        # If you have more digits in front of the decimal place, you may need to extend the minimum width to get things to align
        newline =  "{}     {: 10.6F}   {: 10.6F}   {: 10.6F}\n".format( fields[0], float(fields[1]), float(fields[2]), float(fields[3]) )
        k.write( newline )
        # increment wrote so we know we have written out a line
        wrote+=1
        # Check if we have got to the end
        if wrote >= count:
            # We've written out all the lines we wanted to say that we've finished
            writing=False
            # If you want to stop reading the file here, uncomment the below line - otherwise it will carry on reading the rest
            # of the file
            break
        
m.close()
k.close()

# This section takes the ENERGY file and looks at the specified iteration and splits the loading in to different
# energy regions

f1 = []
f2 = []
f3 = []
f4 = []
header = None

m = open( ENERGY, "r")
reading=False
for line in m:
    fields = line.split()
    label = fields[0]
    label_2 = fields[1]
    count = 0
    
    if label in STEP and reading:
        break
    
    if label in STEP and label_2 in STEP_2:
        header = "{}  {}  {}  {}".format(count, fields[0], fields[1], fields[2], fields[3] )
        # set writing to true so we know that we are in writing mode
        reading=True
        # continue skips to the next line as we don't want to write out the header line - I think
        continue
    
    # This bit will only be executed when we are writing the bits we want
    if reading:
        x = int( float(fields[4]) )
        if x <= -4 :
            f1.append("{}     {: 10.6F}   {: 10.6F}   {: 10.6F}\n".format( fields[0], float(fields[1]), float(fields[2]), float(fields[3]) ))
        elif -4 < x <= 0:
            f2.append("{}     {: 10.6F}   {: 10.6F}   {: 10.6F}\n".format( fields[0], float(fields[1]), float(fields[2]), float(fields[3]) ))    
        elif 0 < x <= 4:
            f3.append("{}     {: 10.6F}   {: 10.6F}   {: 10.6F}\n".format( fields[0], float(fields[1]), float(fields[2]), float(fields[3]) ) )   
        elif x >= 4:   
            f4.append("{}     {: 10.6F}   {: 10.6F}   {: 10.6F}\n".format( fields[0], float(fields[1]), float(fields[2]), float(fields[3]) ) )



l = open( STRUCTURE_E_1, "w")
# Write out number of lines
l.write( "{}\n".format( len(f1) ) )
# header
l.write( "{}\n".format( header ) )
# each coordinate
for line in f1:
    l.write( line )
    
j = open( STRUCTURE_E_2, "w")
# Write out number of lines
j.write( "{}\n".format( len(f2) ) )
# header
j.write( "{}\n".format( header ) )
# each coordinate
for line in f2:
    j.write( line )
    
p = open( STRUCTURE_E_3, "w")
# Write out number of lines
p.write( "{}\n".format( len(f3) ) )
# header
p.write( "{}\n".format( header ) )
# each coordinate
for line in f3:
    p.write( line )
    
q = open( STRUCTURE_E_4, "w")
# Write out number of lines
q.write( "{}\n".format( len(f4) ) )
# header
q.write( "{}\n".format( header ) )
# each coordinate
for line in f4:
    q.write( line )


m.close()
l.close()
j.close()
p.close()
q.close()

     