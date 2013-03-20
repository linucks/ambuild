'''
Created on Feb 22, 2013

@author: abbietrewin
'''
INENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/ENERGY.000"
ENERGY = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/mol_Es/ENERGY_mol_2.out"
STRUCTURE = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/mol_Es/LOADING_mol_2.xyz"
ATOMS_CORE = ['n=','c1','c2', 'cp_1', 'cp_2', 'hc_1', 'hc_2', 'c=1']
STEP = ['ITERATION']
#Change step_2 for which iteration you want out
STEP_2 = ['15000000']
ATOMS_N2 = [ 'NQ', ]
ATOMS_N2_COM = ['COM' ]

def truncateEnergy( INENERGY, ENERGY):
    """
    This takes in the massive ENERGY output file and takes out all the core atoms 
    and writes out to ENERGY_mol.out
    """
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
        elif label in ATOMS_N2_COM:
            label = 'COM'
            newline = label + "   " + fields[2] + "   " + fields[3] + "   " + fields[4] + "   " + fields[5] + "\n"  
            o.write( newline )    
    
    f.close()
    o.close() 


def energy2loading( ENERGY, STRUCTURE ):
    """
    This section takes the ENERGY.out file and checks for the specified iteration 
    and writes this to a LOADING xyz file
    """

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

#I need to loop through N and COM and add the E's together


def splitEnergy( ENERGY ):
    """
    This section takes the ENERGY file and looks at the specified iteration 
    and splits the loading in to different energy regions
    """

    STRUCTURE_E_1 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/mol_Es/EBAND-1_mol.xyz"
    STRUCTURE_E_2 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/mol_Es/EBAND_2_mol.xyz"
    STRUCTURE_E_3 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/mol_Es/EBAND_3_mol.xyz"
    STRUCTURE_E_4 = "/Users/abbietrewin/Dropbox/CC4_private/CC3_run_4_analysis/run_1/mol_Es/EBAND_4_mol.xyz"
    
    f1 = []
    f2 = []
    f3 = []
    f4 = []
    header = None
    
    m = open( ENERGY, "r")
    
    line = m.readline()
    done=False
    while line:
        fields = line.split()
        label = fields[0]
        label_2 = fields[1]
        print "line b4:{}".format(line)
        if label in STEP and label_2 in STEP_2:
            print "line in step:{}".format(line)
            #line = m.readline()
            while not done:
                energies = []  
                for i in range(3):
                    print "line in 3:{}".format(line)
                    line = m.readline()
                    if not line:
                        done=True
                        break
                    fields = line.split()
                    if fields[0] == "ITERATION":
                        done=True
                        break
                    print fields[4] 
    #                energies.append(int( float(fields[4]) ))
                    energies.append( float(fields[4]) )
                    print energies
                    if i == 2:
                        # calculate summed energies
                        x = sum( energies )
                        print x
                        coordline = "{}     {: 10.6F}   {: 10.6F}   {: 10.6F}\n".format( "N", float(fields[1]), float(fields[2]), float(fields[3]), x )
                        if x <= -6 :
                            f1.append( coordline )
                        elif -6 < x <= -5.5:
                            f2.append( coordline )
                        elif -5.5 < x <= -5:
                            f3.append( coordline )
                        elif x > -5:   
                            f4.append( coordline )
            #End while
            if done:
                break
            
        #End if label in STEP and label_2 in STEP_2
        line = m.readline()
    
    
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


truncateEnergy( INENERGY, ENERGY )
print "after truncate"
energy2loading( ENERGY, STRUCTURE )
print "after loading"
splitEnergy( ENERGY )
print "after split"

