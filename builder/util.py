'''
Created on Feb 3, 2013

@author: jmht

Utility functions
'''
import os

def newFilename(filename):
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