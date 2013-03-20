'''
Created on Jan 22, 2013

@author: abbietrewin
'''
import sys
import numpy
import util

a = numpy.array( [1,2,3], dtype=numpy.float64 )
b = numpy.array( [4,5,6], dtype=numpy.float64 )

assert not numpy.array_equal(a,b)


# Calculate normlised cross product to find an axis orthogonal 
#to both that we can rotate about
cross = numpy.cross(a, b)
ncross = cross / numpy.linalg.norm( cross )

if numpy.array_equal( cross, [0,0,0]):
    raise RuntimeError,"ZERO CROSS"

# Find angle
angle = util.vectorAngle(a, b)
print angle

rmat = util.rotation_matrix(ncross, angle)
print rmat
newb = numpy.dot( rmat, b )

print newb
print a / numpy.linalg.norm(a)
print newb / numpy.linalg.norm(newb)




#INCONFIG = "/Applications/DL_MONTE/input_for_DL_MONTE/cage4_n2-DL_MONTE/CONFIG_no_c"
#CONFIG = "/Applications/DL_MONTE/input_for_DL_MONTE/cage4_n2-DL_MONTE/CONFIG.out"
#ATOMS = ['n=','c1','hc_2', 'c2', 'NQ' ]
#ATOMS2 = ['cp_1', 'cp_2', 'hc_1', 'c=1','COM' ]
#
#f = open( INCONFIG, "r" )
#o = open( CONFIG, "w" )
#
#            
## First 7 lines just info - not needed
#for i in range(7):
#    line = f.readline()
#    o.writeXyz( line )
#
#line = f.readline()
#while line:
#    #line = line.strip()
#    fields = line.split()
#    label = fields[0]
#    if label in ATOMS:
#        newline = label + "      c\n"
#        o.writeXyz( newline )
#    else:
#        if label in ATOMS2:
#            newline = label + "    c\n"
#            o.writeXyz( newline )
#        else:
#            o.writeXyz( line )
#        
#    line = f.readline()
#
#f.close()
#o.close() 
#        
#
#        
#        