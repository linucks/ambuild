#!/opt/hoomd-0.11.3/hoomdblue-install/bin/hoomd
#!/Applications/HOOMD-blue.app/Contents/MacOS/hoomd

import math
import os
import sys
import time
import xml.etree.ElementTree as ET

import util
import hoomdblue

class FfieldParameters( object ):
    
    
    def __init__(self):
        
        # REM = bonds are in alphabetical order as in the order of the atoms in a bond
        # e.g. c-h NOT h-c
        self.bonds = { 
                      'cp-cp'   : { 'k' : 1200.0, 'r0' : 1.54 },
                      'cp-ct'   : { 'k' : 1200.0, 'r0' : 1.54 },
                      'c3a-c3a' : { 'k' : 1200.0, 'r0' : 1.54 },
                      'c-c'     : { 'k' : 1200.0, 'r0' : 1.54 },
                      'c-h'     : { 'k' : 1200.0, 'r0' : 1.09 },
                      'cp-hx'   : { 'k' : 1200.0, 'r0' : 1.09 },
                      'ct-hx'   : { 'k' : 1200.0, 'r0' : 1.09 },
                      }

        self.angles = { 
                       #2*math.pi/3 = 120     
                       'cp-cp-cp'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       
                       'c-c-h'       : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'c-c-c'       : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'c3a-c3a-c3a' : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'cp-cp-ct'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'ct-cp-cp'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'cp-cp-ct'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'cp-cp-hx'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'cp-cp-np'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'np-cp-cp'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },
                       'cp-cp-cp'    : { 'k' : 330.0, 't0' : 2*math.pi/3 },


                       #math.pi = 180
                       'ct-ct-cp'    : { 'k' : 330.0, 't0' : math.pi },
                       'cp-ct-ct'    : { 'k' : 330.0, 't0' : math.pi },
                       'ct-ct-hx'    : { 'k' : 330.0, 't0' : math.pi },
                       'ct-ct-cp'    : { 'k' : 330.0, 't0' : math.pi },
                       'cp-ct-nt'    : { 'k' : 330.0, 't0' : math.pi },
                      }

        self.dihedrals = { 
                       'cp-cp-cp-cp'       : { 'k' : 200, 'd' : -1, 'n' : 0 }, # Jens just for testing
                       'cp-cp-cp-np'       : { 'k' : 200, 'd' : -1, 'n' : 2 },
                       'np-cp-cp-cp'       : { 'k' : 200, 'd' : -1, 'n' : 2 },
                       'cp-cp-ct-nt'       : { 'k' : 0, 'd' : -1, 'n' : 2 },
                       'nt-ct-cp-cp'       : { 'k' : 0, 'd' : -1, 'n' : 2 },
                         }
        
        self.impropers = { 
                       'cp-cp-cp-np'       : { 'k' : 200, 'chi' : 0.0 },
                       'np-cp-cp-cp'       : { 'k' : 200, 'chi' : 0.0 },
                       'cp-cp-ct-nt'       : { 'k' : 200, 'chi' : 0.0 },
                       'nt-ct-cp-cp'       : { 'k' : 200, 'chi' : 0.0 },
                         }

        self.pairs = {
        
              # Ones with x in are cap atoms and are ignored - see setPair
              
              #PCFF typed pair potentials
              ('c',  'hc' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              
              ('c2', 'c2' )      : { 'epsilon' : 0.0933,   'sigma' : 3.7736  },
              ('c2', 'ct' )      : { 'epsilon' : 0.0760,   'sigma' : 3.9007  },
              ('c2', 'cp' )      : { 'epsilon' : 0.1016,   'sigma' : 3.7736  },
              ('c2', 'hc' )      : { 'epsilon' : 0.0346,   'sigma' : 3.4893  },
              ('c2', 'nb' )      : { 'epsilon' : 0.0759,   'sigma' : 3.9357  },
              ('c2', 'hx' )      : { 'epsilon' : 0.0346,   'sigma' : 3.4893  },
              
              ('c1', 'c1' )      : { 'epsilon' : 0.0933,   'sigma' : 3.7736  },
              ('c1', 'c2' )      : { 'epsilon' : 0.0933,   'sigma' : 3.7736  },
              
              ('c=1', 'c=1' )    : { 'epsilon' : 0.0640,   'sigma' : 4.0100  },
              ('c=1', 'c1' )     : { 'epsilon' : 0.0760,   'sigma' : 3.9007  },
              ('c=1', 'hc' )     : { 'epsilon' : 0.0254,   'sigma' : 3.6691  },
              ('c=1', 'c2' )     : { 'epsilon' : 0.0760,   'sigma' : 3.9007  },
              ('c=1', 'cp' )     : { 'epsilon' : 0.0828,   'sigma' : 3.9007  },
              ('c=1', 'np' )     : { 'epsilon' : 0.0828,   'sigma' : 3.9007  },#needs changing
              
              ("c3'", 'c4' )     : { 'epsilon' : 0.0933,   'sigma' : 3.7736  },
              
              ('c4', 'c4' )      : { 'epsilon' : 0.0933,   'sigma' : 3.7736  },
              ('c4', 'c4o' )     : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('c4', 'o2e' )     : { 'epsilon' : 0.1598,   'sigma' : 3.6640  },
              ('c4', 'h1' )      : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c4', 'n3a' )     : { 'epsilon' : 0.1187,   'sigma' : 3.9357  },
              ('c4', 'n2a' )     : { 'epsilon' : 0.1187,   'sigma' : 3.9357  },
             # ('c4', 'zn+2' )    : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('c4', 'p' )       : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c4', 'f1p' )     : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c4', 'c3a' )     : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              
              ('c4o', 'c4o' )    : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('c4o', 'c3a' )    : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('c4o', 'o2e' )    : { 'epsilon' : 0.1598,   'sigma' : 3.6640  },
              ('c4o', 'h1' )     : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c4o', 'n3a' )    : { 'epsilon' : 0.1187,   'sigma' : 3.9357  },
              ('c4o', 'n2a' )    : { 'epsilon' : 0.1187,   'sigma' : 3.9357  },
             # ('c4o', 'zn+2' )   : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('c4o', 'p' )      : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c4o', 'f1p' )    : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              
              ('c3a', 'c3a' )    : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('c3a', 'n2a' )    : { 'epsilon' : 0.1187,   'sigma' : 3.9357  },
              ('c3a', 'o2e' )    : { 'epsilon' : 0.1598,   'sigma' : 3.6640  },
              ('c3a', 'h1' )     : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c3a', 'n3a' )    : { 'epsilon' : 0.1187,   'sigma' : 3.9357  },
             # ('c3a', 'zn+2' )   : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('c3a', 'p' )      : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c3a', 'f1p' )    : { 'epsilon' : 0.0376,   'sigma' : 3.4893  },
              ('c3a', 'o2e' )    : { 'epsilon' : 0.1598,   'sigma' : 3.6640  },
              
              ('cp', 'cp' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('cp', 'hc' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('cp', 'ct' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('cp', 'nb' )      : { 'epsilon' : 0.0827,   'sigma' : 3.9357  },
              ('cp', 'nt' )      : { 'epsilon' : 0.0836,   'sigma' : 3.6788  },
              ('cp', 'hx' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('cp', 'np' )      : { 'epsilon' : 0.0664,   'sigma' : 3.6788  },
              ('cp', 'oc' )      : { 'epsilon' : 0.1598,   'sigma' : 3.6640  },
              ('cp', 'c1' )      : { 'epsilon' : 0.1016,   'sigma' : 3.7736  },
              
              ('ct', 'ct' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('ct', 'hc' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('ct', 'nb' )      : { 'epsilon' : 0.0644,   'sigma' : 4.0406  },
              ('ct', 'np' )      : { 'epsilon' : 0.0644,   'sigma' : 4.0406  },
              ('ct', 'nt' )      : { 'epsilon' : 0.0608,   'sigma' : 3.8214  },
              ('ct', 'hx' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              
              ('o2e', 'o2e' )    : { 'epsilon' : 0.2400,   'sigma' : 3.5350  },
              ('o2e', 'h1' )     : { 'epsilon' : 0.0615,   'sigma' : 3.3189  },
              ('o2e', 'n3a' )    : { 'epsilon' : 0.1644,   'sigma' : 3.8484  },
              ('o2e', 'n2a' )    : { 'epsilon' : 0.1644,   'sigma' : 3.8484  },
              #('o2e', 'zn+2' )   : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('o2e', 'p' )      : { 'epsilon' : 0.0615,   'sigma' : 3.3189  },
              ('o2e', 'f1p' )    : { 'epsilon' : 0.0615,   'sigma' : 3.3189  },
              
              ('oc', 'oc' )      : { 'epsilon' : 0.2400,   'sigma' : 3.5350  },
              ('oc', 'np' )      : { 'epsilon' : 0.0992,   'sigma' : 3.5527  },
              ('oc', 'nt' )      : { 'epsilon' : 0.1248,   'sigma' : 3.5527  },
              
              ('h1', 'h1' )      : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
              ('h1', 'n3a' )     : { 'epsilon' : 0.0356,   'sigma' : 3.7161  },
              ('h1', 'n2a' )     : { 'epsilon' : 0.0356,   'sigma' : 3.7161  },
              ('h1', 'p' )       : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
              ('h1', 'f1p' )     : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
              #('h1', 'zn+2' )    : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
              
              ('hc', 'hc' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('hc', 'nb' )      : { 'epsilon' : 0.0248,   'sigma' : 3.7161  },
              ('hc', 'nt' )      : { 'epsilon' : 0.0316,   'sigma' : 3.3431  },
              ('hc', 'hx' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('hc', 'np' )      : { 'epsilon' : 0.0251,   'sigma' : 3.3431  },
              ('hc', 'oc' )      : { 'epsilon' : 0.0615,   'sigma' : 3.3189  },
              ('hc', 'c1' )      : { 'epsilon' : 0.0346,   'sigma' : 3.4893  },
              
              ('hx', 'hx' )      : { 'epsilon' : 0.1106,   'sigma' : 3.7736  },
              ('hx', 'nb' )      : { 'epsilon' : 0.0248,   'sigma' : 3.7161  },
              
              ('n3a', 'n3a' )    : { 'epsilon' : 0.1340,   'sigma' : 4.070   },
              ('n3a', 'n2a' )    : { 'epsilon' : 0.1340,   'sigma' : 4.070   },
             # ('n3a', 'zn+2' )   : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('n3a', 'p' )      : { 'epsilon' : 0.0356,   'sigma' : 3.7161  },
              ('n3a', 'f1p' )    : { 'epsilon' : 0.0356,   'sigma' : 3.7161  },
              
              ('n2a', 'n2a' )    : { 'epsilon' : 0.1340,   'sigma' : 4.070   },
              #('n2a', 'zn+2' )   : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('n2a', 'p' )      : { 'epsilon' : 0.0356,   'sigma' : 3.7161  },
              ('n2a', 'f1p' )    : { 'epsilon' : 0.0356,   'sigma' : 3.7161  },
              
              ('nb', 'nb' )      : { 'epsilon' : 0.0650,   'sigma' : 4.0700  },
              ('nb', 'np' )      : { 'epsilon' : 0.0479,   'sigma' : 3.8600  },
              ('nb', 'ct' )      : { 'epsilon' : 0.0644,   'sigma' : 4.0406  },
              ('nb', 'nt' )      : { 'epsilon' : 0.0603,   'sigma' : 3.8600  },
              
              ('np', 'np' )      : { 'epsilon' : 0.0410,   'sigma' : 3.5700  },
              ('np', 'nt' )      : { 'epsilon' : 0.0516,   'sigma' : 3.5700  },
              ('np', 'c1' )      : { 'epsilon' : 0.1016,   'sigma' : 3.7736  },#needs changing
              ('np', 'c2' )      : { 'epsilon' : 0.1016,   'sigma' : 3.7736  },#needs changing
              
              ('nt', 'nt' )      : { 'epsilon' : 0.0650,   'sigma' : 3.5700  },
              
              ('n=', 'n=' )      : { 'epsilon' : 0.1382,   'sigma' : 3.5759  },
              ('n=', 'cp' )      : { 'epsilon' : 0.1220,   'sigma' : 3.6814  },
              ('n=', 'np' )      : { 'epsilon' : 0.1220,   'sigma' : 3.6814  },#needs changing
              ('n=', 'hc' )      : { 'epsilon' : 0.0459,   'sigma' : 3.3472  },
              ('n=', 'c=1' )     : { 'epsilon' : 0.0888,   'sigma' : 3.8235  },
              ('n=', 'c1' )      : { 'epsilon' : 0.1121,   'sigma' : 3.6814  },
              ('n=', 'c2' )      : { 'epsilon' : 0.1121,   'sigma' : 3.5759  },
              
              ('zn+2', 'zn+2' )  : { 'epsilon' : 0.25106,  'sigma' : 3.296   },
              ('zn+2', 'p' )     : { 'epsilon' : 0.02,     'sigma' : 2.995   },
              ('zn+2', 'f1p' )   : { 'epsilon' : 0.02,     'sigma' : 2.995   },
              
              ('p', 'p' )        : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
              ('p', 'f1p' )      : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
              
              ('f1p', 'f1p' )    : { 'epsilon' : 0.02,     'sigma' : 2.9950  },
        #Zn ghost atom
         #     ('zn+2', 'zn+2' )  : { 'epsilon' : 0.00001,  'sigma' : 0.00001 },
          #    ('zn+2', 'hc' )    : { 'epsilon' : 0.00001,  'sigma' : 0.00001 },
           #   ('zn+2', 'nb' )    : { 'epsilon' : 0.00001,  'sigma' : 0.00001 },
            #  ('zn+2', 'cp' )    : { 'epsilon' : 0.00001,  'sigma' : 0.00001 },
             # ('zn+2', 'c2' )    : { 'epsilon' : 0.00001,  'sigma' : 0.00001 },
              #('zn+2', 'ct' )    : { 'epsilon' : 0.00001,  'sigma' : 0.00001 },
            }
        
        return
    
    def angleParameter( self, angle ):
        return self.angles[ angle ]
    
    def hasAngle(self, angle):
        return angle in self.angles.keys()
    
    def bondParameter( self, bond ):
        return self.bonds[ bond ]
    
    def hasBond(self, bond):
        return bond in self.bonds.keys()
  
    def dihedralParameter( self, dihedral ):
        return self.dihedrals[ dihedral ]
    
    def hasDihedral(self, dihedral):
        return dihedral in self.dihedrals.keys()
    
    def improperParameter( self, improper ):
        return self.impropers[ improper ]
    
    def hasImproper(self, improper):
        return improper in self.impropers.keys()
    
    def pairParameter( self, p1, p2 ):
        """ Return whichever pair is defined """
        if (p1, p2) in self.pairs:
            return self.pairs[ (p1, p2) ]
        if (p2, p1) in self.pairs:
            return self.pairs[ (p2, p1) ]
        
        assert False,"Could not find {0} {1}".format( p1, p2 )
        return
    
    def hasPair(self, p1, p2):
        """ dummy atoms treated specially """
        if p1.lower() == 'x' or p2.lower() == 'x':
            return True
        if self.pairs.has_key( (p1, p2) ) or self.pairs.has_key( (p2, p1) ):
            return True
        return False

class HoomdOptimiser( object ):
    
    def __init__(self):
        
        self.ffield = FfieldParameters()
        self.system = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        self.impropers = None
        self.atomTypes = None
        self.rCut = 5.0
        return

    def checkParameters(self, xmlFilename=None):
        
        assert self.ffield
        self.setAttributesFromFile( xmlFilename )
        #assert self.bonds
        #assert self.angles
        assert self.atomTypes
        
        ok = True
        missingBonds = []
        for bond in self.bonds:
            if not self.ffield.hasBond( bond ):
                ok = False
                missingBonds.append( bond )
        missingAngles = []
        for angle in self.angles:
            if not self.ffield.hasAngle( angle ):
                ok = False
                missingAngles.append( angle )
        missingDihedrals = []
        for dihedral in self.dihedrals:
            if not self.ffield.hasDihedral( dihedral ):
                ok = False
                missingDihedrals.append( dihedral )
        missingImpropers = []
        for improper in self.impropers:
            if not self.ffield.hasImproper( improper ):
                ok = False
                missingImpropers.append( improper )
        missingPairs = []
        for i, atype in enumerate( self.atomTypes ):
            for j, btype in enumerate( self.atomTypes ):
                if j >= i:
                    if not self.ffield.hasPair( atype, btype ):
                        ok = False
                        missingPairs.append( ( atype, btype ) )
        
        if not ok:
            msg = "The following parameters could not be found:\n"
            if missingBonds:
                msg += "Bonds: {0}\n".format( missingBonds )
            if missingAngles:
                msg += "Angles: {0}\n".format( missingAngles )
            if missingDihedrals:
                msg += "Dihedrals: {0}\n".format( missingDihedrals )
            if missingImpropers:
                msg += "Impropers: {0}\n".format( missingImpropers )
            if missingPairs:
                msg += "Pairs: {0}\n".format( missingPairs )
            
            msg += "Please add these to the opt.py file\n"
            
            raise RuntimeError,msg
        
        return

    def writeCar( self, system, filename=None, unwrap=True, pbc=True ):
        """Car File
        """
        
        car = "!BIOSYM archive 3\n"
        car += "PBC=ON\n"
        
        car += "ambuild generated car file\n"
        tstr = time.strftime( "%a %b %d %H:%M:%S %Y", time.gmtime() )
        car += "!DATE {0}\n".format( tstr )
        
        xdim = system.box[0]
        ydim = system.box[1]
        zdim = system.box[2]
        if pbc:
            car += "PBC  {0: < 9.4F} {1: < 9.4F} {2: < 9.4F}  90.0000   90.0000   90.0000 (P1)\n".format( xdim,
                                                                                                      ydim,
                                                                                                      zdim )
        
        for p in system.particles:
            
                label = atype = p.type.strip()
                
                # Treat x-atoms differently
                if label[0].lower() == 'x':
                    symbol = 'x'
                else:
                    symbol = util.label2symbol( label )
                
                x, y, z  = p.position
                ix, iy, iz = p.image
                charge = float( p.charge )
                
                if unwrap:
                    x = util.unWrapCoord( x, ix, xdim, centered=False )
                    y = util.unWrapCoord( y, iy, ydim, centered=False )
                    z = util.unWrapCoord( z, iz, zdim, centered=False )
                else:
                    # Put back with origin at corner
                    x = x + ( xdim / 2 )
                    y = y + ( ydim / 2 )
                    z = z + ( zdim / 2 )
                
                car += "{0: <5} {1: >15.10} {2: >15.10} {3: >15.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format( label, x, y, z, atype, symbol, charge )
        
        car += "end\nend\n\n"
        
        with open( filename, 'w' ) as f:
            f.writelines( car )
            
        return
    
    def writeXyz( self, system, filename=None, unwrap=True ):
        """Car File
        """
        
        xyz = "{0}\n".format( len( system.particles ) )
        xyz += "opt.xyz file\n"
        
        xdim = system.box[0]
        ydim = system.box[1]
        zdim = system.box[2]
        
        for p in system.particles:
            
                label = p.type.strip()
                x, y, z  = p.position
                ix, iy, iz = p.image
                
                if unwrap:
                    x = util.unWrapCoord( x, ix, xdim, centered=False )
                    y = util.unWrapCoord( y, iy, ydim, centered=False )
                    z = util.unWrapCoord( z, iz, zdim, centered=False )
                else:
                    # Put back with origin at corner
                    x = x + ( xdim / 2 )
                    y = y + ( ydim / 2 )
                    z = z + ( zdim / 2 )
                
                xyz += "{0:<5} {1: > 14.10} {2: > 14.10} {3: > 14.10} XXXX 1      {4: <4}    {5: <2} {6: > 2.3f}\n".format( label, x, y, z )
        
        with open( filename, 'w' ) as f:
            f.writelines( xyz )
            
        return
    
    def setAngle( self, anglePotential ):
        for angle in self.angles:
            param = self.ffield.angleParameter( angle )
            anglePotential.set_coeff( angle, k=param['k'], t0=param['t0'] )
        return
    
    def setBond( self, bondPotential ):
        for bond in self.bonds:
            param = self.ffield.bondParameter( bond )
            bondPotential.bond_coeff.set( bond, k=param['k'], r0=param['r0'] )
        return
    
    def setDihedral( self, dihedralPotential ):
        for dihedral in self.dihedrals:
            param = self.ffield.dihedralParameter( dihedral )
            dihedralPotential.set_coeff( dihedral, k=param['k'], d=param['d'], n=param['n']  )
        return
    
    def setImproper( self, improperPotential ):
        for improper in self.impropers:
            param = self.ffield.improperParameter( improper )
            improperPotential.set_coeff( improper, k=param['k'], chi=param['chi']  )
        return
    
    def setPair( self, pairPotential ):
        for i, atype in enumerate( self.atomTypes ):
            for j, btype in enumerate( self.atomTypes ):
                if j >= i:
                    # dummy atoms treated specially
                    if atype.lower() == 'x' or btype.lower() == 'x':
                        pairPotential.pair_coeff.set( atype,
                                                      btype,
                                                      epsilon = 0.0,
                                                      sigma   = 0.0,
                                                      rcut    = 0.0  )
                    else:
                        param = self.ffield.pairParameter( atype, btype )
                        pairPotential.pair_coeff.set( atype, btype, epsilon = param['epsilon'], sigma=param['sigma'] )
        return
    
    def setAttributesFromFile(self, xmlFilename ):
        """Parse the xml file to extract the bonds, angles etc."""
        
        tree = ET.parse( xmlFilename )
        root = tree.getroot()
        
        bonds = []
        x = root.findall(".//bond")
        if len(x):
            btext = x[0].text
            for line in btext.split( os.linesep ):
                line = line.strip()
                if line:
                    bond = line.split()[0]
                    if bond not in bonds:
                        bonds.append( bond )
        self.bonds = bonds
        
        angles = []
        x = root.findall(".//angle")
        if len(x):
            atext = x[0].text
            for line in atext.split( os.linesep ):
                line = line.strip()
                if line:
                    angle = line.split()[0]
                    if angle not in angles:
                        angles.append( angle )
        self.angles = angles
        
        dihedrals = []
        dn = root.findall(".//dihedral")
        if len(dn):
            dtext = dn[0].text
            for line in dtext.split( os.linesep ):
                line = line.strip()
                if line:
                    dihedral = line.split()[0]
                    if dihedral not in dihedrals:
                        dihedrals.append( dihedral )
        self.dihedrals = dihedrals
        
        impropers = []
        dn = root.findall(".//improper")
        if len(dn):
            dtext = dn[0].text
            for line in dtext.split( os.linesep ):
                line = line.strip()
                if line:
                    improper = line.split()[0]
                    if improper not in impropers:
                        impropers.append( improper )
        self.impropers = impropers
        
        atomTypes = []
        atext = root.findall(".//type")[0].text
        for line in atext.split( os.linesep ):
            atomType = line.strip()
            if atomType and atomType not in atomTypes:
                atomTypes.append( atomType )
                
        self.atomTypes = atomTypes
        
        return
    
    def setupSystem(self,
                    xmlFilename,
                    doDihedral=False,
                    doImproper=False,
                    rCut=None,
                    quiet=False ):
        
        # Read parameters from file, check them and set the attributes
        self.checkParameters( xmlFilename=xmlFilename )
        
        if hoomdblue.init.is_initialized():
            hoomdblue.init.reset()
        
        # Init the sytem from the file
        system = hoomdblue.init.read_xml( filename=xmlFilename )
        
        # Below disables pretty much all output
        if quiet:
            print "Disabling HOOMD-Blue output!"
            hoomdblue.globals.msg.setNoticeLevel(0)

        # Set the parameters
        harmonic=None
        if len( self.bonds ):
            harmonic = hoomdblue.bond.harmonic()
            self.setBond( harmonic )
        
        aharmonic=None
        if len( self.angles ):
            aharmonic = hoomdblue.angle.harmonic()
            self.setAngle( aharmonic )
        
        dharmonic = improper = None
        if doDihedral and len( self.dihedrals ):
                dharmonic = hoomdblue.dihedral.harmonic()
                self.setDihedral( dharmonic )
        elif doImproper and len( self.dihedrals ):
            improper = hoomdblue.improper.harmonic()
            self.setImproper( improper )
        
        lj = hoomdblue.pair.lj( r_cut=rCut)
        self.setPair( lj )

        hoomdblue.globals.neighbor_list.reset_exclusions(exclusions = ['1-2', '1-3', '1-4', 'angle', 'body'] )

# Do we need to think about saving the references and deleting them?
#         del harmonic
#         del aharmonic
#         del improper
#         del lj

        return system
    
    def runMD(self,
              xmlFilename,
              doDihedral=False,
              doImproper=False,
              rCut=None,
              quiet=None,
              **kw ):

        if rCut is not None:
            self.rCut = rCut
        
        if doDihedral and doImproper:
            raise RuntimeError,"Cannot have impropers and dihedrals at the same time"
        
        self.system = self.setupSystem( xmlFilename,
                                        doDihedral=doDihedral,
                                        doImproper=doImproper,
                                        rCut=self.rCut,
                                        quiet=quiet )
        
        self._runMD( **kw )
        
        return True
    
    def runMDAndOptimise(self,
                         xmlFilename,
                         doDihedral=False,
                         doImproper=False,
                         rCut=None,
                         quiet=None,
                         **kw ):

        if rCut is not None:
            self.rCut = rCut
        
        if doDihedral and doImproper:
            raise RuntimeError,"Cannot have impropers and dihedrals at the same time"
        
        self.system = self.setupSystem( xmlFilename,
                                        doDihedral=doDihedral,
                                        doImproper=doImproper,
                                        rCut=self.rCut,
                                        quiet=quiet )
        
        # pre-optimise for preopt steps to make sure the sytem is sane - otherwise the MD
        # blows up
        preOptCycles = 1000
        self._optimiseGeometry(optCycles = preOptCycles,
                               maxOptIter=1 )
        # Now run the MD steps
        self._runMD(  **kw )
        
        # Finally do a full optimisation
        optimised = self._optimiseGeometry( **kw )
        
        return optimised
    
    def optimiseGeometry( self,
                          xmlFilename,
                          doDihedral=False,
                          doImproper=False,
                          rCut=None,
                          quiet=None,
                          **kw ):

        if rCut is not None:
            self.rCut = rCut
        
        if doDihedral and doImproper:
            raise RuntimeError,"Cannot have impropers and dihedrals at the same time"
        
        self.system = self.setupSystem( xmlFilename,
                                        doDihedral=doDihedral,
                                        doImproper=doImproper,
                                        rCut=self.rCut,
                                        quiet=quiet )
        
        optimised = self._optimiseGeometry( **kw )
        
        return optimised
    
    def _runMD(self, mdCycles=1000, T=0.1, tau=0.5, dt=0.005, **kw ):
        
        # Added **kw arguments so that we don't get confused by arguments intended for the optimise
        # when MD and optimiser run together
        
        integrator_mode = hoomdblue.integrate.mode_standard( dt=dt )
        nvt_rigid = hoomdblue.integrate.nvt_rigid(group=hoomdblue.group.rigid(), T=T, tau=tau )

        record=False
        if record:
            xmld = hoomdblue.dump.xml(filename="runmd.xml",
                                      vis=True )
            dcdd = hoomdblue.dump.dcd(filename="runmd.dcd",
                                      period=10,
                                      unwrap_full=True,
                                      overwrite=True )

        # run mdCycles time steps
        hoomdblue.run( mdCycles )
        
        nvt_rigid.disable()
        
        del nvt_rigid
        del integrator_mode
        
        return
    
    def _optimiseGeometry(self,
                          carOut="hoomdOpt.car",
                          optCycles = 100000,
                          maxOptIter=100,
                          dt=0.005,
                          **kw ):
        """Optimise the geometry with hoomdblue"""
        
        # Create the integrator with the values specified
        fire = hoomdblue.integrate.mode_minimize_rigid_fire( group=hoomdblue.group.all(),
                                                             dt=dt,
                                                             Nmin=5,
                                                             alpha_start=0.1,
                                                             ftol=1e-2,
                                                             Etol=1e-4,
                                                             finc=1.1,
                                                             fdec=0.5
                                                            )
        
        record=False
        if record:
            # For tracking the optimsation        
            xmld = hoomdblue.dump.xml(filename="runopt.xml", vis=True)
            dcdd = hoomdblue.dump.dcd(filename="runopt.dcd",
                                      period=10, 
                                      unwrap_full=True,
                                      overwrite=True,
                                       )
            
        count = 0
        optimised=True
        while not( fire.has_converged() ):
            hoomdblue.run( optCycles )
            count += 1
            if count >= maxOptIter:
                print "********** DID NOT OPTIMISE!!!! ***************"
                optimised=False
                break
    
        # Delete variables before we return to stop memory leaks
        del fire
#         del harmonic
#         del aharmonic
#         del improper
#         del lj
        
        # Write out a car file so we can see what happened
        #self.writeCar( system=self.system, filename=carOut, unwrap=True )
        
        return optimised
    
def xml2xyz( xmlFilename, xyzFilename ):
    """Convert a hoomdblue xml file to xyz"""


    tree = ET.parse( xmlFilename )
    root = tree.getroot()
    
    atomTypes = []
    atext = root.findall(".//type")[0].text
    atomTypes = [ line.strip() for line in atext.split( os.linesep ) if line.strip() ]
    
    ptext = root.findall(".//position")[0].text
    positions = [ line.strip().split() for line in ptext.split( os.linesep ) if line.strip() ]
        
    assert len(atomTypes) == len(positions)
    
    # Convert atom types to symbols
    symbols = []
    for at in atomTypes:
        if at.lower()[0] == 'x':
            symbols.append( 'x' )
        else:
            symbols.append( util.label2symbol( at ) )
    
    # Now write out xyz
    with open( xyzFilename, 'w') as o:
        o.write("{0}\n".format( len( positions ) ) )
        o.write("XYZ file created from: {0}\n".format( xmlFilename ) )
        for i, symbol in enumerate( symbols ):
            o.write( "{0:5}   {1:0< 15}   {2:0< 15}   {3:0< 15}\n".format( symbol,
                                                                           float( positions[i][0] ),
                                                                           float( positions[i][1] ),
                                                                           float( positions[i][2] )
                                                                           ) )
        
        o.write("\n")
    
    print "Wrote file: {0}".format( xyzFilename )
    
    return
        

if __name__ == "__main__":
    
    optimiser = HoomdOptimiser()
    #optimiser.setAttributesFromFile( '/Users/jmht/Documents/abbie/AMBI/ambuild/builder/hoomdopt.xml' )
    xmlFilename = sys.argv[1]
    xyzFilename = xmlFilename +".xyz"
    xml2xyz( xmlFilename, xyzFilename )
    optimiser.optimiseGeometry( xmlFilename=xmlFilename, doDihedral=True )

