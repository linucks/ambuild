import sys

# can't do below as I think this needs to be set before the python interpreter is started.
#os.environ['DYLD_LIBRARY_PATH'] = "/Applications/HOOMD-blue.app/Contents/MacOS"

sys.path.append("/Applications/HOOMD-blue.app/Contents/lib/hoomd/python-module")
sys.path.append("/Applications/HOOMD-blue.app/Contents/MacOS")

from hoomd_script import *


