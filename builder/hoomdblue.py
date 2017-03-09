import sys
# 
# import platform
# from socket import gethostname
# can't do below as I think this needs to be set before the python interpreter is started.
#os.environ['DYLD_LIBRARY_PATH'] = "/Applications/HOOMD-blue.app/Contents/MacOS"
# if platform.system() == "Darwin":
#     sys.path.append("/Applications/HOOMD-blue.app/Contents/lib/hoomd/python-module")
#     sys.path.append("/Applications/HOOMD-blue.app/Contents/MacOS")
#     
# elif gethostname() == "cytosine":
#     sys.path.append("/home/jmht/Downloads/hoomd-0.11.3/install/lib/hoomd/python-module")
#     sys.path.append("/home/jmht/Downloads/hoomd-0.11.3/install/bin") 
#     
# else:
#     sys.path.append("/opt/hoomd-0.11.3/hoomdblue-install/lib/hoomd/python-module")
#     sys.path.append("/opt/hoomd-0.11.3/hoomdblue-install/bin") 

from hoomd_script import *

# Later versions of hoomd-blue have a context that needs to be updated
if 'context' in locals().keys():
    # Context parses the command-line, so we need to make sure there
    # are no arguments or else we lose the ability to display our own help
    _argv = sys.argv
    sys.argv = [sys.argv[0]]
    context.initialize()
    sys.argv = _argv


