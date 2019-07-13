#!/usr/bin/env python
import os, sys
ambuild_home = "/opt/ambuild.git"
sys.path.insert(0, ambuild_home)

# This imports the builder cell module - this is the only module that should be required
from ambuild import ab_util


mycell = ab_util.cellFromPickle("step_628.pkl")
