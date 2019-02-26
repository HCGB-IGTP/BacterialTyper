#usr/bin/env python
'''
This code generates a virulence and an antibiotic resistance profile.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open

## import my modules
pythonDir = os.path.dirname(os.path.abspath(argv[0]))
sys.path.append(pythonDir)
import functions

configDir = os.path.dirname(os.path.abspath(argv[0])) + '../../config/'
sys.path.append(configDir)
import config
