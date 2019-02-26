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
thisDir = os.path.dirname(os.path.abspath(argv[0]))
sys.path.append(thisDir)

import config
import functions
