#usr/bin/env python
'''
This code generates a bacteriophage identification profile for each sample
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
pythonDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(pythonDir)
import functions

configDir = os.path.dirname(os.path.realpath(__file__)) + '/../../config/'
sys.path.append(configDir)
import config

Phispy_folder = config.EXECUTABLES['Phispy_folder']



