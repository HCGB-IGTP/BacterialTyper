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
pythonDir = os.path.dirname(os.path.abspath(argv[0]))
sys.path.append(pythonDir)
import functions

configDir = os.path.dirname(os.path.abspath(argv[0])) + '../../config/'
sys.path.append(configDir)
import config

RscriptDir = os.path.dirname(os.path.abspath(argv[0])) + '../R/'
MLSTarR_script = RscriptDir + '/MLSTar_call.R'

print (pythonDir)
print (configDir)
print (MLSTarR_script)


## check if profile and sequences are already downloaded
## generate profile folder

## call MLSTar for this sample
## generate sequence folder

## doMLST
## clean makeblastdb files
## clean folders generated
## then
## Rscript MLSTar_call.R --species "saureus" --scheme 1 --dir_seq /home/labs/lslab/jsanchez/DATA/Saures_test/MLSTar/test/seq/ --dir_profile /home/labs/lslab/jsanchez/DATA/Saures_test/MLSTar/test/profile/ --dir ./test2 --file denovo_assembly_scaffolds.fas
