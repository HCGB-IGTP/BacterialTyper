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
thisDir = os.path.dirname(os.path.abspath(argv[0]))
sys.path.append(thisDir)

import config
import functions

MLSTarR_script = thisDir + '/MLSTar_call.R'


## check if profile and sequences are already downloaded
## generate profile folder

## call MLSTar for this sample
## generate sequence folder

## doMLST
## clean makeblastdb files
## clean folders generated
## then
## Rscript MLSTar_call.R --species "saureus" --scheme 1 --dir_seq /home/labs/lslab/jsanchez/DATA/Saures_test/MLSTar/test/seq/ --dir_profile /home/labs/lslab/jsanchez/DATA/Saures_test/MLSTar/test/profile/ --dir ./test2 --file denovo_assembly_scaffolds.fas
