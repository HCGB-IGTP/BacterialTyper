#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Creates sccmec typing
"""
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd
from Bio import SeqIO

## import my modules
from BacterialTyper.config import set_config

## import my HCGB module 
from HCGB.functions import HCGB_files
from HCGB.functions import HCGB_time
from HCGB.functions import HCGB_main

##############
def help_options():
    print ("\nUSAGE: python %s path sample...\n"  %os.path.realpath(__file__))

