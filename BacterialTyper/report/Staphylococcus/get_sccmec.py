#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Creates sccmec typing
"""
## useful imports
import os
from sys import argv
from io import open
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.config import set_config

## import my HCGB module 
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main

##############
def help_options():
    print ("\nUSAGE: python %s path sample...\n"  %os.path.realpath(__file__))

