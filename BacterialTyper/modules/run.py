#!/usr/bin/env python3
###########################################################
## Jose F. Sanchez                                         ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain     ##
###########################################################
"""
This module generates a uniq call of BacterialTyper to 
any other module. It automates the analysis and process.
"""
## useful imports
import time
import io
import os
import re
import sys
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.modules import prep
from BacterialTyper.modules import qc
from BacterialTyper.modules import ident
from BacterialTyper.modules import profile
from BacterialTyper.modules import database
from BacterialTyper.modules import MGE
from BacterialTyper.modules import cluster
from BacterialTyper.modules import assemble
from BacterialTyper.modules import annot
from BacterialTyper.modules import database
from BacterialTyper.modules import trim
from BacterialTyper.modules import report_generation
from BacterialTyper.modules import help_info
from BacterialTyper import __version__ as pipeline_version

from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files

####################################
def run_BacterialTyper(options):

    ## init time
    start_time_total = time.time()

    ## debugging messages
    global Debug
    if (options.debug):
        Debug = True
    else:
        Debug = False

    ##################################
    ### show help messages if desired    
    ##################################
    
    ## if any help_flag provided will print and exit
    help_info.help_info(options)
    
    ### 
    HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
    HCGB_aes.boxymcboxface("BacterialTyper analysis")

    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## set mode: project/detached
    if (options.project):
        outdir = input_dir        
    elif (options.detached):
        outdir = os.path.abspath(options.output_folder)
