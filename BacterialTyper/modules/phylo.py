#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez                                          ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain     ##
##############################################################
"""
Generates a phylogenetic reconstruction
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
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.modules import help_info
from BacterialTyper.scripts import sampleParser

####################################
def run_phylo(options):
    """
    Main function acting as an entry point to the module *phylo*.
    """

    ##################################
    ### show help messages if desired    
    ##################################
    if (options.help_format):
        ## help_format option
        sampleParser.help_format()
        exit()

    elif (options.help_project):
        ## information for project
        help_info.project_help()
        exit()
    
    ## init time
    start_time_total = time.time()

    ## debugging messages
    global Debug
    if (options.debug):
        Debug = True
    else:
        Debug = False
        
    ### set as default paired_end mode
    if (options.single_end):
        options.pair = False
    else:
        options.pair = True

    functions.pipeline_header()
    functions.boxymcboxface("Phylogenetic reconstruction")

    print ("--------- Starting Process ---------")
    functions.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## set mode: project/detached
    global Project
    if (options.project):
        outdir = input_dir        
        Project=True
    elif (options.detached):
        Project=False
        outdir = os.path.abspath(options.output_folder)

    ## get files
    print ("+ Retrieve all files available...")
    pd_samples_retrieved = sampleParser.get_files(options, input_dir, "trim", ['_trim_'])
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        functions.create_folder(outdir)
    
    ## for samples
    outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "phylo")

    
    
    
    
    
    
    
    
    
    
    
    
    