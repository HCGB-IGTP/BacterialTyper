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
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import database_user


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

    ### parse the reference
    print ("+ Retrieve the reference...")
    
    ## get the database 
    options.database = os.path.abspath(options.database)
    db_frame_user_Data = database_generator.getdbs('user_data', options.database, 'user_data', Debug)
    db_frame_ncbi = database_generator.getdbs('user_data', options.database, 'user_data', Debug)
    
    ## return both dataframes
    retrieve_databases = pd.concat([db_frame_user_Data, db_frame_ncbi], sort=True, ignore_index=True)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: retrieve_database **", 'yellow'))
        pd.set_option('display.max_colwidth', None)
        pd.set_option('display.max_columns', None)
        print (retrieve_databases)    
    
    ## Genbank_ID
    if options.Genbank_ID:
        print()
    
    ## user_sample_ID
    elif options.user_sample_ID:
        print()
    
    ## project_sample_ID
    elif options.project_sample_ID:
        print()
    
    ## user_gbk
    elif options.user_gbk:
        print()
    
    else:
        print()







    ## get files to map
    print ("+ Retrieve samples to map available...")
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

    
    
    
    
    
    
    
    
    
    
    
    
    