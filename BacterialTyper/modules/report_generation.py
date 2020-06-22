#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Creates a report in HTML for each step of the process analyzed.
"""
## import useful modules
import os
import sys
import re
import time
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.scripts import functions, database_user
from BacterialTyper.config import set_config
from BacterialTyper.report.Staphylococcus import get_spa_typing

## example report check: https://github.com/tseemann/nullarbor

## for staphyloccocus aureus: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0593-7
## for klebsiella: kleborate
## for micobacterium tuberculosis: MTBSeq

####################################
def run_report(options):
    
    ## init time
    start_time_total = time.time()
    
    ##################################
    ### show help messages if desired    
    ##################################

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

    ## message header
    functions.pipeline_header()
    functions.boxymcboxface("Report generation module")
    print ("--------- Starting Process ---------")
    functions.print_time()
    
    ## call assemble using spades
    start_time_partial = start_time_total
    
    ## absolute path for in & out
    options.database = os.path.abspath(options.database)
    global input_dir
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## set mode: project/detached
    global Project
    if (options.detached):
        options.project = False
        outdir = os.path.abspath(options.output_folder)
        Project=False
    else:
        options.project = True
        outdir = input_dir    
        Project=True

    ##
    print ("\n+ Get project information:")
    
    ## get files: trimm, assembly, annotation
    pd_samples_retrieved = database_user.get_userData_files(options, input_dir)
    
    ## get info: profile, ident, cluster, MGE
    pd_samples_info = database_user.get_userData_info(options, input_dir)
    
    ## get databases to list
    #retrieve_databases = get_options_db(options)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)
        
        print (colored("**DEBUG: pd_samples_info **", 'yellow'))
        print (pd_samples_info)
        
    
    ## generate output folder, if necessary
    print ("\n\n\n+ Generate a report summarizing analysis and sample information")
    if not options.project:
        functions.create_folder(outdir)
        outdir_report = outdir
    else:
        ### report generation
        outdir_report = functions.create_subfolder("report", outdir)
    
    ## create report with all data
    summary_report = functions.create_subfolder("summary_report", outdir_report)
    print ("Folder: ", summary_report)
    
    ## time stamp
    start_time_partial = functions.timestamp(start_time_partial)
    
    ## create species specific report if any
    species_specific_df = pd.DataFrame()
    if options.species_report == "Saures":
        Saureus_specific(pd_samples_retrieved, pd_samples_info, species_specific_df, options)
        
    
    ## time stamp
    start_time_partial = functions.timestamp(start_time_partial)


    print ("\n*************** Finish *******************")
    start_time_partial = functions.timestamp(start_time_total)

    print ("+ Exiting Report generation module.")
    exit()
    
    ## time stamp
    start_time_partial = functions.timestamp(start_time_total)        


def Saureus_specific(samples_df, samples_info, results_df, options):
    """
    Retrieves Saureus specific information.
    
    See additional information in :doc:`../../user_guide/report/Saureus/saureus_report`
    """
    
    ## get European Quality Control genes
    
    
    ## get spatyping
    samples_df.loc[:, "tag"] = "assembly"
    print (samples_df)
    #dictionary_fasta_files = 
    #get_spa_typing.module_call(options.database, dictionary_fasta_files, options.debug)
    
    ## get sccmec
    
    
    