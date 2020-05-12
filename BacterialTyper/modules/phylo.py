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
       
    ####################
    ## Genbank_ID
    ####################
    reference_gbk_file = ""
    if options.Genbank_ID:
        db_frame_ncbi = database_generator.getdbs('NCBI', options.database, 'genbank', Debug)
    
        ## debug message
        if (Debug):
            print (colored("**DEBUG: db_frame_ncbi **", 'yellow'))
            print (db_frame_ncbi) 

        NCBI_folder = os.path.join(options.database, 'NCBI')
        dir_path = os.path.join(NCBI_folder, 'genbank', 'bacteria', options.Genbank_ID)    
        if (options.Genbank_ID in db_frame_ncbi.index): 
            print('+ Reference (%s) available in database provided' %options.Genbank_ID)
        else:
            print ('+ Reference (%s) is not available in database provided' %options.Genbank_ID)
            print ('+ Try to download it.')
            NCBI_folder = os.path.join(options.database, 'NCBI')
            database_generator.ngd_download(dir_path, options.Genbank_ID, NCBI_folder)
    
        ## get files download
        (genome, prot, gff, gbk) = database_generator.get_files_download(dir_path)
        if Debug:
                print ('genome:' + genome)
                print ('prot:' + prot)
                print ('gff:' + gff)
                print ('gbk:' + gbk)
	            
        if functions.is_non_zero_file(gbk):
            print('+ Genbank file format reference available.')
            reference_gbk_file = gbk
        else:
            print(colored('\n+ No genbank file available for the reference specified. Some error occurred while downloading', 'red'))
            exit()
            
    ####################          
    ## user_sample_ID
    ####################
    elif options.user_sample_ID:
        db_frame_user_Data = database_user.get_userData_files(options, os.path.join(options.database, 'user_data'))
        df_data = db_frame_user_Data.groupby('name')

        try:
            this_sample_df = df_data.get_group(options.project_sample_ID)
            print('+ Reference (%s) available in database folder provided' %options.user_sample_ID)
        except:
            print (colored('** WARNING: Reference (%s) not available in database folder provided' %options.user_sample_ID, 'yellow'))
            print ('+ Lets try to update the database first.')
            db_frame_user_dataUpdated = database_user.update_database_user_data(options.database, input_dir, Debug, options)
            df_data = db_frame_user_dataUpdated.groupby('name')
 
            try:
                this_sample_df = df_data.get_group(options.user_sample_ID)
                print('+ Reference (%s) available in database updated' %options.user_sample_ID)
                db_frame_user_Data = db_frame_user_dataUpated

            except:
                print(colored('\n** ERROR: No reference (%s) available in database updated. Some error occurred...' %options.user_sample_ID, 'red'))
                exit()

         ## debug message
        if (Debug):
            print (colored("**DEBUG: db_frame_user_Data **", 'yellow'))
            print (db_frame_user_Data)
            print (colored("**DEBUG: this_sample_df (groupby name)**", 'yellow'))
            print (this_sample_df)


       ## get gbk file
        gbk = this_sample_df.loc[ this_sample_df['ext']=='gbf','sample'].values[0]
        
        ## debug
        if Debug:
            print ("** DEBUG: this_sample_df")
            print (this_sample_df)
            print ('gbk:' + gbk)
  
        ## check if exists
        if functions.is_non_zero_file(gbk):
            print('+ Genbank file format reference available.')
            reference_gbk_file = gbk
        else:
            print(colored('\n+ No genbank file available for the reference specified. Some error occurred while downloading', 'red'))
            exit()
        
    ####################    
    ## project_sample_ID
    ####################
    elif options.project_sample_ID:
        
        db_frame_project_Data = database_user.get_userData_files(options, options.input)
        df_data = db_frame_project_Data.groupby('name')

        try:
            this_sample_df = df_data.get_group(options.project_sample_ID)
            print('+ Reference (%s) available in project folder provided' %options.project_sample_ID)
        except:
            print (colored('** ERROR: Reference (%s) not available in project folder provided' %options.project_sample_ID, 'red'))
            print ('+ Check the spelling or provide a valid ID.')
            exit()
 
        ## debug message
        if (Debug):
            print (colored("**DEBUG: db_frame_project_Data **", 'yellow'))
            print (db_frame_project_Data)
            print (colored("**DEBUG: this_sample_df (groupby name)**", 'yellow'))
            print (this_sample_df)

        ## get gbk file
        gbk = this_sample_df.loc[ this_sample_df['ext']=='gbf','sample'].values[0]

        ## debug
        if Debug:
            print ("** DEBUG: this_sample_df")
            print (this_sample_df)
            print ('gbk:' + gbk)

        ## check if exists
        if functions.is_non_zero_file(gbk):
            print('+ Genbank file format reference available.')
            reference_gbk_file = gbk
        else:
            print(colored('\n+ No genbank file available for the reference specified. Some error occurred while downloading', 'red'))
            exit()

    ####################
    ## user_gbk
    ####################
    elif options.user_gbk:
        options.user_gbk = os.path.abspath(options.user_gbk)
        if functions.is_non_zero_file(options.user_gbk):
            print('+ Reference provided via --user_gbk is available and ready to use.')
        else:
            print('+ Reference provided via --user_gbk not available or accessible.')
            print(colored('\n+ Check the path or integrity of the file. Some error occurred...', 'red'))
            exit()        
    
    else:
        print()




    exit()


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

    
    
    
    
    
    
    
    
    
    
    
    
    
