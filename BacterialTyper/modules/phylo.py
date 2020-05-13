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
from BacterialTyper.scripts import variant_calling


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
    
    ## get the reference
    reference_gbk_file = get_reference_gbk(options)
                 
    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        functions.create_folder(outdir)
    
    ##################################
    ## select samples and map    
    ####################################
    dict_folders = map_samples(options, reference_gbk_file, input_dir, outdir)
    ## time stamp
    start_time_partial = functions.timestamp(start_time_total)

    exit()

    ##################################
    ## Create core alingment
    ##################################
    list_folders = list(dict_folders.values())
    options_string = ""
    variant_calling.snippy_core_call(list_folders, options_string, options.name)
    
    
def map_samples(options, reference_gbk_file, input_dir, outdir):    
    
    pd_samples_retrieved_merge = pd.DataFrame()

    ## all_data // only_project_data
    if (options.all_data or options.only_project_data):
        ## get files to map
        print ("+ Retrieve samples to map available...")
        pd_samples_retrieved = sampleParser.get_files(options, input_dir, "trim", ['_trim_'])
        
        ## discard the sample used as reference if any
        if options.project_sample_ID:
            pd_samples_retrieved = pd_samples_retrieved.drop(index=options.project_sample_ID)
    
        ## create output directories
        outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "phylo")
    
    ####################################
    ## user_data // genbank_data // only_external_data
    if (options.all_data or options.user_data): 
        print ("+ Retrieve samples to map from user database...")
        db_frame_user_Data = database_user.get_userData_files(options, os.path.join(options.database, 'user_data'))

        ## discard the sample used as reference if any
        if options.user_sample_ID:
            db_frame_user_Data = pd_samples_retrieved.drop(index=options.user_sample_ID)
            
        ## create output directories
        outdir_dict2 = functions.outdir_project(os.path.join(options.database, 'user_data'), options.project, db_frame_user_Data, "phylo")
        
        ## merge if both contain data
        if not pd_samples_retrieved.empty:
            pd_samples_retrieved_merge = pd.concat([db_frame_user_Data, pd_samples_retrieved], sort=True, ignore_index=True).drop_duplicates()
            outdir_dict.update(outdir_dict2)
        else:
            outdir_dict = outdir_dict2
            pd_samples_retrieved_merge = db_frame_user_Data   
    
    ## check data
    try:
        if db_frame_user_Data.empty:
            print ()   
    except:
        pd_samples_retrieved_merge = pd_samples_retrieved   
   
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieved_merge **", 'yellow'))
        print (pd_samples_retrieved_merge)
        print (colored("**DEBUG: outdir_dict **", 'yellow'))
        print (outdir_dict)
        
    ####################################
    ## for fastq samples
    ####################################

    # optimize threads
    name_list = set(pd_samples_retrieved_merge["name"].tolist())
    threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    ## call snippy
    print ("\n+ Create mapping of fastq reads for project samples:")
    contig_option = ""
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(snippy_variant_caller,reference_gbk_file, sorted(cluster["sample"].tolist()), threads_job, outdir_dict[name], name, contig_option, options.other_options, Debug): name for name, cluster in sample_frame }
        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    ## TODO: use contig option
    if (options.all_data or options.genbank_data):
        print ()


    return (outdir_dict)

####################
def get_reference_gbk(options):

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
        reference_gbk_file = options.user_gbk

    
    return (reference_gbk_file)

#############################################
def snippy_variant_caller(reference, files, threads, outdir, name, contig_option, other_options, Debug):
    
    ## create subfolder within phylo for this mapping
    subdir = functions.create_subfolder(name, outdir)
       
    ## check if previously process and succeeded
    filename_stamp = subdir + '/.success'

    if os.path.isfile(filename_stamp):
        stamp =    functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
        return()
    else:
         # Call variant calling
        return(variant_calling.snippy_call(reference, files, threads, subdir, name, contig_option, other_options, Debug))
  
