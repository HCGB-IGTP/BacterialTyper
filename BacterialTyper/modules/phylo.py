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
from BacterialTyper.modules import help_info
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import database_user
from BacterialTyper.scripts import variant_calling
from BacterialTyper.scripts import phylo_parser
from BacterialTyper.config import set_config
from BacterialTyper import __version__ as pipeline_version

import HCGB
from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
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

    HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
    HCGB_aes.boxymcboxface("Phylogenetic reconstruction")

    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## set mode: project/detached
    ## Project mode as default
    project_mode=True
    if (options.detached):
        options.project = False
        project_mode=False
        outdir = os.path.abspath(options.output_folder)
    else:
        options.project = True
        outdir = input_dir    
    
    ## get the database 
    options.database = os.path.abspath(options.database)
    
    ### parse the reference
    print ("+ Retrieve the reference...")
    reference_gbk_file = get_reference_gbk(options)
                 
    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        HCGB_files.create_folder(outdir)
    
    ##################################
    ## select samples and map    
    ####################################
    print ("+ Retrieve samples to map available...")
    dict_folders = map_samples(options, reference_gbk_file, input_dir, outdir)
    
    if Debug:
        print (colored("**DEBUG: dict_folders **", 'yellow'))
        print (dict_folders)
    
    ## time stamp
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ##################################
    ## Create core alingment
    ##################################
    outdir_report = HCGB_files.create_subfolder("report", outdir)
    phylo_dir = HCGB_files.create_subfolder("phylo", outdir_report)
    analysis_dir = HCGB_files.create_subfolder(options.name, phylo_dir)
    snippy_dir = HCGB_files.create_subfolder("snippy", analysis_dir)
        
    list_folders = list(dict_folders.values())
    options_string = ""
    variant_calling.snippy_core_call(list_folders, options_string, options.name, 
                                     snippy_dir, options.output_format, Debug)

    ## time stamp
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ## snp distance matrix
    snp_distance_dir = HCGB_files.create_subfolder("snp_distance", analysis_dir)
    name_matrix = os.path.join(snp_distance_dir, "snp_matrix_" + options.name)
    
    countGaps = False
    aln_file = os.path.join(snippy_dir, options.name + '.aln')
    phylo_parser.get_snp_distance(aln_file, options.output_format, countGaps, name_matrix, Debug)
    
    ## time stamp
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ## phylogenetic analysis
    iqtree_output = HCGB_files.create_subfolder("iqtree", analysis_dir)
    phylo_parser.ml_tree(snippy_dir, options.name, options.threads, iqtree_output, Debug)
    
    ## time stamp
    start_time_partial = HCGB_files.timestamp(start_time_total)

    print ("\n*************** Finish *******************")
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ## dump information and parameters
    info_dir = HCGB_files.create_subfolder("info", outdir)
    print("+ Dumping information and parameters")
    runInfo = { "module":"phylo", "time":HCGB_time.timestamp(time.time()),
                "BacterialTyper version":pipeline_version }
    HCGB_info.dump_info_run(info_dir, 'phylo', options, runInfo, options.debug)

    print ("+ Exiting phylo module.")
    return()

#####################
def map_samples(options, reference_gbk_file, input_dir, outdir):    
    """
    """

    ## set it as variable
    contig_option = False
    
    pd_samples_retrieved_merge = pd.DataFrame()
    pd_samples_retrieved = pd.DataFrame()

    ## all_data // only_project_data
    if (options.all_data or options.only_project_data):
        ## get files to map
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
        
        ## discard the sample used as reference if any
        if options.project_sample_ID:
            pd_samples_retrieved = pd_samples_retrieved.drop(index=options.project_sample_ID)
    
        ## create output directories
        outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "phylo", options.debug)
    
    ####################################
    ## user_data // genbank_data // only_external_data
    if (options.all_data or options.user_data): 
        
        ## --------------------------- ##
        ### get user database
        ## --------------------------- ##
        print ("+ Retrieve samples to map from user database...")
        db_frame_user_Data = database_user.get_userData_files(options, os.path.join(options.database, 'user_data'))

        ## discard the sample used as reference if any
        if options.user_sample_ID:
            #db_frame_user_Data = pd_samples_retrieved.drop(index=options.user_sample_ID)
            db_frame_user_Data = db_frame_user_Data.drop(index=options.user_sample_ID) ## Why not this?
            
        ## create output directories in database entries in user_data
        outdir_dict2 = HCGB_files.outdir_subproject(os.path.join(options.database, 'user_data'), db_frame_user_Data, "phylo")
        
        ## If user desires to map contigs, map trimmed as default
        if (contig_option):
            # filter for assemblies retrieved
            db_frame_user_Data = db_frame_user_Data.loc[db_frame_user_Data['tag'] == "assembly",]
        else:
             # filter for raw reads
             db_frame_user_Data = db_frame_user_Data.loc[db_frame_user_Data['tag'] == "reads",]
        
        ## merge if both contain data
        if not pd_samples_retrieved.empty:
            pd_samples_retrieved_merge = pd.concat([db_frame_user_Data, pd_samples_retrieved], sort=True, ignore_index=True).drop_duplicates()
            outdir_dict.update(outdir_dict2)
        else:
            outdir_dict = outdir_dict2
            pd_samples_retrieved_merge = db_frame_user_Data   
    
        ## --------------------------- ##
        ### get genbank database
        ## --------------------------- ##
        ## set contig option for these data only
    
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
        
    ## TODO: use contig option
    if (contig_option or options.genbank_data):
        print ()
        
    ####################################
    ## for fastq samples
    ####################################

    # optimize threads
    name_list = set(pd_samples_retrieved_merge["name"].tolist())
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (options.debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    ## call snippy
    print ("\n+ Create mapping of fastq reads for project samples:")
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved_merge.groupby(["name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(snippy_variant_caller, reference_gbk_file, sorted(cluster["sample"].tolist()), threads_job, outdir_dict[name], options.name, contig_option, options.other_options, name, options.debug): name for name, cluster in sample_frame }
        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    ## subfolder within phylo for this mapping
    new_outdir_dict = {}
    for key,value in outdir_dict.items():
        tag = os.path.join(value, key + '_vs_' + options.name)
        new_outdir_dict[key] = tag 

    return (new_outdir_dict)

####################
def get_reference_gbk(options):

   ####################
    ## Genbank_ID
    ####################
    reference_gbk_file = ""
    if options.Genbank_ID:
        db_frame_ncbi = database_generator.getdbs('NCBI', options.database, 'genbank', options.debug)
    
        ## debug message
        if (options.debug):
            print (colored("**DEBUG: db_frame_ncbi **", 'yellow'))
            print (db_frame_ncbi) 

        NCBI_folder = HCGB_files.create_subfolder('NCBI', options.database)
        dir_path = os.path.join(NCBI_folder, 'genbank', 'bacteria', options.Genbank_ID)    
        if (options.Genbank_ID in db_frame_ncbi.index): 
            print('\t+ Reference (%s) available in database provided' %options.Genbank_ID)
        else:
            print ('\t+ Reference (%s) is not available in database provided' %options.Genbank_ID)
            print ('\t+ Try to download it.')
            database_generator.ngd_download(dir_path, options.Genbank_ID, NCBI_folder)
    
        ## get files download
        (genome, prot, gff, gbk) = database_generator.get_files_download(dir_path)
        if options.debug:
                print (colored("**DEBUG: genome:" + genome, 'yellow'))
                print (colored("**DEBUG: prot:" + prot, 'yellow'))
                print (colored("**DEBUG: gff:" + gff, 'yellow'))
                print (colored("**DEBUG: gbk:" + gbk, 'yellow'))
                
        if HCGB_files.is_non_zero_file(gbk):
            print('\t+ Genbank file format reference available.')
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
            print('\t+ Reference (%s) available in database folder provided' %options.user_sample_ID)
        except:
            print (colored('** WARNING: Reference (%s) not available in database folder provided' %options.user_sample_ID, 'yellow'))
            print ('\t+ Lets try to update the database first.')
            db_frame_user_dataUpdated = database_user.update_database_user_data(options.database, input_dir, options.debug, options)
            df_data = db_frame_user_dataUpdated.groupby('name')
 
            try:
                this_sample_df = df_data.get_group(options.user_sample_ID)
                print('\t+ Reference (%s) available in database updated' %options.user_sample_ID)
                db_frame_user_Data = db_frame_user_dataUpated

            except:
                print(colored('\n** ERROR: No reference (%s) available in database updated. Some error occurred...' %options.user_sample_ID, 'red'))
                exit()

         ## debug message
        if (options.debug):
            print (colored("**DEBUG: db_frame_user_Data **", 'yellow'))
            print (db_frame_user_Data)
            print (colored("**DEBUG: this_sample_df (groupby name)**", 'yellow'))
            print (this_sample_df)


       ## get gbk file
        gbk = this_sample_df.loc[ this_sample_df['ext']=='gbf','sample'].values[0]
        
        ## debug
        if options.debug:
            print ("** DEBUG: this_sample_df")
            print (this_sample_df)
            print ('gbk:' + gbk)
  
        ## check if exists
        if HCGB_files.is_non_zero_file(gbk):
            print('\t+ Genbank file format reference available.')
            reference_gbk_file = gbk
        else:
            print(colored('\n** ERROR: No genbank file available for the reference specified. Some error occurred while downloading', 'red'))
            exit()
        
    ####################    
    ## project_sample_ID
    ####################
    elif options.project_sample_ID:
        
        db_frame_project_Data = database_user.get_userData_files(options, options.input)
        df_data = db_frame_project_Data.groupby('name')

        try:
            this_sample_df = df_data.get_group(options.project_sample_ID)
            print('\t+ Reference (%s) available in project folder provided' %options.project_sample_ID)
        except:
            print (colored('** ERROR: Reference (%s) not available in project folder provided' %options.project_sample_ID, 'red'))
            print ('\t+ Check the spelling or provide a valid ID.')
            exit()
 
        ## debug message
        if (options.debug):
            print (colored("**DEBUG: db_frame_project_Data **", 'yellow'))
            print (db_frame_project_Data)
            print (colored("**DEBUG: this_sample_df (groupby name)**", 'yellow'))
            print (this_sample_df)

        ## get gbk file
        gbk = this_sample_df.loc[ this_sample_df['ext']=='gbf','sample'].values[0]

        ## debug
        if options.debug:
            print ("** DEBUG: this_sample_df")
            print (this_sample_df)
            print ('gbk:' + gbk)

        ## check if exists
        if HCGB_files.is_non_zero_file(gbk):
            print('\t+ Genbank file format reference available.')
            reference_gbk_file = gbk
        else:
            print(colored('\n** ERROR: No genbank file available for the reference specified. Some error occurred while downloading', 'red'))
            exit()

    ####################
    ## user_ref
    ####################
    elif options.user_ref:
        options.user_ref = os.path.abspath(options.user_ref)
        if HCGB_files.is_non_zero_file(options.user_ref):
            print('\t+ Reference provided via --user_ref is available and ready to use.')
        else:
            print('\n** ERROR: Reference provided via --user_ref not available or accessible.')
            print(colored('\n+ Check the path or integrity of the file. Some error occurred...', 'red'))
            exit()
        reference_gbk_file = options.user_ref

    
    return (reference_gbk_file)

#############################################
def snippy_variant_caller(reference, files, threads, outdir, name, contig_option, other_options, sample_name, Debug):
    
    ## create subfolder within phylo for this mapping
    tag = sample_name + '_vs_' + name
    subdir = HCGB_files.create_subfolder(tag, outdir)
       
    ## check if previously process and succeeded
    filename_stamp = subdir + '/.success'
    
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, tag), 'yellow'))
    else:
         # Call variant calling
        code = variant_calling.snippy_call(reference, files, threads, subdir, 
                                           sample_name, contig_option, other_options, Debug)
        if code == 'OK':
            stamp = HCGB_time.print_time_stamp(filename_stamp)

        return(code)    
    
    
