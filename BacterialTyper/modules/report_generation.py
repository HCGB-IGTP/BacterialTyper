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
from BacterialTyper.scripts import functions
from BacterialTyper.scripts import database_user
from BacterialTyper.config import set_config
from BacterialTyper.report.Staphylococcus import get_spa_typing
from BacterialTyper.report import retrieve_genes

## example report check: https://github.com/tseemann/nullarbor

## for staphyloccocus aureus: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0593-7
## for klebsiella: kleborate
## for micobacterium tuberculosis: MTBSeq

## R package: http://bioconductor.org/packages/release/bioc/html/ReportingTools.html

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
    if options.species_report == "Saureus":
        Saureus_specific(pd_samples_retrieved, pd_samples_info, species_specific_df, options)
        
        ## time stamp
        start_time_partial = functions.timestamp(start_time_partial)

    ## create gene specific report if any
    if options.genes_ids:
        ## given a list of genes ids, retrieve sequence for all samples from profile 
        if os.path.isfile(os.path.abspath(options.genes_ids)):
            in_file = os.path.abspath(options.genes_ids)
            gene_names = [line.rstrip('\n') for line in open(in_file)]
            print ('+ Retrieve selected genes from the profile analysis for each sample.')
            print ('+ Searching gene:')
            
            ## get profiles available
            results_geneIDs = pd.DataFrame(columns=('sample', 'gene', 'id', 'sequence'))
            sample_frame = pd_samples_info.groupby(["name"])
            for g in gene_names:
                print ("\t -", g)
                for name, cluster_df in sample_frame:
                    my_list_profiles = cluster_df.loc[cluster_df['tag'] == 'profile']['ext'].to_list()
                    if options.debug:
                        print ("name: ", name)
                        print ("my_list_profiles:")
                        print (my_list_profiles)
                    
                    for p in my_list_profiles:
                        main_profile_folder = cluster_df['ext']['dirname'][0]
                        if p == 'VFDB':
                            p = p.lower() + '_full'
                        
                        profile_folder = os.path.join(main_profile_folder, p)
                        (seq_id, seq_sequence) = retrieve_genes.retrieve_genes(profile_folder, g, Debug)
                        results_geneIDs.loc(length(results_geneIDs)) = (name, g, seq_id, seq_sequence)

        ## debug
        if Debug:
            print (results_geneIDs) 
    
    
    if options.genes_fasta:
        ## given a list of fasta sequences search using blast against proteins annotated or genome
        print()    
    
    
    print ("\n*************** Finish *******************")
    start_time_partial = functions.timestamp(start_time_total)

    print ("+ Exiting Report generation module.")
    exit()
    
    ## time stamp
    start_time_partial = functions.timestamp(start_time_total)        


#######################3
def Saureus_specific(samples_df, samples_info, results_df, options):
    """
    Retrieves Saureus specific information.
    
    See additional information in :doc:`../../user_guide/report/Saureus/saureus_report`
    """
    
    ## get European Quality Control genes
    Staphylococcus_path = os.path.abspath( os.path.join( os.path.realpath(__file__), '..', '..', 'report', 'Staphylococcus'))
    EQC_genes = os.path.join(Staphylococcus_path, "EQC_genes.csv")
    arcA_gene = os.path.join(Staphylococcus_path, "arcA.fasta")
    
    EQC_genes_df = functions.get_data(EQC_genes, ',')


    if options.debug:
        print ("## DEBUG: Saureus_specific")
        print (Staphylococcus_path)
        print (EQC_genes)
        print (arcA_gene)
    
    ## get spatyping
    assembly_files = samples_df.loc[samples_df['tag'] == "assembly", "sample"]
    results_df = get_spa_typing.module_call(options.database, assembly_files.to_dict(), options.debug)
    
    results_df.to_csv("/home/labs/lslab/jsanchez/proc_data/20200405_CPrat_SaureusPaper/analysis/test_spaTyper.csv")    
    ## get sccmec
    
    
    
