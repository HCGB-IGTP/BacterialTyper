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
from Bio import SeqIO
import shutil

## import my modules
from BacterialTyper.scripts import database_user
from BacterialTyper.config import set_config
from BacterialTyper.report import retrieve_genes
from BacterialTyper.report import get_promoter

from BacterialTyper.report.Staphylococcus import get_spa_typing
from BacterialTyper.report.Staphylococcus import agr_typing
from BacterialTyper.report.Staphylococcus import get_sccmec

from BacterialTyper import __version__ as pipeline_version

##
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.info_functions as HCGB_info


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

    if (options.help_spaTyper):
        ## help_format option
        get_spa_typing.help_spaTyper()
        exit()

    elif (options.help_project):
        ## information for project
        help_info.project_help()
        exit()

    ## set default
    options.batch = False

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
    HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
    HCGB_aes.boxymcboxface("Report generation module")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()
    
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
    pd_samples_retrieved['new_name'] = pd_samples_retrieved['name']
    
    ## get info: profile, ident, cluster, MGE
    pd_samples_info = database_user.get_userData_info(options, input_dir)
    
    ## get databases to list
    #retrieve_databases = get_options_db(options)

    ## create output files
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "report", options.debug)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)
        
        print (colored("**DEBUG: pd_samples_info **", 'yellow'))
        print (pd_samples_info)
    
    ## generate output folder, if necessary
    print ("\n\n\n+ Generate a report summarizing analysis and sample information")
    if not options.project:
        HCGB_files.create_folder(outdir)
        outdir_report = outdir
    else:
        ### report generation
        outdir_report = HCGB_files.create_subfolder("report", outdir)
    
    ## create report with all data
    summary_report = HCGB_files.create_subfolder("summary_report", outdir_report)
    print ("Folder: ", summary_report)
    
    ## time stamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)
    
    ## info2dump
    info_dict = {'folder':summary_report}
    
    ########################################
    ## create species specific report if any
    ########################################
    if (options.species_report):
        ## Saureus
        if options.species_report == "Saureus":
            info_returned = Saureus_specific(pd_samples_retrieved, pd_samples_info, options, summary_report, outdir_dict)
            info_dict['Saureus'] = info_returned
        
        ## else
        ## to add accordingly        
            
        ## time stamp
        start_time_partial = HCGB_time.timestamp(start_time_partial)
        
    ###########################################################
    ## create gene fasta sequences retrieval if desired
    ###########################################################
    if options.genes_ids_fasta:
        ## given a list of genes ids, retrieve sequence for all samples from profile 
        if os.path.isfile(os.path.abspath(options.genes_ids_fasta)):
            in_file = os.path.abspath(options.genes_ids_fasta)
            gene_names = [line.rstrip('\n') for line in open(in_file)]
            print ('+ Retrieve selected genes sequences from the profile analysis for each sample.')
            print ('+ Searching gene:')
            
            ## get profiles available
            results_geneIDs = pd.DataFrame(columns=('sample', 'gene', 'id', 'sequence'))
            sample_frame = pd_samples_info.groupby(["name"])
            for g in gene_names:
                print ("\t+", g)
                for name, cluster_df in sample_frame:
                    my_list_profiles = cluster_df.loc[cluster_df['tag'] == 'profile']['ext'].to_list()
                    if options.debug:
                        print ("name: ", name)
                        print ("my_list_profiles:")
                        print (my_list_profiles)
                    
                    for p in my_list_profiles:
                        main_profile_folder = cluster_df.loc[cluster_df['ext'] == p]['dirname'].to_list()[0]
                        p = p.lower()
                        if p == 'vfdb':
                            p = p + '_full'
                        
                        profile_folder = os.path.join(main_profile_folder, p)
                        (seq_id, seq_sequence) = retrieve_genes.retrieve_genes_ids_sequences(profile_folder, g, Debug)
                        if (seq_id):
                             ## save results 
                             results_geneIDs.loc[len(results_geneIDs)] = (name, g, seq_id, seq_sequence)

        ## save for each gene in a separate fasta file
        list_of_genes = set(results_geneIDs['gene'].to_list())

        ## debug
        if Debug:
            print ("** DEBUG **")
            print (results_geneIDs) 
            print (list_of_genes)
        
        ## Save results
        genes_folder = HCGB_files.create_subfolder('genes', summary_report)
        for gene_retrieved in list_of_genes:
            this_frame = results_geneIDs[results_geneIDs['gene'] == gene_retrieved]
        
            gene_retrieved_file = os.path.join(genes_folder, gene_retrieved)
            gene_retrieved_fasta = gene_retrieved_file + ".fasta"
            gene_retrieved_info = gene_retrieved_file + "_info.txt"
            fasta_hd = open(gene_retrieved_fasta, 'w')
            info_hd = open(gene_retrieved_info, 'w')
            
            for item, row in this_frame.iterrows():
                string2write = ">" + row['sample'] + '_' + row['gene'] + '\n' + row['sequence'] + '\n'  
                string2write_info = row['sample'] + '\t' + row['gene'] + '\t' + row['id'] + '\n'
                fasta_hd.write(string2write)
                info_hd.write(string2write_info)
                
            fasta_hd.close()
            info_hd.close()
        
        ## time stamp
        start_time_partial = HCGB_time.timestamp(start_time_partial)
    
        ########################################
        ## create gene promoter fasta sequences retrieval if desired
        ########################################
        if options.promoter_bp:
            ## retrieve as many bp as necessary from genes_ids_fasta
            print("** THIS OPTION IS NOT IMPLEMENTED YET... **")
            #get_promoter.get_promoter(file, geneOfInterest, basePairs, sampleName, option, debug=False):
    
    ########################################
    ## create gene specific report if any
    ########################################
    if options.genes_ids_profile:
        if options.species_report == "Saureus":
            if Debug:
                print ("** options.genes_ids_profile **")
                print ("Analysis already done for Saureus")
        else:
            in_file = os.path.abspath(options.genes_ids_profile)
            gene_names = [line.rstrip('\n') for line in open(in_file)]
            results_Profiles = retrieve_genes.get_genes_profile(pd_samples_info, gene_names, options.debug, "name")
            if options.debug:
                print ("results_Profiles")
                print (results_Profiles)
        
            
            ## open excel writer
            name_excel = summary_report + '/gene_ids_profile.xlsx'
            writer = pd.ExcelWriter(name_excel, engine='xlsxwriter')
            results_Profiles.to_excel(writer, sheet_name="gene_ids") 
            
            ## close
            writer.save()

            ## time stamp
            start_time_partial = HCGB_time.timestamp(start_time_partial)
        
    ###############################################
    ## Search for any additional fasta sequence
    ###############################################
    if options.genes_fasta:
        ## given a list of fasta sequences search using blast against proteins annotated or genome
        print("** THIS OPTION IS NOT IMPLEMENTED YET... **")
    
    print ("\n*************** Finish *******************")
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ## dump
    ## samples information dictionary
    samples_info = {}
    samples_frame = pd_samples_retrieved.groupby('new_name')
    for name, grouped in samples_frame:
        samples_info[name] = grouped['sample'].to_list()

    ## 
    samples_info2 = {'example':'test'}
    print(pd_samples_info)
    
    ## dump information and parameters
    info_dir = HCGB_files.create_subfolder("info", outdir)
    print("+ Dumping information and parameters")
    runInfo = { "module":"report", 
               "time":time.time(),
               "BacterialTyper version":pipeline_version,
               'sample_info': samples_info,
               #'sample_info2': samples_info2,
               'info_analysis': info_dict
                }
    
    HCGB_info.dump_info_run(info_dir, 'report', options, runInfo, options.debug)

    print ("+ Exiting Report generation module.")
    return()

#######################3
def Saureus_specific(samples_df, samples_info, options, folder, outdir_dict):
    """
    Retrieves Saureus specific information.
    
    See additional information in :doc:`../../user_guide/report/Saureus/saureus_report`
    """
    HCGB_aes.boxymcboxface("Staphylococcus aureus report submodule")

    ########################################
    ## get European Quality Control genes
    ########################################
    Staphylococcus_path = os.path.abspath( os.path.join( os.path.realpath(__file__), '..', '..', 'report', 'Staphylococcus'))
    EQC_genes = os.path.join(Staphylococcus_path, "EQC_genes.csv")
    arcA_gene = os.path.join(Staphylococcus_path, "arcA.fasta")
    
    EQC_genes_df = HCGB_main.get_data(EQC_genes, ',', '')
    ## Gene,ID,Source
    ## mecA,ARO:3000617,CARD
    ## mecC,ARO:3001209,CARD
    ## mupA,ARO:3000521,CARD
    
    
    ## info Saureus
    info_Saures = {}
    
    ## debugging messages
    if options.debug:
        HCGB_aes.debug_message("Saureus_specific", 'yellow')
        print (Staphylococcus_path)
        print (arcA_gene)
        
        HCGB_aes.debug_message("EQC_genes", 'yellow')
        print (EQC_genes)
        print (EQC_genes_df)

    ####################
    ## get gene info by unique ID
    ####################
    ## get gene names
    gene_IDs = EQC_genes_df['ID'].to_list()
    
    results_Profiles_ids = pd.DataFrame()
    
    ## dataframe is NOT empty
    if samples_info.shape[0] != 0:
        
        HCGB_aes.print_sepLine('+', 35, 'yellow')
        print("Retrieve genes profiles")
        HCGB_aes.print_sepLine('+', 35, 'yellow')
        
        ## get profiles
        results_Profiles_ids = retrieve_genes.get_genes_profile(samples_info, gene_IDs, options.debug, 'ID')
        if options.debug:
            HCGB_aes.debug_message("results_Profiles_ids", 'yellow')
            print (results_Profiles_ids)
        
    ########################################
    ## add additional genes if required
    ########################################
    if options.genes_ids_profile:
        in_file = os.path.abspath(options.genes_ids_profile)
        gene_names = [line.rstrip('\n') for line in open(in_file)]
        
        if options.debug:
            print ("gene_names")
            print (gene_names)
        
        ## dataframe is NOT empty
        if samples_info.shape[0] != 0:
            ## outdir_dict
            results_Profiles_names = retrieve_genes.get_genes_profile(samples_info, gene_names, options.debug, 'name')
            if options.debug:
                print ("results_Profiles")
                print (results_Profiles)
            
            
    #################################
    ## get blast sequence         ###
    #################################
    # arcA_gene
    
    ####################
    ## get spatyping  ##
    ####################
    HCGB_aes.print_sepLine('+', 35, 'yellow')
    print("Get SPA typing")
    HCGB_aes.print_sepLine('+', 35, 'yellow')

    samples_df = samples_df.set_index('name')
    assembly_files = samples_df.loc[samples_df['tag'] == "assembly", "sample"]
    results_spaType = pd.DataFrame()
    (results_spaType, info_spa) = get_spa_typing.module_call(options.database, assembly_files.to_dict(), outdir_dict, options.debug)
    info_Saures['spa typing'] = info_spa
    
    ####################
    ## get agr typing
    ####################
    print()
    HCGB_aes.print_sepLine('+', 35, 'yellow')
    print("Get agr typing")
    HCGB_aes.print_sepLine('+', 35, 'yellow')

    (agr_results, info_agr) = agr_typing.agrvate_caller(assembly_files.to_dict(), outdir_dict, options.debug)
    info_Saures['agr typing'] = info_agr
    
    ## copy excel file and operon into report folder
    ## remove from dataframe
    
    files2copy = agr_results['operon_fna'].to_list() + agr_results['agr_operon_xlsx'].to_list()
    files2copy = list(filter(None, files2copy))
    arg_results_folder = HCGB_files.create_subfolder('agr_results', folder)
    for f in files2copy:
        shutil.copy(f, arg_results_folder)

    del agr_results['operon_fna']
    del agr_results['agr_operon_xlsx']
    
    ####################
    ## get sccmec
    ####################
    print()
    HCGB_aes.print_sepLine('+', 35, 'yellow')
    print("Get sccmec typing")
    HCGB_aes.print_sepLine('+', 35, 'yellow')

    (sccmec_results, info_sccmec) = get_sccmec.module_call(assembly_files.to_dict(), outdir_dict, options.debug)
    info_Saures['sccmec typing'] = info_sccmec
        
    ####################
    ## save results
    ####################
    print()
    ## open excel writer
    name_excel = os.path.join(folder, 'Saureus_report.xlsx')
    writer = pd.ExcelWriter(name_excel, engine='xlsxwriter')
    
    # results_Profiles ids
    if results_Profiles_ids.shape[0] > 0:
        results_Profiles_ids.to_excel(writer, sheet_name="gene_ids")

    if options.genes_ids_profile:
        # results_Profiles names
        results_Profiles_names.to_excel(writer, sheet_name="gene_names")

    # results_spaType
    results_spaType.to_excel(writer, sheet_name="spaTyper")

    # agr_results
    agr_results.to_excel(writer, sheet_name="agr typing")

    # sccmec_results
    sccmec_results.to_excel(writer, sheet_name="sccmec")

    ## close
    writer.save()
    
    return(info_Saures)

