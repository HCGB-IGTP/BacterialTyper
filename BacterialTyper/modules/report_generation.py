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

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.scripts import database_user
from BacterialTyper.config import set_config
from BacterialTyper.report import retrieve_genes
from BacterialTyper.report.Staphylococcus import get_spa_typing

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
    
    ########################################
    ## create species specific report if any
    ########################################
    if (options.species_report):
        ## Saureus
        if options.species_report == "Saureus":
            Saureus_specific(pd_samples_retrieved, pd_samples_info, options, summary_report)
        
        ## else
        ## to add accordingly        
            
        ## time stamp
        start_time_partial = functions.timestamp(start_time_partial)
        
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
        genes_folder = functions.create_subfolder('genes', summary_report)
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
        start_time_partial = functions.timestamp(start_time_partial)
    
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
            start_time_partial = functions.timestamp(start_time_partial)
        
    ###############################################
    ## Search for any additional fasta sequence
    ###############################################
    if options.genes_fasta:
        ## given a list of fasta sequences search using blast against proteins annotated or genome
        print("** THIS OPTION IS NOT IMPLEMENTED YET... **")
    
    print ("\n*************** Finish *******************")
    start_time_partial = functions.timestamp(start_time_total)

    print ("+ Exiting Report generation module.")
    return()

#######################3
def Saureus_specific(samples_df, samples_info, options, folder):
    """
    Retrieves Saureus specific information.
    
    See additional information in :doc:`../../user_guide/report/Saureus/saureus_report`
    """
    
    ########################################
    ## get European Quality Control genes
    ########################################
    Staphylococcus_path = os.path.abspath( os.path.join( os.path.realpath(__file__), '..', '..', 'report', 'Staphylococcus'))
    EQC_genes = os.path.join(Staphylococcus_path, "EQC_genes.csv")
    arcA_gene = os.path.join(Staphylococcus_path, "arcA.fasta")
    
    EQC_genes_df = functions.get_data(EQC_genes, ',', '')
    ## Gene,ID,Source
    ## mecA,ARO:3000617,CARD
    ## mecC,ARO:3001209,CARD
    ## mupA,ARO:3000521,CARD
    
    ## debugging messages
    if options.debug:
        print ("## DEBUG: Saureus_specific")
        print (Staphylococcus_path)
        print (arcA_gene)
        
        print ("EQC genes")
        print (EQC_genes)
        print (EQC_genes_df)
 

    ####################
    ## get gene info by unique ID
    ####################
    ## get gene names
    gene_IDs = EQC_genes_df['ID'].to_list()
    
    results_Profiles_ids = retrieve_genes.get_genes_profile(samples_info, gene_IDs, options.debug, 'ID')
    if options.debug:
        print ("results_Profiles")
        print (results_Profiles)
    
    ########################################
    ## add additional genes if required
    ########################################
    if options.genes_ids_profile:
        in_file = os.path.abspath(options.genes_ids_profile)
        gene_names = [line.rstrip('\n') for line in open(in_file)]
        
        if options.debug:
            print ("gene_names")
            print (gene_names)
    
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
    samples_df = samples_df.set_index('name')
    assembly_files = samples_df.loc[samples_df['tag'] == "assembly", "sample"]
    results_spaType = pd.DataFrame()
    results_spaType = get_spa_typing.module_call(options.database, assembly_files.to_dict(), options.debug)
    
    ####################
    ## get sccmec
    ####################
    ## todo
    
    ####################
    ## save results
    ####################
    ## open excel writer
    name_excel = folder + '/Saureus_report.xlsx'
    writer = pd.ExcelWriter(name_excel, engine='xlsxwriter')
    
    # results_Profiles ids
    results_Profiles_ids.to_excel(writer, sheet_name="gene_ids")

    if options.genes_ids_profile:
        # results_Profiles names
        results_Profiles_ids.to_excel(writer, sheet_name="gene_names")

    # results_spaType
    results_spaType.to_excel(writer, sheet_name="spaTyper")

    ## close
    writer.save()

