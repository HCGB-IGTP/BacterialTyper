#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Retrieves genes ids from profile analysis generated
"""
## useful imports
import os
import re
import sys
from termcolor import colored
import pandas as pd
from Bio import SeqIO

## import my HCGB module 
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.aesthetics_functions as HCGB_aes

##############
def help_options():
    print ("\nUSAGE: python %s gene_id profile_folder...\n"  %os.path.realpath(__file__))

##############
def help_retrieve_genes():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##############
def retrieve_genes_ids_sequences(profile, gene_ID, debug):
    """    
    """
    ## given a profile folder
    if debug:
        HCGB_aes.debug_message('profile: ', 'yellow')
        print (profile)
        HCGB_aes.debug_message('gene_id: ', 'yellow')
        print (gene_ID)
        
    ##
    assembled_genes_list = []
    assembled_genes_list = HCGB_main.retrieve_matching_files(profile, "assembled_genes.fa", debug)
    assembled_genes_list = [s for s in assembled_genes_list if 'ariba.tmp' not in s]
    
    ## return if empty
    if len(assembled_genes_list) == 0:
        return('','')
    
    if debug:
        HCGB_aes.debug_message('assembled_genes_list: ', 'yellow')
        print(assembled_genes_list)
        
    if os.path.isfile(assembled_genes_list[0]):
        for record in SeqIO.parse(assembled_genes_list[0], "fasta"):
            if debug:
                HCGB_aes.debug_message('record.description: ', 'yellow')
                print(record.description)
 
            search_ID = re.search(gene_ID, record.description)
            if (search_ID):
                return (record.id, str(record.seq))

        return('','')

##############
def get_genes_profile(samples_info, gene_names, debug, option, groupby_id="name"):
    """    
    """
    
    if debug:
        HCGB_aes.debug_message('gene_names: ', 'yellow')
        print(gene_names)
        HCGB_aes.debug_message('option: ' + option, 'yellow')
        
    
    ## search by group id or gene name
    print ('\n+ Retrieve selected genes profile for each sample.')
    results_profileIDs = pd.DataFrame()
    sample_frame = samples_info.groupby([groupby_id])
    for g in gene_names:
        #print ("\t+", g)
        for name_tuple, cluster_df in sample_frame:
            name = name_tuple[0]
            ## get list of files retrieved from profile
            my_list_profiles = cluster_df.loc[cluster_df['tag'] == 'profile']['sample'].to_list() 
	       
            if debug:
                HCGB_aes.debug_message('name: ' + name, 'yellow')
                HCGB_aes.debug_message('my_list_profiles: ', 'yellow')
                print (my_list_profiles)
                HCGB_aes.debug_message('cluster_df: ', 'yellow')
                print (cluster_df)

            ## skip files
            if name == 'report':
                continue

            fill=False
            for profile_csv in my_list_profiles:
                if debug:
                    HCGB_aes.debug_message('profile_csv: ' + profile_csv, 'yellow')

                ## skip files
                if profile_csv.endswith('report_summary.csv'):
                    pass
                
                elif profile_csv.endswith('norm.tsv'): ## amrfinder
                    
                    value = retrieve_genes_ids_profile(profile_csv, g, debug, option, sep_field="\t")
                    ## save results 
                    if (not value.empty):
                        for Name, Data in value.iterrows():
                            results_profileIDs.loc[name,Name] = "Yes"
                        fill=True
                    break
                
                else: 
                    continue
                
                    ## TODO: Control for ariba if necessary
                    ## ARIBA
                    value = retrieve_genes_ids_profile(profile_csv, g, debug, option)
                    ## save results 
                    if (not value.empty):
                        for Name, Data in value.iterrows():
                            results_profileIDs.loc[name,Name] = Data['Status']
                        fill=True
                    break

            if not fill:
                results_profileIDs.loc[name, g] = 'No'

    return (results_profileIDs)

##############
def retrieve_genes_ids_profile(profile, gene_ID, debug, option, sep_field=","):
    """    
    """
    ## read data    
    get_csv_data = HCGB_main.get_data(profile, sep_field, '')
    
    if option == 'name': ## ariba
        list_Genes = get_csv_data['Genes'].to_list()
        get_csv_data.index = get_csv_data['Genes']
    
    elif option == 'ID': ## ariba
        list_Genes = get_csv_data['ID'].to_list()
        get_csv_data.index = get_csv_data['ID']
    
    elif option == 'Gene symbol': ## gene symbol by amrfinderplus
        list_Genes = get_csv_data['Gene symbol'].to_list()
        get_csv_data.index = get_csv_data['Gene symbol']
        

    
    ## debug messages
    if debug:
        HCGB_aes.debug_message('profile: ' + profile, 'yellow')
        HCGB_aes.debug_message('gene_id: ' + str(gene_ID), 'yellow')
        HCGB_aes.debug_message('data: ', 'yellow')
        print(get_csv_data)
        HCGB_aes.debug_message('Option: ' + option, 'yellow')
        HCGB_aes.debug_message('Genes: ', 'yellow')
        print (list_Genes)
        
    ## search accordingly
    if option == 'name':
        regex_search = re.compile("^" + gene_ID + ".*")
        filtered_genes = list(filter(regex_search.match, list_Genes))
        
        ## debug messages
        if debug:
            HCGB_aes.debug_message('filtered_genes: ', 'yellow')
            print (filtered_genes)
            HCGB_aes.debug_message('filtered_genes.loc[filtered_genes]: ', 'yellow')
            print (get_csv_data.loc[filtered_genes])
        
        return (get_csv_data.loc[filtered_genes]) 
        
    else:
        ##
        if gene_ID in list_Genes:
            ## debug messages
            if debug:
                HCGB_aes.debug_message('gene_id: ' + gene_ID, 'yellow')
                print (get_csv_data.loc[gene_ID].to_frame().transpose())
                
            return (get_csv_data.loc[gene_ID].to_frame().transpose())
        else:
            return(pd.DataFrame()) 
    
##############
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()

    ##
    gene_ID = sys.argv[1]
    profile = os.path.abspath( sys.argv[2] )
    
    retrieve_genes_ids_profile(profile, gene_ID, True)

'''******************************************'''
if __name__== "__main__":
    main()
