#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Retrieves genes ids from profile analysis generated
"""
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd
from Bio import SeqIO

## import my modules
from BacterialTyper.config import set_config

## import my HCGB module 
from HCGB.functions import main_functions

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
        print ("** DEBUG **")
        print ("profile: ", profile)
        print ("gene_id: ", gene_ID)
        
    ##
    assembled_genes_list = main_functions.retrieve_matching_files(profile, "assembled_genes.fa", debug)
    assembled_genes_list = [s for s in assembled_genes_list if 'ariba.tmp' not in s]
    
    if debug:
        print ("** DEBUG **")
        print ("assembled_genes_list: ", assembled_genes_list)

    if os.path.isfile(assembled_genes_list[0]):
        for record in SeqIO.parse(assembled_genes_list[0], "fasta"):
            if debug:
                print ("** DEBUG: ", record.description, " **")
                #print (record.seq)
 
            search_ID = re.search(gene_ID, record.description)
            if (search_ID):
                return (record.id, str(record.seq))

        return('','')

##############
def get_genes_profile(samples_info, gene_names, debug, option):
    """    
    """
    ## search by group id or gene name
    print ('\n+ Retrieve selected genes profile for each sample.')
    results_profileIDs = pd.DataFrame()
    sample_frame = samples_info.groupby(["name"])
    for g in gene_names:
        #print ("\t+", g)
        for name, cluster_df in sample_frame:
            my_list_profiles = cluster_df.loc[cluster_df['tag'] == 'profile']['ext'].to_list()
	   
            if debug:
                print ("name: ", name)
                print ("my_list_profiles:")
                print (my_list_profiles)
                print ("cluster_df")
                print (cluster_df)

            fill=False
            for p in my_list_profiles:
                profile_csv = cluster_df.loc[cluster_df['ext'] == p]['sample'].to_list()[0]
                if debug:
                    print ("profile_csv: ", profile_csv)
            
                value = retrieve_genes_ids_profile(profile_csv, g, debug, option)
            
                ## save results 
                if (not value.empty):
                    for Name, Data in value.iterrows():
                        results_profileIDs.loc[name,Name] = Data['Status']
                    fill=True

            if not fill:
                results_profileIDs.loc[name, g] = 'no'

    return (results_profileIDs)

##############
def retrieve_genes_ids_profile(profile, gene_ID, debug, option):
    """    
    """
    ## read data    
    get_csv_data = functions.main_functions.get_data(profile, ',', '')
    
    if option == 'name':
        list_Genes = get_csv_data['Genes'].to_list()
        get_csv_data.index = get_csv_data['Genes']
    elif option == 'ID':
        list_Genes = get_csv_data['ID'].to_list()
        get_csv_data.index = get_csv_data['ID']

    
    ## debug messages
    if debug:
        print ("** DEBUG **")
        print ("profile: ", profile)
        print ("gene_id: ", gene_ID)
        print ("data")
        print (get_csv_data)
        print ("Option: ", option)
        print ("Genes")
        print (list_Genes)
        
    ## search accordingly
    if option == 'name':
        regex_search = re.compile("^" + gene_ID + ".*")
        filtered_genes = list(filter(regex_search.match, list_Genes))
        
        ## debug messages
        if debug:
            print ("** DEBUG **")
            print (filtered_genes)
            print (get_csv_data.loc[filtered_genes])
        
        return (get_csv_data.loc[filtered_genes]) 
        
    else:
        if gene_ID in list_Genes:
            ## debug messages
            if debug:
                print ("** DEBUG **")
                print ("gene_ID: ", gene_ID)
                print (get_csv_data.loc[gene_ID])
            return (get_csv_data.loc[gene_ID])
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
