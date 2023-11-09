#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:34:19 2023

@author: jsanchez
"""

from BacterialTyper.config import set_config
import os
import pandas as pd
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main
from termcolor import colored
import pprint

########################################
def help_MLST(debug=False):

    mlst_bin = set_config.get_exe("mlst",Debug=debug)
    db_path = os.path.abspath(os.path.join( os.path.dirname(mlst_bin), "../db/pubmlst"))
    
    print("\n\n")
    print("MLST help messages:\nAutomatic MLST calling from assembled contigs using mlst software: https://github.com/tseemann/mlst")
    print("\n")
    print("""MLST Sequences and profiles are automatically downloaded from PubMLST REST API (https://rest.pubmlst.org/db/)
into the mlst installation folder within the environment:\n""" + db_path)

    # show species available
    print("Find a list of putative species to use here: ")
    content_file = get_MLST_profiles()

    pprint.pprint(content_file)
    print("\n")


########################################
def get_MLST_profiles(debug=False):
    mlst_bin = set_config.get_exe("mlst",Debug=debug)
    scheme_file = os.path.join( 
        os.path.abspath(os.path.join( os.path.dirname(mlst_bin), "../db")), "scheme_species_map.tab")
    
    print(scheme_file)
    if HCGB_files.is_non_zero_file(scheme_file):
        with open(scheme_file, 'r') as reader:
            content_file = reader.readlines()
        
        my_dict = {}
        for line_file in content_file:
            line_file = line_file.replace("\n","")
            line_file = line_file.replace(" ","")
            line_file_split = line_file.split("\t")
            my_dict[line_file_split[0]] = line_file_split[1] + " " +  line_file_split[2]
                
        return(my_dict)
    
    else:
        print("\n")
        print(colored("ERROR: File not available", 'red'))
        print(colored("ERROR: Check mlst software installation", 'red'))
        raise SystemExit()
        
