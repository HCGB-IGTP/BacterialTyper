#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:34:19 2023

@author: jsanchez
"""

from BacterialTyper.config import set_config
import os
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time


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

########################################
def MLST_call(outfolder, assembly_file, mlst_profile, sample_name, 
              minid=95, mincov=10, minscore=50, debug=False):
    mlst_bin = set_config.get_exe("mlst",Debug=debug)
    
    
    outfolder = HCGB_files.create_subfolder('MLST', outfolder)
    out_csv = os.path.join(outfolder, 'MLST_res.csv')
    out_json = os.path.join(outfolder, 'MLST_res.json')
    out_err  = os.path.join(outfolder, 'MLST_res.log')
    
    filename_stamp = os.path.join(outfolder, '.success_mlst')
    if HCGB_files.is_non_zero_file(filename_stamp):
        return True
    else:
        
        ## init
        cmd_mlst = mlst_bin
        
        ## add label
        cmd_mlst += " --label " + sample_name
        
        ## add options
        cmd_mlst += " --minid %s --mincov %s" %(minid, mincov)
        
        if debug:
            cmd_mlst += " --debug"    
        
        ## add output
        cmd_mlst += " --csv --json " + out_json
        
        cmd_mlst1 = cmd_mlst
        
        if mlst_profile:
            cmd_mlst += " --legacy --scheme %s  --minscore %s" %(mlst_profile, minscore)    
        else:
            ## mlst2use is not available or it was not determined
            pass
        
        ## add input assemble
        ## add output
        cmd_mlst += ' %s > %s 2> %s' %(assembly_file, out_csv, out_err)
        
        if debug:
            HCGB_aes.debug_message("call:", color='yellow')
            print(cmd_mlst)
        
        code = HCGB_sys.system_call(cmd_mlst)
        if (code == 'OK'):
            ## success stamps
            HCGB_time.print_time_stamp(filename_stamp)
            return True
        else:
            ## Sometime the scheme name provided is not 100% correct and creates error: mtuberculosis -> mtuberculosis_2
            ## re-run in auto scheme lineage selection
            cmd_mlst1 += ' %s > %s 2> %s' %(assembly_file, out_csv, out_err)
            code = HCGB_sys.system_call(cmd_mlst1)
            if (code == 'OK'):
                ## success stamps
                HCGB_time.print_time_stamp(filename_stamp)
                return True
            else:
                return False
    
    
    
    
    
    
    