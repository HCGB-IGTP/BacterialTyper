#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Creates sccmec typing
"""
## useful imports
import os
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.config import set_config

## import my HCGB module 
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.aesthetics_functions as HCGB_aes

########################################
def help_options():
    print ("\nUSAGE: python %s name assembly_fasta folder...\n"  %os.path.realpath(__file__))

########################################
def module_call(dict_assemblies, dict_folders, debug=False):
    """Create sccmec call and control for parameters"""
    
    ## info2return
    sccmec_bin = set_config.get_exe('staphopia-sccmec')
    path2database = os.path.abspath( os.path.join( os.path.dirname(sccmec_bin), '..', 'share', "staphopia-sccmec", 'data') )
    info_dict={ 'sccmec database': path2database}
    
    print ("+ Checking SCCmec type and subtypes for each sample retrieved...")
    ## 
    sccmec_results = pd.DataFrame()
    
    ## No need to optimize. There is a problem with the working dir of sccmec and we 
    ## need to change every time.
    for name, assembly_file in dict_assemblies.items():
        report_folder = HCGB_files.create_folder(dict_folders[name])
        sample_folder = HCGB_files.create_subfolder('sccmec_typing', report_folder) 
        ## check if previously done and succeeded
        filename_stamp = sample_folder + '/.success'
        results_csv = os.path.join(sample_folder, 'sccmec_results.csv')

        if os.path.isfile(filename_stamp):
            stamp =  HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
            info_sample = HCGB_main.get_data(results_csv, ',', 'index_col=0')
        else:
            (info_sample, sccmec_type) = sccmec_call(name, assembly_file, sample_folder, debug)
            
            if (info_sample.shape[0] == 0):
                print("+ Some error occurred with sample %s. Please re-run analysis or check log files." %name)
                continue
            else:
                info_sample['sccmec_type'] = sccmec_type
                sccmec_results.to_csv(results_csv)
                
                ## success
                HCGB_time.print_time_stamp(filename_stamp)

        ## merge results
        sccmec_results = pd.concat([sccmec_results, info_sample], join='outer')
                

    print ("+ Jobs finished\n+ Collecting information for all samples...")
    
    ## debug messages
    if debug:
        HCGB_aes.debug_message('sccmec_results', 'yellow')
        HCGB_main.print_all_pandaDF(sccmec_results)
    
    return(sccmec_results, info_dict)


########################################
def sccmec_call(sample, assembly_file, folder, debug=False):
    """staphopia-sccmec call and check results."""
    
    ## prepare call
    sccmec_results = os.path.join(folder, "sccmec_results.json")
    err_call = os.path.join(folder, "sccmec_cmd.err")
    sccmec_bin = set_config.get_exe('staphopia-sccmec')
    
    ## system call
    cmd_call = "%s --assembly %s --json --hamming > %s 2> %s " %(sccmec_bin, 
                                               assembly_file, sccmec_results, err_call) ## use json and hamming distance as options
    status = HCGB_sys.system_call(cmd_call)
    
    ## check results
    ## see https://github.com/staphopia/staphopia-sccmec
    results = pd.DataFrame.from_dict( HCGB_main.read_json_file(sccmec_results, debug)[0], orient='index').transpose()
    res = [col for col in results if (results[col] == 0).any()]
    
    if debug:
        print(results)
        print(res)
    
    ## no results
    if len(res) == 0:
        return(results, "NA")
    
    ## get_summary
    if 'meca' in res:
        sccmec_type = "MRSA:"
        res.remove('meca')
    else:
        sccmec_type = "MSSA:"
    
    sccmec_type += ":".join(res)
    
    ### 
    return(results, sccmec_type)

    
########################################
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
        
    name = argv[1]
    fasta_file = os.path.abspath(argv[2])
    folder = os.path.abspath(argv[3])
    debug=True

    ## path
    folder = HCGB_files.create_folder(folder)
    os.chdir(folder)

    ###
    sccmec_call(name, fasta_file, folder, debug)
    

'''******************************************'''
if __name__== "__main__":
    main()

