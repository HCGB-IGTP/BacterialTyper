#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Creates agr typing
"""
## useful imports
import os
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd
from Bio import SeqIO

## import my modules
from BacterialTyper.config import set_config

## import my HCGB module 
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.aesthetics_functions as HCGB_aes

##############
def help_options():
    print ("\nUSAGE: python %s path sample...\n"  %os.path.realpath(__file__))

##############################
def agrvate_caller(dict_assemblies, dict_folders, debug=False):
    """Create agrvate call and control for parameters"""
    
    ## ATTENTION: agrvate needs to chdir to output folder
    path_here = os.getcwd()
    
    ## info2return
    agrvate_bin = set_config.get_exe('agrvate')
    info_dict={ 'agrvate database': os.path.join(os.path.basename(agrvate_bin), "agrvate_databases")}
    
    
    print ("+ Checking agr genes for each sample retrieved...")
    
    agrvate_results = pd.DataFrame()
    
    ## No need to optimize. There is a problem with the working dir of agrvate and we 
    ## need to change every time.
    for name, assembly_file in dict_assemblies.items():
        report_folder = HCGB_files.create_folder(dict_folders[name])
        sample_folder = HCGB_files.create_subfolder('agr_typing', report_folder) 
        ## check if previously done and succeeded
        filename_stamp = sample_folder + '/.success'
        if os.path.isfile(filename_stamp):
            stamp =  HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
            info_sample = get_results_agrvate(assembly_file, sample_folder, name, debug) 
        else:
            os.chdir(sample_folder)
            info_sample = agrvate_call(name, assembly_file, sample_folder, debug)
        
            if (info_sample.shape[0] == 0):
                print("+ Some error occurred with sample %s. Please re-run analysis or check log files." %name)
            else:
                ## success
                HCGB_time.print_time_stamp(filename_stamp)
        
        ## merge results
        agrvate_results = pd.concat([agrvate_results, info_sample], join='outer')
        
    print ("+ Jobs finished\n+ Collecting information for all samples...")
    os.chdir(path_here)
    
    ## debug messages
    if debug:
        HCGB_aes.debug_message('agrvate_results', 'yellow')
        HCGB_main.print_all_pandaDF(agrvate_results)

    return(agrvate_results, info_dict)

##############################
def agrvate_call(sample, assembly_file, folder, debug=False):
    """agrvate call and check results."""
    
    ## prepare call
    log_call = os.path.join(folder, "agrvate_cmd.log")
    err_call = os.path.join(folder, "agrvate_cmd.err")
    agrvate_bin = set_config.get_exe('agrvate')
    
    ## system call
    cmd_call = "%s -i %s -m -f >  %s 2> %s " %(agrvate_bin, 
                                               assembly_file,
                                               log_call, err_call) ## use mummer (-m) and force results folder (-f)
    status = HCGB_sys.system_call(cmd_call)
    
    if status:
        res = get_results_agrvate(assembly_file, folder, sample, debug)
        return (res)
    else:
        return(False)
    
#########################################
def get_results_agrvate(assembly_file, folder, sample, debug=False):
    ## check results
    ## see https://github.com/VishnuRaghuram94/AgrVATE#results for additional details
    results = pd.DataFrame()
    
    ## check folder is created
    assembly_file_name = os.path.basename(assembly_file).split('.fna')[0]    
    original_results_folder = os.path.join(folder, assembly_file_name + '-results')
    results_folder = os.path.join(folder, 'agrvate_results')
    
    ## rename folder
    if os.path.isdir(original_results_folder):
        print("+ Results folder generated OK")
        print("+ Check results generated:")
        
        ## rename folder
        os.rename(original_results_folder, results_folder)
        os.rename(os.path.join(folder, assembly_file_name + '.fna-error-report.tab'), os.path.join(results_folder, 'error_report.tab'))
    
    ## get results
    if (os.path.isdir(results_folder)):
        ## write to excel1
        file_name_Excel = os.path.join(folder, sample + '_agr_results.xlsx')
        writer_Excel = pd.ExcelWriter(file_name_Excel, engine='xlsxwriter') ## open excel handle
    
        ## get all files
        list_files = HCGB_main.get_fullpath_list(results_folder)
    
        ## summary tab
        summary_tab_file = [s for s in list_files if s.endswith("summary.tab")][0]
        summary_tab =  HCGB_main.get_data(summary_tab_file, '\t', options="")
        summary_tab['sample'] = sample
        
        ## columns
        #agr_group: gp1/gp2/gp3/gp4. 'u' means unknown. 
        ##           If multiple agr groups were found (col 5 = m), 
        ##           the displayed agr group is the majority/highest confidence. 
        # match_score: maximum 15; 0 means untypeable; < 5 means low confidence.
        # canonical_agrD: 1 means canonical; 0 means non-canonical; u means unknown.
        # multiple_agr:  s means single, m means multiple, u means unknown ) 
        ##               Multiple groups are found likely due to multiple S. aureus isolates in sequence
        # frameshifts: Number found in CDS of extracted agr operon ('u' if agr operon not extracted)
        
        ## debug messages
        if debug:
            HCGB_aes.debug_message("agrvate results: Summary tab file", 'yellow')
            print(summary_tab_file)
            print(summary_tab)

        ## add summary results to all results
        del summary_tab['#filename']
        results = summary_tab.copy()

        ## save summary_tab into excel
        ## tab summary
        summary_tab.to_excel(writer_Excel, sheet_name='summary') ## write excel handle

        ## agr_gp tab
        agr_gp_tab_file = [s for s in list_files if s.endswith("agr_gp.tab")][0]
        if HCGB_files.is_non_zero_file(agr_gp_tab_file):
            agr_gp_tab =  HCGB_main.get_data(agr_gp_tab_file, '\t', options='header=None')
            agr_gp_tab.columns = ['contig', 'agr', 'evalue', 'identity', 'start', 'end']
            agr_gp_tab['sample'] = sample
            
            ## columns
            ## Assembly Contig ID
            ## ID of matched agr group kmer
            ## evalue
            ## Percentage identity of match
            ## Start position of kmer alignment on input sequence
            ## End position of kmer alignment on input sequence
    
            ## debug messages
            if debug:
                HCGB_aes.debug_message("agrvate results: agr_gp file", 'yellow')
                print(agr_gp_tab_file)
                print(agr_gp_tab)
            
            ## save agr_gp_tab file into excel
            ## tab operon
            agr_gp_tab.to_excel(writer_Excel, sheet_name='operon') ## write excel handle

        ## agr_operon fna
        try:
            agr_operon_fna_file = [s for s in list_files if s.endswith("agr_operon.fna")][0]
            ## debug messages
            if debug:
                HCGB_aes.debug_message("agrvate results: agr_operon file", 'yellow')
                print(agr_operon_fna_file)
            
            results['operon_fna'] = agr_operon_fna_file
        except:
            results['operon_fna'] = ''

        ## agr_operon fna
        error_report_file = [s for s in list_files if s.endswith("error_report.tab")][0]
        error_report =  HCGB_main.get_data(error_report_file, '\t', options="")
        del error_report['#input_name']

        ## debug messages
        if debug:
            HCGB_aes.debug_message("agrvate results: error_report.tab file", 'yellow')
            print(error_report_file)
            print(error_report)
            
        ## save error_report file into excel
        ## tab steps
        error_report.to_excel(writer_Excel, sheet_name='steps') ## write excel handle
        
        ## merge results
        results = pd.concat([results, error_report], axis=1)

        ## close xlsx file
        writer_Excel.save() ## close excel handle
    
        ## add to pandas dataframe
        results['agr_operon_xlsx'] = file_name_Excel

    ## debug messages
    if debug:
        HCGB_aes.debug_message("agrvate results", 'yellow')
        HCGB_main.print_all_pandaDF(results)
        
    return (results)

##############
def help_options():
    print ("\nUSAGE: python %s name assembly_fasta folder...\n"  %os.path.realpath(__file__))

##############
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
    ## ATTENTION: agrvate needs to chdir to output folder
    os.chdir(folder)

    ###
    agrvate_call(name, fasta_file, folder, debug)
    

'''******************************************'''
if __name__== "__main__":
    main()

