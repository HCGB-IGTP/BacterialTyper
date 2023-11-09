#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 11:52:43 2023

@author: jsanchez
"""
import os
from termcolor import colored
from BacterialTyper.config import set_config
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_file
import pandas as pd


##################################################
def get_dbs(db_name, fold, Debug=False):
    ## TODO
    print(colored("***** TODO: Create a call to retrieve kraken2 dbs ***", 'red'))
    raise SystemExit()
    
##################################################
def help_kraken():
    ## [TODO]
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##################################################
def kraken_caller(sample_name, read_files, threads_num, db_fold, outfolder, hit_groups=3, others="", Debug=False):
    """
    
    Parameters
    ----------
    sample_name : TYPE
        DESCRIPTION.
    read_files : TYPE
        DESCRIPTION.
    threads_num : TYPE
        DESCRIPTION.
    db_fold : TYPE
        DESCRIPTION.
    outfolder : TYPE
        DESCRIPTION.
    hit_groups : TYPE, optional
        DESCRIPTION. The default is 3.
    others : TYPE, optional
        DESCRIPTION. The default is "".
    Debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """    
    
    ##  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9725748/
    ## The --minimum-hit-groups flag specifies the minimum number of ”hit groups” 
    ## needed to make a classification call. Hit groups are overlapping k-mers sharing 
    ## the same minimizer. Kraken 2 uses minimizers to compress the input genomic 
    ## sequences, thereby reducing storage memory needed and run time. In this 
    ## example we increase the minimum number of hit groups from the default 2 
    ## groups to 3 groups for increased accuracy. Lastly, the --report-minimizer-data 
    ## flag reports minimizer and distinct minimizer count information in addition 
    ## to the normal Kraken 2 report.
    
    filename_stamp = os.path.join(outfolder, '.success_kraken')
    if HCGB_file.is_non_zero_file(filename_stamp):
        return True
    else:
     
        if Debug:
            HCGB_aes.debug_message("Calling kraken_caller", color='yellow')
        
        ## output files
        outfolder = HCGB_file.create_folder(outfolder)
        out_file_log = os.path.join(outfolder, sample_name + '.klog')
        outfile_res =  os.path.join(outfolder, sample_name + '.kraken2')
        outfile_report =  os.path.join(outfolder, sample_name + '.kreport')
        
        ## create string call
        kraken_bin = set_config.get_exe("kraken2")
        cmd_kraken = "%s --db %s --memory-mapping --threads %s" %(
            kraken_bin, db_fold, str(threads_num))
    
        ## output    
        cmd_kraken += " --report-minimizer-data --report %s --output %s"%(outfile_report, outfile_res)
        
        ## parameters
        cmd_kraken += " --minimum-hit-groups " + str(hit_groups)
        
        ## add others options if provided
        if others:
            cmd_kraken += " " + others
        else:
            if Debug:
                HCGB_aes.debug_message("No Other options provided:", color='yellow')
                print("Default values for other options will be provided")
        
        ## finish call
        
        ## input files
        read_files_string = " ".join(read_files)
        if len(read_files)>1:
            cmd_kraken += " --paired"
        cmd_kraken += " %s 2> %s " %(read_files_string, out_file_log)
        
        if Debug:
            HCGB_aes.debug_message("call:", color='yellow')
            print(cmd_kraken)
        
        code = HCGB_sys.system_call(cmd_kraken)
        if (code == 'OK'):
            ## success stamps
            HCGB_time.print_time_stamp(filename_stamp)
            return True
        else:
            return False


##################################################
def bracken_caller(sample_name, outfolder, db_fold, level_abundance="S", thres_count = 50, Debug=False):
    """
    

    Parameters
    ----------
    sample_name : TYPE
        DESCRIPTION.
    outfolder : TYPE
        DESCRIPTION.
    db_fold : TYPE
        DESCRIPTION.
    level_abundance : TYPE, optional
        DESCRIPTION. The default is "S".
    thres_count : TYPE, optional
        DESCRIPTION. The default is 50.
    Debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    
    filename_stamp = os.path.join(outfolder, '.success_bracken')
    if HCGB_file.is_non_zero_file(filename_stamp):
       return True
    else:
    
        if Debug:
            HCGB_aes.debug_message("Calling bracken_caller", color='yellow')
        
        ## output files
        outfolder = HCGB_file.create_folder(outfolder)
        out_file_log = os.path.join(outfolder, sample_name + '.blog')
        kraken_res =  os.path.join(outfolder, sample_name + '.kreport')
        
        outfile_report =  os.path.join(outfolder, sample_name + '.breport')
        outfile_res =  os.path.join(outfolder, sample_name + '.bracken')
        
        ## create string call
        bracken_bin = set_config.get_exe("bracken")
        cmd_bracken = "%s -d %s -i %s -w %s -o %s" %(
            bracken_bin, db_fold, kraken_res, outfile_report, outfile_res)
    
        
        cmd_bracken += " -l %s -t %s" %(level_abundance, str(thres_count))
        
        ## finish call
        cmd_bracken += " 2> %s " %(out_file_log)
        
        if Debug:
            HCGB_aes.debug_message("call:", color='yellow')
            print(cmd_bracken)
        
        code = HCGB_sys.system_call(cmd_bracken)
        if (code == 'OK'):
            ## success stamps
            HCGB_time.print_time_stamp(filename_stamp)
            return True
        else:
            return False

##################################################
def kraken_run_all(sample_name, read_files, threads_num, db_fold, outfolder, 
                   level_abundance="S", thres_count = 50, hit_groups=3, others="", debug=False):
    
    
    
    outfolder_here = HCGB_file.create_subfolder('kraken', outfolder)
    filename_stamp = os.path.join(outfolder_here, '.success')
    
    if HCGB_file.is_non_zero_file(filename_stamp):
        stamp =    HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous kraken2 command generated results on: %s [Files: %s]" %(stamp, sample_name), 'green'))
    else:
        if debug:
            print(" *** Call kraken2 for sample %s ***" %sample_name)
        
        codeK = kraken_caller(sample_name, read_files, threads_num, db_fold, outfolder_here, 
                     hit_groups = hit_groups, others=others, Debug=debug)
        if debug:
            print(" *** Call bracken for sample %s ***" %sample_name)
        
        codeB = bracken_caller(sample_name, outfolder_here, db_fold, level_abundance, 
                       thres_count = thres_count, Debug= debug)
        
        if codeB & codeK:
            ## success stamps
            HCGB_time.print_time_stamp(filename_stamp)
            return('OK')
        else:
            return('FAIL')
        
##################################################
def parse_results(bracken_res, folder, cut_off=0.90):
    
    ## read and write in excel
    bracken_df = pd.read_csv(bracken_res, sep="\t")
    bracken_df.to_excel(bracken_res + '.xlsx')
    
    filt_db = bracken_df[bracken_df['fraction_total_reads']>cut_off]
    sp_put = filt_db.iat[0,0]
    
    with open(os.path.join(folder, ".species"), 'w') as opener:
        opener.write(sp_put)
    
    return(sp_put)
    


##################################################
def main():
    
    import sys
    
    if len(sys.argv)<6:
        print("\n## Usage: ")
        print("python " + sys.argv[0] + " sample_name db_fold outfolder read_file(s)\n\n" )
        exit()
    
    print(" *** Call kraken2 ***")
    kraken_caller(sample_name=sys.argv[1], 
                  read_files = sys.argv[4:], 
                  threads_num = 2, 
                  db_fold = sys.argv[2], 
                  outfolder = sys.argv[3], Debug=True)
    
    print("\n\n *** Call bracken ***")
    bracken_caller(sample_name=sys.argv[1],
                   outfolder = sys.argv[3],
                   db_fold = sys.argv[2], Debug=True)

##################################################
if __name__== "__main__":
    main()