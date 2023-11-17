#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:35:24 2023

@author: jsanchez
"""
import os

from termcolor import colored
from BacterialTyper.config import set_config
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_file

###############
def help_amrfinder():
    """
    Print help message information for AMRfinder
    """
    
    print (colored("\n\n***** AMRfinder help message *****\n\n", 'yellow'))
    print("""Identify AMR and virulence genes in proteins and/or contigs and print a report\nDOCUMENTATION: See https://github.com/ncbi/amr/wiki for full documentation""")

###################################
def organisms_amrfinder(db_folder, Debug=False):
    """
    Provide information as generated from AMRfinder for available organisms

    Parameters
    ----------
    db_folder : str
        Absolute path to latest/ AMRfinder folder
    Debug : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    PRints information.

    """    
    
    print (colored("\n\n***** Get the list of organisms available from NCBI AMRFinderplus software *****\n\n", 'yellow'))
    
    amrfinder_bin = set_config.get_exe("amrfinder", Debug=Debug)
    HCGB_sys.system_call(amrfinder_bin + ' --list_organisms --database ' + db_folder + " > tmp.file")
    
    print("")
    with open('tmp.file','r') as reader:
        contents_file = reader.readline()
        while contents_file != "":
            if contents_file == "\n":
                pass

            contents_file = reader.readline()
            contents_file = contents_file.replace(":", ",")
            print(contents_file.replace(",", "\n"))
            

    print ("Read additional information in: https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option")
    
###################
def amrfinder_version_db(db_folder, Debug=False):
    """
    Provide information for AMRfinder version

    Parameters
    ----------
    db_folder : str
        Absolute path to latest/ AMRfinder folder
    Debug : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    Prints version information.

    """
    amrfinder_bin = set_config.get_exe("amrfinder", Debug=Debug)
    HCGB_sys.system_call(amrfinder_bin + '--database ' + db_folder + ' --database_version', Debug=Debug)

####################################
def download_database(db_folder, Debug=False):
    """
    

    Parameters
    ----------
    db_folder : TYPE
        DESCRIPTION.
    Debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    
    ## TODO
    ## check time stamp and rerun if many days since last time
    
    
    time_Stamp = HCGB_time.create_human_timestamp()
    
    amrfinder_bin = set_config.get_exe("amrfinder_update", Debug=Debug)
    log_file = os.path.join(db_folder, time_Stamp + "_update.log")
    ## After updating it might require to index
    cmd_update = "%s --database %s --threads 2 --log %s" %(amrfinder_bin, db_folder, log_file)
    code2return = HCGB_sys.system_call(cmd_update)
    
    return(code2return, time_Stamp)

###################
def amrfinder_caller(sample_name, protein_file, gff_file, nuc_file, threads_num, db_fold, outfolder, others, species_ident = "", Debug=False):
    """
    Call AMRfinder for each sample

    Parameters
    ----------
    sample_name : TYPE
        DESCRIPTION.
    protein_file : TYPE
        DESCRIPTION.
    gff_file : TYPE
        DESCRIPTION.
    nuc_file : TYPE
        DESCRIPTION.
    threads_num : TYPE
        DESCRIPTION.
    db_fold : TYPE
        DESCRIPTION.
    outfolder : TYPE
        DESCRIPTION.
    others : TYPE
        DESCRIPTION.
    Debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    
    if Debug:
        HCGB_aes.debug_message("Calling amrfinder_caller", color='yellow')
    
    
    ## output files
    outfolder = HCGB_file.create_folder(outfolder)
    out_file_log = os.path.join(outfolder, sample_name + '.log')
    outfile_res =  os.path.join(outfolder, sample_name + '.tsv')
    
    ## create string call
    amrfinder_bin = set_config.get_exe("amrfinder", Debug=Debug)
    cmd_amrfinder = "%s --name %s --database %s --protein %s" %(amrfinder_bin, sample_name, db_fold, protein_file)
    
    ## add nucleotide sequence if provided
    if HCGB_file.is_non_zero_file(nuc_file):
        cmd_amrfinder += " --nucleotide " + nuc_file
    
    ## add annotation and format if provided
    if HCGB_file.is_non_zero_file(gff_file):
        cmd_amrfinder += " --gff %s --annotation_format prokka" %gff_file
    
    ## add species name if available
    if species_ident:
        cmd_amrfinder += " --report_common --organism " + species_ident
    else:
        if Debug:
            HCGB_aes.debug_message("No species name option provided:", color='yellow')
    
    ## add others options if provided
    if others:
        cmd_amrfinder += " " + others
    else:
        if Debug:
            HCGB_aes.debug_message("No Other options provided:", color='yellow')
            print("Default values for other options will be provided")
        
    ## finish call
    cmd_amrfinder += " --threads %s --print_node --plus --log %s > %s" %(str(threads_num), out_file_log, outfile_res)
    
    if Debug:
        HCGB_aes.debug_message("call:", color='yellow')
        print(cmd_amrfinder)
    
    code = HCGB_sys.system_call(cmd_amrfinder)

    if (code == 'OK'):
        ## success stamps
        filename_stamp = os.path.join(outfolder, '.success')
        HCGB_time.print_time_stamp(filename_stamp)
        return('OK')
    else:
        return('FAIL')

###################
def parse_organism_provided(org_given, db_folder, Debug=False):
    
    ## string generate using amrfinder -l option
    species_Available = "Acinetobacter_baumannii,Burkholderia_cepacia,Burkholderia_pseudomallei"
    species_Available += ",Citrobacter_freundii,Clostridioides_difficile,Enterobacter_asburiae,Enterobacter_cloacae"
    species_Available += ",Enterococcus_faecalis,Enterococcus_faecium,Klebsiella_oxytoca,Klebsiella_pneumoniae"
    species_Available += ",Neisseria_gonorrhoeae,Neisseria_meningitidis,Pseudomonas_aeruginosa,Serratia_marcescens"
    species_Available += ",Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae"
    species_Available += ",Streptococcus_pneumoniae,Streptococcus_pyogenes,Vibrio_cholerae"

    species_list = species_Available.split(",")
    exception_list = ["Campylobacter" , "Escherichia", "Salmonella"]
    
    if Debug:
        print(species_list)
        print(exception_list)
    
    if org_given.replace(" ","_") in species_list:
        if Debug:
            print("Species provided is available!")
        
        return(org_given.replace(" ","_"))
    
    if org_given.split(" ")[0] in exception_list:
        if Debug:
            print("Species provided is available at the genus level!")
        
        return(org_given.split(" ")[0])

    
    ## None
    if Debug:
        print("Species provided is not available at the species or genus level!")
    return("")    
        
##################################################
def main():
    
    import sys

        
    if len(sys.argv)<7:
        print("Usage: ")
        print("python " + sys.argv[0] + " prot_file gff_file nuc_file db outfile species_ident" )
        exit()
    
    parse_organism_provided(sys.argv[6], sys.argv[4], Debug=True)
    exit()
    
    amrfinder_caller(sample_name="test", 
                     protein_file= sys.argv[1], 
                     gff_file= sys.argv[2],   
                     nuc_file= sys.argv[3], 
                     threads_num=2,  
                     db_fold= sys.argv[4], 
                     outfolder = sys.argv[5], 
                     others="", 
                     species_ident="Klebsiella_pneumoniae",
                     Debug=True)

##################################################
if __name__== "__main__":
    main()
