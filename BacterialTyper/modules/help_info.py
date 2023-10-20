#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
'''
Prints help messages for several modules and options.
'''
## useful imports
from termcolor import colored
    
## import my modules
from BacterialTyper.scripts import bacteriophage
from BacterialTyper.modules import MGE

from BacterialTyper.scripts import variant_calling
from BacterialTyper.scripts import genomic_island
from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.scripts import annotation
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import min_hash_caller
from BacterialTyper.scripts import trimmomatic_call
from BacterialTyper.scripts import KMA_caller
from BacterialTyper.scripts import Kraken2_caller
from BacterialTyper.scripts import MLSTar
from BacterialTyper.report.Staphylococcus import get_spa_typing
from BacterialTyper import __version__ as pipeline_version

import HCGB.functions.aesthetics_functions as HCGB_aes

##########################
def help_info(options):
    """
    Main function to control all help messages requested from any modules.
    
    Given the option provided via options dictionary flags, it calls the specific information in each script.
    
    """
    
    ## project help
    try:
        if (options.help_project):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            project_help()
            return(1)
    except:
        pass

    ## help_format option
    try:
        if (options.help_format):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            help_fastq_format()        
            return(1)
    except:
        pass
    
    ## information for Prokka    
    try:
        if (options.help_Prokka):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            annotation.print_list_prokka()
            return(1)
    except:
        pass
    
    ## information for BUSCO databases    
    try:
        if (options.help_BUSCO):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            BUSCO_caller.print_help_BUSCO()        
            return(1)
    except:
        pass

    ## information for ARIBA databases
    try:
        if (options.help_ARIBA):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            print ("ARIBA databases information:")    
            ariba_caller.help_ARIBA()        
            return(1)
    except:
        pass

    ## information for trimm adapters
    try:
        if (options.help_trimm_adapters):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            trimmomatic_call.print_help_adapters()        
            return(1)
    except:
        pass

    ## information for Multiqc
    try:
        if (options.help_multiqc):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            multiQC_report.multiqc_help()        
            return(1)
    except:
        pass

    ## information for KMA Software
    try:
        if (options.help_KMA):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            KMA_caller.help_kma_database()        
            return(1)
    except:
        pass

    ## information for PhiSpy
    try:
       if (options.help_PhiSpy):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            bacteriophage.help_PhiSpy()        
            return(1)
    except:
        pass
            
    ## information for MGE analysis
    try:
       if (options.help_MGE_analysis):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            MGE.help_MGE_analysis()
            return(1)
    except:
        pass

    ## information for MGE module
    try:
        if (options.help_input_MGE):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            MGE.help_input_MGE()
            return(1)
    except:
        pass
                
    ## information for MLSTar Software
    try:
        if (options.help_MLSTar):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            MLSTar.help_MLSTar()        
            return(1)
    except:
        pass
    
    ## information for Min Hash Software
    try:
        if (options.help_Mash):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            min_hash_caller.helpMash()        
            return(1)
    except:
        pass
    
    ## information for Snippy
    try:
        if (options.help_Snippy):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            variant_calling.help_Snippy()        
            return(1)
    except:
        pass

    ## information for Dimob
    try:
        if (options.help_Dimob):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            genomic_island.help_Dimob()
            return(1)
    except:
        pass

    ## Information for spatyper
    try:
        if (options.help_spaTyper):
            HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
            ## help_format option
            get_spa_typing.help_spaTyper()
            return(1)
    except:
        pass

    return(0)

##########################
def project_help():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def help_fastq_format():
    """
    Explanation of fastq format details.
    """

    HCGB_aes.boxymcboxface("Name format for samples")

    print ("Format for fastq files can be:")
    print ("name.fastq.gz, name_1.fastq.gz, name_R2.fastq.gz, name_L001_R1.fastq.gz, name_L001_R1_001.fastq.gz etc.")
    print ("\nThere are many options and here we provide some guidelines on the name format.")
    print ("\n")

    HCGB_aes.print_sepLine("*",20,"red")
    print ("[1] Length limitation")
    HCGB_aes.print_sepLine("*",20,"red")
    print ("There is a limitation for the sample ID ('name') of 25 characters.")
    print (colored("** BacterialTyper provides an option to rename samples if necessary: module prep option --rename **", 'yellow'))
    print ("\n")

    HCGB_aes.print_sepLine("*",20,"red")
    print ("[2] Single end files")
    HCGB_aes.print_sepLine("*",20,"red")
    print (colored('** Use option --single-end in the different BacterialTyper modules **', 'yellow'))
    print ("name.fastq.gz")
    print ("name.fastq")
    print ("name.fq")
    print ("\n")

    HCGB_aes.print_sepLine("*",20,"red")
    print ("[3] Paired-end files")
    HCGB_aes.print_sepLine("*",20,"red")
    print ("Paired-end files are full supported. The format for these files are:")
    print ("Read1 => name_1.fastq.g or name_R1.fastq.gz")
    print ("Read2 => name_2.fastq.gz or name_R2.fastq.gz")
    print (colored('** See additional details for Lane information **', 'yellow'))
    print ("\n")

    HCGB_aes.print_sepLine("*",55,"red")
    print ("[4] Lane information:")
    HCGB_aes.print_sepLine("*",55,"red")
    print ("In some cases, files might contain lane information (*L00x* and/or *00x*).")
    print ("BacterialTyper supports these names as long as follow these examples:")
    print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_1.fastq.gz\tname_L00x_2.fastq.gz")
    print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
    print ("name_L00x_R1_00x.fastq.gz\tname_L00x_R2_00x.fastq.gz")
    print ("\n")
    print ("Sometimes it might be appropriate to include lane tags (*L00X*) within the name.")
    print (colored("** Use option --include-lane within each module", 'yellow'))
    
    print (colored("\n** If you need to merge fastq files from different lanes, use option within module prep **", 'yellow'))
    print("As an example:")
    print (colored("\n** Option --merge within module prep **", 'yellow'))
    print ("sample1_L001_R1.fastq.gz\tsample1_L001_R2.fastq.gz")
    print ("sample1_L002_R1.fastq.gz\tsample1_L002_R2.fastq.gz")
    print ("sample1_L003_R1.fastq.gz\tsample1_L003_R2.fastq.gz")
    print ("sample1_L004_R1.fastq.gz\tsample1_L004_R2.fastq.gz")
    print ("Result:")
    print ("--------------------------------------------------")
    print ("sample1_R1.fastq.gz\tsample1_R2.fastq.gz")
    print ("\n")
    print (colored("\n** Option --merge-by-lane within module prep **", 'yellow'))
    print ("sample1_L001_R1_001.fastq.gz\tsample1_L001_R2_001.fastq.gz")
    print ("sample1_L001_R1_002.fastq.gz\tsample1_L001_R2_002.fastq.gz")
    print ("sample1_L002_R1_001.fastq.gz\tsample1_L002_R2_001.fastq.gz")
    print ("sample1_L002_R1_002.fastq.gz\tsample1_L002_R2_002.fastq.gz")
    print ("--------------------------------------------------")
    print ("Result:")
    print ("sample1_L001_R1.fastq.gz\tsample1_L001_R2.fastq.gz")
    print ("sample1_L002_R1.fastq.gz\tsample1_L002_R2.fastq.gz")
    print (colored("** Remember to use option --include_lane within each module", 'yellow'))
    print ("\n")
    
    HCGB_aes.print_sepLine("*",55,"red")
    print ("[5] Include all information:")
    HCGB_aes.print_sepLine("*",55,"red")
    print ("In some cases, files might contain other information and it is necessay to " +
           "include it all as a tag nane. See as an example:")
    print ("sample1_L001_XYZ_R1_001.fastq.gz\tsample1_L001_XYZ_R2_001.fastq.gz")
    print (colored("** Remember to use option --include_all within each module", 'yellow'))
    print (colored("** It might be appropiate to change samples names using --rename option under prep module", 'yellow'))
    
    print ("\n")
    HCGB_aes.print_sepLine("*",15,"red")
    print ("[6] Extensions:")
    HCGB_aes.print_sepLine("*",15,"red")
    print ("name_L00x_R2.fastq\tname_L00x_R2.fq\nname_L00x_R2.fastq.gz\tname_L00x_R2.fq.gz")
    print ("\n")

