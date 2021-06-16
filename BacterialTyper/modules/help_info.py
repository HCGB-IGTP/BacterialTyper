#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
import BacterialTyper
'''
Prints help messages for several modules and options.
'''
## useful imports
import time
import io
import os
import re
import sys
from io import open
from termcolor import colored
	
## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.scripts import annotation
from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import bacteriophage
from BacterialTyper.scripts import trimmomatic_call
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.scripts import MLSTar
from BacterialTyper.modules import MGE
from BacterialTyper.scripts import min_hash_caller
from BacterialTyper.scripts import variant_calling
from BacterialTyper.scripts import genomic_island
from BacterialTyper import __version__ as pipeline_version

import HCGB.functions.aesthetics_functions as HCGB_aes

##########################
def run_info(options):

	## project help
	if (options.help_project):
		project_help()
		exit()

	## help_format option
	if (options.help_format):
		help_fastq_format()
		exit()

	## information for Prokka	
	if (options.help_Prokka):
		annotation.print_list_prokka()
		exit()
	
	## information for BUSCO databases	
	if (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
		exit()

	## information for ARIBA databases
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()

	## information for trimm adapters
	if (options.help_trimm_adapters):
		trimmomatic_call.print_help_adapters()
		exit()

	## information for Multiqc
	if (options.help_multiqc):
		multiQC_report.multiqc_help()
		exit()

	## information for KMA Software
	if (options.help_KMA):
		species_identification_KMA.help_kma_database()
		exit()

	## information for PhiSpy
	if (options.help_PhiSpy):
		bacteriophage.help_PhiSpy()
		exit()
		
	## information for MGE analysis
	if (options.help_MGE_analysis):
		MGE.help_MGE_analysis()
		exit()

	## information for MGE module
	if (options.help_input_MGE):
		MGE.help_input_MGE()
		exit()
		
	## information for MLSTar Software
	if (options.help_MLSTar):
		MLSTar.help_MLSTar()
		exit()
	
	## information for Min Hash Software
	if (options.help_Mash):
		min_hash_caller.helpMash()
		exit()
	
	## information for Snippy
	if (options.help_Snippy):
		variant_calling.help_Snippy()
		exit()

	## information for Dimob
	if (options.help_Dimob):
		genomic_island.help_Dimob()
		exit()
	
	
##########################
def project_help():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

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

