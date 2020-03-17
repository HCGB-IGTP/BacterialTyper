#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
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
from BacterialTyper.scripts import functions
from BacterialTyper.scripts import config
from BacterialTyper.scripts import annotation
from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import sampleParser
from BacterialTyper.scripts import bacteriophage
from BacterialTyper.scripts import trimmomatic_call
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.modules import MGE
from BacterialTyper.scripts import min_hash_caller

##########################
def run_info(options):

	## project help
	if (options.help_project):
		project_help()
		exit()

	## help_format option
	if (options.help_format):
		sampleParser.help_format()
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
		
	
##########################
def project_help():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
	
