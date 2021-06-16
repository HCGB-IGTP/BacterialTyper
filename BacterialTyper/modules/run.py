#!/usr/bin/env python3
###########################################################
## Jose F. Sanchez										 ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	 ##
###########################################################
"""
This module generates a uniq call of BacterialTyper to 
any other module. It automates the analysis and process.
"""
## useful imports
import time
import io
import os
import re
import sys
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.modules import prep
from BacterialTyper.modules import qc
from BacterialTyper.modules import ident
from BacterialTyper.modules import profile
from BacterialTyper.modules import database
from BacterialTyper.modules import MGE
from BacterialTyper.modules import cluster
from BacterialTyper.modules import assemble
from BacterialTyper.modules import annot
from BacterialTyper.modules import database
from BacterialTyper.modules import trimm
from BacterialTyper.modules import report_generation
from BacterialTyper.modules import help_info
from BacterialTyper import __version__ as pipeline_version

from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.scripts import annotation
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import min_hash_caller
from BacterialTyper.scripts import trimmomatic_call
from BacterialTyper.scripts import species_identification_KMA
from BacterialTyper.scripts import MLSTar
from BacterialTyper.scripts import bacteriophage

from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files

####################################
def run_BacterialTyper(options):

	## init time
	start_time_total = time.time()

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		help_info.help_fastq_format()
		exit()

	elif (options.help_BUSCO):
		## information for BUSCO
		BUSCO_caller.print_help_BUSCO()
		exit()

	elif (options.help_project):
		## information for project
		help_info.project_help()
		exit()

	elif (options.help_multiqc):
		## information for Multiqc
		multiQC_report.multiqc_help()

	elif (options.help_Prokka):
		## information for Prokka
		annotation.print_list_prokka()
		exit()
	
	elif (options.help_Mash):
		## information for Min Hash Software
		min_hash_caller.helpMash()
		exit()
		
	elif (options.help_ARIBA):
		## information for ARIBA
		ariba_caller.help_ARIBA()
		exit()
        
	elif (options.help_trimm_adapters):
		## help on trimm adapters
		trimmomatic_call.print_help_adapters()
		exit()

	elif (options.help_KMA):
		## information for KMA Software
		species_identification_KMA.help_kma_database()
		exit()

	elif (options.help_MLSTar):
		## information for KMA Software
		MLSTar.help_MLSTar()
		exit()
	
	elif (options.help_PhiSpy):
		## information for PhiSpy software
		bacteriophage.help_PhiSpy()
		exit()

	elif (options.help_MGE_analysis):
		## information for MGE module analysis
		MGE.help_MGE_analysis()
		exit()

	elif (options.help_input_MGE):
		## information for PhiSpy
		MGE.help_input_MGE()
		exit()

	### 
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("BacterialTyper analysis")

	print ("--------- Starting Process ---------")
	HCGB_time.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	## set mode: project/detached
	if (options.project):
		outdir = input_dir		
	elif (options.detached):
		outdir = os.path.abspath(options.output_folder)
