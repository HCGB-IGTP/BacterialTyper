#!/usr/bin/env python3
'''
This code calls several modules and print help messages.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
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
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import annotation
from BacterialTyper import BUSCO_caller
from BacterialTyper import ariba_caller
from BacterialTyper import sampleParser
from BacterialTyper import bacteriophage
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report

##########################
def run(options):

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

	## information for PhiSpy
	if (options.help_PhiSpy):
		bacteriophage.help_PhiSpy()
		exit()

	## information for trimm adapters
	if (options.help_trimm_adapters):
		trimmomatic_call.print_help_adapters()
		exit()

	## information for Multiqc
	if (options.help_multiqc):
		multiQC_report.multiqc_help()
		exit()

	
##########################
def project_help():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
	
