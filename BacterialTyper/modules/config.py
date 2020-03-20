#!/usr/bin/env python3
#########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Sets configuration of the pipeline.
"""
## useful imports
import time
import io
import os
import sys
from termcolor import colored
from distutils.version import LooseVersion

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.config import extern_progs
from BacterialTyper.scripts import functions
from BacterialTyper.config import install_dependencies

## [TODO]
def run(options):

	## init time
	start_time_total = time.time()

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False
	
	functions.pipeline_header()
	functions.boxymcboxface("Pipeline Configuration")
	print ("--------- Starting Process ---------")
	functions.print_time()

	## install or only check
	option_install = False
	if (options.option == 'install'):
		print ("\nCheck dependencies, modules or third party software...")
		print ("\nTry to install all missing dependencies, modules or third party software...")
		option_install = True
	elif (options.option == 'only_check'):
		print ("\nCheck dependencies, modules or third party software and print report...")
	
	functions.print_sepLine("+", 50, False)
	print ('Python:')
	
	## python version
	this_python_version = str(sys.version)
	python_min_version = extern_progs.return_min_version('python')
	if LooseVersion(this_python_version) >= LooseVersion(python_min_version):
		print ("Minimun version (%s) satistied: %s" %( python_min_version, this_python_version))

	functions.print_sepLine("+", 50, False)
	
	## python packages
	print ('\n')
	functions.print_sepLine("+", 50, False)
	print ('Python packages:')
	extern_progs.check_python_packages(Debug, option_install)
	functions.print_sepLine("+", 50, False)
	print ('\n')
	
	functions.print_sepLine("+", 50, False)
	print ('External dependencies:\n')
	functions.print_sepLine("+", 50, False)
	extern_progs.print_dependencies()
	print ('\n')

