#!/usr/bin/env python3
#########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Sets configuration of the pipeline.

This module relies in several scripts and files which are arranged in the 
``BacterialTyper/config`` :doc:`folder <../config/config_index>`:	

	- :doc:`BacterialTyper/config/set_config.py <../config/set_config>`

	- :doc:`BacterialTyper/config/extern_progs.py <../config/extern_progs>`

	- :doc:`BacterialTyper/config/install_dependencies.py <../config/install_dependencies>`

	- Additional information and specifications :doc:`details <../config/additional_info>`.	

.. seealso:: Additional information on BacterialTyper configuration and requirements

	- :doc:`Installation <../../user_guide/installation/installing>`
	
	- :doc:`Requirements <../../user_guide/installation/requirements>` 

	- :doc:`Configuration <../../user_guide/modules/config>`

"""
## useful imports
import time
import io
import os
import sys
from termcolor import colored
from distutils.version import LooseVersion

## import my modules
from BacterialTyper.config import extern_progs
from BacterialTyper.config import install_dependencies
from BacterialTyper.config import set_config
from BacterialTyper import __version__ as pipeline_version

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.files_functions as HCGB_files

###################################################################
def run(options):
	"""
	This is the main function of the module ``config``. It basically checks 
	if the different requirements (python` and third-party software) are
	fulfilled. 

	If any requirement is not available this modules tries to install them or reports to the user to
	manually install them.

	:param option: State whether to check or install missing modules, packages and third party software. Provide: check/install
	:param install_path: Absolute path to install modules or packages missing. Default: ``BacterialTyper`` environment path.
	:param IslandPath: True/False for checking additional perl and software required by this option analysis.
	:param debug: True/false for debugging messages.
	
	:type option: string 
	:type IslandPath: boolean
	:type install_path: string 
	:type debug: boolean	

	.. seealso:: This function depends on several ``BacterialTyper`` functions:

		- :func:`BacterialTyper.config.set_config.check_python_packages`

		- :func:`BacterialTyper.config.set_config.check_perl_packages`

		- :func:`BacterialTyper.config.extern_progs.return_min_version_soft`

		- :func:`BacterialTyper.config.extern_progs.print_dependencies`

	"""

	## init time
	start_time_total = time.time()

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False

	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Pipeline Configuration")
	print ("--------- Starting Process ---------")
	HCGB_time.print_time()

	if (options.install_path):
		if os.path.isdir(options.install_path):
			if (Debug):
				print ("Installation path provided for missing modules, packages, dependencies...")
				print ("Path: " + options.install_path)
		else:
			print (colored("\n*** ERROR ****", 'red'))
			print (colored("Path provided is not a folder", 'red')) 
			print (options.install_path)
			exit()
	else:
		## get python environment path
		env_bin_directory = os.path.dirname(os.environ['_'])

		##os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'templates'))
		options.install_path = os.path.abspath(os.path.join(env_bin_directory, '../software'))

		if (Debug):
			print ("Retrieve environment path as installation path:")
			print ("Path: " + options.install_path)

		HCGB_files.create_folder(options.install_path)

	#######################
	## install or only check
	#######################
	option_install = False
	if (options.option == 'install'):
		print ("\n+ Check dependencies")
		print ("+ Try to install all missing dependencies, modules or third party software...")
		option_install = True

		## check if access and permission
		if os.path.isdir(options.install_path):
			if (set_config.access_check(options.install_path, mode=os.F_OK)):
				print ("Installation path is accessible and has permission for installation if necessary")
			else:
				print (colored("\n*** ERROR ****", 'red'))
				print (colored("No access/permission for this path: %s" %options.install_path, 'red'))
				print (colored("Please provide a valid path with access/permission to install any missing dependencies.", 'red'))
				exit()
		else:
			print (colored("\n*** ERROR ****", 'red'))
			print (colored("Path provided is not a folder", 'red'))
			print (options.install_path)
			exit()

	elif (options.option == 'only_check'):
		print ("\nCheck dependencies, modules or third party software and print report...")

	#######################
	## python version
	#######################
	HCGB_aes.print_sepLine("+", 20, False)
	print ('Python:')
	HCGB_aes.print_sepLine("+", 20, False)

	this_python_version = str(sys.version)
	python_min_version = extern_progs.return_min_version_soft('python')
	if LooseVersion(this_python_version) >= LooseVersion(python_min_version):
		print (colored("Minimum version (%s) satisfied: %s" %( python_min_version, this_python_version), 'green'))
	else:
		print (colored("Minimum version (%s) not satisfied: %s" %(python_min_version, this_python_version), 'red'))
		exit()
		
	#######################
	## perl_version
	#######################
	print ('\n')
	HCGB_aes.print_sepLine("+", 50, False)
	print ('Perl:')
	HCGB_aes.print_sepLine("+", 50, False)
	
	perl_min_version = extern_progs.return_min_version_soft('perl')
	this_perl_path = set_config.get_exe("perl", Debug)
	this_perl_version = set_config.get_version("perl", this_perl_path, Debug)
	if LooseVersion(this_perl_version) >= LooseVersion(perl_min_version):
		print (colored("Minimum version (%s) satisfied: %s" %(perl_min_version, this_perl_version), 'green'))
	else:
		print (colored("Minimum version (%s) not satisfied: %s" %(perl_min_version, this_perl_version), 'red'))
		exit()

	#######################
	## third-party software
	#######################
	print ('\n')
	HCGB_aes.print_sepLine("+", 20, False)
	print ('External dependencies:')
	HCGB_aes.print_sepLine("+", 20, False)
	
	set_config.check_dependencies(option_install, options.install_path, Debug)
	print ('\n')	

	#######################
	## python packages
	#######################
	print ('\n')
	HCGB_aes.print_sepLine("+", 20, False)
	print ('Python packages:')
	HCGB_aes.print_sepLine("+", 20, False)

	set_config.check_python_packages(Debug, option_install, options.install_path)
	HCGB_aes.print_sepLine("+", 20, False)
	print ('\n')

	#######################
	## perl packages
	#######################
	print ('\n')
	HCGB_aes.print_sepLine("+", 20, False)
	print ('Perl packages:')
	HCGB_aes.print_sepLine("+", 20, False)

	set_config.check_perl_packages("perl_dependencies", Debug, option_install, options.install_path)
	HCGB_aes.print_sepLine("+", 20, False)
	print ('\n')

	#######################
	## IslandPath dependencies
	#######################
	if (options.IslandPath):
		print ('\n')
		HCGB_aes.print_sepLine("+", 20, False)
		print ('IslandPath packages and software required:')
		HCGB_aes.print_sepLine("+", 20, False)
	
		set_config.check_IslandPath(Debug, option_install, options.install_path)
		HCGB_aes.print_sepLine("+", 20, False)
		print ('\n')

	#######################
	## R packages
	#######################
	print ('\n')
	HCGB_aes.print_sepLine("+", 20, False)
	print ('R packages:')
	HCGB_aes.print_sepLine("+", 20, False)

	set_config.check_R_packages(option_install, options.install_path, Debug)
	HCGB_aes.print_sepLine("+", 20, False)
	print ('\n')
