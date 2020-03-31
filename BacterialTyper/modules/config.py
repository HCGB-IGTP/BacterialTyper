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
from BacterialTyper.config import set_config
from BacterialTyper.config import extern_progs
from BacterialTyper.scripts import functions
from BacterialTyper.config import install_dependencies

def run(options):
	"""
	This is the main function of the module ``config``. It basically checks 
	if the different requirements (python, perl and third-party software) are
	fulfilled. 

	If any requirement is not available this modules tries to install them or reports to the user to
	manually install them.

	:param option: State wether to check or install missing modules, packages and third party software. Provide: check/install
	:param install_path: Absolute path to install modules or packages missing. Default: ``BacterialTyper`` environment path.
	:param debug: True/false for debugging messages.
	
	:type option: string 
	:param install_path: string 
	:param debug: boolean
	

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

	functions.pipeline_header()
	functions.boxymcboxface("Pipeline Configuration")
	print ("--------- Starting Process ---------")
	functions.print_time()

	if (options.install_path):
		if (Debug):
			print ("Installation path provided for missing modules, packages, dependencies...")
			print ("Path: " + options.install_path)
			
	else:
		## get python environment path
		from distutils.sysconfig import get_python_lib; 
		install_path_tmp = get_python_lib()
		
		##os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'templates'))
		options.install_path = os.path.abspath(os.path.join(install_path_tmp, '../../..'))
			
		if (Debug):
			print ("Retrieve environment path as installation path:")
			print ("Path: " + options.install_path)

	## install or only check
	option_install = False
	if (options.option == 'install'):
		print ("\nCheck dependencies, modules or third party software...")
		print ("\nTry to install all missing dependencies, modules or third party software...")
		option_install = True
	
		## check if access and permission
		if (set_config.access_check(options.install_path)):
			print ("Installation path is accessible and has permission for installation if necessary")
		else:
			print (colored("\n*** ERROR ****", 'red'))
			print (colored("No access/permission for this path: %s" %options.install_path, 'red'))
			print (colored("Please provide a valid path with access/permission to install any missing dependencies.", 'red'))
			exit()
			
	elif (options.option == 'only_check'):
		print ("\nCheck dependencies, modules or third party software and print report...")

	## python version
	functions.print_sepLine("+", 20, False)
	print ('Python:')
	functions.print_sepLine("+", 20, False)

	this_python_version = str(sys.version)
	python_min_version = extern_progs.return_min_version_soft('python')
	if LooseVersion(this_python_version) >= LooseVersion(python_min_version):
		print ("Minimun version (%s) satistied: %s" %( python_min_version, this_python_version))

	## perl_version
	print ('\n')
	functions.print_sepLine("+", 50, False)
	print ('Perl:')
	functions.print_sepLine("+", 50, False)
	
	perl_min_version = extern_progs.return_min_version_soft('perl')
	this_perl_path = set_config.get_exe("perl", Debug)
	this_perl_version = set_config.get_version("perl", this_perl_path, Debug)
	if LooseVersion(this_perl_version) >= LooseVersion(perl_min_version):
		print ("Minimun version (%s) satistied: %s" %(perl_min_version, this_perl_version))

	## third-party software
	print ('\n')
	functions.print_sepLine("+", 20, False)
	print ('External dependencies:')
	functions.print_sepLine("+", 20, False)
	
	set_config.check_dependencies(options.option, options.install_path)
	print ('\n')	

	## python packages
	print ('\n')
	functions.print_sepLine("+", 20, False)
	print ('Python packages:')
	functions.print_sepLine("+", 20, False)

	set_config.check_python_packages(Debug, option_install, options.install_path)
	functions.print_sepLine("+", 20, False)
	print ('\n')

	## perl packages
	print ('\n')
	functions.print_sepLine("+", 20, False)
	print ('Perl packages:')
	functions.print_sepLine("+", 20, False)

	set_config.check_perl_packages(Debug, option_install, options.install_path)
	functions.print_sepLine("+", 20, False)
	print ('\n')

