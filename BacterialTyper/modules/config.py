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
	This is the main function of the module ``config`` that basically checks 
	if the different requirements (python, perl and third-party software) are
	fulfilled. 
	
	If any requirement is not available this modules tries to install them or reports to the user to
	manually install them.
	
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
	python_min_version = extern_progs.return_min_version_soft('python')
	if LooseVersion(this_python_version) >= LooseVersion(python_min_version):
		print ("Minimun version (%s) satistied: %s" %( python_min_version, this_python_version))

	functions.print_sepLine("+", 50, False)
	

	functions.print_sepLine("+", 50, False)
	print ('Perl:')
	
	## perl version
	perl_min_version = extern_progs.return_min_version_soft('perl')
	
	## perl version installed
	this_perl_path = set_config.get_exe("perl", Debug)
	this_perl_version = set_config.get_version("perl", this_perl_path, Debug)
	
	if LooseVersion(this_perl_version) >= LooseVersion(perl_min_version):
		print ("Minimun version (%s) satistied: %s" %(perl_min_version, this_perl_version))

	functions.print_sepLine("+", 50, False)
	
	
	## third-party software
	functions.print_sepLine("+", 50, False)
	print ('External dependencies:\n')
	functions.print_sepLine("+", 50, False)
	#extern_progs.print_dependencies()
	print ('\n')	
	
	## python packages
	print ('\n')
	functions.print_sepLine("+", 50, False)
	print ('Python packages:')
	set_config.check_python_packages(Debug, option_install)
	functions.print_sepLine("+", 50, False)
	print ('\n')
	
	## perl packages
	print ('\n')
	functions.print_sepLine("+", 50, False)
	print ('Perl packages:')
	set_config.check_perl_packages(Debug, option_install)
	functions.print_sepLine("+", 50, False)
	print ('\n')

