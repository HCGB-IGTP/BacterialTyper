#usr/bin/en python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Retrieves files within ``other_tools`` directory and returns path to given script specified
"""
## useful imports
import os
## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import install_dependencies

####################################################################
def perl_scripts(script):
	"""Lists files within ``other_tools/perl`` directory and returns path to given script.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.functions.get_fullpath_list`
	"""
	perlDir = os.path.dirname(os.path.realpath(__file__)) + '/perl/'
	list_perl = functions.get_fullpath_list(perlDir)
	
	dict_perl = {}
	for f in list_perl:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == script):
			return (f)		

####################################################################
def R_scripts(script):
	"""Lists files within ``other_tools/R`` directory and returns path to given script
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.functions.get_fullpath_list`
	"""
	RDir = os.path.dirname(os.path.realpath(__file__)) + '/R/'
	list_R = functions.get_fullpath_list(RDir)
	
	dict_R = {}
	for f in list_R:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == script):
			return (f)
		


		
	