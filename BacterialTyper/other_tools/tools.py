#usr/bin/en python3
'''
This code...
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import os
## import my modules
from BacterialTyper import functions

####################################################################
def perl_scripts(script):
	perlDir = os.path.dirname(os.path.realpath(__file__)) + '/perl/'
	list_perl = functions.get_fullpath_list(perlDir)
	
	dict_perl = {}
	for f in list_perl:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == script):
			return (f)		

####################################################################
def R_scripts(script):
	perlDir = os.path.dirname(os.path.realpath(__file__)) + '/R/'
	list_R = functions.get_fullpath_list(RDir)
	
	dict_R = {}
	for f in list_R:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == script):
			return (f)


