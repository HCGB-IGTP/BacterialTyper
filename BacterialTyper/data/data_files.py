#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Allows access to additional data stored within the repository.
'''
## useful imports
import os
## import my modules
from BacterialTyper.scripts import functions

####################################################################
def data_list(wanted_data):
	data = os.path.dirname(os.path.realpath(__file__))
	list_data = functions.get_fullpath_list(data)
	dict_data = {}
	for f in list_data:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == wanted_data):
			return (f)