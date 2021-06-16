#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Allows access to additional information stored in different files 
under folder ``BacterialTyper/data``.
"""
## useful imports
import os
## import my modules
import HCGB.functions.main_functions as HCGB_main

####################################################################
def data_list(wanted_data):
	"""
	Retrieves information of additional files under folder ``BacterialTyper/data``.
	"""

	data = os.path.dirname(os.path.realpath(__file__))
	list_data = HCGB_main.get_fullpath_list(data)
	for f in list_data:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == wanted_data):
			return (f)