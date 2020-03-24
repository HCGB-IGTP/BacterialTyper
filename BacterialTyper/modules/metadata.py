#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Adds metadata information for each sample analyzed and prepares data for later visualization.
'''
## import useful modules
import os
import sys
import re
import time
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import database_generator
from BacterialTyper import database_user

##############################################
def run(options):

	## init time
	start_time_total = time.time()

	##################################
	### show help messages if desired	
	##################################
	if (options.help_metadata):
		## help_format option
		help_metadata()
		exit()