#!/usr/bin/env python3
'''
This code calls a database generator script for initiating, updating and configure database.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import database_generator
from BacterialTyper import multiQC_report

def init_DB(options):
	database_generator.init_DB(options.path, options.ID_file)


