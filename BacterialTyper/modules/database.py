#!/usr/bin/env python3
'''
This code calls a database generator script for initiating, updating and configure database.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report


def run(options):
	
	print ("Hello world!\n")
