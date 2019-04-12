#!/usr/bin/env python3
'''
This code calls a configuration function to iniate and update the configuration of the pipeline.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report

def run(options):
	print ("This is a configuration script for the pipeline...")
