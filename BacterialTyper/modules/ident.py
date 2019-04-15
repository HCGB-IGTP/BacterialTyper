#!/usr/bin/env python3
'''
This code calls species_identification_KMA and get the most similar taxa then use ariba_caller to check if pubmlst is downloaded and  get MLST profile.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import species_identification_KMA
from BacterialTyper import config

####
def run(options):

	### species_identification_KMA -> most similar taxa
	### ariba_caller -> MLST profile.
	
	print ("Hi there! lets identify the taxa\n")
