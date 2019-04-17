#!/usr/bin/env python3
'''
This code calls 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import sampleParser
from BacterialTyper import functions
from BacterialTyper import config
import os
import sys
from io import open
import shutil

################################
def merge(options):

	## extract files
	print ()

def retrieve(options):
	
	functions.pipeline_header()
	functions.boxymcboxface("Preparing samples")
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get list of input files
	files = []
	
	print ('+ Get input folder(s)')

	if (options.batch):
		if os.path.isfile(input_dir):
			print ('+ Batch file provided exists')
			dir_list = [line.rstrip('\n') for line in open(input_dir)]
			for d in dir_list:
				if os.path.exists(d):
					print ('+ Folder (%s) exists' %d)
					files = files + functions.get_fullpath_list(d)
	else:
		if os.path.exists(input_dir):
			print ('+ Input folder exists')
			## get files in folder
			files = functions.get_fullpath_list(input_dir)
	
	#cprint (files)
	
	## get list of samples
	samples_names = []
	exclude=False

	if (options.in_sample):
		in_file = os.path.abspath(options.in_sample)
		samples_names = [line.rstrip('\n') for line in open(in_file)]
		print ('+ Retrieve selected samples to obtain from the list files available.')		
		exclude=False
	elif (options.ex_sample):
		ex_file = os.path.abspath(options.ex_sample)
		samples_names = [line.rstrip('\n') for line in open(ex_file)]
		print ('+ Retrieve selected samples to exclude from the list files available.')		
		exclude=True
		
		print (samples_names)
	else:
		samples_names = ['.*']
	
	##
	list_samples_retrieved = sampleParser.select_samples(files, samples_names, options.pair, exclude)

	## output: generate symbolic link or copy if desired	
	functions.create_folder(outdir)
	if (options.copy):
		for f in list_samples_retrieved:
			shutil.copy(f, outdir)
			#functions.get_symbolic_link(list_samples_retrieved, outdir)
	else:
		functions.get_symbolic_link(list_samples_retrieved, outdir)
	

