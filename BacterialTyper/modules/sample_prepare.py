#!/usr/bin/env python3
'''
This code calls 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import useful modules
import os
import sys
from io import open
import shutil

## import my modules
from BacterialTyper import sampleParser
from BacterialTyper import functions
from BacterialTyper import config

################################
def merge(options):

	## extract files
	print ('+ Merge files from different sequencing lanes...')

################################
def get_files(options, input_dir):
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
		#print (samples_names)
	else:
		samples_names = ['.*']
	
	##
	pd_samples_retrieved = sampleParser.select_samples(files, samples_names, options.pair, exclude)
	return(pd_samples_retrieved)

################################
def retrieve(options):
	
	functions.pipeline_header()
	functions.boxymcboxface("Preparing samples")
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = get_files(options, input_dir)
	
	## output: generate symbolic link or copy if desired	
	functions.create_folder(outdir)
	list_reads = []
	if (options.copy):
		print ("+ Sample files will be copied...")
	for index, row in pd_samples_retrieved.iterrows():
		if (options.pair):
			#print(row['samples'], row['R1'], row['R2'])
			if (options.copy):
				shutil.copy(row['R1'], outdir)
				shutil.copy(row['R2'], outdir)
			else:
				list_reads.append(row['R1'])
				list_reads.append(row['R2'])
		else:
			if (options.copy):
				shutil.copy(row['read'], outdir)
			else:
				list_reads.append(row['read'])

	if (options.copy):
		print ("+ Sample files have been copied...")
	else:
		print ("+ Sample files will be linked...")
		functions.get_symbolic_link(list_reads, outdir)
	
