#!/usr/bin/env python3
'''
This code calls 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import useful modules
import time
import os
import sys
from io import open
import shutil
import pandas as pd

## import my modules
from BacterialTyper import sampleParser
from BacterialTyper import functions
from BacterialTyper import config

global merge
merge = False

################################
def get_files(options, input_dir, mode, extension):
	## get list of input files
	files = []
	print ('+ Get input folder(s)')
	if (options.batch):
		if os.path.isfile(input_dir):
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
	
	## get information
	if mode == "fastq":
		pd_samples_retrieved = sampleParser.select_samples(files, samples_names, options.pair, exclude, merge)
	else:
		pd_samples_retrieved = sampleParser.select_other_samples(files, samples_names, mode, extension, exclude)
		
	return(pd_samples_retrieved)


################################
def retrieve(options):
	
	## help_format option
	if (options.help_format):
		sampleParser.help_format()
		exit()
		
	functions.pipeline_header()
	functions.boxymcboxface("Preparing samples")
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	### set as defaulta paired_end mode
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True

	## merge option
	if (options.merge):
		global merge
		merge = True

	## get files
	pd_samples_retrieved = get_files(options, input_dir, "fastq", "fastq")

	## output: generate symbolic link or copy if desired	
	functions.create_folder(outdir)

	## merge option
	if (options.merge):
		print ("+ Sample files will be merged...")
		pd_samples_merged = sampleParser.one_file_per_sample(pd_samples_retrieved, outdir, options.threads)
		
		if (options.rename):
			print ("+ Merge files will be renamed...")
			pd_samples_retrieved = pd_samples_merged
		else:
			print ("+ Sample files have been merged...")
			print ("+ Process finished...")
			exit ()

	## Information returned and push into pd_samples_retrieved

	### pd_samples_retrieved from pd_samples_merged
	### name, read_pair, sample, ext, gz

	### pd_samples_retrieved from get_files
	### sample, dirname, name, lane, read_pair, lane_file, ext, gz

	###
	if (options.rename):
		# rename files will be copied into the output folder provided
		names_retrieved = pd.read_csv(options.rename, sep=',', index_col=0, squeeze=True, header=None).to_dict() ## read csv to dictionaru

		## debug
		#print ("")
		#print (names_retrieved)
		#print ("")
		#print (pd_samples_retrieved)
		#print ("")
		
		## rename files 		
		for index, row in pd_samples_retrieved.iterrows():
			renamed = row['dirname'] + '/' + names_retrieved[row['name']] + '_' + row['read_pair'] + row['ext'] + row['gz']
			os.rename(row['sample'], renamed)
			##elif (options.single_end): It should work for both

		print ("+ Sample files have been renamed...")
		
	else:
		## keep original names
		if (options.copy):
			print ("+ Sample files will be copied...")
	
		## LIST READS
		list_reads = []
		for index, row in pd_samples_retrieved.iterrows():
			if (options.copy):
				shutil.copy(row['sample'], outdir)
			else:
				list_reads.append(row['sample'])

		if (options.copy):
			print ("+ Sample files have been copied...")
		else:
			print ("+ Sample files will be linked...")
			functions.get_symbolic_link(list_reads, outdir)
	
