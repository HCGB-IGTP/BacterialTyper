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
	if (options.project):
		### a folder containing a project is provided
		if os.path.exists(input_dir):
			print ('+ Input folder exists')
			## get files in folder
			files_tmp = functions.get_fullpath_list(input_dir)
			files = [s for s in files_tmp if extension in s]
			
	else:
		### provide a single folder or a file with multiple paths (option batch)
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
	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	### set as default paired_end mode
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True

	## merge option
	if (options.merge):
		global merge
		merge = True

	## initial retrieved
	project_mode=False
	if (options.project):
		options.project = False
		project_mode=True
		
	## get files
	pd_samples_retrieved = get_files(options, input_dir, "fastq", "fastq")

	## set option
	if (project_mode):
		options.project = True

	## output: generate symbolic link or copy if desired	
	print ("\n+ Create output folder(s):")
	functions.create_folder(outdir)
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "raw")

	###
	final_dir = ""
	if (options.project):
		print ("+ Generate a directory containinig information within the project folder provided")
		final_dir = functions.create_subfolder("info", outdir)
	else:
		final_dir = outdir
	
	## merge option
	if (options.merge):
		print ("+ Sample files will be merged...")
		pd_samples_merged = sampleParser.one_file_per_sample(pd_samples_retrieved, outdir_dict, options.threads, final_dir)
		
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
	
		## print to a file
		timestamp = functions.create_human_timestamp()
		rename_details = final_dir + '/' + timestamp + '_prep_renameDetails.txt'
		rename_details_hd = open(rename_details, 'w')
	
		## rename files 		
		for index, row in pd_samples_retrieved.iterrows():
			renamed = row['dirname'] + '/' + names_retrieved[row['name']] + '_' + row['read_pair'] + row['ext'] + row['gz']
			os.rename(row['sample'], renamed)
			string = row['sample'] + '\t' + renamed + '\n'
			rename_details_hd.write(string)

		rename_details_hd.close()	

		##elif (options.single_end): It should work for both
		print ("+ Sample files have been renamed...")
		
	else:
		## keep original names
		if (options.copy):
			print ("+ Sample files will be copied...")
			## print to a file
			timestamp = functions.create_human_timestamp()
			copy_details = final_dir + '/' + timestamp + '_prep_copyDetails.txt'
			copy_details_hd = open(copy_details, 'w')

		else:
			print ("+ Sample files will be linked...")
	
		## LIST READS
		list_reads = []
		for index, row in pd_samples_retrieved.iterrows():
			if (options.copy):
				shutil.copy(row['sample'], outdir_dict[row['name']])
				
				string = row['sample'] + '\t' + outdir_dict[row['name']] + '\n'
				copy_details_hd.write(string)
				
			else:
				list_reads.append(row['sample'])
				
				if options.project:
					list_single = [row['sample']] 
					functions.get_symbolic_link(list_single, outdir_dict[row['name']])

		if (options.copy):
			print ("+ Sample files have been copied...")
			copy_details_hd.close()
		else:
			if not options.project:
				functions.get_symbolic_link(list_reads, outdir)
	

			

	
