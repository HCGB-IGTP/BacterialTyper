#!/usr/bin/env python3
'''
This code calls MASH for clustering of fasta files or fastq reads.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## import useful modules
import os
import sys
import re
import time
from io import open
import shutil
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.modules import sample_prepare
from BacterialTyper import sampleParser
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper.modules import info
from BacterialTyper import database_generator
from BacterialTyper import min_hash_caller

##############################################
def run(options):

	## init time
	start_time_total = time.time()

	##################################
	### show help messages if desired	
	##################################
	if (options.help_project):
		## information for project
		info.project_help()
		exit()
	elif (options.help_Mash):
		## information for Min Hash Software
		min_hash_caller.helpMash()
		exit()
		
	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False
		
	### set as default paired_end mode
	#if (options.single_end):
	#	options.pair = False
	#else:
	#	options.pair = True
	
	functions.pipeline_header()
	functions.boxymcboxface("Clustering samples")
	print ("--------- Starting Process ---------")
	functions.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	## set mode: project/detached
	if (options.project):
		outdir = input_dir		
	elif (options.detached):
		outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "assembly", "fna")

	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)

	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		functions.create_folder(outdir)
	
	## for each sample
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "mash")	

	## debug message
	if (Debug):
		print (colored("**DEBUG: outdir_dict **", 'yellow'))
		print (outdir_dict)

	## get databases to check
	retrieve_databases = get_options_db(options)
	
	## time stamp
	start_time_partial = functions.timestamp(start_time_total)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: retrieve_database **", 'yellow'))
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)
		print (retrieve_databases)
		
	### Cluster project samples 
	print ("\n+ Generate mash sketches for each sample analyzed...")
	pd_samples_retrieved = pd_samples_retrieved.set_index('name')
	
	## init dataframe
	colname = ["source", "name", "path", "original"]
	pd_samples_sketched  = pd.DataFrame(columns = colname)
	mash_bin = config.get_exe('mash')
	for index, row in pd_samples_retrieved.iterrows():
		if index in retrieve_databases.index:
			print (colored(' + No need to sketch sample (%s), already available within user data...' %index, 'yellow'))
			continue
		name_out = outdir_dict[index] + '/' + index
		min_hash_caller.sketch_database([row['sample']], mash_bin, name_out, index, '')
		pd_samples_sketched.loc[len(pd_samples_sketched)] = ('project_data', index, name_out + '.msh', row['sample'])

	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_sketched **", 'yellow'))
		print (pd_samples_sketched)

	print ("\n+ Clustering sequences...")
	tmp = retrieve_databases[['source', 'db', 'path', 'original']]
	tmp = tmp.rename({'db': 'name'}).set_index('db')
	pd_samples_sketched = pd_samples_sketched.set_index('name')
	cluster_df = pd.concat([pd_samples_sketched, tmp], join='inner', sort=True)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: cluster_df **", 'yellow'))
		print (cluster_df)



############################################################	
def get_options_db(options):
	##
	## Among all databases available and according to the input options,
	## select the databases to use and set dataframe with this information
	##
	
	print ("\n\n+ Select databases to use for identification:")
	
	### database folder to use
	database2use = os.path.abspath(options.database)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: Database to use: " +  database2use + " **", 'yellow'))
	
	## external file provided: single or batch
	if (options.external_file):
		abs_path_ext_file = os.path.abspath(options.external_file)
		if options.batch_external:
			myList = functions.readList_fromFile(abs_path_ext_file)
			join_str = ','.join(myList)
		else:
			join_str = abs_path_ext_file

	############################################################
	### Options: according to user input: select databases to use
	option_db = ""
	
	############
	## 1) only user data: previously identified and added
	############
	if (options.user_data):
		option_db = "Mash:user_data"
		
	############
	## 2) only genbank data: previously download from NCBI reference genomes
	############
	elif (options.genbank_data):
		option_db = "Mash:genbank"
	
	############
	## 3) only project_data
	############
	elif (options.only_project_data):
		option_db = "Mash:project_data"
		return()

	############
	## 4) only external_data
	############
	elif (options.only_external_data):
		option_db = "Mash_external_data:" + join_str

	#################
	## all databases 
	#################
	else:		
		#############################
		option_db = 'Mash:user_data#Mash:genbank'
		if (options.external_file):
			option_db = option_db + '#Mash_external_data:' + join_str
	
	###############
	### get dbs
	###############
	print ("\n+ Parsing information to retrieve databases")
	print ("+ Reading from database: " + database2use)
	functions.print_sepLine("-",50, False)

	###############
	## debug message
	if (Debug):
		print (colored("**DEBUG: option_db: " +  option_db + " **", 'yellow'))

	pd_MASH = database_generator.getdbs("MASH", database2use, option_db, Debug)
	functions.print_sepLine("-",50, False)

	## return both dataframes
	return (pd_MASH)

