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
	global Project
	if (options.project):
		outdir = input_dir		
		Project=True
	elif (options.detached):
		Project=False
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

	## check if all samples in user_data or genbank are indexed
	siglist_all = []
	for index, row in retrieve_databases.iterrows():
		if not row['path'] == 'NaN':
			#print (row)
			if all([ int(options.kmer_size) == int(row['ksize']), int(options.n_sketch) == int(row['num_sketch']) ]):
				siglist_all.append(min_hash_caller.read_signature(row['path'], options.kmer_size))
				continue

		## index assembly or reads...
		(sigfile, siglist) = generate_sketch(row['folder'], row['original'], index, options.kmer_size, options.n_sketch, Debug)
		retrieve_databases.loc[index]['path'] = sigfile
		retrieve_databases.loc[index]['ksize'] = options.kmer_size
		retrieve_databases.loc[index]['num_sketch'] = options.n_sketch
		siglist_all.append(siglist)
	
	### Cluster project samples
	print (colored("\n+ Collect project data", 'green'))
	print ("+ Generate mash sketches for each sample analyzed...")
	pd_samples_retrieved = pd_samples_retrieved.set_index('name')
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieved **", 'yellow'))
		print (pd_samples_retrieved)

	## init dataframe for project data
	colname = ["source", "name", "path", "original", "ksize", "num_sketch"]
	pd_samples_sketched  = pd.DataFrame(columns = colname)
	for index, row in pd_samples_retrieved.iterrows():
		if index in retrieve_databases.index:
			print (colored('\t+ Sketched signature (%s) available within user data...' %index, 'yellow'))
			continue

		this_sig = outdir_dict[index] + '/' + index + '.sig'
		if os.path.exists(this_sig):
			## File signature might exist

			## read original
			file2print = outdir_dict[index] + '/.original'
			if not os.path.exists(file2print):
				original = ['NaN']
			else:
				original = functions.readList_fromFile(file2print)
				if all([ int(options.kmer_size) == int(original[1]), int(options.n_sketch) == int(original[2])]):
					siglist_all.append(min_hash_caller.read_signature(this_sig, options.kmer_size)) 
					pd_samples_sketched.loc[len(pd_samples_sketched)] = ('project_data', index, this_sig, row['sample'], options.kmer_size, options.n_sketch)
					print (colored('\t+ Sketched signature available (%s) in project folder...' %index, 'green'))
					continue			
					
		print (colored('\t+ Sketched signature to be generated: (%s)...' %index, 'yellow'))
		## index assembly or reads...
		(sigfile, siglist) = generate_sketch(outdir_dict[index], row['sample'], index, options.kmer_size, options.n_sketch, Debug)
		pd_samples_sketched.loc[len(pd_samples_sketched)] = ('project_data', index, sigfile, row['sample'], options.kmer_size, options.n_sketch)
		siglist_all.append(siglist)

	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_sketched **", 'yellow'))
		print (pd_samples_sketched)

	print ("\n+ Clustering sequences...")
	pd_samples_sketched = pd_samples_sketched.set_index('name')
	
	####
	if not retrieve_databases.empty: 
		tmp = retrieve_databases[['source', 'db', 'path', 'original', 'ksize', 'num_sketch']]
		tmp = tmp.rename({'db': 'name'}).set_index('db')
	else:
		tmp = retrieve_databases
	
	## merge both dataframes
	cluster_df = pd.concat([pd_samples_sketched, tmp], join='inner', sort=True)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: cluster_df **", 'yellow'))
		print (cluster_df)
		print (colored("**DEBUG: Signatures **", 'yellow'))
		print (siglist_all)

	## get missing signature files
	print ('\t+ Loading missing signature for comparison...')
	#list_signatures2read = cluster_df['path'].to_list()
	#for l in list_signatures2read:
	#	siglist_all += min_hash_caller.read_signature(l, options.kmer_size)
	
	### siglist_all = set(siglist_all)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: Signatures **", 'yellow'))
		print (siglist_all)
		print (colored("**DEBUG: length siglist_all **", 'yellow'))
		print (len(siglist_all))

	## parse results
	if Project:
		final_dir = outdir + '/report/cluster'
		functions.create_folder(final_dir) 
	else:
		final_dir = outdir

	## compare
	name = 'cluster_' + str(functions.create_human_timestamp())
	tag_cluster_info = final_dir + '/' + name
	print ('+ Saving results in folder: ', final_dir)
	print ('\tFile name: ', name)
	
	pdf = True
	(DataMatrix, labeltext) = min_hash_caller.compare(siglist_all, tag_cluster_info, Debug)
	min_hash_caller.plot(DataMatrix, labeltext, tag_cluster_info, pdf)
	

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
		pd_MASH = pd.DataFrame()
		return(pd_MASH)

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

############################################################	
def generate_sketch(folder, assembly, entry, ksize, n_sketch, Debug):

	(sigfile, siglist) = min_hash_caller.sketch_database({ entry: assembly }, folder, Debug, ksize, n_sketch)
	#functions.print_sepLine("*",50, False)
	
	## print original in file
	file2print = folder + '/.original'
	list_fna = [assembly, str(ksize), str(n_sketch)]
	functions.printList2file(file2print, list_fna)
	return (sigfile, siglist)

	
