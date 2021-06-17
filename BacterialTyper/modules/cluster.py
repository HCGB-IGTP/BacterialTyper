#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Clusters fasta files or fastq reads. Using: project data, genbank entries from database or previous samples. 
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
from HCGB import sampleParser
from BacterialTyper.config import set_config
from BacterialTyper.modules import help_info
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import min_hash_caller
from BacterialTyper import __version__ as pipeline_version

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files

##############################################
def run_cluster(options):

	## init time
	start_time_total = time.time()

	##################################
	### show help messages if desired	
	##################################
	if (options.help_project):
		## information for project
		help_info.project_help()
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
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True
	
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Clustering samples")
	print ("--------- Starting Process ---------")
	HCGB_time.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	## Project mode as default
	project_mode=True
	if (options.detached):
		options.project = False
		project_mode=False
		outdir = os.path.abspath(options.output_folder)
	else:
		options.project = True
		outdir = input_dir	
	
	## get files
	if options.reads:
		if options.noTrim:
			## raw reads
			pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
		else:
			## trimm reads
			pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)

		## keep only R1 reads if paired-end
		if options.pair:
			pd_samples_retrieved = pd_samples_retrieved.loc[pd_samples_retrieved['read_pair']== "R1"]

	else:
		## default
		pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "assembly", ["fna"], options.debug)

	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
	
	# exit if empty
	if pd_samples_retrieved.empty:
		print ("No data has been retrieved from the project folder provided. Exiting now...")
		exit()


	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		HCGB_files.create_folder(outdir)
	
	## for each sample
	outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "mash", options.debug)	

	## debug message
	if (Debug):
		print (colored("**DEBUG: outdir_dict **", 'yellow'))
		print (outdir_dict)

	## get databases to check
	retrieve_databases = get_options_db(options)
	
	## time stamp
	start_time_partial = HCGB_time.timestamp(start_time_total)
	
	## remove samples if specified
	if options.ex_sample:
		ex_samples = HCGB_main.get_info_file(options.ex_sample)
		retrieve_databases = retrieve_databases.loc[~retrieve_databases.index.isin(ex_samples)]
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: retrieve_database **", 'yellow'))
		pd.set_option('display.max_colwidth', None)
		pd.set_option('display.max_columns', None)
		print (retrieve_databases)

	## check if all samples in user_data or genbank are indexed
	siglist_all = []
	for index, row in retrieve_databases.iterrows():
		if not row['path'] == 'NaN':
			if (Debug):
				HCGB_aes.print_sepLine("*",25, False)
				print (row)
				
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
				original = HCGB_main.readList_fromFile(file2print)
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

	print ("\n+ Clustering sequences...")
	pd_samples_sketched = pd_samples_sketched.set_index('name')
	
	####
	if retrieve_databases.empty: 
		cluster_df = pd_samples_sketched
	else:
		tmp = retrieve_databases[['source', 'db', 'path', 'original', 'ksize', 'num_sketch']]
		tmp = tmp.rename(columns={'db': 'name'})
		tmp.set_index('name')
		
		if (Debug):
			print (colored("**DEBUG: tmp **", 'yellow'))
			print (tmp)	
		
		## merge both dataframes
		cluster_df = pd.concat([pd_samples_sketched, tmp], join='inner', sort=True)
	
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_sketched **", 'yellow'))
		print (pd_samples_sketched)
		
		print (colored("**DEBUG: cluster_df **", 'yellow'))
		print (cluster_df)
		
		print (colored("**DEBUG: Signatures **", 'yellow'))
		print (siglist_all)
		
		print (colored("**DEBUG: length siglist_all **", 'yellow'))
		print (len(siglist_all))

	## Assign Colors colorLabels	
	color_df = cluster_df.filter(["source"], axis=1)
	color_df["color"] = "r" ## red::genbank

	## project data
	project_data = list( color_df[ color_df["source"] == "project_data"].index)
	color_df.loc[ color_df.index.isin(project_data), "color"] = "g" ## green::project_data

	## user_data
	user_data = list( color_df[ color_df["source"] == "user_data"].index)
	color_df.loc[color_df.index.isin(user_data), "color"] = "b" ## blue::user_data
	
	colorLabels = color_df['color'].to_dict()
	
	if Debug:
		print(color_df)
		print(colorLabels)

	## parse results
	if options.project:
		outdir_report = HCGB_files.create_subfolder("report", outdir)
		#final_dir = outdir + '/report/cluster'
		final_dir = functions.create_subfolder("cluster", outdir_report) 
	else:
		final_dir = outdir

	## compare
	name = 'cluster_' + str(HCGB_time.create_human_timestamp())
	tag_cluster_info = final_dir + '/' + name
	print ('+ Saving results in folder: ', final_dir)
	print ('\tFile name: ', name)
	(DataMatrix, labeltext) = min_hash_caller.compare(siglist_all, tag_cluster_info, Debug)

	## get colorLabels	
	
	## plot images
	pdf = True
	cluster_returned = min_hash_caller.plot(DataMatrix, labeltext, tag_cluster_info, pdf, colorLabels)
	
	## generate newick tree 
	min_hash_caller.get_Newick_tree(cluster_returned, DataMatrix, labeltext, tag_cluster_info)	
	
	print ("\n*************** Finish *******************")
	start_time_partial = HCGB_time.timestamp(start_time_total)

	## dump information and parameters
	info_dir = HCGB_files.create_subfolder("info", outdir)
	print("+ Dumping information and parameters")
	runInfo = { "module":"cluster", "time":HCGB_time.timestamp(time.time()),
				"BacterialTyper version":pipeline_version }
	HCGB_info.dump_info_run(info_dir, 'cluster', options, runInfo, options.debug)

	print ("+ Exiting cluster module.")
	return()

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
			myList = HCGB_main.readList_fromFile(abs_path_ext_file)
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
	HCGB_aes.print_sepLine("-",50, False)

	###############
	## debug message
	if (Debug):
		print (colored("**DEBUG: option_db: " +  option_db + " **", 'yellow'))

	pd_MASH = database_generator.getdbs("MASH", database2use, option_db, Debug)
	HCGB_aes.print_sepLine("-",50, False)

	## return both dataframes
	return (pd_MASH)

############################################################	
def generate_sketch(folder, assembly, entry, ksize, n_sketch, Debug):

	(sigfile, siglist) = min_hash_caller.sketch_database({ entry: assembly }, folder, Debug, ksize, n_sketch)
	#functions.print_sepLine("*",50, False)
	
	## print original in file
	file2print = folder + '/.original'
	
	## do not write full path because database can move.
	## as long as mantains the folder organization it should work.
	assembly_tmp_path = "../assembly/" + os.path.basename(assembly)
	
	list_fna = [assembly_tmp_path, str(ksize), str(n_sketch)]
	HCGB_main.printList2file(file2print, list_fna)
	return (sigfile[0], siglist[0])

	
