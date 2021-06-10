#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Prepares the database information for further analysis. Several functions are implemented for:
	
	- Project data provided, updates/populates the database of interest

'''
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
import shutil
from sys import argv
from io import open
from termcolor import colored
from Bio import SeqIO
import concurrent.futures

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.scripts import species_identification_KMA
from BacterialTyper.scripts import min_hash_caller
from BacterialTyper.scripts import database_generator

## HCGB module
from HCGB import sampleParser
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

##########################################################################################
def update_database_user_data(database_folder, project_folder, Debug, options):
	"""
	Updates user_data folder within the database folder provided.
	
	It would generate single subfolder for each sample previously analyzed and it would store main information and result files for later interpretation, comparison and/or summarization with new samples analyzed.
	
	:param database_folder:
	:param project_folder:
	:param Debug:
	:param options:
	
	:type database_folder:
	:type project_folder:
	:type Debug:
	:type options:
	
	:returns: Updated database result from :func:`BacterialTyper.scripts.database_generator.update_db_data_file`.
	:rtype: Dataframe
	
	:warnings: Returns **FAIL** if check process failed.
	
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.files_functions.create_subfolder`
		
		- :func:`HCGB.functions.main_functions.functions.get_data`
		
		- :func:`HCGB.functions.main_functions.optimize_threads`
		
		- :func:`BacterialTyper.scripts.database_user.get_userData_files`
		
		- :func:`BacterialTyper.scripts.database_user.update_sample`
		
		- :func:`BacterialTyper.scripts.database_generator.getdbs`
		
		- :func:`BacterialTyper.scripts.database_generator.get_database`
		
		- :func:`BacterialTyper.scripts.database_generator.update_db_data_file`

	"""
	
	print ("\n+ Updating information from user data folder: ", project_folder)
	
	## create folder
	own_data = HCGB_files.create_subfolder("user_data", database_folder)
	
	## Default missing options
	options.project = True
	options.debug = Debug
	if not options.single_end:
		options.pair = True

	####################################
	## get information
	####################################
	
	## get user data files
	project_data_df = get_userData_files(options, project_folder)
	
	## get user data info
	project_info_df = get_userData_info(options, project_folder)
	
	## merge data
	project_all_data = pd.concat([project_data_df, project_info_df], join='outer', sort=True).drop_duplicates()
	#project_all_data.index.name = 'name'

	## debug messages:
	if Debug:
		HCGB_aes.debug_message("project_data_df", 'yellow')
		print(project_data_df)
		
		HCGB_aes.debug_message("project_info_df", 'yellow')
		print(project_info_df)

		HCGB_aes.debug_message("project_all_data", 'yellow')
		print(project_all_data)
	
	print ('\n+ Get database information')
	db_frame = database_generator.getdbs('user_data', database_folder, 'user_data', Debug)
	user_data_db = database_generator.get_database(db_frame, Debug)
	
	## merge dataframe
	sample_frame = project_all_data.groupby("name")
	
	####################################
	## optimize threads
	####################################
	name_list = project_all_data.index.values.tolist()
	threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

	
	print ('\n+ Updating information using %s threads and %s parallel jobs' %(options.threads, max_workers_int))

	####################################
	## loop through frame using multiple threads
	####################################
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		## send for each	
		commandsSent = { executor.submit(update_sample, name, cluster, 
										own_data, user_data_db, Debug): name for name, cluster in sample_frame }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
	
	HCGB_aes.print_sepLine("+", 75, False)
	print ("\n+ Retrieve information ...")

	####################################
	###### populate dataframe
	####################################
	for name, cluster in sample_frame:
		###### dump to file
		info_file = own_data + '/' + name + '/info.txt'
		if os.path.exists(info_file):
			dataGot = HCGB_main.get_data(info_file, ',', 'index_col=0')
			dataGot = dataGot.set_index('ID')
		
			if (options.debug):
				print (colored("**DEBUG: dataGot dataframe **", 'yellow'))
				print (dataGot)
	
			user_data_db = pd.concat([user_data_db, dataGot], join='outer', sort=True).drop_duplicates()
			## concatenating by outer we get all available entries

	if (options.debug):
		print (colored("**DEBUG: user_data_db dataframe **", 'yellow'))
		print (user_data_db)

	HCGB_aes.print_sepLine("+", 75, False)	

	####################################
	## update db		
	####################################
	database_csv = own_data + '/user_database.csv'
	
	dataUpdated = database_generator.update_db_data_file(user_data_db, database_csv)
	print ("+ Database has been generated: \n", database_csv)
	return (dataUpdated)


##########################################################################################
def get_userData_files(options, project_folder):
	## get information regarding files

	## get trimmed ngs files
	print()
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print ("+ Retrieve trimmed reads information:")
	pd_samples_reads = sampleParser.files.get_files(options, project_folder, "trim", ['_trim'], options.debug)
	pd_samples_reads = pd_samples_reads.set_index('name')
	HCGB_aes.print_sepLine("-", 60, 'yellow')

	## get assembly files
	print()
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print ("+ Retrieve assembly information:")
	pd_samples_assembly = sampleParser.files.get_files(options, project_folder, "assembly", ["fna"], options.debug)
	pd_samples_assembly = pd_samples_assembly.set_index('name')
	HCGB_aes.print_sepLine("-", 60, 'yellow')

	## get annotation files
	print()
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print ("+ Retrieve annotation information:")
	pd_samples_annot = sampleParser.files.get_files(options, project_folder, "annot", ['gbf', 'faa', 'gff'], options.debug)
	pd_samples_annot = pd_samples_annot.set_index('name')
	HCGB_aes.print_sepLine("-", 60, 'yellow')

	## debug message
	if (options.debug):
		print (colored("**DEBUG: pd_samples_reads **", 'yellow'))
		print (pd_samples_reads)
		print (colored("**DEBUG: pd_samples_assembly **", 'yellow'))
		print (pd_samples_assembly)
		print (colored("**DEBUG: pd_samples_annot **", 'yellow'))
		print (pd_samples_annot)
	
	## merge
	df = pd.concat([pd_samples_reads, pd_samples_annot, pd_samples_assembly], sort=True, join='inner').drop_duplicates()
	## joining by inner we only get common columns among all

	## debug message
	if (options.debug):
		print (colored("**DEBUG: pd_concat **", 'yellow'))
		print (df)
	
	## set new column with name of samples
	df = df.reset_index()

	## debug message
	if (options.debug):
		print (colored("**DEBUG: pd_concat reset_index**", 'yellow'))
		print (df)
	##
	return(df)
	
##########################################################################################
def get_userData_info(options, project_folder):
	## get information regarding: 
		## genus, species (ident module)
		## card & VFDB (profile module)
		## additional information: MGE, etc

	## get profile information
	print()
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print ("+ Retrieve virulence/resistance profile information:")
	pd_samples_profile = sampleParser.files.get_files(options, project_folder, "profile", ["csv"], options.debug)
	if not pd_samples_profile.empty:
		pd_samples_profile = pd_samples_profile.set_index('name')
	HCGB_aes.print_sepLine("-", 60, 'yellow')

	## get identification information
	print()
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print ("+ Retrieve species identification information:")
	pd_samples_ident = sampleParser.files.get_files(options, project_folder, "ident", ["csv"], options.debug)
	if not pd_samples_ident.empty:
		pd_samples_ident = pd_samples_ident.set_index('name')
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	
	## get mash information
	print()
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print ("+ Retrieve cluster information:")
	pd_samples_mash = sampleParser.files.get_files(options, project_folder, "mash", ["sig"], options.debug)
	if not pd_samples_mash.empty:
		pd_samples_mash = pd_samples_mash.set_index('name')
	HCGB_aes.print_sepLine("-", 60, 'yellow')
	print()

	## add other if necessary

	## debug message	
	if (options.debug):
		print (colored("**DEBUG: pd_samples_profile **", 'yellow'))
		print (pd_samples_profile)
		print (colored("**DEBUG: pd_samples_ident **", 'yellow'))
		print (pd_samples_ident)
		print (colored("**DEBUG: pd_samples_mash **", 'yellow'))
		print (pd_samples_mash)
	
	## merge
	df = pd.concat([pd_samples_profile, pd_samples_ident, pd_samples_mash], join='inner', sort=True).drop_duplicates()
	## joining by inner we only get common columns among all

	## debug message
	if (options.debug):
		print (colored("**DEBUG: pd_concat **", 'yellow'))
		print (df)
	
	## set new column with name of samples
	df = df.reset_index()

	## rename column
	df.rename(columns={'index':'name'}, inplace=True)
	
	## debug message
	if (options.debug):
		print (colored("**DEBUG: pd_concat reset_index**", 'yellow'))
		print (df)

	##
	return(df)

############################################
def update_sample(name, cluster, own_data, user_data_db, Debug):

	## debug message	
	if (Debug):
		print (colored("**DEBUG: sample_frame groupby: name & cluster **", 'yellow'))
		print (name)			
		print (cluster)
	
	if (name == 'report'):
		return()
	
	print ('\t+ Sending command for sample: ', name)

	############################################
	#### check information for this sample
	############################################

	## generate sample
	dir_sample = HCGB_files.create_subfolder(name, own_data)
	
	if name in user_data_db.index:
		print (colored("\t\t+ Data available in database for sample: %s. Checking integrity..." %name, 'yellow'))
		#functions.print_sepLine("+", 75, False)

	## data to generate
	data2dump = pd.DataFrame(columns=('ID','folder','genus','species','name','genome', 'GFF','proteins', 'signature', 'profile', 'ident', 'reads'))
	## iterate over files with different tags: reads, annot, assembly, profile, ident

	##########
	## assembly
	##########
	assembly_dir = HCGB_files.create_subfolder('assembly', dir_sample)
	assembly_file = cluster.loc[cluster['tag'] == 'assembly']['sample'].to_list()
	if assembly_file:
		assembly_file_name = os.path.basename(assembly_file[0])
		genome = assembly_dir + '/' + assembly_file_name
		if not os.path.exists(genome):
			shutil.copy(assembly_file[0], assembly_dir)
	else:
		genome = ""

	##########
	## annot
	##########
	annot_dir = HCGB_files.create_subfolder('annot', dir_sample)
	annot_files = cluster.loc[cluster['tag'] == 'annot']['sample'].to_list()
	prof = ""
	gff = ""
	if annot_files:
		for f in annot_files:
			file_name = os.path.basename(f)
			if f.endswith('faa'):
				prot = annot_dir + '/' + file_name
				if os.path.exists(prot):
					continue
			elif f.endswith('gff'):
				gff = annot_dir + '/' + file_name
				if os.path.exists(gff):
					continue
			shutil.copy(f, annot_dir)
	else:
		gff = ""
		prot = ""

	##########
	## trimm
	##########
	trimm_dir = HCGB_files.create_subfolder('trimm', dir_sample)
	reads_files = cluster.loc[cluster['tag'] == 'reads']['sample'].to_list()
	reads = []
	if reads_files:
		for f in reads_files:
			file_name = os.path.basename(f)
			reads_name = trimm_dir + '/' + file_name
			reads.append(reads_name)
			if not os.path.exists(reads_name):
				shutil.copy(f, trimm_dir)

	##########
	## ident
	##########
	ident_dir = HCGB_files.create_subfolder('ident', dir_sample)
	ident_file = cluster.loc[cluster['tag'] == 'ident']['sample'].to_list()
	if ident_file:
		file_name = os.path.basename(ident_file[0])
		ident_file_name = ident_dir + '/' + file_name
		if not os.path.exists(ident_file_name):
			shutil.copy(ident_file[0], ident_dir)
	else:
		ident_file_name = ""
	
	##########
	## profile
	##########
	profile_dir = HCGB_files.create_subfolder('profile', dir_sample)
	profile_files = cluster.loc[cluster['tag'] == 'profile']['sample'].to_list()
	profile_file = []
	if profile_files:
		for f in profile_files:
			file_name = os.path.basename(f)
			profile_file_name = profile_dir + '/' + file_name
			profile_file.append(profile_file_name)
			if not os.path.exists(profile_file_name):
				shutil.copy(f, profile_dir)
	
	##########
	## mash profile
	##########
	mash_dir = HCGB_files.create_subfolder('mash', dir_sample)
	mash_file = cluster.loc[cluster['tag'] == 'mash']['sample'].to_list()
	if mash_file:
		file_name = os.path.basename(mash_file[0])
		sig_file = mash_dir + '/' + file_name
		if not os.path.exists(sig_file):
			shutil.copy(mash_file[0], mash_dir)
	else:
		sig_file = ""
	
	############################################
	### Dump information

	## TODO: Add species and genus information when parsed from ident csv file
	#####
	data2dump.loc[len(data2dump)] = (name, dir_sample, 'genus', 'species', name, genome, gff, prot, sig_file, '::'.join(sorted(profile_file)), ident_file_name, '::'.join(sorted(reads)) )
	#data2dump = data2dump.set_index('ID')
	
	###### dump to file
	info_file = dir_sample + '/info.txt'
	data2dump.to_csv(info_file)

	###### dump file information to file
	info_file2 = dir_sample + '/info_files.txt'
	cluster.to_csv(info_file2)

	###### timestamp
	filename_stamp = dir_sample + '/.success'
	stamp =	HCGB_time.print_time_stamp(filename_stamp)

	return()

