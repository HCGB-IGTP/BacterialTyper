#!/usr/bin/env python3
'''
This code calls species_identification_KMA and get the most similar taxa then use ariba_caller to check if pubmlst is downloaded and  get MLST profile.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
from io import open
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import species_identification_KMA
from BacterialTyper import ariba_caller
from BacterialTyper.modules import sample_prepare

kma_bin = config.EXECUTABLES["kma"]

####################################
def ARIBA_ident(options, cpu, pd_samples_retrieved, outdir, retrieve_databases):
	functions.boxymcboxface("ARIBA Identification")

####################################
def KMA_ident(options, cpu, pd_samples_retrieved, outdir, retrieve_databases):
	functions.boxymcboxface("KMA Identification")

	## check status
	databases2use = []
	for	index, db2use in retrieve_databases.iterrows():
		## index_name
		if (db2use['source'] == 'KMA_db'):
			print ('+ Check database: ' + db2use['db'])
			index_status = species_identification_KMA.check_db_indexed(db2use['path'] )
			if (index_status == True):
				print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
				databases2use.append(db2use['path'])
			else:
				#databases2use.remove(db2use)
				print (colored("\t**Databases %s is not correctly indexed. Not using it...\n" % db2use['db'], 'red'))
	
	print ("\n+ Send KMA identification jobs...")

	if (options.pair):
		
		cpu_here = int(cpu/len(databases2use))
		if (cpu_here == 0):
			cpu_here = 1
		
		# We can use a with statement to ensure threads are cleaned up promptly
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.threads)) as executor:
			for db2use in databases2use:
				
				## load database on memory
				print ("+ Loading database on memory for faster identification.")
				cmd_load_db = "%s shm -t_db %s -shmLvl 1" %(kma_bin, db2use)
				return_code_load = functions.system_call(cmd_load_db)
				
				## send for each sample
				commandsSent = { executor.submit(species_identification_KMA.kma_ident_module, get_outfile(outdir, row['samples'], db2use), (row['R1'], row['R2']), row['samples'], db2use, cpu_here): index for index, row in pd_samples_retrieved.iterrows() }
	
				for cmd2 in concurrent.futures.as_completed(commandsSent):
					details = commandsSent[cmd2]
					try:
						data = cmd2.result()
					except Exception as exc:
						print ('***ERROR:')
						print (cmd2)
						print('%r generated an exception: %s' % (details, exc))
	
				## remove database from memory
				print ("+ Removing database from memory...")
				cmd_rm_db = "%s shm -t_db %s -shmLvl 1 -destroy" %(kma_bin, db2use)
				return_code_rm = functions.system_call(cmd_rm_db)
				if (return_code_rm == 'FAIL'):
					print (colored("***ERROR: Removing database from memory failed. Please do it! Execute command: %s" %cmd_rm_db,'red'))
	
	else:
		## to do: implement single end mode
		print ('+ No implementation yet. Sorry.')
		exit()
		
	###
	print ("+ KMA identification call finished for all samples...")
	print ("+ Parse results now:")
		
	## parse results
	name_excel = outdir + '/identification_summary.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') ## open excel handle
	
	for db2use in databases2use:
		basename_db = os.path.basename(db2use)
		results_summary = pd.DataFrame()
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)

		for index, row in pd_samples_retrieved.iterrows():
			result = get_outfile(outdir, row['samples'], db2use)
			#print ('\t- File: ' + result + '.spa')
			results = species_identification_KMA.parse_kma_results(row['samples'], result + '.spa')


			### Fix: check if db2use is plasmids as it could be several.
			if (results.index.size > 1):
				print (colored("Sample %s contains multiple strains." %row['samples'], 'yellow'))
				print (colored(results.to_csv, 'yellow'))

			elif (results.index.size == 1):
				results['sample'] = row['samples']
				results_summary = results_summary.append(results)
			else:
				print (colored('\tNo clear strain from database %s has been assigned to sample %s' %(basename_db, row['samples']), 'yellow'))
				
		## subset dataframe	& print result
		results_summary_toPrint = results_summary[['sample','#Template','Query_Coverage','Template_Coverage','Depth']] 
		results_summary_toPrint = results_summary_toPrint.set_index('sample')		
		results_summary_toPrint.to_excel(writer, sheet_name=basename_db) ## write excel handle
	
	writer.save() ## close excel handle	
	return (name_excel, results_summary)

####################################
def get_outfile(output_dir, name, index_name):
	basename_tag = os.path.basename(index_name)
	output_path = functions.create_subfolder(name, output_dir)
	out_file = output_path + '/' + name + '_' + basename_tag	
	return(out_file)
	
####################################
def external_kma(list_files):
	
	external_file_returned = []
	for f in list_files:
		print (f)
		status = species_identification_KMA.check_db_indexed(f)
		if (status): #true
			external_file_returned.append(f)
		else: #false
			dataBase = os.path.abspath(f)
			basename_name = os.path.basename(dataBase)
			status = species_identification_KMA.index_database(dataBase, kma_bin, basename_name, "new")
			if (status): #true
				external_file_returned.append(basename_name)
			else:
				print (colored("***ERROR: Database provided via --kma_external_file option (%s) was not indexed.\n" %f,'orange'))

	return(external_file_returned)

####################################
def get_options_db(options):
	print ("\n\n+ Select databases to use for identification:")
	default_database_folder = config.DATA["database"] ## default: set during configuration
	
	## according to user input: select databases to use
	option_db = ""
	
	## Default
	kma_dbs = ["bacteria", "plasmids"]
	
	############
	## 1) only user data: previously identified and added
	############
	if (options.only_user_data):
		option_db = "user_data"
		
	############
	## 2) only genbank data: previously download from NCBI reference genomes
	############
	elif (options.only_genbank_data):
		option_db = "genbank"
	
	############
	## 3) only kma_db
	############
	elif (options.only_kma_db):
		if (options.kma_dbs):
			options.kma_dbs = set(options.kma_dbs)
			kma_dbs_string = ','.join(options.kma_dbs)
			option_db = "kma:" + kma_dbs_string

		else:
			## rise error & exit
			print (colored("***ERROR: No database provided via --kma_db option.\n",'red'))
			exit()

	############
	## 4) only external kma
	############
	elif (options.only_external_kma):
		if (options.kma_external_file):
			options.kma_external_file = set(options.kma_external_file)		

			## check if indexed and/or index if necessary
			external_kma_dbs_string = external_kma(options.kma_external_files)
			option_db = "kma_external:" + external_kma_dbs_string

		else:
			## rise error & exit
			print (colored("***ERROR: No database provided via --kma_external_file option.\n",'red'))
			exit()

	############
	## 5) all databases 
	############
	else:
	
		############
		## default dbs + user
		############
		
		print ('\t- Selecting kma databases:')
		if (options.kma_dbs):
			options.kma_dbs = options.kma_dbs + kma_dbs
			options.kma_dbs = set(options.kma_dbs)		
		else:
			options.kma_dbs = kma_dbs

		kma_dbs_string = ','.join(options.kma_dbs)
		option_db = "kma:" + kma_dbs_string
		
		for i in options.kma_dbs:
			print (colored('\t\t+ %s' %i, 'green'))

		
		############
		## External file
		############
		if (options.kma_external_files):
			print ('\t- Get additional kma databases:')
			options.kma_external_file = set(options.kma_external_files)		
			
			## check if indexed and/or index if necessary
			external_kma_dbs_list = external_kma(options.kma_external_files)
			external_kma_dbs_string = ','.join(external_kma_dbs_list)

			option_db = option_db + "#kma_external:" + external_kma_dbs_string

			for i in external_kma_dbs_list:
				print (colored('\t\t+ %s' %i, 'green'))

		option_db = option_db + '#user_data'
		option_db = option_db + '#genbank'

	print (option_db)
	exit()
	return (option_db)
	
####################################
def run(options):

	### species_identification_KMA -> most similar taxa
	functions.pipeline_header()

	if (options.fast):
		functions.boxymcboxface("Fast species identification module")
	else:
		functions.boxymcboxface("Species identification")
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir)
	## generate output folder
	functions.create_folder(outdir)

	## optimize threads
	threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size)
	
	print ("+ Generate an species typification for each sample retrieved using two methods.")
	print ("(1) Kmer alignment (KMA) identification")	
	print ("(2) Antimicrobial Resistance Inference By Assembly (ARIBA) identification\n\n")	
	
	## get databases to check
	option_db = get_options_db(options)
	
	### get dbs
	retrieve_databases = species_identification_KMA.getdbs(default_database_folder, option_db)

	########
	(excel_generated, dataFrame) = KMA_ident(options, threads_module, pd_samples_retrieved, outdir, retrieve_databases)
	
	## update database for later usage
	if (options.fast):
		## skip it
		print ("\n+ Check summary of results in file: " + excel_generated)		
	else:
		## update db
		print ("")
		## assembly, annotation, etc...
		## rerun identification with new updated database
	
	########
	ARIBA_ident(options, threads_module, pd_samples_retrieved, outdir, retrieve_databases)

	print ("+ Exiting identification module.")
	exit()
