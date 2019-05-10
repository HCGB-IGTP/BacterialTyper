#!/usr/bin/env python3
'''
This code calls ariba and virulence resistance script and generates a resistance and virulence profile for each sample.
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
import shutil

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import virulence_resistance
from BacterialTyper import database_generator
from BacterialTyper import ariba_caller
from BacterialTyper.modules import sample_prepare

####################################
def run(options):
	
	## init time
	start_time_total = time.time()
	
	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False

	## message header
	functions.pipeline_header()
	functions.boxymcboxface("Virulence & Resistance profile module")
	print ("--------- Starting Process ---------")
	functions.print_time()

	if (options.fast):
		print ("")
	else:
		print ("")

	## print further information for ARIBA databases	
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir)
	## generate output folder
	functions.create_folder(outdir)

	print ("+ Generate a sample profile for virulence and resistance candidate genes for each sample retrieved using:")
	print ("(1) Antimicrobial Resistance Inference By Assembly (ARIBA) software")
	print ("(2) Pre-defined databases by different suppliers or user-defined databases.")
	
	## get databases to check
	retrieve_databases = get_options_db(options)
	
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_total)
		
	########
	ARIBA_ident(options, pd_samples_retrieved, outdir, retrieve_databases)

	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Virulence & Resistance profile module.")
	exit()

####################################
def get_options_db(options):
	print ("\n\n+ Select databases to use for virulence and resistance profile generation:")
	
	### database folder to use
	database2use = options.database
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: Database to use: " +  database2use + " **", 'yellow'))
	
	## according to user input: select databases to use
	option_db = ""
	
	## Default
	ARIBA_dbs = ["CARD", "VFDB"]

	############
	## 1) only provided databases
	############
	if (options.no_def_ARIBA):
		ARIBA_dbs = set(options.ariba_dbs)
		
	############
	## 2) all
	############
	elif (options.ariba_dbs):
		ARIBA_dbs = ARIBA_dbs + options.ariba_dbs
		ARIBA_dbs = set(ARIBA_dbs)

	ARIBA_dbs_string = ','.join(ARIBA_dbs)
	option_db = "ARIBA:" + ARIBA_dbs_string

	## debug message
	if (Debug):
		print (colored("**DEBUG: Database to use: " +  option_db + " **", 'yellow'))
	
	### get dbs	
	return (database_generator.getdbs('ARIBA', database2use, option_db, Debug))


####################################
def ARIBA_ident(options, pd_samples_retrieved, outdir, retrieve_databases):
	functions.boxymcboxface("ARIBA Identification")

	## check status
	databases2use = []
	print ('+ Check databases status: ')
	for	index, db2use in retrieve_databases.iterrows():
		## index_name
		if (db2use['source'] == 'ARIBA'):
			
			index_status = ariba_caller.check_db_indexed(db2use['path'], 'YES')
			if (index_status == True):
				#print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
				databases2use.append(db2use['path'])

	## debug message
	if (Debug):
		print (colored("**DEBUG: databases2use\n" +  "\n".join(databases2use) + "\n**", 'yellow'))

	## Start identification of samples
	print ("\n+ Send ARIBA identification jobs...")

	## get outdir folders
	outdir_samples = pd.DataFrame(columns=('sample', 'db', 'output'))
	for index, row in pd_samples_retrieved.iterrows():
		for db2use in databases2use:
			tmp = get_outfile(outdir, row['samples'], db2use)
			outdir_samples.loc[len(outdir_samples)] = (row['samples'], db2use, tmp)

	## multi-index
	outdir_samples = outdir_samples.set_index(['sample', 'db'])
	
	if (options.pair):
		## message debug
		if (Debug):
			print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		
		#workers = (options.threads/2)
		workers = options.threads
		
		# We can use a with statement to ensure threads are cleaned up promptly
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(workers)) as executor:
			for db2use in databases2use:
				## send for each sample
				commandsSent = { executor.submit(ariba_run_caller, db2use, (row['R1'], row['R2']), outdir_samples.loc[(row['samples'], db2use), 'output'], 1): index for index, row in pd_samples_retrieved.iterrows() }
	
				for cmd2 in concurrent.futures.as_completed(commandsSent):
					details = commandsSent[cmd2]
					try:
						data = cmd2.result()
					except Exception as exc:
						print ('***ERROR:')
						print (cmd2)
						print('%r generated an exception: %s' % (details, exc))
				
				print ("+ Jobs finished for database %s\n+ Collecting information..." %db2use)
				virulence_resistance.check_results(db2use, outdir_samples, outdir)
				print ("")
				
	else:
		## to do: implement single end mode
		print ('+ No implementation yet. Sorry.')
		exit()
		
	## ariba summary results all samples
	print ("\n + Generate a summary file for all samples and one for each database employed...")

	subfolder = functions.create_subfolder("ariba_summary", outdir)
	for database, data in outdir_samples.groupby(level='db'): ## fix
		report_files_databases = {}
		for sample, data2 in data.groupby(level='sample'): ## fix
			report_files_databases[sample] = data2.loc[sample, database]['output'] + '/report.tsv'

		outfile_summary = subfolder + "/"			
		if database.endswith('card_prepareref/'):
			outfile_summary = outfile_summary + 'CARD_summary'
		elif database.endswith('vfdb_full_prepareref/'):
			outfile_summary = outfile_summary + 'VFDB_summary'
		else:
			outfile_summary = outfile_summary + 'Other_summary' ## todo: different databases provided different to VFDB and CARD would collapse file
			
		## call ariba summary
		ariba_caller.ariba_summary_all(outfile_summary, report_files_databases)
		
	print ("\n+ Please check additional summary files generated at folder ", subfolder)
	print ("+ Go to website: https://jameshadfield.github.io/phandango/#/")
	print ("+ For each database upload files *phandango.csv and *phandango.tre and visualize results")
	
	## print additional help?
	
####################################
def ariba_run_caller(db2use, list_files, folder_out, threads):
	## check if already is done
	# generate a stamp when finish parsing each file

	## make stamp time
	filename_stamp = folder_out + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print ("\tA previous command generated results on: ", stamp)
	else:
		if os.path.exists(folder_out):
			shutil.rmtree(folder_out)

		ariba_caller.ariba_run(db2use, list_files, folder_out, threads)
	

####################################
def get_outfile(output_dir, name, index_name):
	
	## message debug
	if (Debug):
		print (colored("**DEBUG: Input names " +  name + '\n' + output_dir + '\n' + index_name + " **\n", 'yellow'))

	basename_tag = index_name.split("_prepareref/")[0]
	basename = os.path.basename(basename_tag)

	output_path = functions.create_subfolder(name, output_dir)
	out_file = output_path + '/' + name + '_' + basename	
	
	## message debug
	if (Debug):
		print (colored("**DEBUG: Output names \n" +  basename + '\n' + basename_tag + '\n' + output_path + '\n' + out_file + " **\n", 'yellow'))

	return(out_file)

