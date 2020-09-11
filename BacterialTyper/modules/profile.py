#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Generates a resistance and virulence profile for each sample.
"""
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
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.scripts import sampleParser
from BacterialTyper.scripts import virulence_resistance
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import card_trick_caller
from BacterialTyper.modules import help_info
from BacterialTyper.scripts import database_user

####################################
def run_profile(options):
	
	## init time
	start_time_total = time.time()
	
	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		sampleParser.help_format()
		exit()

	if (options.help_project):
		## information for project
		help_info.project_help()
		exit()
	
	if (options.help_ARIBA):
		## help_format option
		ariba_caller.help_ARIBA()
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

	## message header
	functions.pipeline_header()
	functions.boxymcboxface("Virulence & Resistance profile module")
	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## absolute path for in & out
	options.database = os.path.abspath(options.database)
	global input_dir
	input_dir = os.path.abspath(options.input)
	outdir=""

	## set mode: project/detached
	global Project
	if (options.detached):
		options.project = False
		outdir = os.path.abspath(options.output_folder)
		Project=False
	else:
		options.project = True
		outdir = input_dir	
		Project=True
	
	## get files
	pd_samples_retrieved = sampleParser.get_files(options, input_dir, "trim", ['_trim'])
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
	
	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		functions.create_folder(outdir)
	## for each sample
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "profile")
	
	###
	print ("+ Generate a sample profile for virulence and resistance candidate genes for each sample retrieved using:")
	print ("(1) Antimicrobial Resistance Inference By Assembly (ARIBA) software")
	print ("(2) Pre-defined databases by different suppliers or user-defined databases.")
	
	## get databases to check
	retrieve_databases = get_options_db(options)
	
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_total)
		
	########
	ARIBA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, start_time_partial)
	
	######################################
	## update database for later usage
	######################################
	if not options.fast:
		## functions.timestamp
		start_time_partial = functions.timestamp(start_time_partial)

		functions.boxymcboxface("Update Sample Database")

		## update db
		print ("+ Update database with samples identified")
		## TODO: check if it works
		dataBase_user = database_user.update_database_user_data(options.database, input_dir, Debug, options)		
		
		## debug message
		if (Debug):
			print (colored("**DEBUG: results obtained **", 'yellow'))
	
	else:
		print ("+ No update of the database has been requested using option --fast")

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Virulence & Resistance profile module.")
	return()

####################################
def get_options_db(options):
	##
	## Among all databases available and according to the input options,
	## select the databases to use and set dataframe with this information
	##
	
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
def ARIBA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, start_time_partial):
	functions.boxymcboxface("ARIBA Identification")

	##################
	## check status	##	
	##################
	databases2use = []
	card_trick_info = ""
	print ('+ Check databases status: ')
	for	index, db2use in retrieve_databases.iterrows():
		## index_name
		if (db2use['source'] == 'ARIBA'):
			index_status = ariba_caller.check_db_indexed(db2use['path'], 'YES')
			if (index_status == True):
				#print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
				databases2use.append(db2use['path'])
			
				## prepare card database ontology for later 		
				if (db2use['db'] == 'card'):
					card_trick_info = card_trick_caller.prepare_card_data(options.database)

		## check status of other databases if any
		# else:		

	## debug message
	if (Debug):
		print (colored("**DEBUG: databases2use\n" +  "\n".join(databases2use) + "\n**", 'yellow'))
		if (card_trick_info):
			print (colored("**DEBUG: card_trick_info: " + card_trick_info + " **", 'yellow'))
		
	######################################################
	## Start identification of samples
	######################################################
	print ("\n+ Send ARIBA identification jobs...")

	## get outdir folders
	outdir_samples = pd.DataFrame(columns=('sample', 'dirname', 'db', 'output'))

	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])

	for name, cluster in sample_frame:
		for db2use in databases2use:
			tmp = get_outfile(outdir_dict[name], name, db2use)
			outdir_samples.loc[len(outdir_samples)] = (name, outdir_dict[name], db2use, tmp)

	## multi-index
	outdir_samples = outdir_samples.set_index(['sample', 'db'])
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: outdir_samples **", 'yellow'))
		print (outdir_samples)	

	######################################################
	## send for each sample
	######################################################
	## ariba assembly cutoff 
	if not (options.ARIBA_cutoff):
		options.ARIBA_cutoff = 0.90

	## optimize threads
	name_list = set(pd_samples_retrieved["name"].tolist())
	threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

	## loop
	results_df = pd.DataFrame()
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		for db2use in databases2use:
			## send for each sample
			commandsSent = { executor.submit(ariba_run_caller, db2use, sorted(cluster["sample"].tolist()), outdir_samples.loc[(name, db2use), 'output'], threads_job, options.ARIBA_cutoff): name for name, cluster in sample_frame }
				
			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))
			
			print ("+ Jobs finished for database %s ..." %db2use)

			## functions.timestamp
			start_time_partial = functions.timestamp(start_time_partial)
			
			print ("+ Collecting information for each sample analyzed:")
			## check results for each database
			results_df_tmp = virulence_resistance.check_results(db2use, outdir_samples, options.ARIBA_cutoff, card_trick_info)
			results_df = pd.concat([results_df, results_df_tmp])
						
			## functions.timestamp
			start_time_partial = functions.timestamp(start_time_partial)

	######################################################
	## Generate final report for all samples
	######################################################
	## ariba summary results all samples
	print ("\n + Generate a summary file for all samples and one for each database employed...")

	## parse results
	if Project:
		final_dir = input_dir + '/report/profile'
		functions.create_folder(final_dir) 
	else:
		final_dir = os.path.abspath(options.output_folder)

	##
	vfdb = False
	subfolder = functions.create_subfolder("ariba_summary", final_dir)
	## subfolder_samples = functions.create_subfolder("samples", final_dir) ## TODO: Copy all xlsx files to a common folder. Is it necessary?

	## open excel writer
	name_excel = final_dir + '/profile_summary.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') 		

	for database, data in outdir_samples.groupby(level='db'): ## fix
		report_files_databases = {}

		for sample, data2 in data.groupby(level='sample'): ## fix
			file_report = data2.loc[sample, database]['output'] + '/report.tsv'
			if os.path.isfile(file_report): ## check if exists
				report_files_databases[sample] = file_report

		outfile_summary = subfolder + "/"			
		if database.endswith('card_prepareref/'):
			outfile_summary = outfile_summary + 'CARD_summary'
			name_db = 'CARD'
		elif database.endswith('vfdb_full_prepareref/'):
			outfile_summary = outfile_summary + 'VFDB_summary'
			name_db = 'VFDB'
			vfdb=True
		else:
			## TODO: check if there are multiple 'other' databases
			## Different databases provided (different to VFDB and CARD) would collapse file
			outfile_summary = outfile_summary + 'Other_summary' 
			name_db = 'other' 
			
		## call ariba summary to summarize results
		csv_all = ariba_caller.ariba_summary_all(outfile_summary, report_files_databases)
		if not csv_all == 'NaN':
			csv2excel = pd.read_csv(csv_all, header=0, sep=',')
			## write excel
			name_tab = name_db + '_found'
			csv2excel.to_excel(writer, sheet_name=name_tab)
	
	## results_df contains excel and csv files for each sample and for each database
	list_databases = set(results_df['database'].to_list())
	for db in list_databases: 
		df_db = results_df[results_df['database'] == db]['csv']
		dict_samples = df_db.to_dict()	
		
		merge_df = pd.DataFrame()
		for sample in dict_samples:

			if os.path.isfile(dict_samples[sample]):
				df = pd.read_csv(dict_samples[sample], header=0, sep=",")
				df = df.set_index('Genes')
				df2 = df.rename(columns={'Status':sample}, inplace=True)
				df2 = df[[sample]]

				## add to a common dataframe
				merge_df = pd.concat([merge_df, df2], axis=1, sort=True)
				merge_df.fillna("NaN", inplace=True)
		
		trans_df = merge_df.transpose()
		## write excel
		name_tab = db + '_all'
		trans_df.to_excel(writer, sheet_name=name_tab)
	
	## close
	writer.save()

	######################################################
	## print additional information for VFDB
	######################################################
	if (vfdb):
		print ("\n\n")
		functions.print_sepLine("*", 50, False)
		print ("+ Check VFDB details in files downloaded from vfdb website:")
		files_VFDB = virulence_resistance.check_VFDB(final_dir + '/VFDB_information')
		functions.print_sepLine("*", 50, False)

	######################################################
	print ("\n+ Please check additional summary files generated at folder ", final_dir)
	print ("+ Go to website: https://jameshadfield.github.io/phandango/#/")
	print ("+ For each database upload files *phandango.csv and *phandango.tre and visualize results")
	## print additional help?
	
####################################
def ariba_run_caller(db2use, list_files, folder_out, threads, cutoff):
	## check if already is done
	# generate a stamp when finish parsing each file

	## make stamp time
	filename_stamp = folder_out + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		files_names = [os.path.basename(s) for s in list_files]
		print (colored("\tA previous command generated results on: %s [Files: %s]" %(stamp, files_names), 'yellow'))
	
	else:
		if os.path.exists(folder_out):
			shutil.rmtree(folder_out) ## delete folder if exists but failed before

		## call
		code = ariba_caller.ariba_run(db2use, list_files, folder_out, threads, cutoff)
		if code == 'FAIL':
			print ("*** ERROR: System call failed for ", folder_out)		

####################################
def get_outfile(output_dir, name, index_name):
	##
	##
	
	basename_tag = index_name.split("_prepareref/")[0]
	basename = os.path.basename(basename_tag)

	if Project:
		out_file = output_dir + '/' + basename	
	else:
		output_path = functions.create_subfolder(name, output_dir)
		out_file = output_path + '/' + name + '_' + basename	
	
	## message debug
	if (Debug):
		print (colored("**DEBUG: Input names " +  name + '\n' + output_dir + '\n' + index_name + "\n", 'yellow'))
		print (colored("**DEBUG: Output names \n" +  basename + '\n' + basename_tag + '\n' +  out_file + " **\n", 'yellow'))

	return(out_file)

