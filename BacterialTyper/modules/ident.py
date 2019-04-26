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
				commandsSent = { executor.submit(species_identification_KMA.kma_ident_module, get_outfile(outdir, row['samples'], db2use), (row['R1'], row['R2']), row['samples'], db2use, cpu_here): index for index, row in pd_samples_retrieved.iterrows() }
	
			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))
	
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
	## according to user input: select databases to use
	print ("\n\n+ Select databases to use for identification:")
	database_folder = config.DATA["database"] ## default: set during configuration
	retrieve_databases = species_identification_KMA.getdbs(database_folder)

	if (options.database):
		## 
		print ("- User provides a database folder. Checking...")
		## check databases integrity
		
		if (options.both_folder):
			## use default and user provided
			print ("- Use both databases: default database and folder user provided.")
		
	
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
