#!/usr/bin/env python3
'''
This code calls species_identification_KMA and get the most similar taxa.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import sampleParser
from BacterialTyper import species_identification_KMA
from BacterialTyper import database_generator
from BacterialTyper import MLSTar
from BacterialTyper import edirect_caller
from BacterialTyper.modules import sample_prepare
from BacterialTyper.modules import database


####################################
def run(options):

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		sampleParser.help_format()
		exit()

	elif (options.help_project):
		## information for project
		help.project_help()
		exit()
	
	elif (options.help_KMA):
		## information for KMA Software
		species_identification_KMA.help_kma_database()
		exit()

	elif (options.help_MLSTar):
		## information for KMA Software
		MLSTar.help_MLSTar()
		exit()

	## init time
	start_time_total = time.time()

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

	### species_identification_KMA -> most similar taxa
	functions.pipeline_header()
	functions.boxymcboxface("Species identification")

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
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "trim", ['_trim_'])
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
	
	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		functions.create_folder(outdir)
	## for each sample
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "ident")	
	
	## let's start the process
	print ("+ Generate an species typification for each sample retrieved using:")
	print ("(1) Kmer alignment (KMA) software.")	
	print ("(2) Pre-defined databases by KMA or user-defined databases.")	
		
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
	
	######## KMA identification
	dataFrame_kma = KMA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, start_time_partial)
	
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)
	
	######## EDirect identification
	dataFrame_edirect = edirect_ident(dataFrame_kma, outdir_dict)
	
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)
	
	######## MLST identification
	MLST_results = MLST_ident(options, dataFrame_kma, outdir_dict, dataFrame_edirect)

	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)
	
	## generate summary for sample: all databases
	## MLST, plasmids, genome, etc
	functions.boxymcboxface("Results Summary")
	##

	#####################################
	## Summary identification results  ##
	#####################################

	## debug message
	if (Debug):
		print (colored("**DEBUG: retrieve results to summarize **", 'yellow'))
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)
		print (dataFrame_kma)
		print (dataFrame_edirect)
		print (MLST_results)

	## parse results
	if Project:
		final_dir = outdir + '/report/ident'
		functions.create_folder(final_dir) 
	else:
		final_dir = outdir

	excel_folder = functions.create_subfolder("samples", final_dir)
	
	# Group dataframe results summary by sample name
	sample_results_summary = dataFrame_kma.groupby(["Sample"])

	## debug message
	if (Debug):
		print (colored("**DEBUG: sample_results_summary **", 'yellow'))
		print (sample_results_summary)
	
	##
	results_summary_KMA = pd.DataFrame()
	MLST_all = pd.DataFrame()
	for name, grouped in sample_results_summary:
		
		## create a excel and txt for sample
		name_sample_excel = excel_folder + '/' + name + '_ident.xlsx'
		name_sample_csv = outdir_dict[name] + '/ident_summary.csv' ## check in detached mode

		writer_sample = pd.ExcelWriter(name_sample_excel, engine='xlsxwriter') ## open excel handle
		
		## subset dataframe	& print result
		results_summary_toPrint_sample = grouped[['Sample','#Template','Query_Coverage','Template_Coverage','Depth', 'Database']] 
		results_summary_toPrint_sample.to_excel(writer_sample, sheet_name="KMA") ## write excel handle
		results_summary_toPrint_sample.to_csv(name_sample_csv) ## write csv for sample
		
		## read MLST
		if MLST_results:
			sample_MLST = pd.read_csv(MLST_results[name], header=0, sep=',')
			print (sample_MLST)

			sample_MLST['genus'] = dataFrame_edirect.loc[dataFrame_edirect['sample'] == name, 'genus'].values[0]
			sample_MLST['species'] = dataFrame_edirect.loc[dataFrame_edirect['sample'] == name, 'species'].values[0]
			sample_MLST.to_excel(writer_sample, sheet_name="MLST") ## write excel handle
		
			## Return information to excel
			MLST_all = pd.concat([MLST_all, sample_MLST]) 
		
		## close excel handle
		writer_sample.save() 		

	##
	name_excel = final_dir + '/identification_summary.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') ## open excel handle

	## KMA dataframe: print result for sources
	results_summary_KMA = dataFrame_kma[['Sample','#Template','Query_Coverage','Template_Coverage','Depth', 'Database']] 
	
	## Sum plasmid and chromosome statistics ##
	## sum coverage
	total_coverage = results_summary_KMA.groupby('Sample')['Query_Coverage'].sum().reset_index()
	
	## debug message
	if (Debug):
		print ("*** Sum: Query_coverage ***")		
		print (total_coverage)
	
	## TODO: FIX SUMMARY REPORT
	results_summary_KMA = results_summary_KMA.set_index('Sample')
	results_summary_KMA = results_summary_KMA.sort_values(by=['Sample', 'Database', 'Query_Coverage'],ascending=[True, True,True])
	results_summary_KMA.to_excel(writer, sheet_name='KMA') ## write excel handle
	
	##
	MLST_all.to_excel(writer, sheet_name='MLST')

	## close excel handle
	writer.save()

	print ("\n+ Check summary of results in file generated" )		
	
	### timestamp
	start_time_partial = functions.timestamp(start_time_partial)					
	
	######################################
	## update database for later usage
	######################################
	if not options.fast:

		functions.boxymcboxface("Update Sample Database")

		## update db
		print ("+ Update database with samples identified")		

		## debug message
		if (Debug):
			print (colored("**DEBUG: dataFrame_edirect **", 'yellow'))
			pd.set_option('display.max_colwidth', -1)
			pd.set_option('display.max_columns', None)
			print (dataFrame_edirect)

		## Download
		print (dataFrame_edirect)
		
		file_toprint = final_dir + '/edirect_info2download.csv'
		dataFrame_edirect.to_csv(file_toprint)

		## dataFrame_edirect
		## assembly, annotation, etc...
		## rerun identification with new updated database
	else:
		print ("+ No update of the database has been requested using option --fast")
		
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting identification module.")
	exit()

####################################
def KMA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, time_partial):
	functions.boxymcboxface("KMA Identification")
	
	## set defaults
	kma_bin = config.get_exe("kma")
	
	## check status
	databases2use = []
	for	index, db2use in retrieve_databases.iterrows():
		## index_name
		if (db2use['source'] == 'KMA'):
			print ('+ Check database: ' + db2use['db'])
			index_status = species_identification_KMA.check_db_indexed(db2use['path'] )
			if (index_status == True):
				print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
				databases2use.append(db2use['path'])
			else:
				#databases2use.remove(db2use)
				print (colored("\t**Databases %s is not correctly indexed. Not using it...\n" % db2use['db'], 'red'))

	## debug message
	if (Debug):
		print (colored("**DEBUG: databases2use\n" +  "\n".join(databases2use) + "\n**", 'yellow'))

	## Start identification of samples
	print ("\n+ Send KMA identification jobs...")

	## optimize threads
	name_list = set(pd_samples_retrieved["name"].tolist())
	threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])
	
	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		for db2use in databases2use:
			
			## load database on memory
			print ("+ Loading database on memory for faster identification.")
			cmd_load_db = "%s shm -t_db %s -shmLvl 1" %(kma_bin, db2use)
			return_code_load = functions.system_call(cmd_load_db)
			
			## send for each sample
			commandsSent = { executor.submit(send_kma_job, outdir_dict[name], sorted(cluster["sample"].tolist()), name, db2use, threads_job, cluster): name for name, cluster in sample_frame }

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
		
			## functions.timestamp
			time_partial = functions.timestamp(time_partial)
		
	###
	print ("+ KMA identification call finished for all samples...")
	print ("+ Parse results now")
		
	## parse results		
	results_summary = pd.DataFrame()
	for db2use in databases2use:
		basename_db = os.path.basename(db2use)
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)

		for name, cluster in sample_frame:
		
			## get result
			## outdir_KMA
			outdir_dict_kma = functions.outdir_subproject(outdir_dict[name], cluster, "kma")
			result = get_outfile(outdir_dict_kma[name], name, db2use)
			#print ('\t- File: ' + result + '.spa')

			## get results using a cutoff value [Defaulta: 80]
			results = species_identification_KMA.parse_kma_results(result + '.spa', options.KMA_cutoff)
			results['Database'] = basename_db

			### check if db2use is plasmids as it could be several.
			if (results.index.size > 1):
				if (basename_db == "plasmids.T"):
					## let it be several entries
					results['Sample'] = name
					results_summary = results_summary.append(results)
				else:
					print (colored("###########################################", 'yellow'))
					print (colored("Sample %s contains multiple strains." %name, 'yellow'))
					print (colored("###########################################", 'yellow'))
					print (colored(results, 'yellow'))
					print ('\n\n')
					
					## add both strains if detected	
					results['Sample'] = name
					results_summary = results_summary.append(results)
					
					## TODO: add multi-isolate flag
		
			elif (results.index.size == 1): ## 1 clear reference
				results['Sample'] = name
				results_summary = results_summary.append(results)
		
			else:
				print (colored('\tNo clear strain from database %s has been assigned to sample %s' %(basename_db, name), 'yellow'))
				## add empty line if no available
				results['Sample'] = name
				results_summary = results_summary.append(results)
	
	print ("+ Finish this step...")
	return (results_summary)

####################################
def send_kma_job(outdir_file, list_files, name, database, threads, dataFrame_sample):

	## outdir_KMA
	outdir_dict_kma = functions.outdir_subproject(outdir_file, dataFrame_sample, "kma")

	## get outfile
	outfile = get_outfile(outdir_dict_kma[name], name, database)

	## check if previously run and succeeded
	basename_tag = os.path.basename(outfile)
	filename_stamp = outdir_dict_kma[name] + '/.success_' + basename_tag
	
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
	else:

		## debug message
		if (Debug):
			print (colored("**DEBUG: species_identification_KMA.kma_ident_module call**", 'yellow'))
			print ("outfile = get_outfile(outdir_dict_kma[name], name, db2use)")
			print ("outfile: ", outfile)
			print ("species_identification_KMA.kma_ident_module(outfile, list_files, name, database, threads) ")
			print ("species_identification_KMA.kma_ident_module" + "\t" + outfile + "\t" + str(list_files) + "\t" + name + "\t" + database + "\t" + str(threads) + "\n") 
	
		# Call KMA
		species_identification_KMA.kma_ident_module(outfile, list_files, name, database, threads) 

####################################
def get_outfile(output_dir, name, index_name):
	basename_tag = os.path.basename(index_name)
	if Project:
		output_path = output_dir
	else:
		output_path = functions.create_subfolder(name, output_dir)
		
	out_file = output_path + '/' + name + '_' + basename_tag	
	return(out_file)	
	
####################################
def edirect_ident(dataFrame, outdir_dict):
	
	## edirect	
	functions.boxymcboxface("EDirect information")
	print ("+ Connect to NCBI to get information from samples identified...")

	# Group dataframe sample name
	sample_results = dataFrame.groupby(["Sample"])

	## create dataframe to return results
	edirect_frame = pd.DataFrame(columns=("sample", "genus", "species", "strain", "BioSample", "genome", "Plasmids"))

	for name, grouped in sample_results:
		## use edirect to get Species_name and entry for later identification
		edirect_folder = functions.create_subfolder('edirect', outdir_dict[name])
		
		################################################
		## TODO: What to do if multi-isolate sample?
		################################################
		
		## chromosome match
		nucc_entry = grouped.loc[grouped['Database'] == 'bacteria.ATG']['#Template'].values[0].split()
		## e.g. NZ_CP029680.1 Staphylococcus aureus strain AR_0215 chromosome, complete genome

		##
		out_docsum_file = edirect_folder + '/nuccore_docsum.txt'
		species_outfile = edirect_folder + '/info.csv'
		filename_stamp = edirect_folder + '/.success_species'
				
		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))

		else: 
			edirect_caller.generate_docsum_call('nuccore', nucc_entry[0], out_docsum_file)
			edirect_caller.generate_xtract_call(out_docsum_file, 'DocumentSummary', 'Organism,BioSample,Strain', species_outfile)
			stamp =	functions.print_time_stamp(filename_stamp)

		## get information from edirect call
		taxa_name_tmp = functions.get_info_file(species_outfile)
		Organism = taxa_name_tmp[0].split(',')[0].split()
		genus = Organism[0] 		## genus
		species = Organism[1] 		## species
		BioSample_name = taxa_name_tmp[0].split(',')[1]	## BioSample
		
		## sometimes strain is missing
		if len(taxa_name_tmp[0].split(',')) > 2:
			strain = taxa_name_tmp[0].split(',')[2] 		## strain
		else:
			strain = 'NaN'

		## plasmid match
		group_plasmid = grouped.loc[grouped['Database'] == 'plasmids.T' ]
		plasmid_entries = group_plasmid['#Template'].tolist()
			## e.g. NZ_CP029083.1 Staphylococcus aureus strain AR464 plasmid unnamed1, complete sequence
		plasmid_entries_str = ",".join([i.split()[0] for i in plasmid_entries])

		#("sample", "taxa", strain, genome "BioSample", "Plasmids"))
		edirect_frame.loc[len(edirect_frame)] = (name, genus, species, strain, BioSample_name, nucc_entry[0], plasmid_entries_str)
	
	return (edirect_frame)

####################################
def MLST_ident(options, dataFrame, outdir_dict, dataFrame_edirect):
	## Add MLST Information

	## debug message
	if (Debug):
		print (colored("**DEBUG: dataFrame_edirect identified**", 'yellow'))
		print (dataFrame_edirect)

	## MLST call	
	functions.boxymcboxface("MLST typing")
	print ("+ Create classical MLST typification of each sample according to species retrieved by kmer...")

	## get assembly files
	input_dir = os.path.abspath(options.input)
	assembly_samples_retrieved = sample_prepare.get_files(options, input_dir, "assembly", "fna")

	## TODO: Samples might not be assembled...to take into account

	## debug message
	if (Debug):
		print (colored("**DEBUG: assembly_samples_retrieved**", 'yellow'))
		print (assembly_samples_retrieved)	
	
	## PubMLST folder
	path_database = os.path.abspath(options.database)
	pubmlst_folder = functions.create_subfolder('PubMLST', path_database)
	
	########################################################################################
	## TODO: Fix and install MLSTar during installation
	rscript = "/soft/general/R-3.5.1-bioc-3.8/bin/Rscript" ##config.get_exe("Rscript") 
	########################################################################################

	# init
	MLST_results = {}

	## Generate MLST call according to species identified for each sample
	## TODO: What to do if multi-isolate sample?
	
	for index, row in dataFrame_edirect.iterrows():
		MLSTar_taxa_name = MLSTar.get_MLSTar_species(row['genus'], row['species'] )
		
		if (MLSTar_taxa_name == 'NaN'):
			continue
		
		## species folder
		species_mlst_folder = functions.create_subfolder(MLSTar_taxa_name, pubmlst_folder)

		## output file
		output_file = species_mlst_folder + '/PubMLST_available_scheme.csv'
		filename_stamp = species_mlst_folder + '/.success_scheme'
				
		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
			# [TODO: Check time passed and download again if >?? days passed]
			
		else: 
			### get scheme available
			MLSTar.getPUBMLST(MLSTar_taxa_name, rscript, output_file)
			stamp =	functions.print_time_stamp(filename_stamp)

		## parse and get scheme for classical MLST
		schemes_MLST = pd.read_csv(output_file, sep=',', header=0)
		
		##
		for item, cluster in schemes_MLST.iterrows():
			if cluster['len'] < 10:
				scheme2use = int(cluster['scheme'])
				continue			
		
		### 
		sample = row['sample']
		MLSTar_folder = functions.create_subfolder('MLST', outdir_dict[sample])
		genome_file = assembly_samples_retrieved.loc[assembly_samples_retrieved['name'] == sample]['sample'].values[0]

		## call MLST
		(results, profile_folder) = MLSTar.run_MLSTar(path_database, rscript, MLSTar_taxa_name, scheme2use, sample, MLSTar_folder, genome_file, options.threads)
		MLST_results[sample] = results

	##	
	print ("+ Finish this step...")
	return (MLST_results)

####################################
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
		print ('\t- Selecting kma databases:')
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
		print ('\t- Get additional kma databases:')
		if (options.kma_external_files):
			options.kma_external_files = set(options.kma_external_files)		

			## check if indexed and/or index if necessary
			external_kma_dbs_list = external_kma(options.kma_external_files) ## fix
			external_kma_dbs_string = ','.join(external_kma_dbs_list)
			option_db = "kma_external:" + external_kma_dbs_string

		else:
			## rise error & exit
			print (colored("***ERROR: No database provided via --kma_external_file option.\n",'red'))
			exit()
			
		## rise attention
		if (options.kma_dbs):
			print (colored("***ATTENTION:\nDefatult databases and databases provided via --kma_dbs option would not be used as --only_external_kma option provided.\n",'red'))

	############
	## 5) all databases 
	############
	elif (options.all_data):
	
		############
		## default dbs + user kma dbs 
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
			print ('\n\n\t- Get additional kma databases:')
			options.kma_external_files = set(options.kma_external_files)		
			
			## check if indexed and/or index if necessary
			external_kma_dbs_list = external_kma(options.kma_external_files) ## fix
			external_kma_dbs_string = ','.join(external_kma_dbs_list)
			option_db = option_db + "#kma_external:" + external_kma_dbs_string

		############
		## Previously identified data
		############
		#if (options.user_data):
		option_db = option_db + '#user_data'
		
		############
		## Genbank reference data
		############
		#if (options.genbank_data):
		option_db = option_db + '#genbank'

	## debug message
	if (Debug):
		print (colored("**DEBUG: option_db: " +  option_db + " **", 'yellow'))
	
	### get dbs	
	return (database.getdbs("KMA", database2use, option_db, Debug))

####################################
def external_kma(list_files):
	
	## set defaults
	kma_bin = config.get_exe("kma")
	
	external_file_returned = []
	for f in list_files:
		f = os.path.abspath(f)
		print (colored('\t+ %s' %f, 'green'))
		status = species_identification_KMA.check_db_indexed(f)
		if (status): #true
			external_file_returned.append(f)
			## debug message
			if (Debug):
				print (colored("**DEBUG: Database (%s) is indexed" %f + " **", 'yellow'))
			
		else: #false
			## debug message
			if (Debug):
				print (colored("**DEBUG: Database (%s) is not indexed" %f + " **", 'yellow'))

			dataBase = os.path.abspath(f)
			basename_name = os.path.basename(dataBase)
			
			## debug message
			if (Debug):
				print (colored("**DEBUG: dataBase " + dataBase + " **", 'yellow'))
				print (colored("**DEBUG: dataBase " + dataBase + " **", 'yellow'))
			
			status = species_identification_KMA.index_database(dataBase, kma_bin, dataBase, "new")
			if (status): #true
				external_file_returned.append(basename_name)
			else:
				print (colored("***ERROR: Database provided via --kma_external_file option (%s) was not indexed.\n" %f,'orange'))

	return(external_file_returned)

