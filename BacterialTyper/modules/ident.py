#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Generates sample identification using KMA software and MLSTar. 
Looks for similar entries on GenBank and retrieves them.  
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
from BacterialTyper.modules import info

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
		info.project_help()
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
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: retrieve results to summarize **", 'yellow'))
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)
		print ("dataframe_kma")
		print (dataFrame_kma)
	
	######## EDirect identification
	dataFrame_edirect = edirect_ident(dataFrame_kma, outdir_dict)
	
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)

	## debug message
	if (Debug):
		print (colored("**DEBUG: retrieve results from NCBI **", 'yellow'))
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)
		print ("dataFrame_edirect")
		print (dataFrame_edirect)
		
	######## MLST identification
	MLST_results = MLST_ident(options, dataFrame_kma, outdir_dict, dataFrame_edirect, retrieve_databases)

	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: retrieve results to summarize **", 'yellow'))
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)
		print ("MLST_results")
		print (MLST_results)

	## generate summary for sample: all databases
	## MLST, plasmids, genome, etc
	functions.boxymcboxface("Results Summary")
	
	#####################################
	## Summary identification results  ##
	#####################################

	## parse results
	if Project:
		final_dir = outdir + '/report/ident'
		functions.create_folder(final_dir) 
	else:
		final_dir = outdir

	###
	excel_folder = functions.create_subfolder("samples", final_dir)
	print ('+ Print summary results in folder: ', final_dir)
	print ('+ Print sample results in folder: ', excel_folder)
	
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
			sample_MLST['genus'] = dataFrame_edirect.loc[dataFrame_edirect['sample'] == name, 'genus'].values[0]
			sample_MLST['species'] = dataFrame_edirect.loc[dataFrame_edirect['sample'] == name, 'species'].values[0]
			sample_MLST.to_excel(writer_sample, sheet_name="MLST") ## write excel handle
		
			## Return information to excel
			MLST_all = pd.concat([MLST_all, sample_MLST]) 
		
		## close excel handle
		writer_sample.save() 		

	##
	name_excel = final_dir + '/identification_summary.xlsx'
	print ('+ Summary information in excel file: ', name_excel)
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
	
	## write excel and close
	MLST_all.to_excel(writer, sheet_name='MLST')
	writer.save() ## close excel handle

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

		## dataFrame_edirect
		file_toprint = final_dir + '/edirect_info2download.csv'
		dataFrame_edirect.to_csv(file_toprint)

		## update database with samples identified
		data2download = dataFrame_edirect.filter(['genus','species', 'strain', 'genome'])
		data2download = data2download.rename(columns={'genome': 'NCBI_assembly_ID', 'strain' : 'name'})
		NCBI_folder = os.path.abspath(options.database) + '/NCBI'
		database_generator.NCBI_DB(data2download, NCBI_folder, Debug)

	else:
		print ("+ No update of the database has been requested using option --fast")
		
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting identification module.")
	exit()

####################################
def KMA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, time_partial):
	"""Kmer identification
	
	Arguments:
		options: options passed to the run main function (threads, KMA_cutoff, etc)
		
		pd_samples_retrieved: pandas dataframe for samples to process.
				
		outdir_dict: dictionary containing  
		
		retrieve_databases: 
		
		time_partial: timestamp for the process when started
	
	Returns:
		results_summary: pandas dataframe containing identification information generated
	
	"""
	
	
	functions.boxymcboxface("KMA Identification")
	## set defaults
	kma_bin = config.get_exe("kma")	

	## check status
	databases2use = []
	for	index, db2use in retrieve_databases.iterrows():
		## index_name
		if (str(db2use['source']).startswith('KMA')):
			print ('+ Check database: ' + db2use['db'])
			fold_name = os.path.dirname(db2use['path'])
			
			index_status = species_identification_KMA.check_db_indexed(db2use['path'], fold_name )
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
		
	## parse results		
	print ("+ KMA identification call finished for all samples...")
	print ("+ Parse results now")
	results_summary = pd.DataFrame()
	for db2use in databases2use:
		### [TODO]: parse data according to database: bacteria, plasmids or user data or genbank data provided
		
		basename_db = os.path.basename(db2use)
		pd.set_option('display.max_colwidth', -1)
		pd.set_option('display.max_columns', None)

		###
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

###################################
def send_kma_job(outdir_file, list_files, name, database, threads, dataFrame_sample):

	## outdir_KMA
	outdir_dict_kma = functions.outdir_subproject(outdir_file, dataFrame_sample, "kma")

	## set defaults
	kma_bin = config.get_exe("kma")

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
	
		## Sparse or not
		#if any(name in basename_tag for name in ['userData_KMA', 'genbank_KMA']):
		if (basename_tag == 'userData_KMA'):
			option = ''
		else:
			option = '-Sparse'
	
		# Call KMA
		species_identification_KMA.kma_ident_call(outfile, list_files, name, database, kma_bin, option, threads)

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
	
	################################################
	## TODO: What to do if multi-isolate sample?
	################################################

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
		
		## chromosome match
		nucc_entry = grouped.loc[grouped['Database'] == 'bacteria.ATG']['#Template'].values[0].split()
		## e.g. NZ_CP029680.1 Staphylococcus aureus strain AR_0215 chromosome, complete genome

		##
		out_docsum_file = edirect_folder + '/nuccore_docsum.txt'
		tmp_species_outfile = edirect_folder + '/info.csv'
		filename_stamp = edirect_folder + '/.success_species'
				
		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))

		else: 
			edirect_caller.generate_docsum_call('nuccore', nucc_entry[0], out_docsum_file)
			edirect_caller.generate_xtract_call(out_docsum_file, 'DocumentSummary', 'Organism,BioSample,AssemblyAcc,Strain', tmp_species_outfile)
			
		########################################	
		## get information from edirect call
		########################################	
		taxa_name_tmp = functions.get_info_file(tmp_species_outfile)
		Organism = taxa_name_tmp[0].split(',')[0].split()
		genus = Organism[0] 							## genus
		species = Organism[1] 							## species
		BioSample_name = taxa_name_tmp[0].split(',')[1]	## BioSample
		AssemblyAcc = taxa_name_tmp[0].split(',')[2] 	## AssemblyAcc
		
		## sometimes strain is missing
		if len(taxa_name_tmp[0].split(',')) > 3:
			strain = taxa_name_tmp[0].split(',')[3] 	## strain
		else:
			strain = 'NaN'
		
		## get GenBank accession ID
		out_docsum_file_assembly = edirect_folder + '/assembly_docsum.txt'
		AssemblyAcc_outfile = edirect_folder + '/AssemblyAcc.csv'
		
		## some error ocurred
		if functions.is_non_zero_file(out_docsum_file_assembly):
			continue
		
		edirect_caller.generate_docsum_call('assembly', AssemblyAcc, out_docsum_file_assembly)
		edirect_caller.generate_xtract_call(out_docsum_file_assembly, 'DocumentSummary', 'Genbank', AssemblyAcc_outfile) 
		
		## Is it better to download Refseq or Genbank?
		## https://www.quora.com/What-is-the-difference-between-Refseq-and-Genbank		
		
		GenbankAcc = functions.get_info_file(AssemblyAcc_outfile)
		if Debug:
			print("Sample: ", name)
			print("Genbank Acc: ", GenbankAcc[0])
		
		## plasmid match
		group_plasmid = grouped.loc[grouped['Database'] == 'plasmids.T' ]
		plasmid_entries = group_plasmid['#Template'].tolist()
			## e.g. NZ_CP029083.1 Staphylococcus aureus strain AR464 plasmid unnamed1, complete sequence
		plasmid_entries_str = ",".join([i.split()[0] for i in plasmid_entries])

		## save edirect_frame
		#("sample", "taxa", strain, genome "BioSample", "Plasmids"))
		edirect_frame.loc[len(edirect_frame)] = (name, genus, species, strain, BioSample_name, GenbankAcc[0], plasmid_entries_str)

		stamp =	functions.print_time_stamp(filename_stamp)

	
	return (edirect_frame)

####################################
def MLST_ident(options, dataFrame, outdir_dict, dataFrame_edirect, retrieve_databases):
	## Add MLST Information

	## TODO: Samples might not be assembled...to take into account and return 0
	## TODO: Fix and install MLSTar during installation
	## TODO: What to do if multi-isolate sample?
	## TODO: Control if a different profile is provided via --MLST_profile
	## TODO: Check time passed and download again if >?? days passed]
	
	return ## skip by now
	
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

	## debug message
	if (Debug):
		print (colored("**DEBUG: assembly_samples_retrieved**", 'yellow'))
		print (assembly_samples_retrieved)	
	
	########################################################################################
	rscript = "/soft/general/R-3.5.1-bioc-3.8/bin/Rscript" ##config.get_exe("Rscript") 
	########################################################################################

	# init
	MLST_results = {}
		
	## get MLST_profile: default or provided
	mlst_profile = retrieve_databases.loc[ retrieve_databases['db'] == 'PubMLST']['path'].item()
	
	## Generate MLST call according to species identified for each sample
	for index, row in dataFrame_edirect.iterrows():
		MLSTar_taxa_name = MLSTar.get_MLSTar_species(row['genus'], row['species'] )
		
		if (MLSTar_taxa_name == 'NaN'):
			continue		
		## species folder
		#species_mlst_folder = functions.create_subfolder(MLSTar_taxa_name, pubmlst_folder)
		species_mlst = mlst_profile.split(',')[0]
		species_mlst_folder = mlst_profile.split(',')[1]
		
		## output file
		output_file = species_mlst_folder + '/PubMLST_available_scheme.csv'
		filename_stamp = species_mlst_folder + '/.success_scheme'
		
		## 
		if MLSTar_taxa_name == species_mlst:		
			if os.path.isfile(filename_stamp):
				stamp =	functions.read_time_stamp(filename_stamp)
				print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))

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
		(results, profile_folder) = MLSTar.run_MLSTar(species_mlst_folder, rscript, MLSTar_taxa_name, scheme2use, sample, MLSTar_folder, genome_file, options.threads)
		MLST_results[sample] = results

	##	
	print ("+ Finish this step...")
	return (MLST_results)

####################################
def get_kma_db(kma_dbs):
	print ('\t- Selecting kma databases:')
	kma_dbs_string = ','.join(kma_dbs)
	option_db = "kma:" + kma_dbs_string
	
	for i in kma_dbs:
		print (colored('\t\t+ %s' %i, 'green'))

	return(option_db)
	

####################################
def get_external_kma(kma_external_files, Debug):
	print ('\t- Get additional kma databases:')
	## external sequences provided are indexed and generated in the same folder provided 
	
	option_db = ""
	if (kma_external_files):
		kma_external_files = set(kma_external_files)		
		kma_external_files = [os.path.abspath(f) for f in kma_external_files]	
		
		## check if indexed and/or index if necessary
		external_kma_dbs_list = []
		
		## set defaults
		kma_bin = config.get_exe("kma")
		for f in kma_external_files:
			file_name = os.path.basename(f)
			fold_name = os.path.dirname(f)
			print (colored('\t\t+ %s' %file_name, 'green'))
			print ()

			## generate db
			databaseKMA = species_identification_KMA.generate_db([f], file_name, fold_name, 'new', 'single', Debug, kma_bin)
			if not databaseKMA:
				print (colored("***ERROR: Database provided is not indexed.\n" %databaseKMA,'orange'))
			else:
				external_kma_dbs_list.append(databaseKMA)
			
		external_kma_dbs_string = ','.join(external_kma_dbs_list)
		option_db = "kma_external:" + external_kma_dbs_string

	else:
		## rise error & exit
		print (colored("***ERROR: No database provided via --kma_external_file option.\n",'red'))
		exit()

	return(option_db)				

####################################
def get_options_db(options):
	'''
	Selects databases to use and set dataframe with this 
	information among all databases available and according 
	to the input options.
	'''
	
	print ("\n\n+ Select databases to use for identification:")
	
	### database folder to use
	database2use = os.path.abspath(options.database)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: Database to use: " +  database2use + " **", 'yellow'))
	
	## according to user input: select databases to use
	option_db = ""
	
	############################################################
	## Default db KMA
	############################################################
	kma_dbs = []
	if not options.only_kma_db: ## exclusive
		kma_dbs = ["bacteria", "plasmids"]
	
	if (options.kma_dbs):
		options.kma_dbs = options.kma_dbs + kma_dbs
		options.kma_dbs = set(options.kma_dbs)		
	else:
		options.kma_dbs = kma_dbs

	## rise error & exit if no dbs provided
	if not (options.kma_dbs):
		print (colored("***ERROR: No database provided via --kma_db option.\n",'red'))
		exit()

	############################################################
	### Options:
	
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
	## 3) only external kma
	############
	elif (options.only_external_kma):
		option_db = get_external_kma(options.kma_external_files, Debug)
		## rise attention
		if (options.kma_dbs):
			print (colored("***ATTENTION:\nDefatult databases and databases provided via --kma_dbs option would not be used as --only_external_kma option provided.\n",'red'))

	#################
	## all databases 
	#################
	else:		
		####################
		## default KMA dbs
		####################
		option_db = get_kma_db(options.kma_dbs)
		
		#################
		## External file
		#################
		if (options.kma_external_files):
			option_db_tmp = get_external_kma(options.kma_external_files, Debug)
			option_db = option_db + '#' + option_db_tmp
			
		#############################
		## Previously identified data
		#############################
		if any([options.user_data, options.all_data]):
			option_db = option_db + '#kma_user_data:user_data'
	
		#############################
		## Genbank reference data
		#############################
		if any([options.genbank_data, options.all_data]):
			option_db = option_db + '#kma_NCBI:genbank'

	###############
	### PubMLST ###
	###############
	print ("\n\t - Select MLST profiles")	
	option_db_PubMLST = 'MLST:PubMLST'
	print (colored("\t\t + Default MLST profile under database provided: PubMLST", 'green'))

	if options.MLST_profile:
		## user provides a PubMLST profile
		options.MLST_profile = os.path.abspath(options.MLST_profile)
		option_db_PubMLST = option_db_PubMLST + '#MLST:'	+ options.MLST_profile
		print (colored("\t\t + User provided MLST profile: %s" %options.MLST_profile, 'green'))
	
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
		print (colored("**DEBUG: option_db_PubMLST : " +  option_db_PubMLST  + " **", 'yellow'))

	pd_KMA = database_generator.getdbs("KMA", database2use, option_db, Debug)
	pd_PubMLST = database_generator.getdbs("MLST", database2use, option_db_PubMLST, Debug)

	functions.print_sepLine("-",50, False)

	## return both dataframes
	pd_Merge = pd.concat([pd_KMA, pd_PubMLST], sort=True, ignore_index=True)
	return (pd_Merge)
