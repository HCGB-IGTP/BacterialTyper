#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Initiates, updates and configures database.
"""
## useful modules
import time
import os
import sys
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import database_user
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.scripts import species_identification_KMA
from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.config import set_config 
from BacterialTyper import __version__ as pipeline_version

import HCGB
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.info_functions as HCGB_info

###############################################################
def run_database(options):

	## print further information if requested
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()

	elif (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
		exit()
		
	elif (options.help_KMA):
		species_identification_KMA.help_kma_database()
		exit()
	
	## init time
	start_time_total = time.time()
	start_time_partial = start_time_total

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
		print ("[Debug mode: ON]")
	else:
		Debug = False

	## message header
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Database")
	print ("--------- Starting Process ---------")
	HCGB_time.print_time()
	
	## set some vars
	kma_bin = set_config.get_exe("kma")

	## create folder absolute
	options.path = os.path.abspath(options.path)	
	HCGB_files.create_folder(options.path)

	## databases
	dict_db_info = {}

	#########
	if Debug:
		print (colored("DEBUG: absolute path folder: " + options.path, 'yellow'))

	##########
	## NCBI	##
	##########
	## if any NCBI options provided
	if any ([options.ID_file, options.descendant]):
		## create folders
		NCBI_folder = HCGB_files.create_subfolder('NCBI', options.path)
		if (options.ID_file):
			## get path and check if it is file
			abs_path_file = os.path.abspath(options.ID_file)
			if os.path.isfile(abs_path_file):
				print ()
				HCGB_aes.print_sepLine("*",50, False)
				print ("--------- Check NCBI ids provided ---------\n")
				HCGB_aes.print_sepLine("*",70, False)
				## get file information
				print ("\t+ Obtaining information from file: %s" %abs_path_file)
				strains2get = HCGB_main.get_data(abs_path_file, ',', '')
				dataBase_NCBI = database_generator.NCBI_DB(strains2get, NCBI_folder, Debug)

				#########
				if Debug:
					print (colored("DEBUG: NCBI data provided: ", 'yellow'))
					print (options.ID_file)

				## functions.timestamp
				start_time_partial = HCGB_time.timestamp(start_time_partial)
				## strains downloaded would be included to a kma index

		## Get all entries belonging to this taxon provided
		if (options.descendant):
			#########
			if Debug:
				print (colored("DEBUG: NCBI descendant option: ON ", 'yellow'))
			
			print ()
			HCGB_aes.print_sepLine("*",70, False)
			print ("--------- Check descendant NCBI taxonomy ids provided ---------\n")
			HCGB_aes.print_sepLine("*",70, False)
			## [TODO]
			dataBase_NCBI = database_generator.NCBI_descendant(options.descendant, NCBI_folder, Debug)
			
		##############################################################
		## update KMA database with NCBI information retrieved
		##############################################################
		print ('\n\n+ Update database for later identification analysis...')
		list_of_files = dataBase_NCBI['genome'].tolist()
		kma_db = HCGB_files.create_subfolder('KMA_db', options.path)	
		genbank_kma_db = HCGB_files.create_subfolder('genbank', kma_db)	
		
		print ('+ Database to update: ', genbank_kma_db)
		species_identification_KMA.generate_db(list_of_files, 'genbank_KMA', genbank_kma_db, 'new', 'batch', Debug, kma_bin)

		## time stamp
		start_time_partial = HCGB_time.timestamp(start_time_total)

	###############
	## user_data ##
	###############	
	if options.project_folder:
		
		##
		dataBase_user = pd.DataFrame()	
		## get absolute path
		abs_project_folder = os.path.abspath(options.project_folder)
		if os.path.exists(abs_project_folder):
			#########
			if Debug:
				print (colored("DEBUG: User provides folder containing project", 'yellow'))

			print ()
			HCGB_aes.print_sepLine("*",70, False)
			print ("--------- Check user provided project folder ---------")
			HCGB_aes.print_sepLine("*",70, False)
			dataBase_user = database_user.update_database_user_data(options.path, abs_project_folder, Debug, options)
		else:
			print (colored("ERROR: Folder provided does not exists: %s" %options.project_folder, 'red'))
			exit()

		##############################################################
		## update KMA database with user_data information retrieved
		##############################################################
		print ('\n\n+ Update database for later identification analysis...')
		list_of_files = dataBase_user['genome'].tolist()
		kma_db = HCGB_files.create_subfolder('KMA_db', options.path)	
		user_kma_db = HCGB_files.create_subfolder('user_data', kma_db)	
		
		print ('+ Database to update: ', user_kma_db)
		species_identification_KMA.generate_db(list_of_files, 'userData_KMA', user_kma_db, 'new', 'batch', Debug, kma_bin)
		
		## time stamp
		start_time_partial = HCGB_time.timestamp(start_time_total)

	##########
	## ARIBA
	##########
	print ()
	HCGB_aes.print_sepLine("*",50, False)
	print ("--------- Check ARIBA parameters provided --------")
	HCGB_aes.print_sepLine("*",50, False)
	if (options.no_ARIBA):
		print ("+ No ARIBA databases would be downloaded...")
		
		#########
		if Debug:
			print (colored("DEBUG: No option ARIBA", 'yellow'))
	
		## set ARIBA information
		dict_db_info['ARIBA'] = {}
	
	else:
		#functions.print_sepLine("*",50, False)
		
		### ariba list databases
		ariba_dbs_list = ['CARD', 'VFDB']
		
		if (options.no_def_ARIBA):
			ariba_dbs_list = options.ariba_dbs
		else:
			if (options.ariba_dbs):
				ariba_dbs_list = ariba_dbs_list + options.ariba_dbs
				ariba_dbs_list = set(ariba_dbs_list)
		
		## set final
		options.ariba_dbs = ariba_dbs_list 

		#########
		if Debug:
			print (colored("DEBUG: Option ARIBA", 'yellow'))
			print (ariba_dbs_list)
		
		dict_ARIBA_info = ariba_caller.download_ariba_databases(ariba_dbs_list, options.path, Debug, options.threads)
	
		#########
		if Debug:
			print (colored("DEBUG: dict_ARIBA_info", 'yellow'))
			print (dict_ARIBA_info)
	
		### ariba list databases		
		if (options.ariba_users_fasta):
			print ("+ Generate ARIBA database for databases provided: prepare fasta and metadata information")

			#########
			if Debug:
				print (colored("DEBUG: Option user ARIBA db", 'yellow'))
				print (ariba_users_fasta)
				print (ariba_users_meta)

			## [TODO]:	
			## ariba prepareref fasta and metadata

		### timestamp
		start_time_partial = HCGB_time.timestamp(start_time_partial)					

		## set ARIBA information
		dict_db_info['ARIBA'] = dict_ARIBA_info

	#########
	## kma ##
	#########
	print ()
	HCGB_aes.print_sepLine("*",50, False)
	print ("--------- Check KMA parameters provided ----------")
	kma_database = options.path + '/KMA_db'	
	HCGB_files.create_folder(kma_database)
	
	## types: bacteria, archaea, protozoa, fungi, plasmids, typestrains
	## downloads all "bacterial" genomes from KMA website

	print ("+ Retrieving information from: KmerFinder CGE databases website")		

	## KMA databases to use	
	## only user dbs	
	if (options.no_def_kma):
		if (options.kma_dbs):
			print("+ Only user databases selected will be indexed...")
		else:
			print ("+ No databases selected.")
			print (colored("ERROR: Please select a kma database.", 'red'))
			exit()
			
	## default dbs + user
	else:
		kma_dbs = ["bacteria"] ## no plasmids available anymore

		## default dbs + user
		if (options.kma_dbs):
			options.kma_dbs = options.kma_dbs + kma_dbs
			options.kma_dbs = set(options.kma_dbs)		
		else:
			options.kma_dbs = kma_dbs

	#########
	if Debug:
		print (colored("DEBUG: options.kma_dbs", 'yellow'))
		print (options.kma_dbs)

	## Get databases
	dict_KMA_db_info = {}
	for db in options.kma_dbs:
		print (colored("\n+ " + db, 'yellow'))
		db_folder = HCGB_files.create_subfolder(db, kma_database)		
		db_dict = species_identification_KMA.download_kma_database(db_folder, db, Debug)
		dict_KMA_db_info[db] = db_dict

	#########
	if Debug:
		print (colored("DEBUG: dict_KMA_db_info", 'yellow'))
		print (dict_KMA_db_info)

	## set dictionary
	dict_db_info['kma'] = dict_KMA_db_info
	
	### timestamp
	start_time_partial = HCGB_time.timestamp(start_time_partial)					

	print ("\n*************** Finish *******************\n")
	start_time_partial = HCGB_time.timestamp(start_time_total)

	## dump information and parameters
	info_dir = HCGB_files.create_subfolder("info", options.path)
	print("+ Dumping information and parameters")
	runInfo = { "module":"database", "time":time.time(),
				"BacterialTyper version":pipeline_version,
				'database_info': dict_db_info}
	
	## dump database information
	HCGB_info.dump_info_run(info_dir, 'database', options, runInfo, options.debug)

	print ("+ Exiting Database module.\n")
	return()
