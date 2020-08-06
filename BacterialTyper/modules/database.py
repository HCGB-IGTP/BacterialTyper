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
from BacterialTyper.config import set_config 
from BacterialTyper.scripts import functions
from BacterialTyper.scripts import BUSCO_caller

###############################################################
def run_database(options):

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
	functions.pipeline_header()
	functions.boxymcboxface("Database")
	print ("--------- Starting Process ---------")
	functions.print_time()

	kma_bin = set_config.get_exe("kma")

	######################################################
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
	######################################################

	## create folder
	## absolute
	options.path = os.path.abspath(options.path)	
	functions.create_folder(options.path)

	#########
	if Debug:
		print (colored("DEBUG: absolute path folder: " + options.path, 'yellow'))

	##########
	## NCBI	##
	##########
	## if any NCBI options provided
	if any ([options.ID_file, options.descendant]):
		## create folders
		NCBI_folder = functions.create_subfolder('NCBI', options.path)
		if (options.ID_file):
			## get path and check if it is file
			abs_path_file = os.path.abspath(options.ID_file)
			if os.path.isfile(abs_path_file):
				print ()
				functions.print_sepLine("*",50, False)
				print ("--------- Check NCBI ids provided ---------\n")
				functions.print_sepLine("*",70, False)
				## get file information
				print ("\t+ Obtaining information from file: %s" %abs_path_file)
				strains2get = functions.get_data(abs_path_file, ',', '')
				dataBase_NCBI = database_generator.NCBI_DB(strains2get, NCBI_folder, Debug)

				#########
				if Debug:
					print (colored("DEBUG: NCBI data provided: ", 'yellow'))
					print (options.ID_file)

				## functions.timestamp
				start_time_partial = functions.timestamp(start_time_partial)
				## strains downloaded would be included to a kma index

		## Get all entries belonging to this taxon provided
		if (options.descendant):
			#########
			if Debug:
				print (colored("DEBUG: NCBI descendant option: ON ", 'yellow'))
			
			print ()
			functions.print_sepLine("*",70, False)
			print ("--------- Check descendant NCBI taxonomy ids provided ---------\n")
			functions.print_sepLine("*",70, False)
			## [TODO]
			dataBase_NCBI = database_generator.NCBI_descendant(options.descendant, NCBI_folder, Debug)
			
		##############################################################
		## update KMA database with NCBI information retrieved
		##############################################################
		print ('\n\n+ Update database for later identification analysis...')
		list_of_files = dataBase_NCBI['genome'].tolist()
		kma_db = functions.create_subfolder('KMA_db', options.path)	
		genbank_kma_db = functions.create_subfolder('genbank', kma_db)	
		
		print ('+ Database to update: ', genbank_kma_db)
		species_identification_KMA.generate_db(list_of_files, 'genbank_KMA', genbank_kma_db, 'new', 'batch', Debug, kma_bin)

		## time stamp
		start_time_partial = functions.timestamp(start_time_total)

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
			functions.print_sepLine("*",70, False)
			print ("--------- Check user provided project folder ---------")
			functions.print_sepLine("*",70, False)
			dataBase_user = database_user.update_database_user_data(options.path, abs_project_folder, Debug, options)
		else:
			print (colored("ERROR: Folder provided does not exists: %s" %options.project_folder, 'red'))
			exit()

		##############################################################
		## update KMA database with user_data information retrieved
		##############################################################
		print ('\n\n+ Update database for later identification analysis...')
		list_of_files = dataBase_user['genome'].tolist()
		kma_db = functions.create_subfolder('KMA_db', options.path)	
		user_kma_db = functions.create_subfolder('user_data', kma_db)	
		
		print ('+ Database to update: ', user_kma_db)
		species_identification_KMA.generate_db(list_of_files, 'userData_KMA', user_kma_db, 'new', 'batch', Debug, kma_bin)
		
		## time stamp
		start_time_partial = functions.timestamp(start_time_total)

	##########
	## ARIBA
	##########
	print ()
	functions.print_sepLine("*",50, False)
	print ("--------- Check ARIBA parameters provided --------")
	functions.print_sepLine("*",50, False)
	if (options.no_ARIBA):
		print ("+ No ARIBA databases would be downloaded...")
		
		#########
		if Debug:
			print (colored("DEBUG: No option ARIBA", 'yellow'))
	
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
			
		#########
		if Debug:
			print (colored("DEBUG: Option ARIBA", 'yellow'))
			print (options.ariba_dbs)
		
		ariba_caller.download_ariba_databases(ariba_dbs_list, options.path, Debug, options.threads)
	
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
		start_time_partial = functions.timestamp(start_time_partial)					

	#########
	## kma ##
	#########
	print ()
	functions.print_sepLine("*",50, False)
	print ("--------- Check KMA parameters provided ----------")
	kma_database = options.path + '/KMA_db'	
	functions.create_folder(kma_database)
	
	## types: bacteria, archaea, protozoa, fungi, plasmids, typestrains
	## downloads all "bacterial" genomes from KMA website
	## kma: ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/

	print ("+ Retrieving information from: ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder website")		

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
		kma_dbs = ["bacteria", "plasmids"]

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
	for db in options.kma_dbs:
		print (colored("\n+ " + db, 'yellow'))
		db_folder = functions.create_subfolder(db, kma_database)		
		species_identification_KMA.download_kma_database(db_folder, db, Debug)

	### timestamp
	start_time_partial = functions.timestamp(start_time_partial)					

	###########	
	## BUSCO ##
	###########
	if (options.BUSCO_dbs):
		print ()
		functions.print_sepLine("*",50, False)
		print ("--------- Check BUSCO datasets provided ---------")
		BUSCO_folder = functions.create_subfolder("BUSCO", options.path)

		#########
		if Debug:
			print (colored("DEBUG: options.BUSCO_dbs", 'yellow'))
			print (options.BUSCO_dbs)

		BUSCO_caller.BUSCO_retrieve_sets(options.BUSCO_dbs, BUSCO_folder)

		### timestamp
		start_time_partial = functions.timestamp(start_time_partial)					

	print ("\n*************** Finish *******************\n")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Database module.\n")
	return()
