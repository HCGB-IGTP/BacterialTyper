#!/usr/bin/env python3
'''
This code calls a database generator script for initiating, updating and configure database.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful modules
import time
import os
import sys
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper import database_generator
from BacterialTyper import ariba_caller
from BacterialTyper import multiQC_report
from BacterialTyper import species_identification_KMA
from BacterialTyper import config 
from BacterialTyper import functions
from BacterialTyper import BUSCO_caller

###############################################################
def run(options):

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
				functions.print_sepLine("*",50, False)
				print ("--------- Check NCBI ids provided ---------\n")
				## get file information
				print ("\t+ Obtaining information from file: %s" %abs_path_file)
				strains2get = functions.get_data(abs_path_file, ',', '')
				dataBase_NCBI = database_generator.NCBI_DB(strains2get, NCBI_folder, Debug)

				#########
				if Debug:
					print (colored("DEBUG: NCBI data provided: ", 'yellow'))
					print (dataFile)

				## functions.timestamp
				start_time_partial = functions.timestamp(start_time_partial)
				## strains downloaded would be included to a kma index

		## Get all entries belonging to this taxon provided
		if (options.descendant):
			#########
			if Debug:
				print (colored("DEBUG: NCBI descendant option: ON ", 'yellow'))
			
			functions.print_sepLine("*",70, False)
			print ("--------- Check descendant NCBI taxonomy ids provided ---------\n")
			## [TODO]
			dataBase_NCBI = database_generator.NCBI_descendant(options.descendant, NCBI_folder, Debug)
			
		##############################################################
		## update KMA database with NCBI information retrieved
		##############################################################
		print ('\n\n+ Update database for later identification analysis...')
		list_of_files = dataBase_NCBI['genome'].tolist()
		kma_bin = config.get_exe("kma")
		kma_db = functions.create_subfolder('KMA_db', options.path)	
		genbank_kma_db = functions.create_subfolder('genbank', kma_db)	
		
		print ('+ Database to update: ', user_kma_db)
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

			functions.print_sepLine("*",70, False)
			print ("--------- Check user provided project folder ---------\n")
			dataBase_user = database_generator.update_database_user_data(options.path, abs_project_folder, Debug, options)
		else:
			print (colored("ERROR: Folder provided does not exists: %s" %options.project_folder, 'red'))
			exit()

		##############################################################
		## update KMA database with user_data information retrieved
		##############################################################
		print ('\n\n+ Update database for later identification analysis...')
		list_of_files = dataBase_user['genome'].tolist()
		kma_bin = config.get_exe("kma")
		kma_db = functions.create_subfolder('KMA_db', options.path)	
		user_kma_db = functions.create_subfolder('user_data', kma_db)	
		
		print ('+ Database to update: ', user_kma_db)
		species_identification_KMA.generate_db(list_of_files, 'userData_KMA', user_kma_db, 'new', 'batch', Debug, kma_bin)
		
		## time stamp
		start_time_partial = functions.timestamp(start_time_total)

	##########
	## ARIBA
	##########
	functions.print_sepLine("*",50, False)
	print ("--------- Check ARIBA parameters provided --------")
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
	exit()

##################################################
def getdbs(source, database_folder, option, debug):
	## option = kma:archaea,plasmids,bacteria#kma_external:/path/to/file1,/path/to/file2#user_data#genbank **
	## read folders within database
	files = os.listdir(database_folder) ## ARIBA/KMA_db/genbank/user_data

	## debug message
	if (debug):
		print (colored("Folders: " + str(files),'yellow'))
		print ()
	
	## init dataframe
	colname = ["source", "db", "path"]
	db_Dataframe  = pd.DataFrame(columns = colname)
	
	## user input
	dbs2use = []
	option_list = option.split("#")
	
	for option_item in option_list:
		
		## debug message
		if (debug):
			print (colored("Option item: " + option_item,'yellow'))
		
		###
		dbs2use_tmp = []
		
		## kma
		if (option_item.startswith('kma')):
			if (option_item.startswith('kma:')):
				dbs2use_tmp = option_item.split(":")[1].split(",")

			elif (option_item.startswith('kma_external:')):
				external = option_item.split(":")[1].split(",")

				## add to dataframe			
				for ext in external:
					name_ext = os.path.basename(ext)
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_External', name_ext, ext]

			elif (option_item.startswith('kma_user_data:')):
				dbs2use_tmp = option_item.split(":")[1].split(",")
			
			elif (option_item.startswith('kma_NCBI:')):
				dbs2use_tmp = option_item.split(":")[1].split(",")
			
		### ARIBA
		elif (option_item.startswith('ARIBA:')):
			dbs2use = option_item.split(":")[1].split(",")
		
		### NCBI: genbank
		elif (option_item.startswith('genbank')):
			dbs2use.append('genbank')
		
		### NCBI: taxonomy ID
		elif (option_item.startswith('tax_id')):
			dbs2use.append('taxonomy_id')
		
		### user_data
		elif (option_item.startswith('user_data')):
			dbs2use.append('user_data')
		
		### PubMLST
		elif (option_item.startswith('MLST')):
			dbs2use_tmp = option_item.split(":")[1].split(",")

		### Other?
		else:
			dbs2use.append(option_item) ## add ARIBA, user_data or genbank option if provided

		## get all		
		dbs2use = dbs2use + dbs2use_tmp
	

	## debug message
	if (debug):
		print (colored("\ndbs2use:\n\t" + "\n\t".join(dbs2use), 'yellow'))
	
	###############
	#### ARIBA ####
	###############
	if (source == 'ARIBA'):
		### Check if folder exists
		functions.create_subfolder('ARIBA', database_folder)
		
		### get information
		ARIBA_dbs = ariba_caller.get_ARIBA_dbs(dbs2use) ## get names
		for ariba_db in ARIBA_dbs:
			this_db = database_folder + '/ARIBA/' + ariba_db + '_prepareref/'
			if os.path.exists(this_db):
				code_check_db = ariba_caller.check_db_indexed(this_db, 'NO')
				if (code_check_db == True):
					db_Dataframe.loc[len(db_Dataframe)] = ['ARIBA', ariba_db, this_db]
					print (colored("\t- ARIBA: including information from database: " + ariba_db, 'green'))
			else:
				print ("+ Database: ", ariba_db, " is not downloaded...")
				print ("+ Download now:")
				folder_db = functions.create_subfolder(ariba_db, database_folder + '/ARIBA')
				code_db = ariba_caller.ariba_getref(ariba_db, folder_db, debug, 2) ## get names 
				if (code_db == 'OK'):
					db_Dataframe.loc[len(db_Dataframe)] = ['ARIBA', ariba_db, this_db]
					print (colored("\t- ARIBA: including information from database: " + ariba_db, 'green'))

	#############
	#### KMA ####
	#############
	elif (source == 'KMA'):
		### Check if folder exists
		KMA_db_abs = functions.create_subfolder('KMA_db', database_folder)
		kma_dbs = os.listdir(KMA_db_abs)

		## debug message
		if (debug):
			print (colored("Folders KMA_db:" + str(kma_dbs) , 'yellow'))

		### get information
		for db in dbs2use:
			this_db = KMA_db_abs + '/' + db
			
			#### genbank	
			if (db == "genbank"):
				## KMA databases exists
				this_db_file = this_db + '/genbank_KMA'
				if os.path.isfile(this_db_file + '.comp.b'):
					print (colored("\t- genbank: including information from different reference strains available.", 'green')) ## include data from NCBI
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_genbank', 'genbank', this_db_file]
		
			#### user_data
			elif (db == "user_data"):
				## KMA databases exists
				this_db_file = this_db + '/userData_KMA'
				if os.path.isfile(this_db_file + '.comp.b'):
					print (colored("\t- user_data: including information from user previously generated results", 'green')) ## include user data
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_user_data', 'user_data', this_db_file]
					
			
			## default KMA databases: bacteria & plasmids
			else:
				##
				if (db == 'plasmids'):
					prefix = '.T'
				else:
					prefix = '.ATG'

				this_db_file = this_db + '/' + db + prefix
				if os.path.isfile(this_db_file + '.comp.b'):
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_db', db, this_db_file]
					print (colored("\t- KMA: including information from database " + db, 'green'))
				else:
					print (colored("\t**KMA: Database %s was not available." %db, 'red'))

					## if missing: call download module
					print ("+ Download missing KMA_db (%s) provided" %db)
					species_identification_KMA.download_kma_database(database_folder + '/KMA_db/' + db, db, debug)

					if os.path.isfile(this_db_file + '.comp.b'):
						db_Dataframe.loc[len(db_Dataframe)] = ['KMA_db', db, this_db_file]
						print (colored("\t- KMA: including information from database " + db, 'green'))
					else:
						print (colored("\t**KMA: Database %s was not available." %db, 'red'))

	##############
	#### NCBI ####
	##############
	elif (source == 'NCBI'):

		### Check if folder exists
		db2use_abs = functions.create_subfolder(dbs2use[0], database_folder)
		
		### genbank entries downloaded
		if dbs2use[0] == 'genbank':
			##
			if os.path.exists(db2use_abs + '/bacteria/'):
				genbank_entries = os.listdir(db2use_abs + '/bacteria/')
				for entry in genbank_entries:
					this_db = db2use_abs + '/bacteria/' + entry
					db_Dataframe.loc[len(db_Dataframe)] = ['NCBI:genbank', entry, this_db]

		elif dbs2use[0] == 'tax_id':		
			tax_id_entries = db2use_abs

	###################
	#### user_data ####
	###################
	elif (source == 'user_data'):
		### Check if folder exists
		db2use_abs = functions.create_subfolder(dbs2use[0], database_folder)

		user_entries = os.listdir(db2use_abs)
		for entry in user_entries:
			this_db = db2use_abs + '/' + entry
			db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, this_db]
	
	#################
	#### PubMLST ####
	#################
	elif (source == 'MLST'):
		### get information
		for db in dbs2use:
			if db == 'PubMLST':
				### Check if folder exists
				db2use_abs = functions.create_subfolder('PubMLST', database_folder)
				list_profiles = os.listdir(db2use_abs)

				for entry in list_profiles:
					this_db = db2use_abs + '/' + entry
					db_Dataframe.loc[len(db_Dataframe)] = ['MLST', 'PubMLST', entry + ',' + this_db]
					print (colored("\t- MLST: including information from profile: " + entry, 'green'))
					
			else:
					db_Dataframe.loc[len(db_Dataframe)] = ['MLST', 'user_profile', db]
					print (colored("\t- MLST: including information from profile provided by user: " + db, 'green'))

	return (db_Dataframe)


