#!/usr/bin/env python3
'''
This code calls a database generator script for initiating, updating and configure database.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import database_generator, ariba_caller, multiQC_report, species_identification_KMA
from BacterialTyper import config, functions
from termcolor import colored
import os, sys

###############################################################
def initial_run(options):

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
	functions.boxymcboxface("Initiate Database")
	print ("--------- Starting Process ---------")
	functions.print_time()

	## print further information for ARIBA databases	
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()


	## create folder and call modules:	
	functions.create_folder(os.path.abspath(options.path))

	## NCBI
	functions.print_sepLine("*",50)
	dataFile = database_generator.init_DB(options.ID_file, options.path)

	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_total)

	## ARIBA
	if (options.no_ARIBA):
		print ()
	else:
		functions.print_sepLine("*",50)
		ariba_caller.download_ariba_databases(list_dbs, options.path, Debug)
		### timestamp
		start_time_partial = functions.timestamp(start_time_partial)					

	### kma databases
	functions.print_sepLine("*",50)
	kma_database = options.path + '/KMA_db'	
	if (options.index_KMA):
		index_db_kma(kma_database, '/KMA_user')		
	else:
		kma_download(options, kma_database)		

	### timestamp
	start_time_partial = functions.timestamp(start_time_partial)					
	
	#### plasmid_data
	##
	
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Initiate Database module.")
	exit()

###############################################################
def updateDB_NCBI(options):
	
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
	functions.boxymcboxface("Update Database")
	print ("--------- Starting Process ---------")
	functions.print_time()

	options.path = os.path.abspath(options.path)
	
	## genbank
	if (options.ID_file):
		dataFile = database_generator.update_database(options.ID_file, options.path)

	## print further information for ARIBA databases	
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()

	
	## ARIBA
	if (options.ARIBA_db):
		functions.print_sepLine("*",50)
		ariba_caller.download_ariba_databases(options.path, Debug)
	
	### kma databases
	functions.print_sepLine("*",50)
	kma_database = options.path + '/KMA_db'	
	if (options.index_kma):
		index_db_kma(kma_database, '/KMA_user')		
	else:
		kma_download(options, kma_database)		

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Update Database module.")
	exit()

###############################################################
def kma_download(options, database_folder):

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
	
	for db in options.kma_dbs:
		print (colored("\n+ " + db, 'yellow'))
		db_folder = functions.create_subfolder(db, database_folder)		
		species_identification_KMA.download_kma_database(db_folder, db, Debug)


###############################################################
def index_db_kma():
	## ToDo implement self index of database
	print()


