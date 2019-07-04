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

## import my modules
from BacterialTyper import database_generator
from BacterialTyper import ariba_caller
from BacterialTyper import multiQC_report
from BacterialTyper import species_identification_KMA
from BacterialTyper import config 
from BacterialTyper import functions
from BacterialTyper import BUSCO_caller

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

	## print further information for ARIBA databases or BUSCO and exit
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()
	elif (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
		exit()

	## create folder and call modules:	
	functions.create_folder(os.path.abspath(options.path))

	## NCBI
	functions.print_sepLine("*",50, False)
	dataFile = database_generator.init_DB(options.ID_file, options.path)

	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_total)

	## ARIBA
	if (options.no_ARIBA):
		print ()
	else:
		functions.print_sepLine("*",50, False)
		ariba_caller.download_ariba_databases(options.ariba_dbs, options.path, Debug)
		### timestamp
		start_time_partial = functions.timestamp(start_time_partial)					

	### kma databases
	functions.print_sepLine("*",50, False)
	kma_database = options.path + '/KMA_db'	
	if (options.index_KMA):
		index_db_kma(kma_database, '/KMA_user')		
	else:
		kma_download(options, kma_database)		

	### timestamp
	start_time_partial = functions.timestamp(start_time_partial)					
	
	#### plasmid_data
	##

	## BUSCO datasets
	if (options.BUSCO_dbs):
		BUSCO_folder = functions.create_subfolder("BUSCO", options.path)
		BUSCO_caller.BUSCO_retrieve_sets(options.BUSCO_dbs, BUSCO_folder)
	
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

	## print further information for ARIBA databases or BUSCO and exit
	if (options.help_ARIBA):
		print ("ARIBA databases information:")	
		ariba_caller.help_ARIBA()
		exit()
	elif (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
		exit()

	## absolute
	options.path = os.path.abspath(options.path)
	
	## time
	start_time_partial = start_time_total
	
	## genbank
	if (options.ID_file):
		print ("--------- Check ID file provided ---------")
		dataFile = database_generator.update_database(options.ID_file, options.path)
		### timestamp
		start_time_partial = functions.timestamp(start_time_partial)					
	
	## ARIBA
	if (options.ariba_dbs):
		print ("--------- Check ARIBA databases provided ---------")
		ariba_caller.download_ariba_databases(options.ariba_dbs, options.path, Debug)
		### timestamp
		start_time_partial = functions.timestamp(start_time_partial)					
		functions.print_sepLine("*",50, False)
	
	### kma databases
	print ("--------- Check KMA databases provided ---------")
	kma_database = options.path + '/KMA_db'	
	if (options.index_kma):
		index_db_kma(kma_database, '/KMA_user')		
	else:
		kma_download(options, kma_database)		
	
	### timestamp
	print ("")
	start_time_partial = functions.timestamp(start_time_partial)					
	print ("")

	## BUSCO datasets
	if (options.BUSCO_dbs):
		print ("--------- Check BUSCO datasets provided ---------")
		options.BUSCO_dbs = set(options.BUSCO_dbs)
		BUSCO_folder = functions.create_subfolder("BUSCO", options.path)
		BUSCO_caller.BUSCO_retrieve_sets(options.BUSCO_dbs, BUSCO_folder)

	print ("\n\n*************** Finish *******************")
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
	## [TODO]
	## ToDo implement self index of database
	print()


