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
def index_db_kma():
	print()

############################################################### 
def download_ariba_databases(main_folder):
	
	print("\n\n+ Download databases for Antimicrobial Resistance Identification By Assembly (ARIBA).")
	ariba_folder = functions.create_subfolder("ARIBA", main_folder)
	
	## where database is one of: 
	print ("+ Available databases:")

	out_info = main_folder + '/ARIBA_information.txt'
	hd = open(out_info, 'w')

	dbs = ["argannot", "card", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"]
	for db_set in dbs:
		print (colored("+ " + db_set,'yellow'))
		folder_set = functions.create_subfolder(db_set, ariba_folder)
		ariba_caller.ariba_getref(db_set, folder_set + "/" + db_set)

		## print citation for each database
		functions.print_sepLine("'", 75)

		hd.write(db_set)
		hd.write('\n')
	
	hd.close()
	

###############################################################
def updateDB_NCBI(options):
	## genbank
	dataFile = database_generator.update_database(options.ID_file, options.path)
	
	## ARIBA
	if (options.ARIBA_db):
		download_ariba_databases(options.path)
	
	## KMA_index
	if (options.index_KMA):
		index_db_kma(kma_database, 'user')		


###############################################################
def initial_run(options):

	## create folder and call modules:	
	functions.pipeline_header()
	functions.create_folder(os.path.abspath(options.path))

	## NCBI
	dataFile = database_generator.init_DB(options.ID_file, options.path)

	## ARIBA
	if (options.no_ARIBA):
		print ()
	else:
		download_ariba_databases(options.path)

	### kma databases
	kma_database = options.path + '/kma_db'	
	if (options.index_KMA):
		index_db_kma(kma_database, 'user')		
	else:
	
		## types: bacteria, archaea, protozoa, fungi, plasmids, typestrains
	
		## downloads all "bacterial" genomes from KMA website
		## kma: ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/
		print ("+ Retrieving information from: ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/...")		
		species_identification_KMA.download_kma_database(kma_database, options.kma_db)	

	#### plasmid_data


