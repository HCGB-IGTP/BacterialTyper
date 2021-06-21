#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
This code prepares the database information for further analysis.
Several functions are implemented for:
	
	- Manipulate data entries, update information
	
	- Download NCBI assembly IDs provided
"""
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
import ncbi_genome_download as ngd
import shutil
from sys import argv
from io import open
from termcolor import colored
from Bio import SeqIO
import concurrent.futures

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import species_identification_KMA
from BacterialTyper.scripts import min_hash_caller

import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

##########################################################################################
def NCBI_DB(strains2get, data_folder, Debug):
	"""Donwloads given taxa from NCBI if not available and updates database information.
	
	This function checks in the given folder if strain of interest is available. If not it would connect to NCBI using python module ncbi_genome_download and downloads some information.
	
	:param strains2get: dataframe containing genus, species and NCBI assembly columns among others. See example below.
	:param data_folder: Absolute path to database NCBI folder.
	:param Debug: Print messages for debugging purposes if desired. 
	:type strains2get: dataframe
	:type data_folder: string
	:type Debug: bool
	:return: Dataframe of genbank database updated for all available entries.

	Columns for the dataframe :file:`strains2get` consist of:
	
	sample,genus,species,strain,BioSample,genome,Plasmids
 
	See and example in file: :file:`/devel/results/strains2get_NCBI_DB.csv` and shown here:
	
	.. include:: ../../devel/results/strains2get_NCBI_DB.csv
		:literal:
		
	See example of the return dataframe, containing database information updated in file: :file:`/devel/results/genbank_database.csv` here:
	
	.. include:: ../../devel/results/genbank_database.csv
		:literal:
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.file_funtcions.create_folder`
	
		- :func:`HCGB.functions.main_functions.get_data`
	
		- :func:`BacterialTyper.scripts.database_generator.get_dbs`
	
		- :func:`BacterialTyper.scripts.database_generator.get_database`
		
		- :func:`BacterialTyper.scripts.database_generator.NCBIdownload`
		
		- :func:`BacterialTyper.scripts.database_generator.update_db_data_file`
		
	.. include:: ../../links.inc	 	
	
	"""
	
	## set index
	strains2get = strains2get.set_index('NCBI_assembly_ID', drop=False) ## set new index but keep column
	strains2get.index.names = ['ID'] ## rename index
	strains2get = strains2get.drop_duplicates()
	
	#########
	if Debug:
		print (colored("DEBUG: NCBI data provided: ", 'yellow'))
		print (strains2get)

	## get data existing database
	print ("+ Create the database in folder: \n", data_folder)
	HCGB_files.create_folder(data_folder)
	
	## read database 
	db_frame = getdbs('NCBI', data_folder, 'genbank', Debug)
	database_df = get_database(db_frame, Debug)
	
	#########
	if Debug:
		print (colored("DEBUG: NCBI genbank database retrieved: ", 'yellow'))
		print("db_frame")
		print(db_frame)
		print()
		
		print ("database_df")
		print (database_df)
	
	## loop and download
	for index, row in strains2get.iterrows():
		HCGB_aes.print_sepLine("+", 75, False)
		acc_ID = index #strains2get.loc[index]['NCBI_assembly_ID']
		info = "Genus: " + strains2get.loc[index]['genus'] + '\n' + "Species: " +  strains2get.loc[index]['species'] + '\n' + "Strain: " +  strains2get.loc[index]['name'] + '\n' + "ID accession: " +  acc_ID + '\n'
		dir_path = data_folder + '/genbank/bacteria/' + acc_ID ## module ngd requires to download data in bacteria subfolder under genbank folder

		## check if already exists
		if acc_ID in database_df.index:
			print ("\n+ Data is already available in database for: ")
			print (colored(info, 'green'))

		else:
			## download
			print ("\n+ Downloading data for:")	
			print (colored(info, 'green'))
			data_accID = NCBIdownload(acc_ID, strains2get, data_folder)
			this_db = HCGB_main.get_data(data_accID, ',', 'index_col=0')
			this_db = this_db.set_index('ID')
			database_df = database_df.append(this_db)

	## Generate/Update database
	database_csv = data_folder + '/genbank_database.csv'
	db_updated = update_db_data_file(database_df, database_csv)
	print ("+ Database has been generated in file: ", database_csv)
	return (db_updated)

##########################################################################################
def NCBI_descendant(tax_ID, NCBI_folder, Debug):
	print ()
	
	## [TODO]	
	## get information from NCBItaxonomy from module ete3
	## and generate a dataframe containing as header: genus,species,name,NCBI_assembly_ID
	## and call NCBI_DB to download and update database, then returned dataframe.
	
	# Use BacDup implementation for this purpose
	


def ngd_download(dir_path, acc_ID, data_folder):
	download = False
	print ('+ Check data for ID: ', acc_ID)
	if os.path.exists(dir_path):
		print ('+ Folder already exists: ', dir_path)
		## get files download
		(genome, prot, gff, gbk) = get_files_download(dir_path)
		if all([genome, prot, gff, gbk]):
			download = False
		else:
			print ('+ Not all necessary data is available. Download it again.')
			download = True
	else:
		download = True
	
	if download:
		print ('+ Downloading:')
		## download in data folder provided
		ngd.download(section='genbank', file_formats='fasta,gff,protein-fasta,genbank', assembly_accessions=acc_ID, output=data_folder, groups='bacteria')

		## check if files are gunzip
		files = os.listdir(dir_path)
		files_list = []		
		for f in files:
			if f.endswith('gz'):
				files_list.append(f)
				print ("\t- Extracting files: ", f)
				HCGB_files.extract(dir_path + '/' + f, dir_path)
				#os.remove(dir_path + '/' + f)
	else:
		print ('+ Data is already available, no need to download it again')


##########################################################################################
def NCBIdownload(acc_ID, data, data_folder):	
	
	## module ngd requires to download data in bacteria subfolder under genbank folder
	dir_path = os.path.join(data_folder, 'genbank', 'bacteria', acc_ID) 
	ngd_download(dir_path, acc_ID, data_folder)
	
	## get files download
	(genome, prot, gff, gbk) = get_files_download(dir_path)

	## check if any plasmids downloaded
	plasmid_count = 0
	plasmid_id = []
	contig_out_file = dir_path + '/' + acc_ID + '_chromosome.fna'
	plasmid_out_file = dir_path + '/' + acc_ID + '_plasmid.fna' 
	
	## open
	contig_out_file_handle = open(contig_out_file, 'w')
	for seq_record in SeqIO.parse(genome, "fasta"):
		plasmid_search = re.search(r".*plasmid.*", seq_record.description)
		if plasmid_search:
			## count and get names for plasmids
			plasmid_count += 1
			name = str( seq_record.id )
			plasmid_id.append(name)
		
			### Separate plasmids from main sequence
			plasmid_out_file_handle = open(plasmid_out_file, 'a')
			plasmid_out_file_handle.write(seq_record.format("fasta"))
			plasmid_out_file_handle.write('\n')
			plasmid_out_file_handle.close()
		else:
			contig_out_file_handle.write(seq_record.format("fasta"))
			contig_out_file_handle.write('\n')
			contig_out_file_handle.close()

	## no plasmids found
	if plasmid_count == 0:
		plasmid_out_file = ""
		
	data2download=pd.DataFrame(columns=('ID','folder','genus','species','name','genome', 'chr', 'GFF','GBK', 'proteins','plasmids_number','plasmids_ID','plasmids'))
	data2download.loc[len(data2download)] = (acc_ID, dir_path, data.loc[acc_ID]['genus'], 
											data.loc[acc_ID]['species'], data.loc[acc_ID]['name'], genome, 
											contig_out_file, gff, prot, gbk,
											plasmid_count, "::".join(plasmid_id), plasmid_out_file)

	## dump to file
	info_file = dir_path + '/info.txt'
	data2download.to_csv(info_file)
	
	## timestamp
	filename_stamp = dir_path + '/.success'
	stamp =	HCGB_time.print_time_stamp(filename_stamp)

	## return data
	return(info_file)		
	
##########################################################################################
def get_files_download(folder):
	## check if files are gunzip
	files = os.listdir(folder)
	genome=""
	prot=""
	gff=""
	gbk=""
	for f in files:
		if f.endswith('genomic.fna'):
			genome = os.path.join(folder, f)
		elif f.endswith('genomic.gff'):
			gff = os.path.join(folder, f)
		elif f.endswith('genomic.gbk'):
			gbk = os.path.join(folder, f)
		elif f.endswith('genomic.gbff'):
			gbk = os.path.join(folder, f)
		elif f.endswith('protein.faa'):
			prot = os.path.join(folder, f)

	return(genome, prot, gff, gbk)			

##########################################################################################
def get_database(db_frame, Debug):
	data4db = pd.DataFrame()
	for index, row in db_frame.iterrows():
		## information
		this_file = db_frame.loc[index]['path'] + '/info.txt'
		if os.path.isfile(this_file):
			print ('+ Reading information for sample: ', db_frame.loc[index]['db'])
			print (colored("\t+ Obtaining information from file: %s" %this_file, 'yellow'))
			this_db = HCGB_main.get_data(this_file, ',', 'index_col=0')
			data4db = data4db.append(this_db)
			timestamp = db_frame.loc[index]['path'] + '/.success'
			if os.path.isfile(timestamp):
				stamp =	HCGB_time.read_time_stamp(timestamp)
				print (colored("\t+ Data generated on: %s" %stamp, 'yellow'))

			HCGB_aes.print_sepLine("*",25, False)

	## index by ID
	if not data4db.empty:
		data4db = data4db.set_index('ID')

	return(data4db)

##########################################################################################
def update_db_data_file(data, csv):
	if os.path.isfile(csv):
		print ("\n+ Updating database")
		print ("+ Obtaining information from database file: %s" %csv)
		db2update = HCGB_main.get_data(csv, ',', 'index_col=0')
		
		## TODO: provide preference to db2update
		df = pd.concat([db2update, data], join='inner', sort=True).drop_duplicates()
		df.to_csv(csv)
		return (df)
	else:
		data.to_csv(csv)
		return (data)
		
##################################################
def getdbs(source, database_folder, option, debug):
	"""Get databases available within the folder provided.
	
	:param source: Type of database to search: ARIBA, KMA, NCBI, MLST, user_data
	:param database_folder: Absolute path to database folder.
	:param option: String containing multiple entries separated by '#' that indicate the type of database entries to search within each source type.
	:param debug: True/False for debugging messages.
	
	:type source: string
	:type database_folder: string
	:type option: string
	:type debug: bool
	
	:returns: Dataframe containing absolute paths to the available databases for each type requested. It contains columns for: "source", "db", "path"
		
	e.g.: 	source = KMA
			option = kma:archaea,plasmids,bacteria#kma_external:/path/to/file1,/path/to/file2#user_data#genbank **
			
	e.g.: 	source = NCBI
			option = genbank
	
	"""
	
	## init dataframe
	colname = ["source", "db", "path", "timestamp"]
	db_Dataframe  = pd.DataFrame(columns = colname)

	## read folders within database
	if os.path.isdir(database_folder):
		files = os.listdir(database_folder) ## ARIBA/KMA_db/genbank/user_data
	else:
		return db_Dataframe

	## debug message
	if (debug):
		print (colored("Folders: " + str(files),'yellow'))
		print ()
	
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
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_External', name_ext, ext, '']

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
		
		### MLST
		elif (option_item.startswith('MLST')):
			dbs2use_tmp = option_item.split(":")[1].split(",")

		### Mash
		elif (option_item.startswith('Mash')):
			if (option_item.startswith('Mash_external_data:')):
				external = option_item.split(":")[1].split(",")
				## add to dataframe			
				for ext in external:
					name_ext = os.path.basename(ext)
					name_ext_ = name_ext.split('.fna')[0]
					db_Dataframe.loc[len(db_Dataframe)] = ['Mash_external', name_ext_, ext, '']
			else:
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
		ARIBA_folder = HCGB_files.create_subfolder('ARIBA', database_folder)
		
		### get information
		ARIBA_dbs = ariba_caller.get_ARIBA_dbs(dbs2use) ## get names
		for ariba_db in ARIBA_dbs:
			this_db = os.path.join(ARIBA_folder, ariba_db + '_prepareref')
			## debug message
			if (debug):
				print (colored("Checking: " + this_db, 'yellow'))

			if os.path.exists(this_db):
				(code_check_db, stamp) = ariba_caller.check_db_indexed(this_db, 'NO')
				if (code_check_db == True):
					db_Dataframe.loc[len(db_Dataframe)] = ['ARIBA', ariba_db, this_db, stamp]
					print (colored("\t- ARIBA: including information from database: " + ariba_db, 'green'))
				else:
					print("Fail...")
			else:
				print ("+ Database: ", ariba_db, " is not downloaded...")
				print ("+ Download now:")
				folder_db = HCGB_files.create_subfolder(ariba_db, ARIBA_folder)
				code_db = ariba_caller.ariba_getref(ariba_db, folder_db, debug, 2) ## get names 
				if (code_db == 'OK'):
					db_Dataframe.loc[len(db_Dataframe)] = ['ARIBA', ariba_db, this_db, time.time()]
					print (colored("\t- ARIBA: including information from database: " + ariba_db, 'green'))

	#############
	#### KMA ####
	#############
	elif (source == 'KMA'):
		### Check if folder exists
		KMA_db_abs = HCGB_files.create_subfolder('KMA_db', database_folder)
		kma_dbs = os.listdir(KMA_db_abs)

		## debug message
		if (debug):
			print (colored("Folders KMA_db:" + str(kma_dbs) , 'yellow'))

		### get information
		for db in dbs2use:
			this_db = KMA_db_abs + '/' + db

			## debug message
			if (debug):
				print (colored("this_db:" + this_db , 'yellow'))
			
			#### genbank	
			if (db == "genbank"):
				## KMA databases exists
				this_db_file = this_db + '/genbank_KMA'
				if os.path.isfile(this_db_file + '.comp.b'):
					print (colored("\t- genbank: including information from different reference strains available.", 'green')) ## include data from NCBI
					stamp = time.time() ## to do add time
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_genbank', 'genbank', this_db_file, stamp]
		
			#### user_data
			elif (db == "user_data"):
				## KMA databases exists
				this_db_file = this_db + '/userData_KMA'
				if os.path.isfile(this_db_file + '.comp.b'):
					print (colored("\t- user_data: including information from user previously generated results", 'green')) ## include user data
					stamp = time.time() ## to do add time
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_user_data', 'user_data', this_db_file, stamp]
					
			
			## default KMA databases: bacteria & plasmids
			else:
				##
				if (db == 'plasmids'):
					prefix = '.T'
				elif (db == 'viral'):
					prefix = '.TG'
				else:
					prefix = '.ATG'

				this_db_file =os.path.join(this_db, db, db + prefix)
				## debug message
				if (debug):
					print (colored("this_db_file:" + this_db_file , 'yellow'))

				if os.path.isfile(this_db_file + '.comp.b'):
					stamp = os.path.join(this_db, '.success')
					db_Dataframe.loc[len(db_Dataframe)] = ['KMA_db', db, this_db_file, HCGB_main.get_info_file(stamp)[0]]
					print (colored("\t- KMA: including information from database " + db, 'green'))
				else:
					print (colored("\t**KMA: Database %s was not available." %db, 'red'))

					## if missing: call download module
					print ("+ Download missing KMA_db (%s) provided" %db)
					species_identification_KMA.download_kma_database(
						os.path.join(database_folder, 'KMA_db', db), db, debug)

					if os.path.isfile(this_db_file + '.comp.b'):
						stamp = time.time() ## to do add time
						db_Dataframe.loc[len(db_Dataframe)] = ['KMA_db', db, this_db_file, stamp]
						print (colored("\t- KMA: including information from database " + db, 'green'))
					else:
						print (colored("\t**KMA: Database %s was not available." %db, 'red'))

	##############
	#### NCBI ####
	##############
	elif (source == 'NCBI'):
		
		## TODO: get additional information from 
		## info_file = dir_path + '/info.txt'

		### Check if folder exists
		path_genbank = os.path.join(database_folder, source, 'genbank')
		db2use_abs = HCGB_files.create_subfolder(dbs2use[0], database_folder)
		
		### genbank entries downloaded
		if dbs2use[0] == 'genbank':
			##
			if os.path.exists(path_genbank + '/bacteria'):
				genbank_entries = os.listdir(os.path.join(path_genbank, 'bacteria'))
				for entry in genbank_entries:
					this_db = os.path.join(path_genbank,'bacteria', entry)
					stamp = time.time() ## to do add time
					db_Dataframe.loc[len(db_Dataframe)] = ['NCBI:genbank', entry, this_db, stamp]

		elif dbs2use[0] == 'tax_id':		
			tax_id_entries = db2use_abs

	###################
	#### user_data ####
	###################
	elif (source == 'user_data'):
		### Check if folder exists
		db2use_abs = HCGB_files.create_subfolder(dbs2use[0], database_folder)

		user_entries = os.listdir(db2use_abs)
		for entry in user_entries:
			this_db = db2use_abs + '/' + entry
			stamp = time.time() ## to do add time
			db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, this_db, stamp]
	
	#################
	#### PubMLST ####
	#################
	elif (source == 'MLST'):
		### get information
		for db in dbs2use:
			if db == 'PubMLST':
				### Check if folder exists
				db2use_abs = HCGB_files.create_subfolder('PubMLST', database_folder)
				list_profiles = os.listdir(db2use_abs)

				for entry in list_profiles:
					this_db = db2use_abs + '/' + entry
					stamp = time.time() ## to do add time
					db_Dataframe.loc[len(db_Dataframe)] = ['MLST', 'PubMLST', entry + ',' + this_db, stamp]
					print (colored("\t- MLST: including information from profile: " + entry, 'green'))
					
			else:
				stamp = time.time() ## to do add time
				db_Dataframe.loc[len(db_Dataframe)] = ['MLST', 'user_profile', db, stamp]
				print (colored("\t- MLST: including information from profile provided by user: " + db, 'green'))

	##################
	#### Min Hash ####
	##################
	elif (source == 'MASH'):

		db_Dataframe['original'] = ''
		db_Dataframe['ksize'] = ''
		db_Dataframe['num_sketch'] = ''
		db_Dataframe['folder'] = ''
	
		### get information
		for db in dbs2use:
			#### genbank	
			if (db == "genbank"):
				
				### Check if folder exists
				db2use_abs = database_folder + '/NCBI/genbank/bacteria'
				if os.path.exists(db2use_abs):
					print (colored("\n\t- genbank: including information from different reference strains available.", 'green')) ## include data from NCBI
					genbank_entries = os.listdir(db2use_abs)
					for entry in genbank_entries:
						print ('\t+ Reading information from sample: ', entry)
						this_db = db2use_abs + '/' + entry
						
						## get additional information from 
						info_file = this_db + '/info.txt'
						info_data = pd.read_csv(info_file).set_index('ID')
						
						info_data.fillna("NaN", inplace=True)
						
						## get readable name for each strain
						entry_strain = str(info_data.loc[entry]['name'])
						
						if entry_strain == 'NaN': ## TODO: debug if it works
							entry_strain = entry
							print()
						else:
							print ('\t\t+ Rename into: ', entry_strain)
													
																	 
						list_msh = HCGB_main.retrieve_matching_files(this_db, '.sig', debug)
						if (list_msh):
							## print original in file
							file2print = this_db + '/.original'
							if not os.path.exists(file2print):
								original = ['NaN']
							else:
								original = HCGB_main.readList_fromFile(file2print)
							stamp = time.time() ## to do add time
							db_Dataframe.loc[len(db_Dataframe)] = ['genbank', entry_strain, list_msh[0], stamp, this_db + '/mash/' + original[0], original[1], original[2], this_db]
						else:
							## index assembly or reads...
							list_fna = HCGB_main.retrieve_matching_files(this_db, 'genomic.fna', debug)

							## not available
							stamp = time.time() ## to do add time
							db_Dataframe.loc[len(db_Dataframe)] = ['genbank', entry_strain, 'NaN', stamp, list_fna[0], 'NaN', 'NaN', this_db]

			#### user_data
			elif (db == "user_data"):
				print (colored("\n\t- user_data: including information from user previously generated results", 'green')) ## include user data
				db2use_abs = HCGB_files.create_subfolder('user_data', database_folder)
				user_entries = os.listdir(db2use_abs)
				for entry in user_entries:
					if entry == 'user_database.csv':
						continue
					
					print ('\t+ Reading information from sample: ', entry)
					this_db = db2use_abs + '/' + entry
					this_mash_db = this_db + '/mash/' + entry + '.sig'
					if os.path.exists(this_mash_db):
						## print original in file
						file2print = this_db + '/mash/.original'
						if not os.path.exists(file2print):
							original = ['NaN', 'NaN', 'NaN']
						else:
							original = HCGB_main.readList_fromFile(file2print)

						##
						db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, this_mash_db, this_db + '/mash/' + original[0], original[1], original[2], this_db + '/mash']
					else:
						## not available
						list_fna = HCGB_main.retrieve_matching_files(this_db + '/assembly', '.fna', debug)
						db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, 'NaN', list_fna[0], 'NaN', 'NaN', this_db + '/mash']

	#### external_data
	### TODO: Fix this
	mash_bin = "" #set_config.get_exe('mash')
	if any(name in 'Mash_external' for name in db_Dataframe['source'].to_list()):
		print (colored("\t- external_data: including information from external data provided by user", 'green')) ## include user data
		db_Dataframe = db_Dataframe.set_index("db", drop = False)
		frame = db_Dataframe[ db_Dataframe['source'] == 'Mash_external' ]
		for index, row in frame.iterrows():
			print ('\t+ Reading information for file: ', row['db'])
			outfile = row['path'] + '.msh'
			if not os.path.exists(outfile):
				path_file = os.path.dirname(row['path'])
				this_db_file = min_hash_caller.sketch_database([row['path']], mash_bin, row['path'], row['db'], path_file)
				HCGB_aes.print_sepLine("*",50, False)

			db_Dataframe.loc[row['db']] = ['Mash_external', row['db'], outfile, row['path']]
	
	## index by id	
	db_Dataframe = db_Dataframe.set_index("db", drop = False)			
	return (db_Dataframe)
