#usr/bin/env python
'''
This code prepares the database information for further analysis.
Several functions are implemented for:
	- Manipulate data entries, update information
	- Download NCBI assembly IDs provided
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
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
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import ariba_caller
from BacterialTyper import species_identification_KMA
from BacterialTyper.modules import sample_prepare
from BacterialTyper import min_hash_caller

##########################################################################################
def NCBI_DB(strains2get, data_folder, Debug):
	## set index
	strains2get = strains2get.set_index('NCBI_assembly_ID', drop=False) ## set new index but keep column
	strains2get.index.names = ['ID'] ## rename index

	#########
	if Debug:
		print (colored("DEBUG: NCBI data provided: ", 'yellow'))
		print (strains2get)

	## get data existing database
	print ("+ Create the database in folder: \n", data_folder)
	## read database 
	db_frame = getdbs('NCBI', data_folder, 'genbank', Debug)
	database_df = get_database(db_frame, Debug)
	
	#########
	if Debug:
		print (colored("DEBUG: NCBI genbank database retrieved: ", 'yellow'))
		print (database_df)
	
	## loop and download
	for index, row in strains2get.iterrows():
		functions.print_sepLine("+", 75, False)
		acc_ID = strains2get.loc[index]['NCBI_assembly_ID']
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
			this_db = functions.get_data(data_accID, ',', 'index_col=0')
			this_db = this_db.set_index('ID')
			database_df = database_df.append(this_db)

	## Generate/Update database
	database_csv = data_folder + '/genbank_database.csv'
	db_updated = update_db_data_file(database_df, database_csv)
	print ("+ Database has been generated: \n", db_updated)
	return (db_updated)

##########################################################################################
def NCBI_descendant(tax_ID, NCBI_folder, Debug):
	print ()
	
	## [TODO]	
	## get information from NCBItaxonomy from module ete3
	## and generate a dataframe containing as header: genus,species,name,NCBI_assembly_ID
	## and call NCBI_DB to download and update database, then returned dataframe.


##########################################################################################
def NCBIdownload(acc_ID, data, data_folder):	
	
	## module ngd requires to download data in bacteria subfolder under genbank folder
	dir_path = data_folder + '/genbank/bacteria/' + acc_ID 
	download = False
	if os.path.exists(dir_path):
		print ('+ Folder already exists: ', dir_path)
		## get files download
		(genome, prot, gff) = get_files_download(dir_path)
		if all([genome, prot, gff]):
			download = False
		else:
			download = True
	else:
		download = True
	
	if download:
		## download in data folder provided
		ngd.download(section='genbank', file_format='fasta,gff,protein-fasta', assembly_accessions=acc_ID, output=data_folder, group='bacteria')

		## check if files are gunzip
		files = os.listdir(dir_path)
		files_list = []		
		for f in files:
			if f.endswith('gz'):
				files_list.append(f)
				print ("\t- Extracting files: ", f)
				functions.extract(dir_path + '/' + f, dir_path)
				#os.remove(dir_path + '/' + f)
	else:
		print ('+ Data is already available, no need to download it again')

	## get files download
	(genome, prot, gff) = get_files_download(dir_path)

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
		
	data2download=pd.DataFrame(columns=('ID','folder','genus','species','name','genome', 'chr', 'GFF','proteins','plasmids_number','plasmids_ID','plasmids'))
	data2download.loc[len(data2download)] = (acc_ID, dir_path, data.loc[acc_ID]['genus'], data.loc[acc_ID]['species'], data.loc[acc_ID]['name'], genome, contig_out_file, gff, prot, plasmid_count, "::".join(plasmid_id), plasmid_out_file)

	## dump to file
	info_file = dir_path + '/info.txt'
	data2download.to_csv(info_file)
	
	## timestamp
	filename_stamp = dir_path + '/.success'
	stamp =	functions.print_time_stamp(filename_stamp)

	## return data
	return(info_file)		
	
##########################################################################################
def get_files_download(folder):
	## check if files are gunzip
	files = os.listdir(folder)
	genome=()
	prot=()
	gff=()

	for f in files:
		if f.endswith('genomic.fna'):
			genome = folder + '/' + f
		elif f.endswith('genomic.gff'):
			gff = folder + '/' + f
		elif f.endswith('protein.faa'):
			prot = folder + '/' + f

	return(genome, prot, gff)			

##########################################################################################
def get_database(db_frame, Debug):
	data4db = pd.DataFrame()
	for index, row in db_frame.iterrows():
		## information
		this_file = db_frame.loc[index]['path'] + '/info.txt'
		if os.path.isfile(this_file):
			print ('+ Reading information for sample: ', db_frame.loc[index]['db'])
			print (colored("\t+ Obtaining information from file: %s" %this_file, 'yellow'))
			this_db = functions.get_data(this_file, ',', 'index_col=0')
			data4db = data4db.append(this_db)
			timestamp = db_frame.loc[index]['path'] + '/.success'
			if os.path.isfile(timestamp):
				stamp =	functions.read_time_stamp(timestamp)
				print (colored("\t+ Data generated on: %s" %stamp, 'yellow'))

			functions.print_sepLine("*",25, False)

	## index by ID
	if not data4db.empty:
		data4db = data4db.set_index('ID')

	return(data4db)

##########################################################################################
def update_db_data_file(data, csv):
	if os.path.isfile(csv):
		print ("\n+ Updating database")
		print ("+ Obtaining information from database file: %s" %csv)
		db2update = functions.get_data(csv, ',', 'index_col=0')
		
		## TODO: provide preference to db2update
		df = pd.concat([db2update, data], join='inner', sort=True).drop_duplicates()
		df.to_csv(csv)
		return (df)
	else:
		data.to_csv(csv)
		return (data)
		
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
					db_Dataframe.loc[len(db_Dataframe)] = ['Mash_external', name_ext_, ext]
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

	####
	return(getdbs_df(source, dbs2use, database_folder, debug, db_Dataframe))

##################################################
def getdbs_df(source, dbs2use, database_folder, Debug, db_Dataframe):
	
	## init dataframe
	#colname = ["source", "db", "path"]
	#db_Dataframe  = pd.DataFrame(columns = colname)
	
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

	##################
	#### Min Hash ####
	##################
	elif (source == 'MASH'):

		db_Dataframe['original'] = ''
		db_Dataframe['ksize'] = ''
		db_Dataframe['num_sketch'] = ''
		db_Dataframe['folder'] = ''
	
		### Check if folder exists
		Mash_db_abs = functions.create_subfolder('Mash_db', database_folder)

		### get information
		for db in dbs2use:
			this_db = Mash_db_abs + '/' + db

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
						list_msh = functions.retrieve_matching_files(this_db, '.sig')
						if (list_msh):
							## print original in file
							file2print = this_db + '/.original'
							if not os.path.exists(file2print):
								original = ['NaN']
							else:
								original = functions.readList_fromFile(file2print)

							db_Dataframe.loc[len(db_Dataframe)] = ['genbank', entry, list_msh[0], original[0], original[1], original[2], this_db]
						else:
							## index assembly or reads...
							list_fna = functions.retrieve_matching_files(this_db, 'genomic.fna')

							## not available
							db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, 'NaN', list_fna[0], 'NaN', 'NaN', this_db]

			#### user_data
			elif (db == "user_data"):
				print (colored("\n\t- user_data: including information from user previously generated results", 'green')) ## include user data
				db2use_abs = functions.create_subfolder('user_data', database_folder)
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
							original = functions.readList_fromFile(file2print)

						##
						db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, this_mash_db, original[0], original[1], original[2], this_db + '/mash']
					else:
						## not available
						list_fna = functions.retrieve_matching_files(this_db + '/assembly', '.fna')
						db_Dataframe.loc[len(db_Dataframe)] = ['user_data', entry, 'NaN', list_fna[0], 'NaN', 'NaN', this_db + '/mash']

	#### external_data
	### TODO: Fix this
	mash_bin = config.get_exe('mash')
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
				functions.print_sepLine("*",50, False)

			db_Dataframe.loc[row['db']] = ['Mash_external', row['db'], outfile, row['path']]
	
	## index by id	
	db_Dataframe = db_Dataframe.set_index("db", drop = False)			
	return (db_Dataframe)
