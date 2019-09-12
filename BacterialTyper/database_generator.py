#usr/bin/env python
'''
This code prepares the database information for further analysis.
Several functions are implemented for:
- If NCBI assembly IDs provided, downloaded them updates/populates the database of interest
- If own assembly data provided, updates/populates the database of interest
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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import ariba_caller
from BacterialTyper import species_identification_KMA
from BacterialTyper.modules import database
from BacterialTyper.modules import sample_prepare

## import data
dataDir = os.path.dirname(os.path.realpath(__file__)) + '/../../data/'
plasmid_groups = dataDir + '/available_plasmids_data.txt'

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
	db_frame = database.getdbs('NCBI', folder, 'genbank', Debug)
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
	
	## download in data folder provided
	ngd.download(section='genbank', file_format='fasta,gff,protein-fasta', assembly_accessions=acc_ID, output=data_folder, group='bacteria')

	## module ngd requires to download data in bacteria subfolder under genbank folder
	dir_path = data_folder + '/genbank/bacteria/' + acc_ID 

	## check if files are gunzip
	files = os.listdir(dir_path)
	files_list = []		
	for f in files:
		if f.endswith('gz'):
			files_list.append(f)
			print ("\t- Extracting files: ", f)
			functions.extract(dir_path + '/' + f, dir_path)
			#os.remove(dir_path + '/' + f)

	## get files download
	(genome, prot, gff) = get_files_download(data_folder)

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
		if f.endswith('.fna'):
			genome = folder + '/' + f
		elif f.endswith('.gff'):
			gff = folder + '/' + f
		elif f.endswith('.faa'):
			prot = folder + '/' + f
	return(genome, prot, gff)			

##########################################################################################
def get_database(db_frame, Debug):
	data4db = pd.DataFrame()
	for index, row in db_frame.iterrows():
		print ('+ Reading information for sample: ', db_frame.loc[index]['db'])
		timestamp = db_frame.loc[index]['path'] + '/.success'
		if os.path.isfile(timestamp):
			stamp =	functions.read_time_stamp(timestamp)
			print (colored("\t+ Data generated on: %s" %stamp, 'yellow'))

		## information
		this_file = db_frame.loc[index]['path'] + '/info.txt'
		if os.path.isfile(this_file):
			print (colored("\t+ Obtaining information from file: %s" %this_file, 'yellow'))
			this_db = functions.get_data(this_file, ',', 'index_col=0')
			data4db = data4db.append(this_db)

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
		
		## TODO: join and provide preference to db2update
		df = pd.concat([db2update, data], join='inner', sort=True).drop_duplicates()
		df.to_csv(csv)
		return (df)
	else:
		data.to_csv(csv)
		return (data)

##########################################################################################
def get_userData_files(options, project_folder):
	## get information regarding files

	## get assembly files
	pd_samples_assembly = sample_prepare.get_files(options, project_folder, "assembly", "fna")
	pd_samples_assembly = pd_samples_assembly.set_index('name')

	## get annotation files
	pd_samples_annot = sample_prepare.get_files(options, project_folder, "annot", ['gbf', 'faa', 'gff'])
	pd_samples_annot = pd_samples_annot.set_index('name')

	## get trimmed ngs files
	pd_samples_reads = sample_prepare.get_files(options, project_folder, "trim", ['_trim_'])
	pd_samples_reads = pd_samples_reads.set_index('name')

	## debug message
	if (options.debug):
		print (colored("**DEBUG: pd_samples_reads **", 'yellow'))
		print (pd_samples_reads)
		print (colored("**DEBUG: pd_samples_assembly **", 'yellow'))
		print (pd_samples_assembly)
		print (colored("**DEBUG: pd_samples_annot **", 'yellow'))
		print (pd_samples_annot)
	
	## merge
	df = pd.concat([pd_samples_reads, pd_samples_annot, pd_samples_assembly], join='inner').drop_duplicates()
	## joining by inner we only get common columns among all

	##
	return(df)
	
##########################################################################################
def get_userData_info(options, project_folder):
	## get information regarding: 
		## genus, species (ident module)
		## card & VFDB (profile module)
		## additional information: MGE, etc

	## get profile information
	pd_samples_profile = sample_prepare.get_files(options, project_folder, "profile", ["csv"])
	pd_samples_profile = pd_samples_profile.set_index('name')

	## get identification information
	pd_samples_ident = sample_prepare.get_files(options, project_folder, "ident", ["csv"])
	pd_samples_ident = pd_samples_ident.set_index('name')

	## add other if necessary

	## debug message	
	if (options.debug):
		print (colored("**DEBUG: pd_samples_profile **", 'yellow'))
		print (pd_samples_profile)
		print (colored("**DEBUG: pd_samples_ident **", 'yellow'))
		print (pd_samples_ident)
		## add other if necessary
	
	## merge
	df = pd.concat([pd_samples_profile, pd_samples_ident], join='inner', sort=True).drop_duplicates()
	## joining by inner we only get common columns among all

	return(df)

##########################################################################################
def update_database_user_data(database_folder, project_folder, Debug, options):

	print ("\n+ Updating information from user data folder:")
	print (project_folder)
	
	## create folder
	own_data = functions.create_subfolder("user_data", database_folder)
	
	## Default missing options
	options.project = True
	options.debug = Debug
	if not options.single_end:
		options.pair = True

	## get user data files
	project_data_df = get_userData_files(options, project_folder)
	
	## get user data info
	project_info_df = get_userData_info(options, project_folder)
	
	## merge data
	project_all_data = pd.concat([project_data_df, project_info_df], join='inner').drop_duplicates()

	## read database 
	db_frame = database.getdbs('user_data', database_folder, 'user_data', Debug)
	user_data_db = get_database(db_frame, Debug)
	
	## merge dataframe
	sample_frame = project_all_data.groupby("name")
	for name, cluster in sample_frame:

		## debug message	
		if (options.debug):
			print (colored("**DEBUG: sample_frame groupby: name & cluster **", 'yellow'))
			print (name)			
			print (cluster)
	
		############################################
		#### check information for this sample
		############################################

		## generate sample
		dir_sample = functions.create_subfolder(name, own_data)

		if name in user_data_db.index:
			print (colored("+ Data is already available in database for: %s" %name, 'green'))
		else:
		
			## data to generate
			data2dump = pd.DataFrame(columns=('ID','folder','genus','species','name','genome', 'GFF','proteins'))

			## iterate over files with different tags: reads, annot, assembly, profile, ident

			##########
			## assembly
			##########
			assembly_dir = functions.create_subfolder('assembly', dir_sample)
			assembly_file = cluster.loc[cluster['tag'] == 'assembly']['sample'].to_list()
			shutil.copy(assembly_file[0], assembly_dir)
			genome = assembly_file[0]
		
			##########
			## annot
			##########
			annot_dir = functions.create_subfolder('annot', dir_sample)
			annot_files = cluster.loc[cluster['tag'] == 'annot']['sample'].to_list()
			prof = ""
			gff = ""
			for f in annot_files:
				shutil.copy(f, annot_dir)
				if f.endswith('faa'):
					prot = f
				elif f.endswith('gff'):
					gff = f
		
			##########
			## trimm
			##########
			trimm_dir = functions.create_subfolder('trimm', dir_sample)
			reads_files = cluster.loc[cluster['tag'] == 'reads']['sample'].to_list()
			for f in reads_files:
				shutil.copy(f, trimm_dir)

			## profile and ident information would place in folders accordingly but would no be incorporated
			## in the database file
			
			##########
			## ident
			##########
			ident_dir = functions.create_subfolder('ident', dir_sample)
			ident_file = cluster.loc[cluster['tag'] == 'ident']['sample'].to_list()
			shutil.copy(ident_file[0], ident_dir)
			
			##########
			## profile
			##########
			profile_dir = functions.create_subfolder('profile', dir_sample)
			profile_files = cluster.loc[cluster['tag'] == 'profile']['sample'].to_list()
			for f in profile_files:
				shutil.copy(f, profile_dir)
			
			############################################
			### Dump information
		
			#####
			data2dump.loc[len(data2dump)] = (name, dir_sample, 'genus', 'species', name, genome, gff, prot)

			###### dump to file
			info_file = dir_sample + '/info.txt'
			data2dump.to_csv(info_file)
	
			###### dump file information to file
			info_file2 = dir_sample + '/info_files.txt'
			cluster.to_csv(info_file2)
		
			###### timestamp
			filename_stamp = dir_sample + '/.success'
			stamp =	functions.print_time_stamp(filename_stamp)

			###### populate dataframe
			user_data_db.append(data2dump)
			functions.print_sepLine("+", 75, False)

	## update db		
	if (options.debug):
		print (colored("**DEBUG: user_data_db dataframe **", 'yellow'))
		print (user_data_db)
	
	database_csv = own_data + '/user_database.csv'
	dataUpdated = update_db_data_file(user_data_db, database_csv)
	print ("+ Database has been generated: \n", database_csv)
	return (dataUpdated)

