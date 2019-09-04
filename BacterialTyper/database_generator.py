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

## import data
dataDir = os.path.dirname(os.path.realpath(__file__)) + '/../../data/'
plasmid_groups = dataDir + '/available_plasmids_data.txt'

##########################################################################################
def NCBI_DB(ID_file, data_folder, Debug):
	## get file information
	strains2get = functions.get_data(ID_file, ',', '')	

	## set index
	strains2get = strains2get.set_index('NCBI_assembly_ID', drop=False) ## set new index but keep column
	strains2get.index.names = ['ID'] ## rename index

	#########
	if Debug:
		print (colored("DEBUG: NCBI data provided: ", 'yellow'))
		print (strains2get)

	print ("+ Create the database in folder: \n", data_folder)	
	data2download=pd.DataFrame(columns=('ID','folder','genus','species','name','genome', 'chr', 'GFF','proteins','plasmids_number','plasmids_ID','plasmids'))
	#data2download = data2download.set_index('ID')
	print ()
	
	## get data existing database
	database_df = get_database(data_folder, Debug)
	
	#########
	if Debug:
		print (colored("DEBUG: NCBI genbank database retrieved: ", 'yellow'))
		print (database_df)
	
	## loop and download
	for index, row in strains2get.iterrows():
		functions.print_sepLine("+", 75, False)
		acc_ID = strains2get.loc[index]['NCBI_assembly_ID']
		info = "Genus: " + strains2get.loc[index]['##genus'] + '\n' + "Species: " +  strains2get.loc[index]['species'] + '\n' + "Strain: " +  strains2get.loc[index]['name'] + '\n' + "ID accession: " +  acc_ID + '\n'
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
	data2download.loc[len(data2download)] = (acc_ID, dir_path, data.loc[acc_ID]['##genus'], data.loc[acc_ID]['species'], data.loc[acc_ID]['name'], genome, contig_out_file, gff, prot, plasmid_count, "::".join(plasmid_id), plasmid_out_file)

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
def get_database(folder, Debug):
	## read database 
	db_frame = database.getdbs('NCBI', folder, 'genbank', Debug)
	data4db = pd.DataFrame()
	for index, row in db_frame.iterrows():
		print ('+ Reading information for sample: ', db_frame.loc[index]['db'])
		timestamp = db_frame.loc[index]['path'] + '/.success'
		if os.path.isfile(timestamp):
			stamp =	functions.read_time_stamp(timestamp)
			print (colored("\t+ Data downloaded on: %s" %stamp, 'yellow'))

		this_file = db_frame.loc[index]['path'] + '/info.txt'
		print (colored("\t+ Obtaining information from file: %s" %this_file, 'yellow'))
		this_db = functions.get_data(this_file, ',', 'index_col=0')
		data4db = data4db.append(this_db)

	data4db = data4db.set_index('ID')
	return(data4db)

##########################################################################################
def update_db_data_file(data, csv):
	if os.path.isfile(csv):
		print ("+ Updating database")
		print ("+ Obtaining information from database file: %s" %csv)
		db2update = functions.get_data(csv, ',', 'index_col=0')
		df = pd.concat([db2update, data], join='inner').drop_duplicates()
		df.to_csv(csv)
		return (df)
	else:
		data.to_csv(csv)
		return (data)


##########################################################################################
def update_database_user_data(data, folder):

	print ("+ Updating information from user data...")
	data2update = functions.get_data(folder + '/database.csv', '\t')
	data2update = data2update.set_index('ID')
	
	## create folder
	own_data = functions.create_subfolder("user_data", folder)
	
	## get info to update
	data2be_updated = functions.get_data(data)
	data2be_updated = data2be_updated.set_index('ID')
	
	for index, row in data2be_updated.iterrows():
		
		if index in data2update.index:
			print ("\n+ Data is already available in database for:", index)
		else:
			folder_ID = data2be_updated.loc[index]['folder']
			data2update.loc[index] = "." ## init row
			dir_sample = functions.create_subfolder(index, own_data)
		
			## get files
			(genome, prot, gff, plasmids, number, plasm_ids) = get_files_folder(folder_ID)

			## populate dataframe
			data2update.loc[index]['folder'] = data2be_updated.loc[index]['folder']
			data2update.loc[index]['genus'] = data2be_updated.loc[index]['##genus']
			data2update.loc[index]['species'] = data2be_updated.loc[index]['species']
			data2update.loc[index]['name'] = data2be_updated.loc[index]['name']
			data2update.loc[index]['ID'] = index
		
			# copy file
			genome_name =  dir_sample + '/' + os.path.basename(genome)
			shutil.copy(genome, genome_name)
			data2update.loc[index]['genome'] = genome_name
		
			# copy file
			gff_name =  dir_sample + '/' + os.path.basename(gff)
			shutil.copy(gff, gff_name)
			data2update.loc[index]['GFF'] = gff_name

			# copy file
			prot_name =  dir_sample + '/' + os.path.basename(prot)
			shutil.copy(prot, prot_name)
			data2update.loc[index]['proteins'] = prot_name

			data2update.loc[index]['plasmids_number'] = number
			data2update.loc[index]['plasmids_ID'] = plasm_ids		


			# copy file
			if (os.path.isfile(plasmids)):
				plasm_name =  dir_sample + '/' + os.path.basename(plasmids)
				shutil.copy(plasmids, plasm_name)
				data2update.loc[index]['plasmids'] = plasm_name
			else:
				data2update.loc[index]['plasmids'] = ""

			functions.print_sepLine("+", 75, False)

		
	dataUpdated = update_db_data_file(data2update, folder)
	dataUpdated.to_csv(folder + "/database.csv")
	return(folder + "/database.csv")

##########################################################################################
def get_files_folder(folder):
	## check if files are gunzip
	files = os.listdir(folder)
	genome=""
	plasmids=""
	number=0
	plasm_ids=[]
	prot=()
	gff=()

	for f in files:
		if f.endswith('plasmid.fna'):
			plasmids = folder + '/' + f
			for seq_record in SeqIO.parse(plasmids, "fasta"):
				number += 1
				name = str( seq_record.id )
				plasm_ids.append(name)			
		elif f.endswith('chromosome.fna'):
			genome = folder + '/' + f
		elif f.endswith('.gff'):
			gff = folder + '/' + f
		elif f.endswith('.faa'):
			prot = folder + '/' + f
			
	return (genome, prot, gff, plasmids, number, "::".join(plasm_ids))

## Download all genomes from a taxa and descendent
## https://www.biostars.org/p/302533/	 ## > ncbi.get_descendant_taxa() 


##########################################################################################
def plasmidID_user_data(folder):
	print ("+ Adding available plasmid sequences to the plasmid database:")

	## to do....
	##

##########################################################################################
def plasmidID_db_NCBI(path, name):

	print ("+ Preparing folder for downloading data:")
	plasmid_data = functions.create_subfolder("plasmid_data", path)

	## name = Firmicutes,Alphaproteobacteria or all
	seq_folder, info = download_plasmid_NCBI(plasmid_data, name)
	print ("+ Information from plasmids downloaded is in file: " + info)
	
	#
	plasmids_fna = plasmid_data + '/plasmids_database.fna'
	functions.concat_fasta(seq_folder, plasmids_fna)
	print ("+ Plasmids are availabe in file: " + plasmids_fna)

	return(plasmids_fna)


##########################################################################################
def download_plasmid_NCBI(folder, name):

	print ("+ Downloading information from available plasmids on NCBI...")
	path = functions.create_subfolder("download", folder)

	## NCBI Genome Plasmids ftp site.
	functions.wget_download("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt", path)
	
	plasmids_info = path + '/plasmids.txt'
	plasmids = functions.get_data(plasmids_info, '\t')
	sequence_folder = functions.create_subfolder("seqs", path)
	
	## get desirable group
	names = []
	if (name == 'all'):
		plasmids_filter = plasmids
	else:
		names = name.split(",")
		plasmids_filter = plasmids.loc[plasmids['SubGroup'].isin(names)]

	## 
	ids_download = path + "/plasmid_ids_Downloaded.txt"
	out_handle = open(ids_download, "a")

	## loop the table and download nucleotide sequences
	count = 0
	for index, entry in plasmids_filter.iterrows():
		#print (plasmids_filter.iloc[index])
		sequence = plasmids_filter.loc[index]['RefSeq']
		
		### get id to download	
		if (sequence == '-'):
			continue
		else:
			## file
			seq = sequence_folder + '/' + sequence + '.fna'
		
			## download if not available
			if os.path.isfile(seq):
				print ("%s is already available in folder...." %sequence)
			else:
				## eutils NCBI
				cmd='curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=' + sequence + '&amp;retmode=text&amp;rettype=fasta" > ' + seq
				functions.system_call(cmd)
				string = sequence + '\t' + plasmids_filter.loc[index]['Kingdom'] + '\t'+  plasmids_filter.loc[index]['#Organism/Name'] + '\t' + plasmids_filter.loc[index]['Group'] +'\t' + plasmids_filter.loc[index]['SubGroup'] + '\t' + plasmids_filter.loc[index]['Plasmid Name']
				out_handle.write(string)
				out_handle.write("\n")
				count += 1
	
	out_handle.close()
	
	print ("\n\n+ %s sequences have been downloaded from NCBI belonging to: %s" %(count, names))
	return(sequence_folder, ids_download)


##########################################################################################
def update_database(strains2get_file, folder):
	## get file information from database
	database = functions.get_data(folder + '/database.csv', ',')	
	strains2get = functions.get_data(strains2get_file, ',')	
	
	## download
	data = NCBIdownload(strains2get, database, folder)
	print ("+ Database has been updated: \n", data)
	return (data)


