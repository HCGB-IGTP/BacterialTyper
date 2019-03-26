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
from sys import argv
from io import open
import pandas as pd
import ncbi_genome_download as ngd
import functions
from Bio import SeqIO
import shutil

######
def NCBIdownload(data, data2download, folder):	
	
	## set index
	data = data.set_index('NCBI_assembly_ID', drop=False) ## set new index but keep column
	data.index.names = ['ID'] ## rename index
	data2download = data2download.set_index('ID')
	
	for index, row in data.iterrows():
		print ("###########################################################")
		acc_ID = data.loc[index]['NCBI_assembly_ID']
		info = "Genus: " + data.loc[index]['genus'] + '\n' + "Species: " +  data.loc[index]['species'] + '\n' + "Strain: " +  data.loc[index]['name'] + '\n' + "ID accession: " +  acc_ID + '\n'
		dir_path = folder + '/genbank/bacteria/' + acc_ID
		#print (info)

		if acc_ID in data2download.index:
			print ("\n+ Data is already available in database for: ")
			print (info)
		else:
			if os.path.exists(dir_path):
				## path exist
				print ("\n+ Data is already downloaded for: ")
				print (info)
			else:
				print ("\n+ Downloading data for:")	
				print (info)
				ngd.download(section='genbank', file_format='fasta,gff,protein-fasta', assembly_accessions=acc_ID, output=folder,  group='bacteria')


			## set database
			data2download.loc[acc_ID] = '.'
			data2download.loc[acc_ID]['folder'] = dir_path
			data2download.loc[acc_ID]['genus'] = data.loc[acc_ID]['genus']
			data2download.loc[acc_ID]['species'] = data.loc[acc_ID]['species']
			data2download.loc[acc_ID]['name'] = data.loc[acc_ID]['name']

			## check if files are gunzip
			files = os.listdir(dir_path)
			files_list = []		
			for f in files:
				if f.endswith('gz'):
					files_list.append(f)
					print ("\t- Extracting files: ", f)
					functions.extract(dir_path + '/' + f)
					os.remove(dir_path + '/' + f)

			## get files download
			(genome, prot, gff) = get_files_download(dir_path)
	
			## check if any plasmids downloaded
			plasmid_count = 0
			plasmid_id = []
		
			contig_out_file = dir_path + '/' + acc_ID + '_chromosome.fna'
			plasmid_out_file = dir_path + '/' + acc_ID + '_plasmid.fna' 

			contig_out_file_handle = open(contig_out_file, 'w')
			plasmid_out_file_handle = open(plasmid_out_file, 'w')
		
			for seq_record in SeqIO.parse(genome, "fasta"):
				plasmid_search = re.search(r".*plasmid.*", seq_record.description)
				if plasmid_search:
					plasmid_count += 1
					name = str( seq_record.id )
					plasmid_id.append(name)
				
					plasmid_out_file_handle.write(seq_record.format("fasta"))
					plasmid_out_file_handle.write('\n')
				else:
					contig_out_file_handle.write(seq_record.format("fasta"))
					contig_out_file_handle.write('\n')
				
			##
			contig_out_file_handle.close()
			plasmid_out_file_handle.close()
		
			## populate DB
			data2download.loc[acc_ID]['genome'] = contig_out_file
			data2download.loc[acc_ID]['GFF'] = gff
			data2download.loc[acc_ID]['proteins'] = prot	
			data2download.loc[acc_ID]['plasmids_number'] = plasmid_count
			data2download.loc[acc_ID]['plasmids_ID'] = "::".join(plasmid_id)		
			if plasmid_count > 0:
				data2download.loc[acc_ID]['plasmids'] = plasmid_out_file
			else:
				data2download.loc[acc_ID]['plasmids'] = ""

			print ("###########################################################\n\n")
			
	db_updated = update_db_data_file(data2download, folder)
	db_updated.to_csv(folder + "/database.csv")
	return(folder + "/database.csv")

######
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

######
def update_database(strains2get_file, folder):
	## get file information from database
	database = get_data(folder + '/database.csv')	
	strains2get = get_data(strains2get_file)	
	
	## download
	data = NCBIdownload(strains2get, database, folder)
	print ("+ Database has been updated: \n", data)
	return (data)

######
def update_db_data_file(data, folder):
	if os.path.isfile(folder + "/database.csv"):
		db2update = pd.read_csv(folder + "/database.csv", header=0)
		db2update = db2update.set_index('ID')
		df = pd.concat([db2update, data], join='inner').drop_duplicates()
		return (df)
	else:
		return (data)

######
def update_database_user_data(data, folder):

	print ("+ Updating information from assembly data...")
	data2update = get_data(folder + '/database.csv')
	data2update = data2update.set_index('ID')
	
	## create folder
	own_data = functions.create_subfolder("user_data", folder)
	
	## get info to update
	data2be_updated = get_data(data)
	data2be_updated = data2be_updated.set_index('ID')
	
	for index, row in data2be_updated.iterrows():
		
		if index in data2update.index:
			print ("\n+ Data is already available in database for:", index)
		else:
	
			print ("###########################################################")
			folder_ID = data2be_updated.loc[index]['folder']
			data2update.loc[index] = "." ## init row
			dir_sample = functions.create_subfolder(index, own_data)
		
			## get files
			(genome, prot, gff, plasmids, number, plasm_ids) = get_files_folder(folder_ID)

			## populate dataframe
			data2update.loc[index]['folder'] = data2be_updated.loc[index]['folder']
			data2update.loc[index]['genus'] = data2be_updated.loc[index]['genus']
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
		
	dataUpdated = update_db_data_file(data2update, folder)
	dataUpdated.to_csv(folder + "/database.csv")
	return(folder + "/database.csv")
	
######
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

######
def get_data(ID_file):	
	print ("+ Obtaining information from file: ", ID_file)
	data = pd.read_csv(ID_file, header=0)
	print ("\n+ Data:")
	print (data)
	print ("\n\n")	
	return(data)

######
def init_DB(ID_file, folder):
	## get file information
	strains2get = get_data(ID_file)	
	print ("+ Create the database in folder: \n", folder)	
	data2download=pd.DataFrame(columns=('ID','folder','genus','species','name','genome','GFF','proteins','plasmids_number','plasmids_ID','plasmids'))
	
	## download
	data = NCBIdownload(strains2get, data2download, folder)
	print ("+ Database has been updated: \n", data)
	return (data)

#####
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	## get arguments provided
	option = argv[1]
	folder = os.path.abspath(argv[2])
	ID_file = os.path.abspath(argv[3])
	
	if (option == "init_db"):
		init_DB(ID_file, folder)
	elif (option == "update_NCBI"):
		update_database(ID_file, folder)
	elif (option == "update_user_data"):
		update_database_user_data(ID_file, folder)

######
def help_options():
	print ("\nUSAGE: python %s option database_folder ID_file\n"  %os.path.realpath(__file__))
	print ("---------------")
	print ("+ Option: init_db | update_NCBI | update_user_data")
	print ("+ database_folder: path to the database [if does not exist will be created]")
	print ("+ ID_file is a CSV file containing different values according to the option:")
	print ("\n#####################")
	print ("Option 1: init_db & update_NCBI: genus,species,name,NCBI_assembly_ID")	
	print ("---- Example 1 -------")
	print ("genus,species,name,NCBI_assembly_ID")
	print ("Staphylococcus,aureus,MRSA252,GCx_0000001")
	print ("Staphylococcus,epidermis,strain2,GCx_0000002")
	print ("Staphylococcus,aureus,strain3,GCx_0000003")
	print ("----------------------")
	print ("#####################")
	print ("\nOption 2: update_own_data: ID,folder,genus,species,name")	
	print ("** If genus,species or names is not known just leave it blank")
	print ("---- Example 2 -------")	
	print ("ID,folder,genus,species,name")
	print ("sample1,/path/to/spades_assembly_folder1/,,,")
	print ("sample2,/path/to/spades_assembly_folder2/,,,")
	print ("----------------------")
	print ("#####################\n")
		
######
'''******************************************'''
if __name__== "__main__":
	main()
	
	
## Download all genomes from a taxa and descendent
## https://www.biostars.org/p/302533/	

