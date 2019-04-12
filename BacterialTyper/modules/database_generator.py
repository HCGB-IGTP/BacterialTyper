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
from Bio import SeqIO
import shutil

## functions
pythonDir = os.path.dirname(os.path.realpath(__file__)) + '/../tools/python'
sys.path.append(pythonDir)
import functions

## import data
dataDir = os.path.dirname(os.path.realpath(__file__)) + '/../data/'
plasmid_groups = dataDir + '/available_plasmids_data.txt'

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
			
			## open
			contig_out_file_handle = open(contig_out_file, 'w')
		
			for seq_record in SeqIO.parse(genome, "fasta"):
				plasmid_search = re.search(r".*plasmid.*", seq_record.description)
				if plasmid_search:
					plasmid_count += 1
					name = str( seq_record.id )
					plasmid_id.append(name)
				
					plasmid_out_file_handle = open(plasmid_out_file, 'a')
					plasmid_out_file_handle.write(seq_record.format("fasta"))
					plasmid_out_file_handle.write('\n')
					plasmid_out_file_handle.close()
				else:
					contig_out_file_handle.write(seq_record.format("fasta"))
					contig_out_file_handle.write('\n')
				
			##
			contig_out_file_handle.close()

		
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
	database = get_data(folder + '/database.csv', '\t')	
	strains2get = get_data(strains2get_file, '\t')	
	
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

	print ("+ Updating information from user data...")
	data2update = get_data(folder + '/database.csv', '\t')
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

###
def plasmidID_user_data(folder):
	print ("+ Adding available plasmid sequences to the plasmid database:")


####
def concat_fasta(dirFasta, Fasta):
	print ("+ Concatenating all information into one file...")	
	cmd = 'cat ' + dirFasta + '/*fna > ' + Fasta
	return(functions.system_call(cmd))

####
def plasmidID_db_NCBI(path, name):

	print ("+ Preparing folder for downloading data:")
	plasmid_data = functions.create_subfolder("plasmid_data", path)

	## name = Firmicutes,Alphaproteobacteria or all
	seq_folder, info = download_plasmid_NCBI(plasmid_data, name)
	print ("+ Information from plasmids downloaded is in file: " + info)
	
	#
	plasmids_fna = plasmid_data + '/plasmids_database.fna'
	concat_fasta(seq_folder, plasmids_fna)
	print ("+ Plasmids are availabe in file: " + plasmids_fna)

	return(plasmids_fna)

####
def download_plasmid_NCBI(folder, name):

	print ("+ Downloading information from available plasmids on NCBI...")

	path = functions.create_subfolder("download", folder)


	## NCBI Genome Plasmids ftp site.
	functions.wget_download("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt", path)
	
	plasmids_info = path + '/plasmids.txt'
	plasmids = get_data(plasmids_info, '\t')
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
def get_data(ID_file, SEP):	
	print ("+ Obtaining information from file: ", ID_file)
	data = pd.read_csv(ID_file, header=0, sep=SEP)
	print ("\n+ Data:")
	print (data)
	print ("\n\n")	
	return(data)

######
def init_DB(ID_file, folder):
	## get file information
	strains2get = get_data(ID_file, ',')	
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
	
	if (option == "init_db"):
		ID_file = os.path.abspath(argv[3])
		init_DB(ID_file, folder)
	elif (option == "update_NCBI"):
		ID_file = os.path.abspath(argv[3])
		update_database(ID_file, folder)
	elif (option == "update_user_data"):
		ID_file = os.path.abspath(argv[3])
		update_database_user_data(ID_file, folder)
	elif (option == "plasmidID_db_NCBI"):
		plasmidID_db_NCBI(folder, argv[3])
	elif (option == "plasmidID_user_data"):
		plasmidID_user_data(folder)

######
def help_options():

	print ("\nUSAGE: python %s option database_folder [ID_file]\n"  %os.path.realpath(__file__))
	print ("+ Option:\n\t|-init_db\n\t|-update_NCBI\n\t|-update_user_data\n\t|-plasmidID_db_NCBI\n\t|-plasmidID_user_data\n")
	print ("+ database_folder: path to the database [if does not exist will be created]")
	print ("\n+ ID_file is a CSV file containing different values according to the option:")
	print ("\n#####################")
	print (" Option 1: init_db & update_NCBI:")
	print ("#####################")
	print ("ID_file format: genus,species,name,NCBI_assembly_ID")	
	print ("---- Example 1 -------")
	print ("genus,species,name,NCBI_assembly_ID")
	print ("Staphylococcus,aureus,MRSA252,GCx_0000001")
	print ("Staphylococcus,epidermis,strain2,GCx_0000002")
	print ("Staphylococcus,aureus,strain3,GCx_0000003")
	print ("----------------------")
	print ("\n#####################")
	print (" Option 2: update_own_data:")
	print ("#####################")
	print ("ID_file format: ID,folder,genus,species,name")	
	print ("** If genus,species or names is not known just leave it blank")
	print ("---- Example 2 -------")	
	print ("ID,folder,genus,species,name")
	print ("sample1,/path/to/spades_assembly_folder1/,,,")
	print ("sample2,/path/to/spades_assembly_folder2/,,,")
	print ("----------------------")
	print ("\n#####################")
	print (" Option 3: plasmidID_db_NCBI")
	print ("#####################")
	print ("+ No ID_file provided. ")
	print ("+ Provide a comma separated string instead with any of the groups below or 'all' to download them all:")
	get_data(plasmid_groups, ',')
	print ("\n\n+ Example: python %s plasmidID_db_NCBI database_folder Firmicutes,Alphaproteobacteria" %os.path.realpath(__file__))
	print ("----------------------")
	print ("\n#####################")
	print (" Option 4: plasmidID_user_data")
	print ("#####################")
	print ("+ No ID_file provided. ")
	print ("+ Provide word: user. Plasmids previously identified from user samples would be added to the dabase generated for plasmidID identification.")
	print ("----------------------")
	
		
######
'''******************************************'''
if __name__== "__main__":
	main()
	
	
## Download all genomes from a taxa and descendent
## https://www.biostars.org/p/302533/	

