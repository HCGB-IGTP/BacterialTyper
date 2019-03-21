#usr/bin/env python
'''
This code downloads information for NCBI assembly IDs provided and updates/populates de database of interest
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

######
def NCBIdownload(data2download, folder):	
	
	## Init columns
	data2download['Genome'] = "."
	data2download['GFF'] = "."
	data2download['Proteins'] = "."
	data2download['Plasmids_number'] = "."
	data2download['Plasmids_NCBI_ID'] = "."

	for index, row in data2download.iterrows():
		print ("###########################################################")
		acc_ID = data2download.loc[index]['NCBI_assembly_ID']
		info = "Genus: " + data2download.loc[index]['Genus'] + '\n' + "Species: " +  data2download.loc[index]['species'] + '\n' + "Strain: " +  data2download.loc[index]['name'] + '\n' + "ID accession: " +  acc_ID + '\n'
		dir_path = folder + '/genbank/bacteria/' + acc_ID
		
		if os.path.exists(dir_path):
			## path exist
			print ("\n+ Data is already available for: ")
			print (info)

		else:
			print ("\n+ Downloading data for:")	
			print (info)
			ngd.download(section='genbank', file_format='fasta,gff,protein-fasta', assembly_accessions=acc_ID, output=folder,  group='bacteria')

		## check if files are gunzip
		files = os.listdir(dir_path)
		files_list = []		
		for f in files:
			if f.endswith('gz'):
				files_list.append(f)
				print ("\t- Extracting files: ", f)
				functions.extract(dir_path + '/' + f)
				os.remove(dir_path + '/' + f)

		(genome, prot, gff) = get_files_folder(dir_path, acc_ID)
		data2download.loc[index]['Genome'] = genome
		data2download.loc[index]['GFF'] = gff
		data2download.loc[index]['Proteins'] = prot
		
		plasmid_count=0
		plasmid_id = []
		
		for seq_record in SeqIO.parse(genome, "fasta"):
			plasmid_search = re.search(r".*plasmid.*", seq_record.description)
			if plasmid_search:
				plasmid_count += 1
				name = str( seq_record.id )
				plasmid_id.append(name)
			
		data2download.loc[index]['Plasmids_number'] = plasmid_count
		data2download.loc[index]['Plasmids_NCBI_ID'] = "::".join(plasmid_id)		

		print ("###########################################################\n\n")

	db_updated = update_db_data(data2download, folder)
	db_updated.to_csv(folder + "/database.csv")
	return(folder + "/database.csv")

######
def get_files_folder(folder, ident):
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

######
def update_db_data(data, folder):
	if os.path.isfile(folder + "/database.csv"):
		db2update = pd.read_csv(folder + "/database.csv", header=0)
		df = pd.concat([db2update, data], join='inner').drop_duplicates()
		return (df)
	else:
		return (data)

######
def get_data(ID_file):	
	print ("+ Obtaining information from file:")
	print (ID_file)
	data = pd.read_csv(ID_file, header=0)
	print ("\n+ Data:")
	print (data)
	print ("\n\n")	
	return(data)

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	## get arguments provided
	ID_file = os.path.abspath(argv[1])
	folder = os.path.abspath(argv[2])
	
	## get file information
	strains2get = get_data(ID_file)
	
	## download
	data = NCBIdownload(strains2get, folder)
	
	print ("+ Database has been updated: \n", data)



######
def help_options():
	print ("\nUSAGE: python %s ID_file database_folder\n"  %os.path.realpath(__file__))
	print ("---------------")
	print ("ID_file example")
	print ("CSV file containing: genus,species name,strain name,Genbank or Refseq ID")	
	print ("---------------")
	print ("Genus,species,name,NCBI_assembly_ID")
	print ("genus,species,strain1,GCx_0000001")
	print ("genus,species,strain2,GCx_0000002")
	print ("genus,species,strain3,GCx_0000003")
	print ("genus,species,strain4,GCx_0000004")
	print ("genus,species,strain5,GCx_0000005")
	print ("*******")
	print ("______________")
		
######
'''******************************************'''
if __name__== "__main__":
	main()
	
	
## Download all genomes from a taxa and descendent
## https://www.biostars.org/p/302533/	

