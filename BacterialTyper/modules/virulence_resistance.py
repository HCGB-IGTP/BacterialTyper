#usr/bin/env python
'''
This code generates a virulence and an antibiotic resistance profile.
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
from sys import argv
from io import open
from datetime import datetime

## import modules
pythonDir = os.path.dirname(os.path.realpath(__file__)) + '/../tools/python'
sys.path.append(pythonDir)
import functions

## config
configDir = os.path.dirname(os.path.realpath(__file__)) + '/../config/'
sys.path.append(configDir)
import config

############
def printNote():
	print ("\n******\nNOTE:\nCore dataset (dataset_A) includes genes associated with experimentally verified VFs only.")
	print ("Full dataset (dataset_B) covers all genes related to known and predicted VFs in the database.\n******\n")

############
def subset_VFs_species(species, subset_xlsx_file, VF_xlsx_file):
	
	dfs = pd.read_excel(VF_xlsx_file, sheet_name="VFs", index_cols=0, skiprows=1)
	df2 = dfs.set_index("Bacteria", drop = False) 				## order
	df_species = df2.loc[species:species] 						## filter by species
	df_species2 = df_species.set_index("VFID", drop = False) 	## order
	ls_items = list(df_species2)								## get column headers
	df_species2.to_excel(name, sheet_name=species, header=ls_items)
	
######
def check_VFDB(VFDB_path, links_VFDB):
	if os.path.exists(VFDB_path):
		code = check_download_VFDB(VFDB_path)
	else:
		functions.create_folder(VFDB_path)
		code = 'FAIL'
	
	## download file if necessary
	if code == 'FAIL':
		download_VFDB_files(links_VFDB, VFDB_path)
	
	## get files
	files = os.listdir(VFDB_path)
	VFDB_files = {}
	print ('+ Several files are available with multiple information:')
	for item in files:
		if (item == 'VFs.xls'):
			VFDB_files['VFs-xls'] = VFDB_path + '/' + item
			print ('\tVFs.xls: Virulence Factor description file')

		if (item == 'Comparative_tables_from_VFDB'):
			VFDB_files['Comparative_tables'] = VFDB_path + '/' + item
			print ('\tComparative_tables_from_VFDB: Intra-genera comparison tables')
			
		if (item == 'VFDB_setA_pro.fas'):
			VFDB_files['VFDB_setA'] = VFDB_path + '/' + item
			print ('\tVFDB_setA_pro.fas: Protein sequences of core dataset')

		if (item == 'VFDB_setB_pro.fas'):
			VFDB_files['VFDB_setB'] = VFDB_path + '/' + item
			print ('\tVFDB_setB_pro.fas: Protein sequences of full dataset')
	
	printNote()
	print ('\nFor further details please visit: <http://www.mgc.ac.cn/VFs/main.htm>')
	
	return (VFDB_files)
	
######
def download_VFDB_files(links_file, path):
	links_file_handle = open(links_file)
	text = links_file_handle.read()
	lines = text.splitlines()
	
	## Open file and readlines
	print ('+ Downloading files:\n')
	for line in lines:
		if not line.startswith('#'):  
			functions.wget_download(line, path)
	
	print ('+ Decompressing gzip files\n')
	## decompress files			
	files = os.listdir(path)
	for item in files:
		functions.extract(path + '/' + item)

	## make stamp time
	filename = path + '/README.txt' 
	timefile = open(filename, 'w')    
	string2write = str(time.time())
	timefile.write(string2write)

######
def check_download_VFDB(path):
	filename = path + '/README.txt'	
	if os.path.isfile(filename):
		time_handle = open(filename)
		text = time_handle.read()	
		date_time = datetime.fromtimestamp(float(text))
		d = date_time.strftime("%d %B, %Y")
		print("Files were downloaded on:", d)
		return ('OK')
	else:
		return ('FAIL')

######
def subsetting_species(species, path, VFDB_files):
	print ('+ Subsetting information for species of interest')
	subset_xlsx_file = path + '/VFs_' + species + '.xls'
	subset_VFs_species(species, subset_xlsx_file, VFDB_files['VFs-xls'])
	
	## subset A
	subsetA_out = path + '/' + species + '_subsetA_prot.fasta'
	functions.subset_fasta(species, VFDB_files['VFDB_setA'], subsetA_out)

	## subset B
	subsetB_out = path + '/' + species + '_subsetB_prot.fasta'
	functions.subset_fasta(species, VFDB_files['VFDB_setB'], subsetB_out)

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
		
	species = argv[1]
	path = argv[2]
	VFDB_path = argv[3]
	links_VFDB = argv[4]
	
	## Start
	## check if files are downloaded first	
	VFDB_files = check_VFDB(VFDB_path, links_VFDB)
	
	## subset for genera of interest
	if os.path.exists(path):
		print ('')
	else:
		functions.create_folder(path)
	
	# subset for species of interest	
	subsetting_species(species, path, VFDB_files)


######
def help_options():
	print ("\nUSAGE: python %s species path VFDB_path VFs_links_file\n"  %os.path.realpath(__file__))
		
######
'''******************************************'''
if __name__== "__main__":
	main()
		
