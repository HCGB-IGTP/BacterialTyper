 #!/usr/bin/env python3
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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import ariba_caller

######
def parse_vfdb(fileResults):
	print ("")

######
def parse_card(fileResults):
	print ("+ Parsing CARD result file")
	
	data = pd.read_csv(fileResults, header=0, sep='\t')
	
	print (data)
	
	cluster_len = len(data.groupby(['cluster']))
	print ("\n%s clusters generated...\n" %cluster_len)
	
	print (data.groupby(['cluster', 'var_type'])['reads'].count())
	
	#cluster = data.groupby(['cluster'])
	#for cl in cluster:
	#	print (cl)
	#	var_type_len = len(cl.groupby(['var_type']))
	
	
	## ariba_ref_name	ref_name	gene	var_only	flag	reads	cluster	ref_len	ref_base_assembled	pc_ident	ctg	ctg_len	ctg_cov	known_var	var_type	var_seq_type	known_var_change	has_known_var	ref_ctg_change	ref_ctg_effect	ref_start	ref_end	ref_nt	ctg_start	ctg_end	ctg_nt	smtls_total_depth	smtls_nts	smtls_nts_depth	var_description	free_text
	## murA.3003776.BX571856.1.2257045_2258311.4633	murA.3003776.BX571856.1.2257045_2258311.4633	1	1	27	2998	murA_2	1266	1266	97.24	murA_2.l15.c4.ctg.1	2087	212.7	1	SNP	p	V65L	0	.	.	193	195	GTT	542	544	GTT	304;306;307	G;T;T	285;276;275	murA.3003776.BX571856.1.2257045_2258311.4633:1:1:V65L:.:b'murA or UDP-N-acetylglucosamine enolpyruvyl transferase catalyses the initial step in peptidoglycan biosynthesis and is inhibited by fosfomycin. Overexpression of murA through mutations confers fosfomycin resistance.'	Staphylococcus aureus murA with mutation conferring resistance to fosfomycin

	### CARD
	### murA.3003776.BX571856.1.2257045_2258311.4633
	# murA --> gene
	# 3003776 --> ARO class
	# BX571856.1 --> Reference genome id Genbank
	# 2257045_2258311 --> Start & End coordinates
	# 4633
	

######
def parse_results(fileResults):
	print ("")

######
def results_parser(database, fileResults):	
	###
	if (database == 'vfdb_full'):
		parse_vfdb(fileResults)
	elif (database == 'card'):
		parse_card(fileResults)
	else:
		parse_results(fileResults)				

######
def check_VFDB(VFDB_path):
	download_VFDB_files(VFDB_path)
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
	print ('\nFor further details please visit: <http://www.mgc.ac.cn/VFs/main.htm>')
	return (VFDB_files)
	
######
def download_VFDB_files(folder):
	links = ("http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz","http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz")	
	filename_stamp = folder + '/download_timestamp.txt' 
	## check how old is the data and if it is necessary to download again
	if os.path.exists(folder):
		if os.path.isfile(filename_stamp):
			st_hd = open(filename_stamp, 'r')
			st = st_hd.read()
			st_hd.close()
			stamp = datetime.fromtimestamp(float(st)).strftime('%Y-%m-%d %H:%M:%S')
			print ("+ A previous download generated results on: ", stamp)
			time_today = time.time()
			elapsed = time_today - float(st)
			days_passed = int((elapsed/3600)/24)		
			print ("+ %s days ago" %days_passed)		
			if (days_passed > 30):
				print ("+ Downloading information again just to be sure...")
			else:
				print ("+ No need to download data again.")
				return()
	else:
		functions.create_folder(folder)
	
	## Open file and readlines
	print ('+ Downloading files:\n')
	for line in links:
		if not line.startswith('#'):  
			functions.wget_download(line, folder)
	
	print ('+ Decompressing gzip files\n')
	## decompress files			
	files = os.listdir(folder)
	for item in files:
		#print (folder)
		functions.extract(folder + '/' + item, folder)

	## make stamp time
	timefile = open(filename_stamp, 'w')    
	string2write = str(time.time())
	timefile.write(string2write)
	return()

######
def subsetting_species(species, folder, excel_file):
	print ('+ Subsetting information for species of interest')
	subset_xlsx_file = folder + '/VFs_' + species + '.xls'
	
	## subset
	dfs = pd.read_excel(excel_file, sheet_name="VFs", index_cols=0, skiprows=1)
	df2 = dfs.set_index("Bacteria", drop = False) 				## order
	df_species = df2.loc[species:species] 						## filter by species
	df_species2 = df_species.set_index("VFID", drop = False) 	## order
	ls_items = list(df_species2)								## get column headers
	df_species2.to_excel(subset_xlsx_file, sheet_name=species, header=ls_items)

######
def check_CARD(folder):
	##
	# to do
	# parse download file
	# parse html file
	# link = "https://card.mcmaster.ca/download"
	# get timestamp for different entries
	# get newest
	# tmp_folder = functions.create_subfolder("tmp", folder)
	# functions.wget(link, tmp_folder)
	
	link = "https://card.mcmaster.ca/download/0/broadstreet-v3.0.2.tar.gz"
	
	## folder would be in ARIBA databases, CARD folder additional data
	functions.wget(link, folder)
	
	## extract files
	files = os.listdir(folder)
	for item in files:
		functions.extract(folder + '/' + item, folder)
	
		
######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
		
	#species = argv[1]
	#folder = os.path.abspath(argv[2])
	#VFDB_path = os.path.abspath(argv[3])
	
	## Start
	## check if files are downloaded first	
	#VFDB_files = check_VFDB(VFDB_path)
	
	## subset for genera of interest
	#functions.create_folder(folder)
	
	# subset for species of interest	
	#subsetting_species(species, folder, VFDB_files['VFs-xls'])

	parse_card(argv[1])

######
def help_options():
	print ("\nUSAGE: python %s species path VFDB_path\n"  %os.path.realpath(__file__))
		
######
'''******************************************'''
if __name__== "__main__":
	main()
		
