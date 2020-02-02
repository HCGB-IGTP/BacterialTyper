#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Parses information and generates a virulence and an antibiotic resistance profile using ARIBA results.
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
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import ariba_caller
from BacterialTyper import card_trick_caller

#############################################################
def parse_vfdb(folder, sampleName, fileResults, fileFlags, summary, assembly_cutoff):
	## 
	## Parses results from VFDB database.
	## Input is a folder for output results, sample name, and
	## fileResults is report.tsv generated by ariba
	## fileFlags was generated by results_parser calling ariba expand flags, 
	## summary is the summary file generated by ariba summary
	##

	## get data
	summary_data = pd.read_csv(summary, header=0, sep=',') 			## report_summary.csv :: parse information from ARIBA 
	fileFlags_data = pd.read_csv(fileFlags, header=0, sep='\t')		## flags_explain.tsv :: ariba expand flag: explained flags
	original_data = pd.read_csv(fileResults, header=0, sep='\t')	## report.tsv :: ariba report generated

	## summary data
	summary_data = summary_data.set_index('name')
	list_found_genes = summary_data.columns
	cluster_len = len(list_found_genes)

	# print info
	print ("\tCheck VFDB result: ", sampleName)
	print ("\t%s genes putatively involved in the strain virulence..." %cluster_len)

	############################################################################	
	## analyze each cluster confering virulence
	############################################################################	
	data = original_data.loc[original_data['cluster'].isin(list_found_genes)]
	colnames = ['Reference', 'ID', 'Protein-coding','Presence/Abscence', 'Species', 'Variants', 'Description', 'Additional information']
	df_results = found_results(colnames, data, list_found_genes, 'VFDB')
	
	############################################################################	
	## analyze each cluster identified: might or might not confer virulence
	############################################################################	
	## get results: found, identified, partial
	df_identified = identified_results(original_data, "VFDB", list_found_genes, assembly_cutoff)

	##########################
	## generate excel sheet
	##########################
	name_excel = folder + '/' + sampleName + '_VFDB_results.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') 	## open excel handle

	## write excel handle
	df_results.to_excel(writer, sheet_name='results')				## write results
	original_data.to_excel(writer, sheet_name='VFDB_ARIBA_report') 	## Original data from ARIBA
	fileFlags_data.to_excel(writer, sheet_name='flags') 			## ARIBA flags explained
	summary_data.to_excel(writer, sheet_name='ARIBA_summary') 		## ARIBA summary generated
	df_identified.to_excel(writer, sheet_name='identified') 		## Identified genes: ARIBA flags explained
	
	name_csv = folder + '/' + sampleName + '_VFDB_summary.csv'
	df_identified.to_csv(name_csv)

	## close excel handle		
	writer.save()											

	return (name_excel, name_csv)

#############################################################
def parse_card(folder, sampleName, fileResults, fileFlags, summary, assembly_cutoff, card_trick_info):
	## 
	## Parses results from CARD database.
	## Input is a folder for output results, sample name, and
	## fileResults is report.tsv generated by ariba
	## fileFlags was generated by results_parser calling ariba expand flags, 
	## summary is the summary file generated by ariba summary
	##
	
	## get data
	summary_data = pd.read_csv(summary, header=0, sep=',') 			## report_summary.csv :: parse information from ARIBA 
	fileFlags_data = pd.read_csv(fileFlags, header=0, sep='\t')		## flags_explain.tsv :: ariba expand flag: explained flags
	original_data = pd.read_csv(fileResults, header=0, sep='\t')	## report.tsv :: ariba report generated
	card_ontology = functions.get_data(card_trick_info + '/aro.obo.csv', ',', 'index_col=0') 		## read card_info generated for card_trick parse
	
	## summary data
	summary_data = summary_data.set_index('name')
	list_found_genes = summary_data.columns
	cluster_len = len(list_found_genes)
	
	## print info
	print ("\tCheck CARD result: ", sampleName)
	print ("\t%s genes putatively involved in resistance to some antibiotics..." %cluster_len)

	## subset
	data = original_data.loc[original_data['cluster'].isin(summary_data.columns)]
	
	############################################################################	
	## analyze each cluster confering resistance
	############################################################################	
	colnames = ['Reference', 'ID', 'Protein-coding','Presence/Abscence', 'Variants', 'Description', 'Additional information']
	
	## get results: confering resistance
	df_results = found_results(colnames, data, list_found_genes, 'CARD')

	## get results: found, identified, partial
	df_identified = identified_results(original_data, "CARD", list_found_genes, assembly_cutoff)

	############################################################################
	## use card-trick python package to get ontology for each term
	AROS_identified = list(df_identified['ID'])
	information_ontology = card_trick_caller.get_info_CARD(AROS_identified, 'ARO', card_ontology)

	##########################
	## generate excel sheet
	##########################

	## open excel handle
	name_excel = folder + '/' + sampleName + '_CARD_results.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') 		
	
	## write excel handle
	df_results.to_excel(writer, sheet_name='results') 					## write results
	df_identified.to_excel(writer, sheet_name='identified') 			## Identified genes: ARIBA flags explained
	information_ontology.to_excel(writer, sheet_name='CARD_ontology') 	## CARD ontology
	original_data.to_excel(writer, sheet_name='ARIBA_report') 			## Original data from ARIBA
	summary_data.to_excel(writer, sheet_name='ARIBA_summary') 			## ARIBA summary generated
	fileFlags_data.to_excel(writer, sheet_name='flags') 				## ARIBA flags explained
	
	name_csv = folder + '/' + sampleName + '_CARD_summary.csv'
	df_identified.to_csv(name_csv)
	
	## close excel handle
	writer.save()														

	return (name_excel, name_csv)

########################################
def found_results(colnames, dataFrame, list_found_genes, db2use_name):
	##
	## For the information provided by ARIBA summary, get only genes that confer resistance/virulence
	## according to the database due to presence or identified variants
	##
	
	## create dataframe for parsing results and later printing
	df_results = pd.DataFrame(columns=colnames, index=list_found_genes)
	
	## loop
	cluster = dataFrame.groupby(['cluster'])
	for name, group in cluster:
		
		## get ID according to database
		info = get_id(db2use_name, group)
		
		if (db2use_name == 'CARD'):
			df_results.loc[name]['ID'] = info['ID']		
			df_results.loc[name]['Reference'] = info['reference']				

		elif (db2use_name == 'VFDB'):
			df_results.loc[name]['ID'] = info['ID']
			df_results.loc[name]['Reference'] = info['protein_id']
			df_results.loc[name]['Description'] = info['name_id']
			df_results.loc[name]['Species'] = info['species_id']
		else:
			df_results.loc[name]['ID'] = info['ID']		
		
		## coding or non-coding
		sum_gene = int(group['gene'].sum())
		if (sum_gene >= 1):
			df_results.loc[name]['Protein-coding'] = 'yes'
		else:
			df_results.loc[name]['Protein-coding'] = 'no'
			
		## presence_absence
		summary_var_Only=0
		for var_Only in group['var_only']:
			if (var_Only == '.'):
				var_Only = 0
			summary_var_Only = summary_var_Only + int(var_Only)			

		if (summary_var_Only >= 1):
			df_results.loc[name]['Presence/Abscence'] = 'no'
		else:
			df_results.loc[name]['Presence/Abscence'] = 'yes'

		## sum types: has_known_var
		summary = 0
		for var_name in group['has_known_var']:
			if (var_name == '.'):
				var_name = 0
			summary = summary + int(var_name)		

		if (summary >= 1):
			hash_known_var_data = group.loc[group['has_known_var'] == '1']
			df_results.loc[name]['Variants'] = "; ".join(list(hash_known_var_data['known_var_change']))
			
			var_description = list(hash_known_var_data['var_description'])
			desc_list = []
			for desc in var_description:
				text = desc.split(":.:")[1]
				desc_list.append(text)
			df_results.loc[name]['Additional information'] = "; ".join(set(desc_list))
			
		else:
			df_results.loc[name]['Variants'] = summary
			df_results.loc[name]['Additional information'] = '-'

		## other
		description = "; ".join(set(list(group['free_text'])))
		description = description.replace('b\'', '')
		df_results.loc[name]['Description'] = description

	## debugging
	#pd.set_option('display.max_colwidth', -1)
	#pd.set_option('display.max_columns', None)
	#print (df_results)

	return (df_results)

########################################
def get_id(db2use_name, group):
	## 
	## Retrieve gene ID of a given CARD, VFDB or any other database
	## 
	info = {}

	### Get ID according to Database
	if db2use_name == 'CARD':
		ariba_ref_name = str(group['ref_name']).split('.')
		info['ID'] = 'ARO:' + ariba_ref_name[1]
		info['reference'] = ariba_ref_name[2]

	elif db2use_name == 'VFDB':
		## Get VFDB id
		text_free = "".join(set(list(group['free_text'])))
		text_search = re.search(r"Original name:\s(VF.*\d+)\((.*\d+)\)\s(\(.*)\[(.*)\]", text_free)
		if text_search:
			info['ID'] = text_search.group(1)
			info['protein_id'] = text_search.group(2)
			info['name_id'] = text_search.group(3)
			info['species_id'] = text_search.group(4)
		else:
			info['ID'] = '-'
			info['protein_id'] = '-'
			info['name_id'] = text_free
			info['species_id'] = '-'

	else: 
		## different database
		info['ID'] = '-'
	
	return (info)

########################################
def identified_results(original_data, db2use_name, list_found_genes, assembly_threshold):

	######################################################################################	
	## analyze each cluster identified: might or might not confer resistance/virulence
	######################################################################################
	##
	## For the information provided by ARIBA original data, get all genes identified
	## (might produce or not resistance/virulence due to presence or identified variants)
	##
	## Generate tag:
	##	identified: identified gene as confering resistance/virulence
	##	found: found gene and assembly but not confering resistance/virulence
	##	partial: not fulfilling > assembly_threshold % 		
	## 
	found_genes_cluster = original_data.groupby(['cluster'])
	index_names = set(original_data['cluster'].to_list())

	## create dataframe for parsing results and later printing
	colnames_identified = ['Status', 'ID', 'Protein-coding','Presence_Abscence', 'Variants', 'pc_ident', 'pc_len']
	df_identified = pd.DataFrame(columns=colnames_identified, index=index_names)
	df_identified.index.names = ['Genes']

	## loop
	for name, group in found_genes_cluster:

		## get ID according to database
		info = get_id(db2use_name, group)
		df_identified.loc[name]['ID'] = info['ID']		
	
		## coding or non-coding
		sum_gene = int(group['gene'].sum())
		if (sum_gene >= 1):
			df_identified.loc[name]['Protein-coding'] = 'yes'
		else:
			df_identified.loc[name]['Protein-coding'] = 'no'
			
		## presence_absence
		summary_var_Only=0
		for var_Only in group['var_only']:
			if (var_Only == '.'):
				var_Only = 0
			summary_var_Only = summary_var_Only + int(var_Only)
			
		if (summary_var_Only >= 1):
			df_identified.loc[name]['Presence_Abscence'] = 'no'
		else:
			df_identified.loc[name]['Presence_Abscence'] = 'yes'

		## sum types: has_known_var
		summary = 0
		for var_name in group['has_known_var']:
			if (var_name == '.'):
				var_name = 0
			summary = summary + int(var_name)		

		if (summary >= 1):
			df_identified.loc[name]['Variants'] = 'yes'
		else:
			df_identified.loc[name]['Variants'] = 'no'
			
		## pc_ident
		pc_ident = group['pc_ident'].to_list()
		df_identified.loc[name]['pc_ident'] = pc_ident[0]

		## pc_len = (ref_base_assembled / ref_len)*100
		ref_base = group['ref_base_assembled'].to_list()
		ref_len = group['ref_len'].to_list()
		pc_float = float(ref_base[0]/ref_len[0])
		pc_len = pc_float*100		
		df_identified.loc[name]['pc_len'] = pc_len
		
		### status
		if any(name in s for s in list_found_genes):
			#print ("Gene: ", name, "-- Identify to confer resistance")
			df_identified.loc[name]['Status'] = 'Identified'
		else:
			if pc_float > float(assembly_threshold):
				#print ("Gene: ", name, "-- Found, but not confering resistance.")
				df_identified.loc[name]['Status'] = 'Found'
			else:
				#print ("Gene: ", name, "-- Found, but not confering resistance -- Partial.")
				df_identified.loc[name]['Status'] = 'Partial'
	
	### debugging
	##print (df_identified)
	return (df_identified)
	
#############################################################
def parse_results(folder, sampleName, fileResults, fileFlags, summary):
	## [TODO]
	print (colored("\n\n***** TODO: Implement if a different database provided (!= CARD, VFDB) *****\n\n", 'red'))
	#return(name_excel, name_csv)

#############################################################
def check_results(db2use, outdir_sample, assembly_cutoff, card_trick_info):
	
	"""
	.. seealso:: Additional information to ARIBA results generated.
	
		- :doc:`ARIBA results explained <../../../data/ARIBA_explained>`
	
	"""
	
	## 
	## outdir_sample is a dataframe containing information of the output folder generated by ariba. 
	## It is index for each database and for each sample.
	## This function iterates for each sample and generates call to specific function to parse results.
	## 

	## iterate multi-index dataframe
	dataFrame_results = pd.DataFrame(columns=("csv", "excel", "database"))
	for sample, data in outdir_sample.groupby(level='sample'):
		for database, data2 in data.groupby(level='db'): 
			if (database != db2use):
				continue

			folderResults = data2.loc[sample, db2use]['output']
			outfolder = data2.loc[sample, db2use]['dirname']
			if db2use.endswith('card_prepareref/'):
				database = 'card'
				name_db = 'CARD'
			elif db2use.endswith('vfdb_full_prepareref/'):
				database = 'vfdb_full'
				name_db = 'VFDB'
			else:
				database = 'other'
				name_db = 'other'

			## might generate conflicts if several other databases provided
			## TODO: check
			filename_stamp = outfolder + '/.success_' + database
			if os.path.isfile(filename_stamp):
				stamp =	functions.read_time_stamp(filename_stamp)
				print (colored("\tA previous command generated results on: %s [%s]" %(stamp, sample), 'yellow'))
				name_excel = outfolder + '/' + sample + '_' + name_db + '_results.xlsx'
				name_csv = outfolder + '/' + sample + '_' + name_db + '_summary.csv'

			else:
				(name_excel, name_csv) = results_parser(database, folderResults, sample, outfolder, assembly_cutoff, card_trick_info)

			dataFrame_results.loc[sample] = (name_csv, name_excel, name_db) ## to return

	return (dataFrame_results)

##################################################################
def results_parser(database, folderResults, sampleName, outfolder, assembly_cutoff, card_trick_info):	
	"""
	.. seealso:: Additional information to ARIBA results generated.
	
		- :doc:`ARIBA results explained <../../../data/ARIBA_explained>`
	
	"""
	
	## 
	## Parse ARIBA results generated for all databases. There is a lot of information generated
	## that is common for all databases and listed as an example in folder results here.
	##
	# folder results
	# |-- assembled_seqs.fa.gz: reference sequences identified
	# |-- assemblies.fa.gz: query sequences retrieved from sample
	# |-- assembled_genes.fa.gz: encoding genes from assemblies.fa
	# |-- debug.report.tsv: initial report.tsv before filtering
	# |-- log.clusters.gz: log details for each cluster
	# |-- report.tsv: report for each sample
	# |-- version_info.txt: version and additional information
	##
	## This function basically extracts files and generated additionally information for later
	## parse according to type of database provided.
	## Input is the type of database, the folder of results, name and output folder
	##
	
	if not os.path.exists(folderResults):
		print ("+ Finish parsing information for sample [%s]. Results folder does not exist." %sampleName)
		return('NaN','NaN')
		
	## get files
	list_files = os.listdir(folderResults)

	## init
	assemblies=""
	assemled_genes=""
	fileResults=""

	print ("\n+ Parsing result file for sample: ", sampleName)

	## extract files
	print ("\n+ Extracting files if necessary:")
	for f in list_files:
		if f.endswith('.gz'):
			functions.extract(folderResults + '/' + f, folderResults)
		if (f=='report.tsv'):
			fileResults=folderResults + '/' + f
		elif (f=='assemblies.fa.gz'):
			assemblies=folderResults + '/assemblies.fa'
		elif (f=='assembled_genes.fa.gz'):
			assemled_genes=folderResults + '/assembled_genes.fa'
	print ("\n")
	
	### expand flags
	flagResults = folderResults + '/flags_explain.tsv'
	fileFlags = ariba_caller.ariba_expandflag(fileResults, flagResults)
	
	######################
	## generate summary
	######################
	##
	## ariba has function that generates a summary for samples
	##
	summary_results_tmp = folderResults + '/report_summary_tmp'
	summary_results = folderResults + '/report_summary.csv'
	options = "--no_tree"
		## Info
		## https://github.com/sanger-pathogens/ariba/wiki/The-assembled-column-from-ariba-summary
	
	ariba_caller.ariba_summary(summary_results_tmp, [fileResults], options)
	
	## fix names: just for aesthetics
	fake_dict = {sampleName : fileResults}
	ariba_caller.fix_ariba_summary(summary_results_tmp + '.csv', summary_results, fake_dict)
	os.remove(summary_results_tmp + '.csv')
	
	############################################
	### check results according to database
	############################################
	if (database == 'vfdb_full'):
		(name_excel, name_csv) = parse_vfdb(outfolder, sampleName, fileResults, fileFlags, summary_results, assembly_cutoff)
	elif (database == 'card'):
		(name_excel, name_csv) = parse_card(outfolder, sampleName, fileResults, fileFlags, summary_results, assembly_cutoff, card_trick_info)
	else: 
		## [TODO] check results according to databases different than CARD/VFDB
		(name_excel, name_csv) = parse_results(outfolder, sampleName, fileResults, fileFlags, summary_results)				
	
	print ('\tCheck additional information on ', name_excel)
	
	## print success timestamp
	filename_stamp = outfolder + '/.success_' + database
	stamp =	functions.print_time_stamp(filename_stamp)

	return (name_excel, name_csv)

#############################################################
def check_VFDB(VFDB_path):
	## 
	## Given a folder, check/download VFDB information
	## and print information.
	## Returns files downloaded or available
	## 
	
	## download
	download_VFDB_files(VFDB_path)

	## Print info files
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
	
#############################################################
def download_VFDB_files(folder):
	## 
	## Given a folder, check if it contains VFDB information
	## or download it from website: http://www.mgc.ac.cn
	## 
	links = ("http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz","http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz")	

	## check if data is downloaded, how old is the data and if it is necessary to download again
	## consider >30 days long enough to be updated again
	
	## time stamp
	filename_stamp = folder + '/download_timestamp.txt' 
	if os.path.exists(folder):
		if os.path.isfile(filename_stamp):
			stamp = functions.read_time_stamp(filename_stamp)
			print ("+ A previous download generated results on: ", stamp)
			days_passed = functions.get_diff_time(stamp)
			print ("+ %s days ago" %days_passed)		
			if (days_passed > 30): ## download again
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
	
	## decompress files			
	print ('+ Decompressing gzip files\n')
	files = os.listdir(folder)
	for item in files:
		#print (folder)
		if item.endswith('.gz'):
			functions.extract(folder + '/' + item, folder)

	## make stamp time
	functions.print_time_stamp(filename_stamp)
	
	return()

#############################################################
def subsetting_species(species, folder, excel_file):
	## 
	## Given a Species name, this function subsets information from VFDB
	## and generates excel file from input name provided (excel_file)
	## 	

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

	#sampleName = "test_sample"
	#folder = argv[4]
	#parse_vfdb(folder, sampleName, argv[1], argv[2], argv[3])
	
######
def help_options():
	print ("\nUSAGE: python %s species path VFDB_path\n"  %os.path.realpath(__file__))
		
######
'''******************************************'''
if __name__== "__main__":
	main()
		
