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
from BacterialTyper.scripts import ariba_caller
from BacterialTyper.scripts import card_trick_caller
from BacterialTyper.config import set_config

## Import HCGB functions
import HCGB
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys


#############################################################
def parse_vfdb(folder, sampleName, fileResults, fileFlags, summary, assembly_cutoff):
	"""Parses results from VFDB database.
	
	:param folder: ARIBA results generated using VFDB database
	:param sampleName: the sample name id
	:param fileResults: report.tsv file generated
	:param fileFlags: contains information for ``ariba expand flags`` (results_parser)
	:param summary: Summary file create by ``ariba summary`` (ariba_summary)
	:param assembly_cutoff: ARIBA assembly threshold cutoff [0-1].
	
	:type folder: string
	:type sampleName: string
	:type fileResults: string
	:type fileFlags: string
	:type summary: string
	:type assembly_cutoff: float
	 
	"""
	## 
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
	colnames = ['Reference', 'ID', 'Protein-coding','Presence/Absence', 'Species', 'Variants', 'Description', 'Additional information']
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
	card_ontology = HCGB_main.get_data(card_trick_info + '/aro.obo.csv', ',', 'index_col=0') 		## read card_info generated for card_trick parse
	
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
	colnames = ['Reference', 'ID', 'Protein-coding','Presence/Absence', 'Variants', 'Description', 'Additional information']
	
	## get results: conferring resistance
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
			df_results.loc[name]['Presence/Absence'] = 'no'
		else:
			df_results.loc[name]['Presence/Absence'] = 'yes'

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
	colnames_identified = ['Status', 'ID', 'Protein-coding','Presence_Absence', 'Variants', 'pc_ident', 'pc_len']
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
			df_identified.loc[name]['Presence_Absence'] = 'no'
		else:
			df_identified.loc[name]['Presence_Absence'] = 'yes'

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
	return("", "")

#############################################################
def check_results(db2use, outdir_sample, assembly_cutoff, card_trick_info):
	
	"""
	.. seealso:: Additional information to ARIBA results generated.
	
		- :ref:`ARIBA-explained`
	
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
			if db2use == 'card':
				database = 'card'
				name_db = 'CARD'
			elif db2use == 'vfdb_full':
				database = 'vfdb_full'
				name_db = 'VFDB'
			else:
				## might generate conflicts if several other databases provided
				database = 'other'
				name_db = 'other'

			## Timestamp created in profile/ folder for each db analyzed
			filename_stamp = outfolder + '/.success_' + database
			if os.path.isfile(filename_stamp):
				stamp =	HCGB_time.read_time_stamp(filename_stamp)
				print (colored("\tA previous command generated results on: %s [%s]" %(stamp, sample), 'yellow'))
				name_excel = outfolder + '/' + sample + '_' + name_db + '_results.xlsx'
				name_csv = outfolder + '/' + sample + '_' + name_db + '_summary.csv'

			else:
				(name_excel, name_csv) = results_parser(database, folderResults, sample, outfolder, assembly_cutoff, card_trick_info)

			dataFrame_results.loc[sample] = (name_csv, name_excel, name_db) ## to return

	return (dataFrame_results)

##################################################################
def results_parser(database, folderResults, sampleName, outfolder, assembly_cutoff, card_trick_info):	
	"""Parse ARIBA results
	
	This function basically extracts files and generated additionally information for later
	parse according to type of database provided.
	
	.. seealso:: Additional information to ARIBA results generated.
	
		- :ref:`ARIBA-explained`
	"""	
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
		filePath = os.path.join(folderResults, f)
		if f.endswith('.gz'):
			HCGB_files.extract(filePath, folderResults)
		if (f=='report.tsv'):
			fileResults=filePath
		elif (f=='assemblies.fa.gz'):
			assemblies=os.path.join(folderResults, 'assemblies.fa')
		elif (f=='assembled_genes.fa.gz'):
			assemled_genes=os.path.join(folderResults, 'assembled_genes.fa')
	print ("\n")
	
	## no results generated
	if not HCGB_files.is_non_zero_file(fileResults):
		print('+ No results generated for sample: ', sampleName)
		return('','')
		
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
	stamp =	HCGB_time.print_time_stamp(filename_stamp)

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
			stamp = HCGB_time.read_time_stamp(filename_stamp)
			print ("+ A previous download generated results on: ", stamp)
			days_passed = HCGB_time.get_diff_time(filename_stamp)
			print ("\t\t** %s days ago" %days_passed)		
			if (days_passed > 30): ## download again
				print ("\t\t** Downloading information again just to be sure...")
			else:
				print ("\t\t** No need to download data again.")
				return()
	else:
		HCGB_files.create_folder(folder)
	
	## Open file and readlines
	print ('+ Downloading files:\n')
	for line in links:
		if not line.startswith('#'):  
			HCGB_sys.wget_download(line, folder)
	
	## decompress files			
	print ('+ Decompressing gzip files\n')
	files = os.listdir(folder)
	for item in files:
		#print (folder)
		if item.endswith('.gz'):
			HCGB_files.extract(folder + '/' + item, folder)

	## make stamp time
	HCGB_time.print_time_stamp(filename_stamp)
	
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

## http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID=VFG049753
## http://www.mgc.ac.cn/cgi-bin/VFs/vfs.cgi?VFID=VF0403#VF0403
