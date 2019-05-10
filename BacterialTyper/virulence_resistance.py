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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import ariba_caller

#############################################################
def parse_vfdb(folder, sampleName, fileResults, fileFlags, summary):
	print ("\tCheck VFDB result: ", sampleName)
	summary_data = pd.read_csv(summary, header=0, sep=',')
	fileFlags_data = pd.read_csv(fileFlags, header=0, sep='\t')
	original_data = pd.read_csv(fileResults, header=0, sep='\t')

	## summary data
	summary_data = summary_data.set_index('name')
	cluster_len = len(summary_data.columns)
	print ("\t%s genes putatively involved in the strain virulence..." %cluster_len)

	## subset
	data = original_data.loc[original_data['cluster'].isin(summary_data.columns)]
	
	## create dataframe
	colnames = ['Reference', 'VFDB ID', 'Protein-coding','Presence/Abscence', 'Species', 'Variants', 'Description', 'Additional information']
	df_results = pd.DataFrame(columns=colnames, index=summary_data.columns)
	cluster = data.groupby(['cluster'])
	
	for name, group in cluster:
		#print ("#####################\n")
		#print ("Group:", name)
		#print ("Dataframe:\n", group)
		#print ("#####################\n")		
		## Example
		##Original name: VFG002420(gb|NP_644838) (adsA) Adenosine synthase A [AdsA (VF0422)] [Staphylococcus aureus subsp. aureus MW2]
		text_free = "".join(set(list(group['free_text'])))
		text_search = re.search(r"Original name:\s(VF.*\d+)\((.*\d+)\)\s(\(.*)\[(.*)\]", text_free)
		if text_search:
			VF_id = text_search.group(1)
			protein_id = text_search.group(2)
			name_id = text_search.group(3)
			species_id = text_search.group(4)
		
		else:
			VF_id = '-'
			protein_id = '-'
			name_id = text_free
			species_id = '-'
		
		
		df_results.loc[name]['VFDB ID'] = VF_id
		df_results.loc[name]['Reference'] = protein_id
		df_results.loc[name]['Species'] = species_id
		df_results.loc[name]['Description'] = name_id
		
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

	#pd.set_option('display.max_colwidth', -1)
	#pd.set_option('display.max_columns', None)
	#print (df_results)

	## generate excel sheet
	name_excel = folder + '/' + sampleName + '_VFDB_summary.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') 	## open excel handle
	df_results.to_excel(writer, sheet_name='summary') 		## write excel handle
	original_data.to_excel(writer, sheet_name='VFDB_ARIBA_report') 		## write excel handle
	fileFlags_data.to_excel(writer, sheet_name='flags') 		## write excel handle
	writer.save()											## close excel handle		
	return (name_excel)


#############################################################
def parse_card(folder, sampleName, fileResults, fileFlags, summary):
	
	print ("\tCheck CARD result: ", sampleName)
	summary_data = pd.read_csv(summary, header=0, sep=',')
	fileFlags_data = pd.read_csv(fileFlags, header=0, sep='\t')
	original_data = pd.read_csv(fileResults, header=0, sep='\t')
	
	## summary data
	summary_data = summary_data.set_index('name')
	cluster_len = len(summary_data.columns)
	print ("\t%s genes putatively involved in resistance to some antibiotics..." %cluster_len)

	## subset
	data = original_data.loc[original_data['cluster'].isin(summary_data.columns)]
	
	resistance=[]
	freetext = set(list(data['free_text']))
	for text in freetext:
		search_freetext = re.search(r".*conferring resistance to (.*)", text)
		if (search_freetext):
			resistance.append(search_freetext.group(1))

	resistance_String = "; ".join(set(resistance))
	print ('\tPutative resistance(s) to: ' + resistance_String + ', ...')
	
	## create dataframe
	colnames = ['Reference', 'ARO id', 'Protein-coding','Presence/Abscence', 'Variants', 'Description', 'Additional information']
	df_results = pd.DataFrame(columns=colnames, index=summary_data.columns)
	
	## analyze each cluster
	cluster = data.groupby(['cluster'])
	for name, group in cluster:
		#print ("#####################\n")
		#print ("Group:", name)
		#print ("Dataframe:\n", group)
		#print ("#####################\n")
		ariba_ref_name = str(group['ref_name']).split('.')
		df_results.loc[name]['ARO id'] = 'ARO:' + ariba_ref_name[1]
		df_results.loc[name]['Reference'] = ariba_ref_name[2]
		
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

	#pd.set_option('display.max_colwidth', -1)
	#pd.set_option('display.max_columns', None)
	#print (df_results)

	## generate excel sheet
	name_excel = folder + '/' + sampleName + '_CARD_summary.xlsx'
	writer = pd.ExcelWriter(name_excel, engine='xlsxwriter') 	## open excel handle
	df_results.to_excel(writer, sheet_name='results') 		## write excel handle
	original_data.to_excel(writer, sheet_name='CARD_ARIBA_report') 		## write excel handle
	summary_data.to_excel(writer, sheet_name='ARIBA_summary') 		## write excel handle
	fileFlags_data.to_excel(writer, sheet_name='flags') 		## write excel handle
	writer.save()											## close excel handle	
	return (name_excel)
	
#############################################################
def parse_results(folder, sampleName, fileResults, fileFlags, summary):
	print ("")

#############################################################
def check_results(db2use, outdir_sample, outfolder):
	vfdb = False
	## iterate multi-index dataframe
	for sample, data in outdir_sample.groupby(level='sample'): ## fix
		for database, data2 in data.groupby(level='db'): ## fix
			if (database != db2use):
				continue

			folderResults = data2.loc[sample, db2use]['output']			
			if db2use.endswith('card_prepareref/'):
				results_parser('card', folderResults, sample, outfolder)
			elif db2use.endswith('vfdb_full_prepareref/'):
				vfdb = True
				results_parser('vfdb_full', folderResults, sample, outfolder)
			else:
				results_parser('.', folderResults, sample, outfolder)

	if (vfdb):
		print ("\n\n")
		functions.print_sepLine("*", 50, False)
		print ("+ Check VFDB details in files downloaded from vfdb website:")
		files_VFDB = check_VFDB(outfolder + '/VFDB_information')
		functions.print_sepLine("*", 50, False)
	return()

##################################################################
def results_parser(database, folderResults, sampleName, outfolder):	

	### folder results
	# assembled_seqs.fa.gz: reference sequences identified
	# assemblies.fa.gz: query sequences retrieved from sample
	# assembled_genes.fa.gz: encoding genes from assemblies.fa
	# debug.report.tsv: initial report.tsv before filtering
	# log.clusters.gz: log details for each cluster
	# report.tsv: report for each sample
	# version_info.txt: version and additional information
	
	##
	list_files = os.listdir(folderResults)
	
	## init
	assemblies=""
	assemled_genes=""
	fileResults=""
	print ("\n+ Parsing result file for sample: ", sampleName)
	print ("\n+ Extracting files if necessary:")
	## extract files
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
	
	## generate summary
	summary_results_tmp = folderResults + '/report_summary_tmp'
	summary_results = folderResults + '/report_summary.csv'
	ariba_caller.ariba_summary(summary_results_tmp, [fileResults])
	
	## fix names
	fake_dict = {sampleName : fileResults}
	ariba_caller.fix_ariba_summary(summary_results_tmp + '.csv', summary_results, fake_dict)
	os.remove(summary_results_tmp + '.csv')
	
	###
	if (database == 'vfdb_full'):
		name_excel = parse_vfdb(outfolder, sampleName, fileResults, fileFlags, summary_results)
	elif (database == 'card'):
		name_excel = parse_card(outfolder, sampleName, fileResults, fileFlags, summary_results)
	else: ## to do
		name_excel= parse_results(outfolder, sampleName, fileResults, fileFlags, summary_results)				
	
	print ('\tCheck additional information on ', name_excel)

#############################################################
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
	
#############################################################
def download_VFDB_files(folder):
	links = ("http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz","http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz")	
	filename_stamp = folder + '/download_timestamp.txt' 
	## check how old is the data and if it is necessary to download again
	if os.path.exists(folder):
		if os.path.isfile(filename_stamp):
			stamp = functions.read_time_stamp(filename_stamp)
			print ("+ A previous download generated results on: ", stamp)
			days_passed = functions.get_diff_time(stamp)
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
		if item.endswith('.gz'):
			functions.extract(folder + '/' + item, folder)

	## make stamp time
	functions.print_time_stamp(filename_stamp)
	return()

#############################################################
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
		
