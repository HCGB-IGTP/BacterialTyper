#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Generates a multi locus sequence typing (MLST) using MLSTar R package.
- Downloads profiles
- Searches assembly fasta sequences
- Plots MLST profiles
"""

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
import pandas as pd
from termcolor import colored
import shutil

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table
from matplotlib.backends.backend_pdf import PdfPages

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.config import install_dependencies
from BacterialTyper.other_tools import tools
from BacterialTyper.data import data_files

## R scripts
MLSTarR_script = tools.R_scripts('MLSTar_call')
MLSTarR_plot = tools.R_scripts('MLSTar_plot')
MLSTarR_download_seq = tools.R_scripts('MLSTar_downloadPubMLST_seq')
MLSTarR_download_prf = tools.R_scripts('MLSTar_downloadPubMLST_profile')
MLSTarR_getpubmlst = tools.R_scripts('MLSTar_getpubmlst')

##########################################
def help_MLSTar():
	'''
	Provides help when script call as a single script.
	'''
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
	get_MLSTar_species()
	print ("")

##########################################
def run_MLSTar(species_mlst_folder, rscript, species, scheme, name, path, fileGiven, threads):
	
	## fileGiven = fasta file

	## PubMLST folder
	#pubmlst_folder = functions.create_subfolder('PubMLST', database_folder)
	
	## species folder
	#species_mlst_folder = functions.create_subfolder(species, pubmlst_folder)
	
	## scheme folder
	scheme_name = 'scheme_' + str(scheme)
	scheme_folder = functions.create_subfolder(scheme_name, species_mlst_folder)

	## seq/profile folder
	seq_folder = functions.create_subfolder('seq', scheme_folder)
	profile_folder = functions.create_subfolder('prf', scheme_folder)

	print ("##################################################")
	print ("+ MLST profiling for sample: %s" %name)
	print ("##################################################")

	## check if profile and sequences are already downloaded
	download_PubMLST(profile_folder, scheme, seq_folder, rscript, species)

	## call MLSTar for this sample
	results = run_doMLST(profile_folder, seq_folder, name, rscript, path, fileGiven, threads)
	
	return (results, profile_folder)

##########################################
def run_doMLST(profile_folder, seq_folder, name, rscript, path, fileGiven, threads):

	print ('+ Generating profile for sample...')

	folder_results = os.path.join(path, name  + '_alleles') 
	
	## success timestamp
	filename_stamp = path + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
		res_file = path + '/' + name + "_MLST_results.csv"
		return(res_file)
		
	else:
		if os.path.exists(folder_results):
			shutil.rmtree(folder_results)
			
		logFile = path + '_logFile.txt'
		cmd_profiler = "%s %s --dir_profile %s --dir_seq %s --file %s --dir %s --name %s --threads %s --lib.loc %s 2> %s" %(rscript, MLSTarR_script, profile_folder, seq_folder, fileGiven, path, name, threads, get_MLSTar_package_installed(), logFile)
		callCode = functions.system_call(cmd_profiler)

	if callCode == 'OK':
		res_file = path + '/' + name + "_MLST_results.csv"
		stamp =	functions.print_time_stamp(filename_stamp)
		return (res_file)
	else:
		return ('FAIL')

##########################################
def update_MLSTar_profile_alleles():
	return("")
	# [TODO: update_MLSTar_profile_alleles():
	## check folder *_alleles/ and add information *_MLST_results.csv
	
##########################################
def get_MLSTar_species(genus, species):

	"""
	Retrieve the correct name within `PubMLST databases`_ for the given species and/or genus specified.
	
	:param genus: Genus name.
	:param species: Species name.
	:type genus: string
	:type species: string

	:returns: Name in the PubMLST database for the taxa of interest. Returns `NaN` if not available. 
	
	Available datasets are stored in file :file:`/data/PubMLST_datasets.csv`. See available :doc:`PubMLST datasets<../../../data/PubMLST_datasets>`.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.data.data_files.data_list`
	
	"""
	
	## MLSTar available data
	MLSTar_species = data_files.data_list("PubMLST_datasets")
	
	# pandas from csv file
	data = pd.read_csv(MLSTar_species, names=["Kingdom","name","Species"], sep=",")
	
	## check if name matches	
	taxa_name = genus + ' ' + species
	sp_exists = data.loc[data["Species"] == taxa_name]['name']
	if not sp_exists.empty: ## check if exists
		return (sp_exists.values[0])
	else:
		## TODO: Check if it works
		## Check if there is a genus entry in database
		genus_name = genus + ' spp.'
		genus_exist = data.loc[data["Species"] == genus_name]['name']

		if not genus_exist.empty: ## check if exists
			return (genus_exist.values[0])
			
	## we have a problem...
	## neither a taxa or genus exists...
	return ('NaN')
		
##########################################
def getPUBMLST(species, rscript, out_name):
	"""
	Using `MLSTar software`_ retrieve for the given `species` the available schemes in PubMLST_.  

	It generates information in file `out_name` in csv format. 
	
	See example in :file:`/devel/results/getPubMLST_example.csv`
	
	.. include:: ../../devel/results/getPubMLST_example.csv
		:literal:
	   
	:param species: Sample name or tag to identify sample
	:param rscript: Path to Rscript to generate system call.
	:param out_name: Output file to generate profile information.

	:type species: string
	:type rscript: string
	:type out_name: string

	:return: OK/FAIL
	:warnings: Returns **FAIL** if process stopped.
	
	.. include:: ../../links.inc

	"""
	
	MLSTarR_getpubmlst = tools.R_scripts('MLSTar_getpubmlst')

	## species is a comma separated string
	path_package = get_MLSTar_package_installed()
	#MLSTar_getpubmlst
	cmd_getPUBMLST = "%s %s --species %s --output %s --lib.loc %s 2> /dev/null" %(rscript, MLSTarR_getpubmlst, species, out_name, path_package)
	return(functions.system_call(cmd_getPUBMLST))

##########################################
def plot_MLST(results, profile, rscript):
	path_package = get_MLSTar_package_installed()
	path_folder = os.path.dirname(results)
	cmd_plotter = "%s %s --output %s --folder_profile %s --file_result %s --lib.loc %s" %(rscript, MLSTarR_plot, path_folder, profile, results, path_package)
	return(functions.system_call(cmd_plotter))
	
##########################################
def download_PubMLST(profile_folder, scheme, seq_folder, rscript, species):

	## Check if profile exists
	file_prof = profile_folder + '/profile_scheme' + str(scheme) + '.tab'

	######################
	## download profile ##
	######################
	if os.path.exists(file_prof):
		print ('+ Profile file for scheme exists...')

		## check if previously download profile and succeeded
		filename_stamp = profile_folder + '/.success'

		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, species, 'profile'), 'yellow'))
		#
		#else
		# [TODO: Check time passed and download again if >?? days passed]
		
	else:
		print ('+ Profile file for scheme will be downloaded...')
		print ('+ Downloading profile...')
		logFile = profile_folder + '_download.log'
		cmd_profile = "%s %s --species %s --scheme %s --dir_profile %s --lib.loc %s 2> %s" %(rscript, MLSTarR_download_prf, species, scheme, profile_folder, get_MLSTar_package_installed(), logFile)
		callCode = functions.system_call(cmd_profile)
		
		if (callCode == 'OK'):
			## success timestamp
			filename_stamp = profile_folder + '/.success'
			stamp =	functions.print_time_stamp(filename_stamp)
	
	#######################
	## download sequence ##
	#######################
	seq_bool = 0
	if os.path.exists(seq_folder):
		print ('+ Sequence folder exists...')
		
		## check if previously download sequence and succeeded
		filename_stamp = seq_folder + '/.success'
		
		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, species,'sequence'), 'yellow'))

			########################################################################
			#else
			# [TODO: Check time passed and download again if >?? days passed]
			########################################################################
			#files = os.listdir(seq_folder)
			#count_fas = 0		
			#for f in files:
			#	if f.endswith('.fas'):
			#		count_fas += 1
			#	else:				
			#		os.remove(seq_folder + '/' + f)
			#
			#if count_fas > 6:
			#	print ("+ Assuming sequences are previously downloaded...")
			#else:
			#	seq_bool = 1		
			#	os.rmdir(seq_folder)
			#	seq_folder_path = functions.create_folder(seq_folder)
			
		else:
			print ('+ Sequence files for scheme will be downloaded...')
			print ('+ Downloading sequences...')
			logFile = seq_folder + '_download.log'
			cmd_seq = "%s %s --species %s --scheme %s --dir_seq %s --lib.loc %s 2> %s" %(rscript, MLSTarR_download_seq, species, scheme, seq_folder, get_MLSTar_package_installed(), logFile)
			callCode = functions.system_call(cmd_seq)
			
			if (callCode == 'OK'):
				## success timestamp
				filename_stamp = seq_folder + '/.success'
				stamp =	functions.print_time_stamp(filename_stamp)

def get_MLSTar_package_installed(debug=False):

	install_path = set_config.R_package_path_installed()
	(check_install_system, check_install_path) = set_config.get_check_R_files()
	R_script_exe = set_config.get_exe('Rscript')
	
	if debug:
		print ('\n+ Check package: MLSTar')
		
	## ATTENTION: optparse library missing, no installation withn conda
	
	## first try to check if package available in system
	cmd_check = R_script_exe + ' ' + check_install_system + ' -l MLSTar'
	code = functions.system_call(cmd_check, message=False, returned=False)
	if (code=='OK'):
		return('system')
	else:
		## check if installed in path
		cmd_check_path = R_script_exe + ' ' + check_install_path + ' -l MLSTar -p ' + install_path
		code2 = functions.system_call(cmd_check_path, message=False, returned=False)

		if (code2=='OK'):
			return(install_path)
		else:
			(install_R, install_github_package) = install_dependencies.get_install_R_files()
			cmd_R = '%s %s -l iferres/MLSTar -p %s' %(R_script_exe, install_github_package, install_path)
			code3 = functions.system_call(cmd_R, message=False, returned=False)
			if (code3):
				return (install_path)
			else:
				print ('ERROR')
				exit()	

##########################################
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	#fileGiven = os.path.abspath(argv[1])
	genome = os.path.abspath(argv[1]) ## assembly
	species = argv[2]
	rscript = argv[3]
	database = os.path.abspath(argv[4])
	name = argv[5]
	path = os.path.abspath(argv[6])
	
	##
	threads = 1
	functions.create_folder(path)
	
	## PubMLST folder
	pubmlst_folder = functions.create_subfolder('PubMLST', database)
	
	## species folder
	species_mlst_folder = functions.create_subfolder(species, pubmlst_folder)

	## output file
	output_file = species_mlst_folder + '/PubMLST_available_scheme.csv'
	
	### get scheme available
	getPUBMLST(species, rscript, output_file)
	
	## parse PubMLST results	
	scheme = 1
		
	## call MLST
	(results, profile_folder) = run_MLSTar(species_mlst_folder, rscript, species, scheme, name, path, genome, threads)
		
	## plot MLST
	plot_MLST(results, profile_folder, rscript)

##########################################
def help_options():
	print ("\nUSAGE: python %s genome species rscript_bin database name path\n"  %os.path.realpath(__file__))

		
##########################################
'''******************************************'''
if __name__== "__main__":
	main()


#	fileName = "MLSTar_species"
#	data.to_csv(fileName + '.txt', sep='\t')
#
#	## plot table
#	pdf_name = fileName + '.pdf'
#	pp = PdfPages(pdf_name)
#
#	fig,ax = plt.subplots()
#	fig.patch.set_visible(False)
#	
#	ax.axis('off')
#	ax.axis('tight')
#	
#	## color according to kingdom
#	colors = data.applymap(lambda x: 'palegreen' if x== 'bacteria' else ('lightyellow' if x== 'eukarya' else ('lightsalmon' if x=='other' else 'white' )))
#	## https://www.rapidtables.com/web/color/html-color-codes.html
#	
#	tab = ax.table(
#		cellText=data.values,
#		colLabels=data.columns,
#		cellColours =colors.values,
#		colWidths=[0.26 for x in data.columns],
#		loc='center', colLoc = 'center', rowLoc='left', cellLoc='center')
#	tab.auto_set_font_size(False)
#	tab.set_fontsize(5)
#
#	fig.tight_layout()	
#	pp.savefig()
#	plt.close()
#	pp.close()	
