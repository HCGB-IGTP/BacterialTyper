#!/usr/bin/env python3
'''
This code generates a bacteriophage identification profile for each sample
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
from termcolor import colored

#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
#from pandas.plotting import table
#from matplotlib.backends.backend_pdf import PdfPages

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper.other_tools import tools
from BacterialTyper.data import data_files

## R scripts
MLSTarR_script = tools.R_scripts('MLSTar_call')
MLSTarR_plot =tools.R_scripts('MLSTar_plot')
MLSTarR_download_seq = tools.R_scripts('MLSTar_downloadPubMLST_seq')
MLSTarR_download_prf = tools.R_scripts('MLSTar_downloadPubMLST_profile')

######
def help_MLSTar():

	print ("")
	get_MLSTar_species()
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

######
def run_MLSTar(species, scheme, name, path, fileGiven, threads):

	profile_folder = config.MLSTar['profile_folder']
	seq_folder = config.MLSTar['sequence_folder']
	rscript = config.EXECUTABLES['Rscript']
	
	## check if profile and sequences are already downloaded
	download_PubMLST(profile_folder, scheme, seq_folder, name, rscript, species)

	## call MLSTar for this sample
	results = run_doMLST(profile_folder, seq_folder, name, rscript, path, fileGiven, threads)
	return (results)

######
def run_doMLST(profile_folder, seq_folder, name, rscript, path, fileGiven, threads):
	print ('+ Create sample folder...')
	folder_path = functions.create_subfolder(name, path)

	print ('+ Generating profile for sample...')
	cmd_profiler = "%s %s --dir_profile %s --dir_seq %s --file %s --dir %s --name %s --threads %s" %(rscript, MLSTarR_script, profile_folder, seq_folder, fileGiven, folder_path, name, threads)
	callCode = functions.system_call(cmd_profiler)

	if callCode == 'OK':
		res_file = folder_path + '/' + name + "_results.txt"
		return (res_file)
	else:
		exit()

######
def update_MLSTar_profile_alleles():
	return("")
	
######
def get_MLSTar_species():

	## MLSTar available data
	MLSTar_species = data_files.data_list("MLSTar_species")

	# pandas from csv file
	print (MLSTar_species)
	data = pd.read_csv(MLSTar_species, header=0, sep=",")
	print (data)
	return(data)

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

######
def call_plot(results):
	profile_folder = config.MLSTar['profile_folder']
	rscript = config.EXECUTABLES['Rscript']
	plot_MLST(results, profile_folder, rscript)
	
######
def plot_MLST(results, profile, rscript):
	path_folder = os.path.dirname(results)
	cmd_plotter = "%s %s --output %s --folder_profile %s --file_result %s" %(rscript, MLSTarR_plot, path_folder, profile, results)
	return(functions.system_call(cmd_plotter))
	
######
def download_PubMLST(profile_folder, scheme, seq_folder, name, rscript, species):

	print ("##################################################")
	print ("+ MLST profiling for sample: %s" %name)
	print ("##################################################")
	
	## Check if profile exist
	file_prof = profile_folder + '/profile_scheme' + str(scheme) + '.tab'
	
	#print (file_prof)
	
	if os.path.exists(file_prof):
		print ('+ Profile file for scheme exists...')
	else:
		print ('+ Profile file for scheme will be downloaded...')
		
		if os.path.exists(profile_folder):
			print ('+ Profile folder already exists...')
		else:
			print ('+ Create profile folder...')
			profile_folder_path = functions.create_folder(profile_folder)
		
		print ('+ Downloading profile...')
		cmd_profile = "%s %s --species %s --scheme %s --dir_profile %s" %(rscript, MLSTarR_download_prf, species, scheme, profile_folder)
		callCode = functions.system_call(cmd_profile)
	
	## Check if seqs exist
	seq_bool = 0
	if os.path.exists(seq_folder):
		print ('+ Sequence folder exists...')
		files = os.listdir(seq_folder)
		count_fas = 0		
		for f in files:
			if f.endswith('.fas'):
				count_fas += 1
			else:				
				os.remove(seq_folder + '/' + f)

		if count_fas > 6:
			print ("+ Assuming sequences are previously downloaded...")
		else:
			seq_bool = 1		
			os.rmdir(seq_folder)
			seq_folder_path = functions.create_folder(seq_folder)
			
	else:
		print ('+ Create sequence folder...')
		seq_folder_path = functions.create_folder(seq_folder)
		seq_bool = 1
		

	if seq_bool == 1:
		print ('+ Sequence files for scheme will be downloaded...')
		print ('+ Downloading sequences...')
		cmd_seq = "%s %s --species %s --scheme %s --dir_seq %s" %(rscript, MLSTarR_download_seq, species, scheme, seq_folder)
		callCode = functions.system_call(cmd_seq)	


######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	fileGiven = os.path.abspath(argv[1])
	species = argv[2]
	scheme = int(argv[3])
	
	name = argv[4]
	rscript = argv[5]
	threads = argv[6]
	
	profile_folder = os.path.abspath(argv[7])
	seq_folder = os.path.abspath(argv[8])
	path = os.path.abspath(argv[9])
	
	## start
	download_PubMLST(profile_folder, scheme, seq_folder, name, rscript, species)
	
	## call MLSTar for this sample
	results = run_doMLST(profile_folder, seq_folder, name, rscript, path, fileGiven, threads)

	plot_MLST(results, profile_folder, rscript)

######
def help_options():
	print ("\nUSAGE: python %s genome species scheme name rscript_bin threads profile_folder seq_folder path\n"  %os.path.realpath(__file__))

		
######
'''******************************************'''
if __name__== "__main__":
	main()
