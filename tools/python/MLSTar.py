#usr/bin/env python
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

## import my modules
pythonDir = os.path.dirname(os.path.abspath(argv[0]))
sys.path.append(pythonDir)
import functions

configDir = os.path.dirname(os.path.abspath(argv[0])) + '/../../config/'
sys.path.append(configDir)
import config

RscriptDir = os.path.dirname(os.path.abspath(argv[0])) + '/../R'
MLSTarR_script = RscriptDir + '/MLSTar_call.R'
MLSTarR_plot = RscriptDir + '/MLSTar_plot.R'
MLSTarR_download_seq = RscriptDir + '/MLSTar_downloadPubMLST_seq.R'
MLSTarR_download_prf = RscriptDir + '/MLSTar_downloadPubMLST_profile.R'

######
def run_MLSTar(species, scheme, name, path, fileGiven, threads):

	profile_folder = config.MLSTar['profile_folder']
	seq_folder = config.MLSTar['sequence_folder']
	rscript = config.EXECUTABLES['Rscript']
	
	## check if profile and sequences are already downloaded
	download_PubMLST(profile_folder, scheme, seq_folder, name, rscript, species)

	## call MLSTar for this sample
	run_doMLST(profile_folder, scheme, seq_folder, name, rscript, species, path, fileGiven, threads)

######
def help_options():
	print ("\nUSAGE: python %s genome species scheme name rscript_bin threads profile_folder seq_folder path\n"  %os.path.abspath(argv[0]))

######
def run_doMLST(profile_folder, scheme, seq_folder, name, rscript, species, path, fileGiven, threads):
	print ('+ Create sample folder...')
	folder_path = functions.create_subfolder(name, path)

	print ('+ Generating profile for sample...')
	cmd_profiler = "%s %s --species %s --scheme %s --dir_profile %s --dir_seq %s --file %s --dir %s --name %s --threads %s" %(rscript, MLSTarR_script, species, scheme, profile_folder, seq_folder, fileGiven, folder_path, name, threads)
	callCode = functions.system_call(cmd_profiler)

######
def update_MLSTar_profile_alleles():
	return("")

######
def plot_MLST(results, profile):
	return("")
	
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
	run_doMLST(profile_folder, scheme, seq_folder, name, rscript, species, path, fileGiven, threads)

	
		
######
'''******************************************'''
if __name__== "__main__":
	main()
