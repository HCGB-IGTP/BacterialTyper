#usr/bin/en python3
'''
This code...
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab"," IGTP"," Spain
'''
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

###############
def busco_datasets():
	busco_data = pd.DataFrame(columns=("Taxonomic range","Dataset","ftp_site"))
	busco_data.loc[len(busco_data)] = ("All bacteria dataset","bacteria","http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Proteobacteria all","proteobacteria","http://busco.ezlab.org/v2/datasets/proteobacteria_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Proteobacteria – Alphaproteobacteria – Order Rhizobiales","rhizobiales","http://busco.ezlab.org/v2/datasets/rhizobiales_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Proteobacteria – Betaproteobacteria","betaproteobacteria","http://busco.ezlab.org/v2/datasets/betaproteobacteria_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Proteobacteria – Gammaproteobacteria all","gammaproteobacteria","http://busco.ezlab.org/v2/datasets/gammaproteobacteria_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Proteobacteria – Gammaproteobacteria – Order Enterobacteriales","enterobacteriales","http://busco.ezlab.org/v2/datasets/enterobacteriales_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Proteobacteria – Delta and Epsilon proteobacteria","deltaepsilonsub","http://busco.ezlab.org/v2/datasets/deltaepsilonsub_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Actinobacteria","actinobacteria","http://busco.ezlab.org/v2/datasets/actinobacteria_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Cyanobacteria","cyanobacteria","http://busco.ezlab.org/v2/datasets/cyanobacteria_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Firmcutes all","firmicutes","http://busco.ezlab.org/v2/datasets/firmicutes_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Firmcutes – Clostridia","clostridia","http://busco.ezlab.org/v2/datasets/clostridia_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Firmcutes – Bacilli – order Lactobacillales","lactobacillales","http://busco.ezlab.org/v2/datasets/lactobacillales_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Firmcutes – Bacilli – order Bacillales","bacillales","http://busco.ezlab.org/v2/datasets/bacillales_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Bacteroidetes","bacteroidetes","http://busco.ezlab.org/v2/datasets/bacteroidetes_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Spirochaetes","spirochaetes","http://busco.ezlab.org/v2/datasets/spirochaetes_odb9.tar.gz")
	busco_data.loc[len(busco_data)] = ("phylum Tenericutes","tenericutes","http://busco.ezlab.org/v2/datasets/tenericutes_odb9.tar.gz")
	busco_data = busco_data.set_index('Dataset')
	return(busco_data)

###############
def print_help_BUSCO():
	functions.print_sepLine("*", 50, 'yellow')
	print ('BUSCO Help')
	functions.print_sepLine("*", 50, 'yellow')
	print ("Benchmarking of Universal Single Copy Orhtologs (BUSCO)\n")	
	
	print ("BUSCO v3 provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness.")
	print ("It is based on evolutionarily-informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB v9.\n")
	
	print ("BUSCO assessments are implemented in open-source software, with a large selection of lineage-specific sets of Single-Copy Orthologs.")
	print ("These conserved orthologs are ideal candidates for large-scale phylogenomics studies, and the annotated BUSCO gene models built during ")
	print ("genome assessments provide a comprehensive gene predictor training set for use as part of genome annotation pipelines.\n")
	
	print ("For more information visit the website: https://busco.ezlab.org/\n")
	functions.print_sepLine("*", 50, False)

	print ("\nPlease select among the available datasets according to your sample expected taxonomic range.") 
	print ("Bear in mind that several datasets could be provided as some represent a broad taxonomic range: ")
	print ("E.g. Bacteria > Proteobacteria > Gammaproteobacteria > Enterobacteriales")
	print ("E.g. Bacteria > Firmicutes > Bacillales")
	print ("Datasets:\n")	
	print_available_BUSCO()

###############
def print_available_BUSCO():
	print_df = busco_datasets()
	print_df = print_df[['Taxonomic range']]
	pd.set_option('display.max_colwidth', -1)
	functions.print_sepLine("-", 100, False)
	print (print_df)
	functions.print_sepLine("-", 100, False)
	print ("\n")

###############
def BUSCO_download(name, ftp, folder):
	
	print (colored("\n+ BUSCO dataset: " + name + " - v9 OrthoDB", 'yellow'))
	subfolder = folder + '/' + name
	file_name = os.path.basename(ftp)
	path_file = subfolder + '/' + file_name
	folderName = subfolder + '/' + file_name.split('.tar.gz')[0]

	if os.path.exists(subfolder):
		print ('Already available in the path provided...') 
	else:
		## does not exists
		functions.create_folder(subfolder)

		functions.wget_download(ftp, subfolder)

		## extract
		print ("+ Extract file...")
		functions.extract(path_file, subfolder)

		## timestamp
		filename_stamp = subfolder + '/.success'
		functions.print_time_stamp(filename_stamp)

	## check if it is allright
	code = BUSCO_check_dataset(folderName)
	functions.print_sepLine("-", 50, False)
	
	if (code == 'FAIL'):
		print (colored('*** Dataset failed. Try to download it again...','red'))
		shutil.rmtree(folder)
		BUSCO_download(name, ftp, folder)
		
	return (folderName)

###############
def BUSCO_check_dataset(folder):
	config_file = folder + '/dataset.cfg'
	##	name=bacteria_odb9
	##	species=E_coli_K12
	##	domain=prokaryota
	##	creation_date=2016-11-01
	##	number_of_BUSCOs=148
	##	number_of_species=3663

	#print ("+ Checking the integrity of BUSCO dataset in folder: ", folder)
	functions.print_sepLine("+", 10, False)
	print ("Statistics")
	functions.print_sepLine("+", 10, False)
	if os.path.isfile(config_file):
		list_config = functions.readList_fromFile(config_file)
		for elem in list_config:
			line = elem.split("=")
			line[0] = line[0].replace("_", " ")
			print ("\t".join(line))
		
		if os.path.isfile(folder + '/../.success'):
			timestamp = functions.read_time_stamp(folder + '/../.success')
			print ("download date: ", timestamp)
		
		print ("Available in folder: ", folder)
		print (colored("Dataset....[ OK ]\n", 'green'))	
	else:
		print (colored("Dataset....[ FAIL ]\n", 'red'))	
		return ('FAIL')

###############
def BUSCO_retrieve_sets(list_datasets, folder):

	## check if name matches
	data_df = busco_datasets()
	for elem in list_datasets:
		if not elem in data_df.index:
			print (colored("***ATTENTION: Please check the datasets provided: " + elem, 'red'))
			print ("\n\nAvailable datasets:")
			print_available_BUSCO()
			exit()
	
	## database: folder
	dataset = {}
	for elem in list_datasets:
		ftp = data_df.loc[elem]['ftp_site']
		dataset[elem] = BUSCO_download(elem, ftp, folder)
	
	return (dataset)

###############
def BUSCO_run(dataset, fasta, threads, output_name, dataset_name, mode):

	my_out_folder = output_name + '/run_' + dataset_name
	## timestamp
	filename_stamp =  my_out_folder + '/.success'

	print (colored("\tBUSCO Dataset [%s]; Sample [%s]" %(dataset_name, fasta), 'yellow'))
		
	## check previous run
	if os.path.isfile(filename_stamp):
		timestamp = functions.read_time_stamp(filename_stamp)
		print (colored("\tSuccessfully run on date: %s"  %timestamp, 'green'))
	else:
	
		busco_bin = config.get_exe('busco')
		os.chdir(output_name)
		logFile = dataset_name + '.log'
		cmd = '%s -i %s -f -c %s --blast_single_core --mode %s -l %s -o %s > %s' %(busco_bin, fasta, threads, mode, dataset, dataset_name, logFile)
		functions.system_call(cmd)
	
		if os.path.isfile(my_out_folder + '/short_summary_' + dataset_name + '.txt'):
			## timestamp
			functions.print_time_stamp(filename_stamp)
		else:
			print (colored("BUSCO failed: Dataset [%s]; Sample [%s]" %(dataset_name, fasta), 'red'))
			return ('FAIL')

	return()

###############
def BUSCO_plot(outfolder):
	busco_plot_bin = config.get_exe('busco_plot')
	#logFile = dataset_name + '.log'
	cmd = '%s -wd %s' %(busco_plot_bin, outfolder)
	functions.system_call(cmd)


###############
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		#help_options()
		exit()
	
	BUSCO_check_dataset(sys.argv[1])
	
######
'''******************************************'''
if __name__== "__main__":
	main()
		

