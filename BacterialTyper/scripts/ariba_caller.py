#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Calls ARIBA_ software for several databases and custom dataset.
"""
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd
from ete3 import Tree

## import my modules
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

from BacterialTyper.config import set_config
from BacterialTyper.modules import citation

############################################################### 
def get_ARIBA_dbs(list_dbs):
	
	## General
	dbs = [
		"argannot", "card", "megares", 
		"plasmidfinder", "resfinder", 
		"srst2_argannot", "vfdb_core", 
		"vfdb_full", "virulencefinder"]
	
	conversion_dbs = {
		"ARG-ANNOT":"argannot", 
		"CARD":"card", 
		"MEGARes":"megares", 
		"PlasmidFinder":"plasmidfinder", 
		"ResFinder":"resfinder", 
		"srst2":"srst2_argannot", 
		"VFDB":"vfdb_full", 
		"VirulenceFinder":"virulencefinder"
	}
	
	db2returns = []
	for i in list_dbs:
		if (conversion_dbs[i] in dbs):
			db2returns.append(conversion_dbs[i])
	return (db2returns)
	
############################################################### 
def help_ARIBA():

	dict_ariba = citation.ariba_citation()
	## to do: finish filling information for different databases
	print ("")
	HCGB_aes.print_sepLine("*", 50, False)
	print ("CARD:")
	print ("The Comprehensive Antibiotic Resistance Database (CARD) is a rigorously curated collection of characterized, peer-reviewed  resistance determinants and associated antibiotics, organized by the Antibiotic Resistance Ontology (ARO) and AMR gene detection models.")
	print ('Citation:', dict_ariba['CARD'])
	print ("")	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("VFDB:")
	print ("The virulence factor database (VFDB) is an integrated and comprehensive online resource for curating information about virulence factors of bacterial pathogens. Since its inception in 2004, VFDB has been dedicated to providing up-to-date knowledge of VFs from various medically significant bacterial pathogens.")
	print ('Citation:', dict_ariba['VFDB'])
	print ("")	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("ARG-ANNOT:\n")
	print ("...")
	print ('Citation:', dict_ariba['ARG-ANNOT'])
	print ("")	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("MEGARes:")
	print ("The MEGARes database contains sequence data for approximately 4,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing.")
	print ('Citation:', dict_ariba['MEGARes'])
	print ("")	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("PlasmidFinder:\n")
	print ("...")
	print ('Citation:', dict_ariba['PlasmidFinder'])
	print ("")	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("ResFinder:\n")
	print ("The ResFinder database is a curated database of acquired resistance genes.")
	print ('Citation:', dict_ariba['ResFinder'])
	print ("")	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("srst2:")
	print ("...")
	print ('Citation:', dict_ariba['srst2'])
	print ("")
	HCGB_aes.print_sepLine("*", 50, False)
	print ("")
	
############################################################### 
def download_ariba_databases(list_dbs, main_folder, Debug, threads):

	"""Download ARIBA_ databases.
	
	Using ARIBA software this function retrieves desired databases and prepare them for later analysis.
	
	:param list_dbs: List of databases to download.
	:param main_folder: Absolute path to database folder.
	:param Debug: True/false for printing developer messages
	:param threads: Number of CPUs to use.
	
	:type list_dbs: string 
	:type main_folder: string
	:type Debug: Boolean
	:type threads: integer
	
	 .. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.file_functions.create_subfolder`
		
		- :func:`HCGB.functions.time_functions.read_time_stamp`
		
		- :func:`BacterialTyper.scripts.ariba_caller.get_ARIBA_dbs`
	
		- :func:`BacterialTyper.scripts.ariba_caller.ariba_getref`		
		
	 
	.. include:: ../../links.inc
	"""

	print("\n\n+ Download databases for Antimicrobial Resistance Identification By Assembly (ARIBA).")
	ariba_folder = HCGB_files.create_subfolder("ARIBA", main_folder)

	## information
	dict_info = {}

	## print ARIBA databases: 
	print ("+ Available databases:")
	dbs = get_ARIBA_dbs(list_dbs)
	
	for db_set in dbs:

		HCGB_aes.print_sepLine("-",30, False)
		print (colored("+ " + db_set,'yellow'))
		
		## prepare folders
		folder_set = HCGB_files.create_subfolder(db_set, ariba_folder)
		outdir_prepare_ref = folder_set + '_prepareref'

		## stamp time file
		filename_stamp_prepare = outdir_prepare_ref + '/.success'
		stamp_ariba_getref  = False
		
		## check if previously done
		if os.path.isfile(filename_stamp_prepare):
			stamp =	HCGB_time.read_time_stamp(filename_stamp_prepare)
			print ("\t+ Database is downloaded in folder: ", folder_set)
			print ("\t+ Data is available and indexed in folder: ", outdir_prepare_ref)
			print (colored("\tDatabase was previously downloaded and prepared on: %s" %stamp, 'yellow'))
		
			## Check if necessary to download again after several months/days
			days_passed = HCGB_time.get_diff_time(filename_stamp_prepare)
			print ("\t\t** %s days ago" %days_passed)		
			if (days_passed > 30): ## download again
				print ("\t\t** Downloading information again just to be sure...")

		(stamp_ariba_getref, fasta, metadata) = ariba_getref(db_set, folder_set, Debug, threads)
		
		#else:
		#	(stamp_ariba_getref, fasta, metadata) = ariba_getref(db_set, folder_set, Debug, threads)
		
		if not stamp_ariba_getref:
			print (colored("** ARIBA getref failed or generated a warning for " + db_set, 'red'))
			stamp_ariba_getref='ERROR'
			fasta  = ""
			metadata = ""
			
		## fill information
		dict_info[db_set] = {'folder': outdir_prepare_ref,
							'stamp': stamp_ariba_getref,
							'fasta':fasta,
							'metadata': metadata }
			
	return (dict_info)

##########
def check_db_indexed(folder, option):
	"""Check if ARIBA_ database is indexed.
	
	In the given folder it looks for '00.info.txt' file.
	
	:param folder: Absolute path to database folder.
	:param option: Whether to print more information messages or not [Yes/No]. 
	
	:type folder: string 
	:type option: string
	:returns: Boolean True/False.
	 
	 
	 .. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.time_functions.read_time_stamp`
	 
	.. include:: ../../links.inc
	"""
	
	path_basename = folder.split('/')
	if os.path.isfile(os.path.join(folder,'00.info.txt')):
		time_stamp = os.path.join(folder, '.success')
		if os.path.isfile(time_stamp):
			stamp =	HCGB_main.get_info_file(time_stamp)
			print (colored("\tA previous command generated results on: %s [%s]" %(stamp, path_basename[-2]), 'yellow'))
			return (True, stamp[0])
	else:
		if (option == 'YES'):
			print (colored("\t- ARIBA database: " + path_basename + " [ ERROR ]", 'red'))
		return (False, "")
	
##########
def ariba_summary(out, infileList, options):
	"""Create ARIBA summary
	
	This function calls ARIBA_ summary and generates a summary of multiple ARIBA report files. It also creates Phandango_ files.
	
	:param out: Prefix of output files
	:param infileList: Files to be summarised
	:param options: Additional options for ariba summary. See details below. Provide them within quotes.
	
	:type out: string
	:type infileList: list
	:type options: string
	
	:returns: functions.system_call message
	
	.. note:: Additional options for ariba_summary (as specified by ARIBA_)
	
		--cluster_cols col_ids			Comma separated list of cluster columns to include. Choose from: assembled, match, ref_seq, pct_id, ctg_cov, known_var, novel_var [match]

		--no_tree             Do not make phandango tree

		--min_id FLOAT        Minimum percent identity cutoff to count as assembled [90]

		--only_clusters Cluster_names			Only report data for the given comma-separated list of cluster names, eg: cluster1,cluster2,cluster42
		
		--v_groups      Show a group column for each group of variants
		
		--known_variants      Report all known variants
		
		--novel_variants      Report all novel variants
	
		--col_filter string			Choose whether columns where all values are "no" or "NA" are removed [yes/no]

		--row_filter string			Choose whether rows where all values are "no" or "NA" are removed [yes/no]


	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.system_call_functions.system_call`
	
	"""
	
	logFile = out + '.log'
	infileList_string = " ".join(infileList)
	cmd_summary = 'ariba summary %s %s %s 2> %s' %(options, out, infileList_string, logFile)
	
	return(HCGB_sys.system_call(cmd_summary))

##########
def ariba_getref(database, outdir, Debug, threads):
	######################################################################################
	## usage: ariba getref [options] <db> <outprefix>
	######################################################################################
	## Download reference data from one of a few supported public resources
	## positional arguments:
	##	DB name            Database to download. Must be one of: argannot card megares plasmidfinder resfinder srst2_argannot vfdb_core vfdb_full virulencefinder
	##  outprefix          Prefix of output filenames
	######################################################################################

	## where database is one of: 
	##	argannot, card, megares, plasmidfinder, resfinder,
	##	srst2_argannot, vfdb_core, vfdb_full, virulencefinder.

	## folders
	outdir_name = outdir + '/' + database
	outdir_prepare_ref = outdir + '_prepareref'

	## download information in database folder provided by config
	print ("\t+ Retrieve information from database: " + database)

	## check if previously downloaded and succeeded
	filename_stamp = outdir + '/.success'

	if os.path.isfile(filename_stamp):
		stamp =	HCGB_time.read_time_stamp(filename_stamp)
		download_ariba_cmd = 'OK'

		## debug message
		if (Debug):
			print (colored("**DEBUG: ariba getref %s succeed " %database + "**", 'yellow'))

	else:
		cmd_getref = 'ariba getref %s %s' %(database, outdir_name)
		download_ariba_cmd = HCGB_sys.system_call(cmd_getref)
	
	if not download_ariba_cmd:
		## rise error & exit
		print (colored("***ERROR: ariba getref %s failed " %database + " **",'red'))
		return(False, "", "")

	## debug message
	if (Debug):
		print (colored("**DEBUG: Run ariba prepareref %s " %database + "**", 'yellow'))

	## check if previously prepareref and succeeded
	## get information
	list_files = os.listdir(outdir)
	fasta = ""
	metadata = ""
	for f in list_files:
		if f.endswith('tsv'):
			metadata = outdir + '/' + f
		elif f.endswith('fa'):
			fasta = outdir + '/' + f

	filename_stamp_prepare = outdir_prepare_ref + '/.success'
	if os.path.isfile(filename_stamp_prepare):
		stamp =	HCGB_time.read_time_stamp(filename_stamp_prepare)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, database), 'yellow'))
		return(stamp, fasta, metadata)
	else:
		code = ariba_prepareref(fasta, metadata, outdir_prepare_ref, threads)
		
		if (code == 'OK'):
			filename_stamp = outdir_prepare_ref + '/.success'
			HCGB_time.print_time_stamp(filename_stamp_prepare)
			return(time.time(), fasta, metadata)
		else:
			return(False, fasta, metadata)
				
##########
def ariba_prepareref(fasta, metadata, outfolder, threads):
	## prepareref
	cmd_prepareref = 'ariba prepareref -f %s -m %s %s --threads %s' %(fasta, metadata, outfolder, threads)
	return(HCGB_sys.system_call(cmd_prepareref))
	
	######################################################################################
	## usage: ariba prepareref [options] <outdir>
	######################################################################################
	## Prepare reference data for running the pipeline with "ariba run"
	## positional arguments:
	##   outdir                Output directory (must not already exist)
	## optional arguments:
	## -h, --help            show this help message and exit
	## input files options:
	##   -f FILENAME, --fasta FILENAME
	## 		REQUIRED. Name of fasta file. Can be used more than once if your sequences are spread over more than on file
	##   -m FILENAME, --metadata FILENAME
	##      Name of tsv file of metadata about the input sequences. Can be used more than once if your metadata is spread over more than one file. Incompatible with --all_coding
	##   --all_coding {yes,no}
	##      Use this if you only have a fasta of presence absence sequences as input, and no metadata. Use "yes" if all sequences are coding, or "no" if they are all non- coding. Incompatible with -m/--metadata
	## 
	## cd-hit options:
	##   --no_cdhit Do not run cd-hit. Each input sequence is put into its own "cluster". Incompatible with --cdhit_clusters.
	##   --cdhit_clusters FILENAME --> File specifying how the sequences should be clustered. Will be used instead of running cdhit. Format is one cluster per line. Sequence names separated by whitespace. Incompatible with --no_cdhit
	##   --cdhit_min_id FLOAT  Sequence identity threshold (cd-hit option -c) [0.9]
  	## 	 --cdhit_min_length FLOAT Length difference cutoff (cd-hit option -s) [0.0]
	##   --cdhit_max_memory INT Memory limit in MB (cd-hit option -M) [None]. Use 0 for unlimited.
	##
	## other options:
	##   --min_gene_length INT -->  Minimum allowed length in nucleotides of reference genes [6]
	##   --max_gene_length INT -->  Maximum allowed length in nucleotides of reference genes [10000]
	##   --genetic_code INT   -->   Number of genetic code to use. Currently supported 1,4,11 [11]
	##   --force               Overwrite output directory, if it already exists
	##   --threads INT         Number of threads (currently only applies to cdhit)[1]
	##  --verbose             Be verbose
	######################################################################################

##########
def ariba_expandflag(input_file, output_file):
	######################################################################################
	## usage: ariba expandflag infile.tsv outfile.tsv
	######################################################################################	
	#	Expands the flag column in a report file from number to comma-separated list
	#	of flag bits
	#
	#	positional arguments:
	#	  infile      Name of input report TSV file
	#	  outfile     Name of output report TSV file
	######################################################################################
	
	## download information in database folder provided by config
	#print ("+ Call ariba module 'expandflag' to add additional information to each entry.")
	
	## Somehow when I do expandflag of report.tsv file I get spaces convert to tabs and everything is 
	## distorted. I had to make this 
	
	data = pd.read_csv(input_file, header=0, sep='\t')
	list_data2 = set(list(data['flag']))
	tmp_name = os.path.splitext(input_file)[0]
	
	## print only flag
	flag_file = tmp_name + '-tmp.tsv'
	flag_file_hd = open(flag_file, "w")	
	flag_file_hd.write('#flag_name\tflag')
	flag_file_hd.write('\n')
	
	for line in zip(list_data2, list_data2):
		string = ('flag_{}\t{}'.format(*line))
		flag_file_hd.write(string)
		flag_file_hd.write('\n')
	flag_file_hd.close()
	
	## generate description
	flag_file_out = tmp_name + '-description.tsv'
	cmd = 'ariba expandflag "%s" %s' %(flag_file, flag_file_out)
	HCGB_sys.system_call(cmd)
	
	os.remove(flag_file)	
	return(flag_file_out)

##########
def ariba_pubmlstget(species, outdir):
	######################################################################################
	## usage: ariba pubmlstget [options] <"species in quotes"> <output_directory>
	######################################################################################	
	## Download typing scheme for a given species from PubMLST, and make an ARIBA db
	## positional arguments:
	##  species     Species to download. Put it in quotes
	##	outdir      Name of output directory to be made (must not already exist)
	######################################################################################
	
	## download information in database folder provided by config
	print ("+ Call ariba module 'pubmlstget' to retrieve MLST information.")
	HCGB_files.create_folder(outdir)
	cmd = 'ariba pubmlstget "%s" %s' %(species, outdir)
	return(HCGB_sys.system_call(cmd))
	
##########
def ariba_run(database, files, outdir, threads, threshold):
	"""ARIBA search call
	
	Given a database, it generates an ARIBA search of the reads provided. Results would be later processed using :func:`BacterialTyper.scripts.virulence_resistance.check_results` and :func:`BacterialTyper.scripts.virulence_resistance.results_parser`
	
	:param database: Folder containing ARIBA database previously indexed.
	:param files: List of files with absolute path to fastq reads (R1 & R2).
	:param outdir: Absolute path to folder, must not exist.
	:param threads: Number of CPUs to use.
	:param threshold: If proportion of gene assembled (regardless of into how many contigs) is at least this value then the flag gene_assembled is set [Default: 0.90].
	
	:type database: string
	:type files: list
	:type outdir: string
	:type threads: integer
	:type threshold: float [0-1]
	
	:return: Message: OK/FAIL.
	:rtype: string
	
	:warnings: No single mode available, only paired-end.
	
	:warnings: Output directory must not exist.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.system_call_functions.system_call`
		
		- :func:`HCGB.functions.time_functions.print_time_stamp`
		
		- :func:`BacterialTyper.scripts.virulence_resistance.check_results` 
		
		- :func:`BacterialTyper.scripts.virulence_resistance.results_parser`
	
	
	"""
	
	######################################################################################
	##	usage: ariba run [options] <prepareref_dir> <reads1.fq> <reads2.fq> <outdir>
	##	Runs the local assembly pipeline. Input is dir made by prepareref, and paired reads
	##
	## positional arguments:
	##  prepareref_dir        Name of output directory when "ariba prepareref" was run
	##  reads_1               Name of fwd reads fastq file
	##  reads_2               Name of rev reads fastq file
	##  outdir                Output directory (must not already exist)
	## optional arguments:
	##  -h, --help            show this help message and exit
	## nucmer options:
	##  --nucmer_min_id INT   Minimum alignment identity (delta-filter -i) [90]
	##  --nucmer_min_len INT  Minimum alignment length (delta-filter -i) [20]
	##  --nucmer_breaklen INT Value to use for -breaklen when running nucmer [200]
	## Assembly options:
	##  --assembler {fermilite,spades} Assembler to use
	##  --assembly_cov INT    Target read coverage when sampling reads for assembly [50]
	##  --min_scaff_depth INT Minimum number of read pairs needed as evidence for scaffold link between two contigs [10]
	##  --spades_mode {wgs,sc,rna}  If using Spades assembler, either use default WGS mode, Single Cell mode (`spades.py --sc`) or RNA mode (`spades.py --rna`). Use SC or RNA mode if your input is from a viral sequencing with very uneven and deep                        coverage. Set `--assembly_cov` to some high value if using SC or RNA mode
	##  --spades_options SPADES_OPTIONS: Extra options to pass to Spades assembler. Sensible default options will be picked based on `--spades_mode` argument.  Anything set here will replace the defaults completely
	## Other options:
	##  --threads INT         Experimental. Number of threads. Will run clusters in parallel, but not minimap (yet) [1]
	##  --assembled_threshold FLOAT (between 0 and 1)  If proportion of gene assembled (regardless of into how many contigs) is at least this value then the flag gene_assembled is set [0.95]
	##  --gene_nt_extend INT  Max number of nucleotides to extend ends of gene matches to look for start/stop codons [30]
	##  --unique_threshold FLOAT (between 0 and 1) If proportion of bases in gene assembled more than once is <= this value, then the flag unique_contig is set [0.03]
	##  --force               Overwrite output directory, if it already exists
	##  --noclean             Do not clean up intermediate files
	##  --tmp_dir TMP_DIR     Existing directory in which to create a temporary directory used for local assemblies
	##  --verbose             Be verbose
	######################################################################################

	# outdir must not exist.
	logFile = outdir + '_run.log'
	
	######################################################
	## default options: might be interesting to change?
	##  --assembled_threshold 0.95
	##  --nucmer_min_id INT	90
	##  --nucmer_min_len INT 20
	######################################################
	
	if (len(files) == 2):
		cmd = 'ariba run --assembled_threshold %s %s %s %s %s --threads %s 2> %s' %(threshold, database, files[0], files[1], outdir, threads, logFile)
	else:
	
		## [TODO]
		print (colored("\n\n***** TODO: Check what to do with single end when calling ARIBA *****\n\n", 'red'))
		print (colored("\n\n***** No implementation yet *****\n\n", 'red'))
		exit()

	##
	code = HCGB_sys.system_call(cmd)

	## make stamp time
	if (code == 'OK'): 
		filename_stamp = outdir + '/.success'
		return ('OK')
	else:
		return ('FAIL')
	
#############################################################
def ariba_summary_all(outfile, dict_files):
	"""Create final ARIBA summary 
	
	Generates a temporal summary of all the reports generated by ARIBA (using :func:`BacterialTyper.scripts.ariba_caller.ariba_summary`). The result is a csv file and newick tree for later Phandango_ visualization.
	
	These files contained for each entry a ".match" string appended. Using :func:`BacterialTyper.scripts.ariba_caller.fix_ariba_summary` and :func:`BacterialTyper.scripts.ariba_caller.fix_phandando_tree`, we would removed this ".match" and prepare files 
	
	Final result files are a csv file (containing code colors) and newick tree for later Phandango_ visualization.
	
	:param outfile:
	:param dict_files: 
	
	:type outfile:
	:type dict_files: 
	
	:return: Output file
	:rtype: string
		
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.ariba_caller.ariba_summary`
		
		- :func:`BacterialTyper.scripts.ariba_caller.fix_ariba_summary` 
		
		- :func:`BacterialTyper.scripts.ariba_caller.fix_phandango_tree`
	
	"""
	
	## Create tmp list of files 
	file_list = []
	for files in dict_files:
		file_list.append(dict_files[files])
	
	## call ariba_summary
	outfile_tmp = outfile + '_tmp'
	options = ''
	ariba_summary(outfile_tmp, file_list, options)

	## fix output
	fix_ariba_summary(outfile_tmp + '.csv', outfile + '.csv', dict_files)
	fix_ariba_summary(outfile_tmp + '.phandango.csv', outfile + '.phandango.csv', dict_files)
	fix_phangando_tree(outfile_tmp + '.phandango.tre', outfile + '.phandango.tree', dict_files)
	
	if os.path.isfile(outfile_tmp + '.csv'):
		os.remove(outfile_tmp + '.csv')
	else:
		return('NaN')

	if os.path.isfile(outfile_tmp + '.phandango.csv'):
		os.remove(outfile_tmp + '.phandango.csv')

	if os.path.isfile(outfile_tmp + '.phandango.tree'):
		os.remove(outfile_tmp + '.phandango.tree')


	return(outfile + '.csv')
	
#############################################################
def fix_ariba_summary(csv_file, outfile, dict_files):

	if not os.path.isfile(csv_file):
		## summary file not generated
		print (colored("*** ERROR: summary file (%s) not generated ***" %csv_file, 'red'))
		return()

	summary_data = pd.read_csv(csv_file, header=0, sep=',')
	## parse results
	cluster = ['name']
	for column in summary_data:
		if (column == 'name'):
			continue
		cluster_name = column.replace(".match", "")
		cluster.append(cluster_name)

	## format summary data
	summary_data.columns = cluster 
	for index, row in summary_data.iterrows():
		for key, value in dict_files.items():
			if (value == row['name'] ):
				summary_data.loc[index]['name'] = key
	
	summary_data = summary_data.set_index('name')
	summary_data.to_csv(outfile)
	return()
	
#############################################################
def fix_phangando_tree(tree_file, outfile_name, dict_files): 
	
	if not os.path.isfile(tree_file):
		## summary file not generated
		print (colored("*** ERROR: phandango tree file (%s) not generated ***" %tree_file, 'red'))
		return()

	t = Tree(tree_file)
	for leaf in t:
		for key, value in dict_files.items():
			if (value == leaf.name):
				leaf.name = key

	# We can also write into a file
	t.write(format=1, outfile=outfile_name)
	return()


######
def	help_options():
	print ("\nUSAGE: python %s file name xx threads path\n"  %os.path.realpath(__file__))
	print ("*** if paired-end data provide a csv string for file argument ***\n\n")

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	help_ARIBA()
		
######
if __name__== "__main__":
	main()

