#usr/bin/env python
'''
This code calls ARIBA software 
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
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

############################################################### 
def get_ARIBA_dbs():
	dbs = ["argannot", "card", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"]
	return (dbs)

############################################################### 
def download_ariba_databases(main_folder, Debug):

	## ToDo check if already download	
	print("\n\n+ Download databases for Antimicrobial Resistance Identification By Assembly (ARIBA).")
	ariba_folder = functions.create_subfolder("ARIBA", main_folder)
	## where database is one of: 
	print ("+ Available databases:")
	out_info = main_folder + '/ARIBA_information.txt'
	hd = open(out_info, 'w')

	dbs = get_ARIBA_dbs()
	for db_set in dbs:
		functions.print_sepLine("-",30)
		print (colored("+ " + db_set,'yellow'))
		
		folder_set = functions.create_subfolder(db_set, ariba_folder)
		return_ariba_getref = ariba_getref(db_set, folder_set, Debug)
		
		if (return_ariba_getref == 'OK'):
			## print citation for each database
			functions.print_sepLine("'", 75)
			hd.write(db_set)
			hd.write('\n')
		else:
			print (colored("** ARIBA getref failed for " + db_set, 'red'))
	
	hd.close()


##########
def ariba_summary():
	######################################################################################
	## usage: ariba summary [options] <outprefix> [report1.tsv report2.tsv ...]
	######################################################################################
	## Make a summary of multiple ARIBA report files, and also make Phandango files
	## positional arguments:
	##  outprefix             Prefix of output files
	##  infiles               Files to be summarised
	##
	## optional arguments:
	##  -h, --help            show this help message and exit
	##  -f FILENAME, --fofn FILENAME File of filenames of ariba reports to be summarised.
	##                  Must be used if no input files listed after the 
    ##                  outfile. The first column should be the filename. An
    ##                  optional second column can be used to specify a sample
    ##                  name for that file, which will be used instead of the
    ##                  filename in output files. Columns separated by whitespace.
	##  --preset minimal|cluster_small|cluster_all|cluster_var_groups|all|all_no_filter
	##                  Shorthand for setting --cluster_cols,--col_filter,--
    ##                  row_filter,--v_groups,--variants. Using this overrides those options
	##  --cluster_cols col1,col2,...
    ##                  Comma separated list of cluster columns to include.
    ##                  Choose from: assembled, match, ref_seq, pct_id,
    ##                  ctg_cov, known_var, novel_var [match]
	##  --col_filter y|n      Choose whether columns where all values are "no" or "NA" are removed [y]
	##  --no_tree             Do not make phandango tree
	##  --row_filter y|n      Choose whether rows where all values are "no" or "NA" are removed [y]
	##  --min_id FLOAT        Minimum percent identity cutoff to count as assembled [90]
	##  --only_clusters Cluster_names
    ##                  Only report data for the given comma-separated list of
    ##                  cluster names, eg: cluster1,cluster2,cluster42
	##  --v_groups      Show a group column for each group of variants
	##  --known_variants      Report all known variants
	##  --novel_variants      Report all novel variants
	##  --verbose             Be verbose
	##Files must be listed after the output file and/or the option --fofn must be
	##used. If both used, all files in the filename specified by --fofn AND the
	##files listed after the output file will be used as input.
	######################################################################################

	print ()
	
##########
def ariba_getref(database, outdir, Debug):
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
	cmd_getref = 'ariba getref %s %s' %(database, outdir_name)
	download_ariba_cmd = functions.system_call(cmd_getref)
	
	if (download_ariba_cmd == 'OK'):
		## debug message
		if (Debug):
			print (colored("**DEBUG: ariba getref %s succeed " %database + "**", 'yellow'))
	
	else: 
		## rise error & exit
		print (colored("***ERROR: ariba getref %s failed " %database + " **",'red'))
		return('FAIL')

	## debug message
	if (Debug):
		print (colored("**DEBUG: Run ariba prepareref %s " %database + "**", 'yellow'))

	## get information
	list_files = os.listdir(outdir)
	fasta = ""
	metadata = ""
	for f in list_files:
		if f.endswith('tsv'):
			metadata = outdir + '/' + f
		elif f.endswith('fa'):
			fasta = outdir + '/' + f
	
	return(ariba_prepareref(fasta, metadata, outdir_prepare_ref))	
	
##########
def ariba_prepareref(fasta, metadata, outfolder):
	## prepareref
	cmd_prepareref = 'ariba prepareref -f %s -m %s %s' %(fasta, metadata, outfolder)
	return(functions.system_call(cmd_prepareref))
	
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
	functions.create_folder(outdir)
	cmd = 'ariba pubmlstget "%s" %s' %(species, outdir)
	return(functions.system_call(cmd))
	
##########
def ariba_run():
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

	print ()
	
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

		
######
if __name__== "__main__":
	main()

