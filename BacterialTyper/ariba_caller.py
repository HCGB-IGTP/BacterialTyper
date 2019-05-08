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
import pandas as pd

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper.modules import citation

## General
dbs = ["argannot", "card", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"]
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

############################################################### 
def get_ARIBA_dbs(list_dbs):
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
	functions.print_sepLine("*", 50)
	print ("CARD: https://card.mcmaster.ca/")
	print ("The Comprehensive Antibiotic Resistance Database (CARD) is a rigorously curated collection of characterized, peer-reviewed  resistance determinants and associated antibiotics, organized by the Antibiotic Resistance Ontology (ARO) and AMR gene detection models.")
	print ('Citation:', dict_ariba['CARD'])
	print ("")	
	functions.print_sepLine("*", 50)
	print ("VFDB: http://www.mgc.ac.cn/VFs/main.htm")
	print ("The virulence factor database (VFDB) is an integrated and comprehensive online resource for curating information about virulence factors of bacterial pathogens. Since its inception in 2004, VFDB has been dedicated to providing up-to-date knowledge of VFs from various medically significant bacterial pathogens.")
	print ('Citation:', dict_ariba['VFDB'])
	print ("")	
	functions.print_sepLine("*", 50)
	print ("ARG-ANNOT:\n")
	print ("...")
	print ('Citation:', dict_ariba['ARG-ANNOT'])
	print ("")	
	functions.print_sepLine("*", 50)
	print ("MEGARes: http://megares.meglab.org/")
	print ("The MEGARes database contains sequence data for approximately 4,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing.")
	print ('Citation:', dict_ariba['MEGARes'])
	print ("")	
	functions.print_sepLine("*", 50)
	print ("PlasmidFinder:\n")
	print ("...")
	print ('Citation:', dict_ariba['PlasmidFinder'])
	print ("")	
	functions.print_sepLine("*", 50)
	print ("ResFinder:\n")
	print ("The ResFinder database is a curated database of acquired resistance genes.")
	print ('Citation:', dict_ariba['ResFinder'])
	print ("")	
	functions.print_sepLine("*", 50)
	print ("srst2: https://github.com/katholt/srst2")
	print ("...")
	print ('Citation:', dict_ariba['srst2'])
	print ("")
	functions.print_sepLine("*", 50)
	print ("")
	
############################################################### 
def download_ariba_databases(list_dbs, main_folder, Debug):

	## ToDo check if already download	
	print("\n\n+ Download databases for Antimicrobial Resistance Identification By Assembly (ARIBA).")
	ariba_folder = functions.create_subfolder("ARIBA", main_folder)
	## where database is one of: 
	print ("+ Available databases:")
	out_info = main_folder + '/ARIBA_information.txt'
	hd = open(out_info, 'a')

	dbs = get_ARIBA_dbs(list_dbs)
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
			return ('FAIL')
	
	hd.close()
	return('OK')

##########
def check_db_indexed(folder, option):
	# get databases downloaded
	lineList = [line.rstrip('\n') for line in open(folder + '../../' + 'ARIBA_information.txt')]
	if os.path.isfile(folder + '00.info.txt'):
		name = folder.split("_prepareref/")[0]
		path_basename = os.path.basename(name)
		if path_basename in lineList: ## double check
			if (option == 'YES'):
				print (colored("\t- ARIBA database: " + path_basename + " [ OK ]", 'green'))
			return True
		else:
			if (option == 'YES'):
				print (colored("\t- ARIBA database: " + path_basename + " [ ERROR ]", 'red'))
			return False
	else:
		if (option == 'YES'):
			print (colored("\t- ARIBA database: " + path_basename + " [ ERROR ]", 'red'))
		return False
	
##########
def ariba_summary(out, infileList):
	######################################################################################
	## usage: ariba summary [options] <outprefix> [report1.tsv report2.tsv ...]
	######################################################################################
	## Make a summary of multiple ARIBA report files, and also make Phandango files
	## positional arguments:
	##  outprefix             Prefix of output files
	##  infiles               Files to be summarised
	##
	## optional arguments:
	##  -h, --help      show this help message and exit
	##  -f FILENAME		File of filenames of ariba reports to be summarised. Must be used if no input files listed after the outfile. The first column should be the filename. An optional second column can be used to specify a sample name for that file, which will be used instead of the filename in output files. Columns separated by whitespace.
	##  --preset minimal|cluster_small|cluster_all|cluster_var_groups|all|all_no_filter
	##                  Shorthand for setting --cluster_cols,--col_filter,--row_filter,--v_groups,--variants. Using this overrides those options
	##  --cluster_cols col1,col2,...  Comma separated list of cluster columns to include. Choose from: assembled, match, ref_seq, pct_id, ctg_cov, known_var, novel_var [match]
	##  --col_filter y|n      Choose whether columns where all values are "no" or "NA" are removed [y]
	##  --no_tree             Do not make phandango tree
	##  --row_filter y|n      Choose whether rows where all values are "no" or "NA" are removed [y]
	##  --min_id FLOAT        Minimum percent identity cutoff to count as assembled [90]
	##  --only_clusters Cluster_names. Only report data for the given comma-separated list of cluster names, eg: cluster1,cluster2,cluster42
	##  --v_groups      Show a group column for each group of variants
	##  --known_variants      Report all known variants
	##  --novel_variants      Report all novel variants
	##  --verbose             Be verbose
	######################################################################################
	
	logFile = out + '.log'
	infileList_string = " ".join(infileList)
	cmd_summary = 'ariba summary %s %s 2> %s' %(out, infileList_string, logFile)
	return(functions.system_call(cmd_summary))

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
	functions.system_call(cmd)
	
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
	functions.create_folder(outdir)
	cmd = 'ariba pubmlstget "%s" %s' %(species, outdir)
	return(functions.system_call(cmd))
	
##########
def ariba_run(database, files, outdir, threads):
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

	if (len(files) > 1):
		cmd = 'ariba run %s %s %s %s --threads %s' %(database, files[0], files[1], outdir, threads)
	else:
		print ("") ## todo
	##
	#return(functions.system_call(cmd))
	return()
	
	#########################################################
	# 	Column	Description
	# 	1. ariba_ref_name	ariba name of reference sequence chosen from cluster (needs to rename to stop some tools breaking)
	# 	2. ref_name	original name of reference sequence chosen from cluster, before renaming
	# 	3. gene	1=gene, 0=non-coding (same as metadata column 2)
	# 	4. var_only	1=variant only, 0=presence/absence (same as metadata column 3)
	# 	5. flag	cluster flag
	# 	6. reads	number of reads in this cluster
	# 	7. cluster	name of cluster
	# 	8. ref_len	 length of reference sequence
	# 	9. ref_base_assembled	number of reference nucleotides assembled by this contig
	# 	10. pc_ident	%identity between reference sequence and contig
	# 	11. ctg	name of contig matching reference
	# 	12. ctg_len	length of contig
	# 	13. ctg_cov	mean mapped read depth of this contig
	# 	14. known_var	is this a known SNP from reference metadata? 1 or 0
	# 	15. var_type	The type of variant. Currently only SNP supported
	# 	16. var_seq_type	Variant sequence type. if known_var=1, n or p for nucleotide or protein
	# 	17. known_var_change	if known_var=1, the wild/variant change, eg I42L
	# 	18. has_known_var	if known_var=1, 1 or 0 for whether or not the assembly has the variant
	# 	19. ref_ctg_change	amino acid or nucleotide change between reference and contig, eg I42L
	# 	20. ref_ctg_effect	effect of change between reference and contig, eg SYS, NONSYN (amino acid changes only)
	# 	21. ref_start	start position of variant in reference
	# 	22. ref_end	end position of variant in reference
	# 	23. ref_nt	nucleotide(s) in reference at variant position
	# 	24. ctg_start	start position of variant in contig
	# 	25. ctg_end	end position of variant in contig
	# 	26. ctg_nt	nucleotide(s) in contig at variant position
	# 	27. smtls_total_depth	 total read depth at variant start position in contig, reported by mpileup
	# 	28. smtls_nts	nucleotides on contig, as reported by mpileup. The first is the contig nucleotide
	# 	29. smtls_nts_depth	depths on contig, as reported by mpileup. One number per nucleotide in the previous column
	# 	30. var_description	description of variant from reference metdata
	# 	31. free_text	other free text about reference sequence, from reference metadata
	#########################################################
	
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

