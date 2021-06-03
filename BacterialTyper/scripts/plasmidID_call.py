#usr/bin/env python
'''
This code calls plasmidID and identifies the plasmids given
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
from Bio import SeqIO
import shutil
import subprocess
import concurrent.futures

## import modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
import HCGB
import HCGB.functions.blast_functions as HCGB_blast

## perl scripts
perlDir = os.path.dirname(os.path.realpath(__file__)) + '/../tools/perl'
retrieve_seqs_script = perlDir + '/get-seq_ids.pl'
fasta_split_script  = perlDir + '/fasta_splitter.pl'

######
def plasmidID_module(file1, file2, plasmid_contigs, plasmid_database_path, outfolder):
	
	clustering_file = filter_and_cluster_database(plasmid_database_path, threads)
	plasmidID_bin = config.EXECUTABLES['plasmidID_bin'] + '/plasmidID.sh'
	plasmidID_call(plasmidID_bin, file1, file2, plasmid_contigs, clustering_file, threads, outfolder)

######
def filter_and_cluster_database(database_path, threads):

	original_file = database_path + '/plasmid_database.fna'
	clustering_file = database_path + '/plasmids_clustered.fna'

	if os.path.isfile(clustering_file):
		print ("+ Plasmid database was previously clustered...")
		return(clustering_file)
	else:
		print ('+ Cluster putative duplicated entries on database...')

	folder = functions.create_subfolder('blast_search', database_path)
	
	## makeblastDB
	dbName = folder + '/plasmids_DB'
	functions.makeblastdb(dbName, original_file)
	
	## split plasmid seqs ids
	tmp_folder = functions.create_subfolder("tmp", folder)
	cmd_split= 'perl %s %s %s %s ' %(fasta_split_script, original_file, threads, tmp_folder)

	print ("\n+ Splitting database file to speed up the computation...")
	functions.system_call(cmd_split)

	## get files in folder
	split_files = os.listdir(tmp_folder)	
	num_threads = 1
	print ("\n+ Sending blastn commands: ")

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
		# Start the load operations and mark each future with its URL
		commandsSent = { executor.submit(functions.blastn, tmp_folder + '/' + fil + '-blast-out.txt', dbName, tmp_folder + '/' + fil, num_threads): fil for fil in split_files }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print('%r generated an exception: %s' % (details, exc))
	
	########################
	## parseBlast results
	########################
	print ("+ Concatenating all results into one file...")	
	outFile_merged = folder + '/blastn_output_concat.txt'
	cmd_cat = 'cat ' + tmp_folder + '/*-blast-out.txt > ' + outFile_merged
	functions.system_call(cmd_cat)
	
	## thresholds
	eval_thresh_float = float(1e-20)
	aln_thresh_given = 90
	min_length = 2000
	
	outFile_parsed = folder + '/blastn_output_parsed.txt'
	output_file = open(outFile_parsed, 'w')	

	print ('+ Parsing BLAST results generated...\n')
	my_cluster_list = []
	sequences2discard = []

	## get number lines
	lines = functions.get_number_lines(outFile_merged)
	n = 0
	
	## get results
	fh = open(outFile_merged)
	for blast_record in HCGB_blast.parse(fh, eval_thresh=eval_thresh_float, aln_thresh=aln_thresh_given, length_thresh=min_length):
		### show progress bar
		n +=10
		functions.progbar(n, lines, 100)

		for hit in blast_record.hits:
			for hsp in hit:
				output_file.write('****Alignment****\n')
				output_file.write('query id: {}\n'.format(blast_record.qid))
				sequences2discard.append(hsp.sid) ## add subject
				output_file.write('sequence: %s\n' %hsp.sid)
				output_file.write('e value: %s\n' %hsp.evalue)
				output_file.write('aln_length: %s\n' %hsp.length)
				output_file.write('qlen: %s [>%s]\n' %(hsp.qlen, min_length))
				aln_perc = (int(hsp.length)/int(hsp.qlen))*100
				output_file.write('aln_perc: %s [> %s]\n\n' %(aln_perc, aln_thresh_given))

				#print('query id: {}'.format(blast_record.qid))
				#print('sequence: ', hsp.sid)
				#print('e value:', hsp.evalue)
				#print('aln_length:', hsp.length)
				#print('qlen:', hsp.qlen, ' [>', min_length, ' bp]')
				#print ('aln_perc:', aln_perc, ' [>', aln_thresh_given ,'%]')
		
			## add query:
			qseq = format(blast_record.qid)
			if (qseq in sequences2discard):
				continue ## avoid adding query if already reported as subject
			else:
				my_cluster_list.append(qseq)

	print ("\n\n")

	fh.close()
	output_file.close()
	
	## get unique ids & print in file
	print ("+ Obtaining clustered sequences...")
	my_cluster_list = set(my_cluster_list)
	list_cluster_file = folder + '/plasmids_clustered_ids.txt'
	functions.printList2file(list_cluster_file, my_cluster_list)
	
	## retrieve ids
	cluster_fasta_file = folder + '/plasmids_clustered.fna'
	cmd= 'perl %s %s %s > %s' %(retrieve_seqs_script, list_cluster_file, original_file, cluster_fasta_file)
	functions.system_call(cmd)
	
	final_cluster_fasta_file = database_path + '/plasmids_clustered.fna'	
	shutil.copy(cluster_fasta_file, final_cluster_fasta_file)
	
	## maybe add a step to build bowtie index from database to speed up later process
	
	
	return(final_cluster_fasta_file)
	
######
def plasmidID_call(plasmidID_bin, file1, file2, plasmid_contigs, database, threads, outfolder):
	name_split = os.path.basename(file1).split("_trim_R1.fastq")
	name = name_split[0]
	cmd = "%s --R1 %s --R2 %s --database %s --sample %s --no-trim --contigs %s --threads %s --group %s" %(plasmidID_bin, file1, file2, database, name, plasmid_contigs, threads, outfolder)
	functions.system_call(cmd)

######
def	help_options():
	print ("\nUSAGE: python %s fileR1 fileR2 database_path plasmidID.sh output_folder threads plasmids_assembled\n"  %os.path.realpath(__file__))

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	file1 = os.path.abspath(argv[1])
	file2 = os.path.abspath(argv[2])
	database_path = os.path.abspath(argv[3])
	plasmidID_bin = os.path.abspath(argv[4])
	path = os.path.abspath(argv[5])
	threads = int(argv[6])
	plasmid_contigs = os.path.abspath(argv[7])

	##
	cluster_database = filter_and_cluster_database(database_path, threads)
	plasmidID_call(plasmidID_bin, file1, file2, plasmid_contigs, cluster_database, threads, path)

######
'''******************************************'''
if __name__== "__main__":
	main()


###############
# version 1.4.1
# plasmidID is a computational pipeline tha reconstruct and annotate the most likely plasmids present in one sample
#
# usage : /imppc/labs/lslab/jsanchez/git_repo/plasmidID/plasmidID.sh <-1 R1> <-2 R2> <-d database(fasta)> <-s sample_name> [-g group_name] [options]
#
#	Mandatory input data:
#	-1 | --R1	<filename>	reads corresponding to paired-end R1 (mandatory)
#	-2 | --R2	<filename>	reads corresponding to paired-end R2 (mandatory)
#	-d | --database	<filename>	database to map and reconstruct (mandatory)
#	-s | --sample	<string>	sample name (mandatory), less than 37 characters
#
#	Optional input data:
#	-g | --group	<string>	group name (optional). If unset, samples will be gathered in NO_GROUP group
#	-c | --contigs	<filename>	file with contigs. If supplied, plasmidID will not assembly reads
#	-a | --annotate <filename>	file with configuration file for specific annotation
#	-o 		<output_dir>	output directory, by default is the current directory
#
#	Pipeline options:
#	--explore	Relaxes default parameters to find less reliable relationships within data supplied and database
#	--only-reconstruct	Database supplied will not be filtered and all sequences will be used as scaffold
#						This option does not require R1 and R2, instead a contig file can be supplied
#
#	Trimming:
#	--trimmomatic-directory Indicate directory holding trimmomatic .jar executable
#	--no-trim	Reads supplied will not be quality trimmed
#
#	Coverage and Clustering:
#	-C | --coverage-cutoff	<int>	minimun coverage percentage to select a plasmid as scafold (0-100), default 80
#	-S | --coverage-summary	<int>	minimun coverage percentage to include plasmids in summary image (0-100), default 90
#	-f | --cluster		<int>	identity percentage to cluster plasmids into the same representative sequence (0-100), default 80
#
#	Contig local alignment
#	-i | --alignment-identity <int>	minimun identity percentage aligned for a contig to annotate, default 90
#	-l | --alignment-percentage <int>	minimun length percentage aligned for a contig to annotate, default 30
#	-L | --length-total	<int>	minimun alignment length to filter blast analysis
#	--extend-annotation <int>	look for annotation over regions with no homology found (base pairs), default 500bp
#
#	Draw images:
#	--config-directory <dir>	directory holding config files, default config_files/
#	--config-file-individual <file-name> file name of the individual file used to reconstruct
#	Additional options:
#
#	-M | --memory	<int>	max memory allowed to use
#	-T | --threads	<int>	number of threads
#	-v | --version		version
#	-h | --help		display usage message
#
#example: ./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d database.fasta -s ECO_553 -G ENTERO
#		./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d PacBio_sample.fasta -c scaffolds.fasta -C 60 -s ECO_60 -G ENTERO --no-trim
