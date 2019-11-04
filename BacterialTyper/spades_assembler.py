#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
This code calls SPADES assembler and plasmidSPADES.
Retrieves plasmids if any and discard them from the main assembly.
"""
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
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper.blast_parser import parse
from BacterialTyper.other_tools import tools

## perl scripts
contig_stats_script = tools.perl_scripts('contig_stats')
rename_seqs_script = tools.perl_scripts('rename_FASTA_seqs')

################################################
def run_SPADES_plasmid_assembly(path, file1, file2, sample, SPADES_bin, threads):
	"""Generates plasmid assembly: call run_SPADES with plasmid option
	
	Arguments:
		path : absolute path to store results. It must exists.
		file1 : FASTQ R1 reads.
		file2 : FASTQ R2 reads.
		sample :	Sample name or tag to identify.
		SPADES_bin : Binary executable for SPADES assembly software.
	 	threads : Number of CPUs to use.
	
	Returns:
		plasmids assembled, if succeded
		FAIL 
	"""
	
	print ('+ Running plasmid assembly...')
	name = sample + '_plasmid'
	options = '--plasmid '
	message_return = run_SPADES(path, file1, file2, name, SPADES_bin, options, threads)

	if 	message_return == 'FAIL':	
		print ("\n\n***ERROR: plasmidSPADES failed for sample " + sample)	
		exit()

	scaffolds_retrieved = functions.retrieve_matching_files(path + '/' + name, "scaffolds.fasta")
	if scaffolds_retrieved == '':	
		print ('\n\n***ATTENTION: No plasmids assembly...')

	return (scaffolds_retrieved[0])

################################################
def run_SPADES_assembly(path, file1, file2, sample, SPADES_bin, threads):
	"""Generate main assembly: call run_SPADES and rename contigs
	
	Arguments:
		path : absolute path to store results. It must exists.
		file1 : FASTQ R1 reads.
		file2 : FASTQ R2 reads.
		sample :	Sample name or tag to identify.
		SPADES_bin : Binary executable for SPADES assembly software.
	 	threads : Number of CPUs to use.
	
	Returns:
		contigs assembled renamed, if succeded
		FAIL 
	"""
	##print ('+ Running main assembly...')
	options = ''
	message_return = run_SPADES(path, file1, file2, sample, SPADES_bin, options, threads)
	if 	message_return == 'FAIL':	
		print ("\n\n***ERROR: SPADES failed for sample " + sample)
		return ('FAIL')

	scaffolds_retrieved = functions.retrieve_matching_files(path, "scaffolds.fasta")
	new_contigs = path + '/' + sample + '_assembly.fna'
	rename_contigs(scaffolds_retrieved[0], "scaffolds_" + sample, new_contigs)
		
	if scaffolds_retrieved == '':	
		print ('\n\n***ERROR: No scaffolds assembly...')
		return ('FAIL')
			
	return (new_contigs)

################################################
def run_SPADES(sample_folder, file1, file2, name, SPADES_bin, options, threads):
	"""Generate SPADES system call.
	
	Arguments:
		sample_folder : absolute path to store results. It must exists.
		file1 : FASTQ R1 reads.
		file2 : FASTQ R2 reads.
		name :	Sample name or tag to identify.
		SPADES_bin : Binary executable for SPADES assembly software.
		options : Plasmid assembly is possible if specificed via options (--plasmid). 
	 	threads : Number of CPUs to use.
	 
	 Returns:
	 	OK: If assembly succeeded generates timestamp file (.success_assembly) within sample_folder provided.
	 	FAIL: If some error ocurred during the assembly.
	 	
	"""
	
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success_assembly'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
		return('OK')

	## call system for SPADES sample given
	logFile = sample_folder + '/' + name + '.log'
	
	## command	
	cmd_SPADES = '%s %s-t %s -o %s -1 %s -2 %s > %s 2> %s' %(SPADES_bin, options, threads, sample_folder, file1, file2, logFile, logFile)
	code = functions.system_call(cmd_SPADES)
	
	if (code == 'OK'):
		## success stamps
		filename_stamp = sample_folder + '/.success_assembly'
		stamp =	functions.print_time_stamp(filename_stamp)
		return('OK')

	return "FAIL"
	
################################################
def run_module_SPADES(name, folder, file1, file2, threads):
	"""Prepare SPADES system call.
	
	Arguments:
		name :	Sample name or tag to identify.
		folder : absolute path to store results. It must exists.
		file1 : FASTQ R1 reads.
		file2 : FASTQ R2 reads.
		threads : Number of CPUs to use.
	 
	 Returns:
	 	Assembly statistics file, if succeeded.
	 	FAIL: If some error ocurred during the assembly.
	 	
	"""
	
	print ("+ Calling spades assembly for sample...", name)	
	
	## get configuration
	SPADES_bin = config.get_exe('spades')
	
	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, name, SPADES_bin, threads)
	
	if path_to_contigs == 'FAIL':
		return ('FAIL')
	else:
		## contig stats
		#print ('+ Get assembly statistics:...\n')
		contig_out = contig_stats(path_to_contigs)	
	
		## check statistics in file
		print ("+ Check statistics for sample %s in file:\n%s" %(name, contig_out))
		return(contig_out)

################################################
def discardPlasmids(contigs, plasmids, path, sample):
	
	## check if any plasmids
	if (plasmids == 'FAIL'):
		#print ('+ No plasmids assembled.')
		#print ('+ No need to discard any plasmids from the main assembly')
		
		contig_out_file = os.path.dirname(path) + '/' + sample + '/' + sample + '_chromosome.fna.tmp'
		shutil.copy(contigs, contig_out_file)
		return (contig_out_file, plasmids)
	
	## discard 
	print ('+ Check if any plasmids are also reported in main assembly...')

	folder = functions.create_subfolder('blast_search', path)	
	
	## makeblastDB
	dbName = folder + '/mainAssembly'
	functions.makeblastdb(dbName, contigs)
	
	## blastn command
	outFile = folder + '/blastn_output.txt'
	threads = 1
	functions.blastn(outFile, dbName, plasmids, threads)
	
	########################
	## parseBlast results
	########################
	
	## thresholds
	eval_thresh_float = float(1e-20)
	aln_thresh_given = 90
	min_length = 1000
	
	outFile_parsed = folder + '/blastn_output_parsed.txt'
	output_file = open(outFile_parsed, 'w')	
	sequences2discard = []

	print ('+ Parsing BLAST results generated...\n')
	## get results
	fh = open(outFile)
	for blast_record in parse(fh, eval_thresh=eval_thresh_float, aln_thresh=aln_thresh_given, length_thresh=min_length):
		for hit in blast_record.hits:
			for hsp in hit:
				output_file.write('****Alignment****')
				output_file.write('\n')
				
				output_file.write('query id: {}'.format(blast_record.qid))
				output_file.write('\n')
				
				sequences2discard.append(hsp.sid)
				output_file.write('sequence: %s' %hsp.sid)
				output_file.write('\n')
				
				output_file.write('e value: %s' %hsp.evalue)
				output_file.write('\n')
				
				output_file.write('aln: %s' %hsp.length)
				output_file.write('\n')
				
				output_file.write('qlen: %s [>%s]' %(hsp.qlen, min_length))
				output_file.write('\n')
				
				aln = (int(hsp.qlen)/int(hsp.slen))*100
				output_file.write('aln/slen: %s [> %s]' %(aln, aln_thresh_given))
				output_file.write('\n\n')

	fh.close()
	output_file.close()
	
	items = len(sequences2discard)
	print ('There are %s sequences to discard from main assembly identified as plasmids' %items)

	## print filtered contigs
	contig_out_file = os.path.dirname(path) + '/' + sample + '/' + sample + '_chromosome.fna.tmp'
	plasmid_out_file = os.path.dirname(path) + '/' + sample + '/' + sample + '_plasmid.fna.tmp'
		
	contig_out_file_handle = open(contig_out_file, 'w')
	plasmid_out_file_handle = open(plasmid_out_file, 'w')
	
	contig_items = SeqIO.parse(contigs, 'fasta')
	for seq in contig_items:
		if seq.id in sequences2discard:
			plasmid_out_file_handle.write(seq.format("fasta"))
			plasmid_out_file_handle.write('\n')
		else:
			contig_out_file_handle.write(seq.format("fasta"))
			contig_out_file_handle.write('\n')

	contig_out_file_handle.close()
	plasmid_out_file_handle.close()	
	
	return (contig_out_file, plasmid_out_file)

################################################
def contig_stats(sequences):
	"""Generate assembly statistics
	
	Calls additional perl script to generate contig statistics:
	Usage:	
	perl BacterialTyper/other_tools/perl/contig_stats.pl fasta_file 
	- Provide a single fasta file for Contig Statistics...
	- Default splitting sets: 0, 150, 500, 1000, 5000, 10000
	- Provide new parts using a csv argument for the script: perl BacterialTyper/other_tools/perl/contig_stats.pl fasta_file 1000,10000
	"""
	file_out = sequences + '_stats.txt'
	cmd_stats = 'perl %s %s 1000,10000 > %s' %(contig_stats_script, sequences, file_out) ## [TODO] Generate this code in python
	code_chr = functions.system_call(cmd_stats)
	return (file_out)

################################################
def stats(new_contigs, new_plasmids):
	"""Generate assembly sequence statistics.
	
	For contigs and plasmids (if any) prints statistics in screen and in csv/txt files 
	"""
	## generate contig statistics
	print ('+ Get assembly statistics:...\n')

	## get contig statistics	
	contig_out = contig_stats(new_contigs)	
	contig_out_file = open(contig_out, 'r')
	contig_out_file_read = contig_out_file.read()
	contig_out_file.close()
	
	## dump in screen
	print (contig_out_file_read)
	print ()	

	if (new_plasmids == 'FAIL'):
		print ('+ No plasmids identified...\n')
	else:
		print ('+ Plasmids assembly')
		plasmid_out = contig_stats(new_plasmids)	

		## dump in screen
		plasmid_out_file = open(plasmid_out, 'r')
		plasmid_file_read = plasmid_out_file.read()
		plasmid_out_file.close()
		print(plasmid_file_read)	
	
################################################
def rename_contigs(fasta_file, name, new_fasta):
	"""
	Calls additional perl script to rename contigs:

	Usage:
	Please provide the next arguments:
	perl BacterialTyper/other_tools/perl/rename_FASTA_seqs.pl fasta_file name_file name2add ADD|REPLACE|BEGIN|ADD_BEGIN
	
	ADD: will add the id plus a counter
	REPLACE: will discard the previous name and add this unique id and counter
	BEGIN: will keep the first split of the id and add at the beginning the given name
	"""
	perl_call = "perl %s %s %s %s REPLACE" %(rename_seqs_script, fasta_file, new_fasta, name) ## [TODO] Generate this code in python
	return (functions.system_call(perl_call))
	
################################################
def	help_options():
	print ("\nUSAGE: python %s file1 file2 name SPADES_bin threads path\n"  %os.path.realpath(__file__))

################################################
def main():
	  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	file1 = os.path.abspath(argv[1])
	file2 = os.path.abspath(argv[2])
	sample = argv[3]
	SPADES_bin = argv[4]
	threads = int(argv[5])
	path = 	argv[6]

	folder = functions.create_subfolder(sample, path)

	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, sample, SPADES_bin, threads)
	
	## assembly plasmids
	path_to_plasmids = run_SPADES_plasmid_assembly(folder, file1, file2, sample, SPADES_bin, threads)

	## discard plasmids from main
	tmp_contigs, tmp_plasmids = discardPlasmids(path_to_contigs, path_to_plasmids, folder, sample)
	
	## rename fasta sequences
	new_contigs_list = tmp_contigs.split(".tmp")
	new_contigs = new_contigs_list[0]
	rename_contigs(tmp_contigs, "scaffolds_chr", new_contigs)
	
	new_plasmids=""
	if os.path.isfile(tmp_plasmids):
		new_plasmids_list = tmp_plasmids.split(".tmp")
		new_plasmids = new_plasmids_list[0]
		rename_contigs(tmp_plasmids, "scaffolds_plasmids", new_plasmids)
	
	## contig stats
	stats(new_contigs, new_plasmids)	
		
######

'''******************************************'''
if __name__== "__main__":
	main()

################################################
def run_module_SPADES_old(name, folder, file1, file2, threads):

	print ("+ Calling spades assembly for sample...", name)	

	## folder create
	functions.create_folder(folder)
	
	## get configuration
	SPADES_bin = config.get_exe('spades')
	
	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, name, SPADES_bin, threads)

	## assembly plasmids
	path_to_plasmids = run_SPADES_plasmid_assembly(folder, file1, file2, name, SPADES_bin, threads)
	
	## discard plasmids from main
	(tmp_contigs, tmp_plasmids) = discardPlasmids(path_to_contigs, path_to_plasmids, folder, name)
	
	## rename fasta sequences
	new_contigs = tmp_contigs.split(".fna.tmp")[0] + '.fna'	
	rename_contigs(tmp_contigs, "scaffolds_chr", new_contigs)
	
	new_plasmids=""
	if os.path.isfile(tmp_plasmids):
		new_plasmids = tmp_plasmids.split(".fna.tmp")[0] + '.fna'	
		rename_contigs(tmp_plasmids, "scaffolds_plasmids", new_plasmids)
	
	## contig stats
	stats(new_contigs, new_plasmids)
	
	## success stamps
	filename_stamp = folder + '/.success'
	stamp =	functions.print_time_stamp(filename_stamp)

