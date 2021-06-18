#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
This script calls SPADES assembler and plasmidSPADES. It also generates descriptive statistics of the assembly process.

This script can be called as a single script if desired or called from the BacterialTyper module assemble or as a python api.
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
from BacterialTyper.config import set_config
from BacterialTyper.scripts import assembly_stats_caller

import HCGB
from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.fasta_functions as HCGB_fasta
import HCGB.functions.blast_functions as HCGB_blast
import HCGB.functions.system_call_functions as HCGB_sys

################################################
def run_SPADES_assembly(path, file1, file2, sample, SPADES_bin, threads, debug=False):
	"""Generate main assembly using SPADES
	
	- Calls SPADES to assemble reads (using :func:`BacterialTyper.scripts.spades_assembler.SPADES_systemCall`) 
	
	- SPADES generates a file named as *scaffolds.fasta* within the directory provided. This function retrieves path to contigs/scaffolds assembled (using :func:`HCGB_main.retrieve_matching_files`).
	
	- Renames contigs retrieved using sample name (using :func:`BacterialTyper.scripts.spades_assembler.rename_contigs`).
		
	:param path: Absolute path to folder.
	:param file1: Absolute path to fastq reads (R1).
	:param file2: Absolute path to fastq reads (R2).
	:param sample: Sample name or tag to identify sample
	:param SPADES_bin: Binary executable for SPADES assembly software.
	:param threads: Number of CPUs to use.
	:type path: string
	:type file1: string
	:type file2: string
	:type name: string
	:type threads: integer
	:return: Contigs/scaffolds assembled renamed.
	:rtype: string : Path to assembly fasta file.
	:warnings: Returns **FAIL** if assembly process stopped.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.spades_assembler.SPADES_systemCall`
	
		- :func:`HCGB.functions.main_functions.retrieve_matching_files`
	
		- :func:`HCGB.functions.fasta_functions.rename_fasta_seqs`
	"""
	##print ('+ Running main assembly...')
	options = ''
	message_return = SPADES_systemCall(path, file1, file2, sample, SPADES_bin, options, threads, debug)
	if 	message_return == 'FAIL':	
		print ("\n\n***ERROR: SPADES failed for sample " + sample)
		return ('FAIL')

	scaffolds_retrieved = HCGB_main.retrieve_matching_files(path, "scaffolds.fasta", debug)
	if scaffolds_retrieved == '':	
		print ('\n\n***ERROR: No scaffolds assembly...')
		return ('FAIL')
	
	### Due to limiations with Genbank format, no more thatn 37 characters are supported for 
	### locus tag identification. This might affect later annotation process and subsequent analysis
	### https://github.com/tseemann/prokka/issues/337 
	new_contigs = path + '/' + sample + '_assembly.fna'
	id_conversion_file = HCGB_fasta.rename_fasta_seqs(scaffolds_retrieved[0], sample, new_contigs)
		
	if	id_conversion_file == 'FAIL':	
		print ("\n\n***ERROR: Rename contigs failed for sample " + sample)
		return ('FAIL')
	else:
		print ("+ Name conversion details saved in file " + id_conversion_file)
	
	return (new_contigs)

################################################
def SPADES_systemCall(sample_folder, file1, file2, name, SPADES_bin, options, threads, debug=False):
	"""Generate SPADES system call.
	
	It calls system for SPADES and generates time stamp file in the folder provided (sample_folder + '/.success_assembly') for later analysis.
	
	Steps:
	
	- It generates system call for SPADES assembly. 
	
	- It generates timestamp file.
	
	:param sample_folder: Absolute path to store results. It must exists.
	:param file1: Absolute path to fastq reads (R1).
	:param file2: Absolute path to fastq reads (R2).
	:param name: Sample name or tag to identify sample.
	:param SPADES_bin: Binary executable for SPADES assembly software.
	:param options: Plasmid assembly is possible if specificed via options (--plasmid).
	:param threads: Number of CPUs to use.
	
	:type name: string
	:type sample_folder: string
	:type file1: string
	:type file2: string
	:type SPADES_bin: string
	:type options: string
	:type threads: integer
	
	:return: Returns **OK** if assembly process succeeded and fasta file is generated.
	:rtype: string.
	:warnings: Returns **FAIL** if assembly process stopped.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.main_functions.system_call`
	
		- :func:`HCGB.functions.time_functions.print_time_stamp`
	"""
	
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success_assembly'
	if os.path.isfile(filename_stamp):
		stamp =	HCGB_time.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
		return('OK')

	## call system for SPADES sample given
	logFile = sample_folder + '/' + name + '.log'
	
	## command	
	cmd_SPADES = '%s %s-t %s -o %s -1 %s -2 %s > %s 2> %s' %(SPADES_bin, options, threads, 
															sample_folder, file1, file2, logFile, logFile)
	code = HCGB_sys.system_call(cmd_SPADES)
	
	if (code == 'OK'):
		## success stamps
		filename_stamp = sample_folder + '/.success_assembly'
		stamp =	HCGB_time.print_time_stamp(filename_stamp)
		return('OK')

	return "FAIL"
	
################################################
def run_module_assembly(name, folder, file1, file2, threads, debug=False):
	"""Assembly main module call.
	
	It calls assembly function to process data provided and returns genome statistics. Steps: 
	
	- Retrieves SPADES_ executable (See details :func:`BacterialTyper.scripts.set_config.get_exe`) using the minimun version required (See :func:`BacterialTyper.scripts.set_config.min_version_programs` for details)
	
	- It generates a call to SPADES_ assembler (See :func:`BacterialTyper.scripts.spades_assembler.run_SPADES_assembly`). 
		
	- If assembly succeeds and fasta file is generated under the directory provided, contig statistics are generated (:func:`BacterialTyper.scripts.spades_assembler.contig_stats`).
	
	- It retrieves spades executable using 
	
	:param name: Sample name or tag to identify sample
	:param folder: Absolute path to folder.
	:param file1: Absolute path to fastq reads (R1).
	:param file2: Absolute path to fastq reads (R2).
	:param threads: Number of CPUs to use.
	:type name: string
	:type folder: string
	:type file1: string
	:type file2: string
	:type threads: integer
	:return: Assembly statistics file.
	:rtype: string : Path to file assembly statistics file.
	:warnings: Returns **FAIL** if assembly process stopped.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.set_config.get_exe`
	
		- :func:`BacterialTyper.scripts.spades_assembler.run_SPADES_assembly`
	
		- :func:`BacterialTyper.scripts.set_config.min_version_programs`
	
		- :func:`BacterialTyper.scripts.spades_assembler.contig_stats`
	
	.. include:: ../../links.inc	 	
	"""
	
	print ("+ Calling spades assembly for sample...", name)	
	
	## get configuration
	SPADES_bin = set_config.get_exe('spades')
	
	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, name, SPADES_bin, threads, debug)
	
	if path_to_contigs == 'FAIL':
		return ('FAIL')
	else:
		## contig stats
		#print ('+ Get assembly statistics:...\n')
		(stats_dict, excel_file) = contig_stats(path_to_contigs, debug)
	
		## check statistics in file
		print ("+ Check statistics for sample %s in file:\n\t%s" %(name, excel_file))
		return([stats_dict, excel_file])


################################################

def contig_stats(assembly_file, debug):
	"""Generate assembly statistics
	
	Create assembly statistics using the script assembly_stats_caller and pip module assembly_stats
	
	:param assembly_file: Absolute path to assembly fasta file.
	:type assembly_file: string
	:param debug: Boolean for debugging messages
	:type debug: bool
	:return: Text file containing statistics for later analysis.
	:rtype: string
	
	"""
	name_file_list = assembly_file.split(".fna")
	info = assembly_stats_caller.assembly_stats_caller(assembly_file, name_file_list[0], debug)
	return (info)	

################################################
def	help_options():
	"""Help options when run spades_assembler.py as a single script
	
	Parameters reflected here refer to spades_assembler.py main function.
	
	:param file1: Absolute path to fastq reads (R1).
	:param file2: Absolute path to fastq reads (R2).
	:param sample: Sample name or tag to identify sample
	:param SPADES_bin: Binary executable for SPADES assembly software.
	:param threads: Number of CPUs to use.
	:param path: Absolute path to folder.
	
	:type file1: string
	:type file2: string
	:type sample: string
	:type SPADES_bin: string
	:type threads: integer
	:type path: string	
	"""
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

	folder = HCGB_files.create_subfolder(sample, path)

	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, sample, SPADES_bin, threads, debug=True)
	
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
	
	
	## generate contig statistics
	print ('+ Get assembly statistics:...\n')

	## get contig statistics	
	contig_out = contig_stats(new_contigs, True)	
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
		plasmid_out = contig_stats(new_plasmids, True)	

		## dump in screen
		plasmid_out_file = open(plasmid_out, 'r')
		plasmid_file_read = plasmid_out_file.read()
		plasmid_out_file.close()
		print(plasmid_file_read)	
	
######

'''******************************************'''
if __name__== "__main__":
	main()

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

	folder = HCGB_files.create_subfolder('blast_search', path)	
	
	## makeblastDB
	dbName = folder + '/mainAssembly'
	HCGB_.makeblastdb(dbName, contigs)
	
	## blastn command
	outFile = folder + '/blastn_output.txt'
	threads = 1
	HCGB_blast.blastn(outFile, dbName, plasmids, threads)
	
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
	for blast_record in HCGB_blast.parse(fh, eval_thresh=eval_thresh_float, aln_thresh=aln_thresh_given, length_thresh=min_length):
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
def run_module_SPADES_old(name, folder, file1, file2, threads):

	print ("+ Calling spades assembly for sample...", name)	

	## folder create
	HCGB_files.create_folder(folder)
	
	## get configuration
	SPADES_bin = set_config.get_exe('spades')
	
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
	stamp =	HCGB_time.print_time_stamp(filename_stamp)

################################################
def run_SPADES_plasmid_assembly(path, file1, file2, sample, SPADES_bin, threads, debug=False):
	"""Generate plasmid assembly using SPADES
	
	- Calls SPADES to assemble plasmids using --plasmid option (using :func:`BacterialTyper.scripts.spades_assembler.SPADES_systemCall`) 
	
	- SPADES generates a file named as *scaffolds.fasta* within the directory provided. This function retrieves path to contigs/scaffolds assembled.
	
	:param path: Absolute path to folder.
	:param file1: Absolute path to fastq reads (R1).
	:param file2: Absolute path to fastq reads (R2).
	:param sample: Sample name or tag to identify sample
	:param SPADES_bin: Binary executable for SPADES assembly software.
	:param threads: Number of CPUs to use.
	:type path: string
	:type file1: string
	:type file2: string
	:type name: string
	:type threads: integer
	:return: Plasmid contigs/scaffolds assembled.
	:rtype: string : Path to assembly fasta file.
	:warnings: Returns **FAIL** if assembly process stopped.
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.spades_assembler.SPADES_systemCall`
	
		- :func:`HCGB_main.retrieve_matching_files`
	"""
	print ('+ Running plasmid assembly...')
	name = sample + '_plasmid'
	options = '--plasmid '
	message_return = SPADES_systemCall(path, file1, file2, name, SPADES_bin, options, threads)

	if 	message_return == 'FAIL':	
		print ("\n\n***ERROR: plasmidSPADES failed for sample " + sample)	
		exit()

	scaffolds_retrieved = HCGB_main.retrieve_matching_files(path + '/' + name, "scaffolds.fasta", debug)
	if scaffolds_retrieved == '':	
		print ('\n\n***ATTENTION: No plasmids assembly...')

	return (scaffolds_retrieved[0])