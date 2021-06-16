#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Calls Snippy for an haploid variant calling and core genome alignment
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

## import my modules
from BacterialTyper.config import set_config
import HCGB.functions.system_call_functions as HCGB_sys

## https://pyvcf.readthedocs.io/en/latest/

##############################
def	help_options():
	print ("\nUSAGE: python %s arguments...\n"  %os.path.realpath(__file__))

##############################
def help_Snippy():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##############################
def parse_snippy_files(folder_results):
	"""
	Parse Snippy output files 
	
	Check details of results generated in :ref:`snippy-output-files`
	"""
	
	return ()

##############################
def snippy_call(reference_fasta, list_files, threads, outdir, name, contig_option, other_options, Debug):
	"""
	Creates variant calling for a sample vs. a reference.
	
	By default, it uses ``rgid`` option with argument ``name`` provided. Argument ``list_files``
	contains files to map that could be a single end file, two paired-end fastq files or a contig file. If 
	fasta contig files provided, set ``contig_option`` True. 
	
	All output files within ``outdir`` folder would containg tag ``snps``
	
	:param reference_fasta: Absolute path to reference fasta file.
	:param list_files: List of absolute path to fastq files (.fq / .fq.gz / fastq / fastq.gz)
	:param threads: Number of CPU cores to use.
	:param outdir: Output folder.
	:param name: Name of the sample
	:param contig_option: True/false to map contigs provided instead of files. Contigs provided via list_files.
	:param other_options: String of options to include in snippy call
	:param Debug: True/false for debugging messages
	
	:type reference_fasta: string
	:type list_files: list
	:type threads: int
	:type outdir: string 
	:type name: string 
	:type contig_options: bool
	:type other_options: string
	:type Debug: bool
	"""
	
	## create snippy call
	snippy_exe = set_config.get_exe('snippy', Debug)
	
	## start snippy_cmd 
	log_file = os.path.join(outdir, "snippy_cmd.log")
	snippy_cmd = '%s --cpus %s --reference %s --force --unmapped --outdir %s --rgid %s' %(
		snippy_exe, threads, reference_fasta, outdir, name)
	
	## force option: prevent finish early if folder exists
	## unmapped option: keep unmapped reads
	
	## add files to map
	if contig_option:
		snippy_cmd = snippy_cmd + ' --ctgs ' + list_files[0]
	else:
		if (len(list_files) == 1):
			snippy_cmd = snippy_cmd + ' --se ' + list_files[0]
		elif (len(list_files) == 2):
			snippy_cmd = snippy_cmd + ' --pe1 ' + list_files[0] + ' --pe2 ' + list_files[1]
		else:
			print(colored("** ERROR: No reads or contigs provided...", "red"))
			return(False)		
	
	## add log
	snippy_cmd = snippy_cmd + ' 2> ' + log_file
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: snippy_cmd **", 'yellow'))	
		print (snippy_cmd)
	
	## create system call
	return(HCGB_sys.system_call(snippy_cmd, returned=False, message=True))

###############################
def snippy_core_call(list_folder, options, name, output_dir, output_format, Debug):
	"""
	Create core alignment for samples align to the same reference
	
	ATTENTION: Requires sample names to be different within the first 10 characters.
	
	:param list_folder:
	:param options:
	:param name:
	:param output_dir:
	:param output_format:
	:param Debug:
	
	:type list_folder: list
	:type options: string
	:type name: string
	:type output_dir:
	:type output_format:
	:type Debug:
	 
	"""
	
	## create snippy-core call
	snippy_core_exe = set_config.get_exe('snippy_core', Debug)
	
	## start snippy_cmd 
	list_folder_string = " ".join(list_folder)
	log_file = os.path.join(output_dir, "snippy_cmd.log")
	name_outdir =  os.path.join(output_dir, name)
	
	## use one reference: must be the same for all comparisons
	reference_fasta = list_folder[1] + "/ref.fa"	
	snippy_core_cmd = '%s -aformat %s --ref %s --prefix %s %s 2> %s' %(snippy_core_exe, output_format, reference_fasta, name_outdir, list_folder_string, log_file)
	
	return (HCGB_sys.system_call(snippy_core_cmd))

##############################
def main():

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	
	
	
####################################
'''******************************************'''
if __name__== "__main__":
	main()



