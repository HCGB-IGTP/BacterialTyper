#usr/bin/env python
'''
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import subprocess
import sys
import wget
from datetime import datetime
from xtract import xtract
from Bio import SeqIO
from termcolor import colored
import pandas as pd

## import configuration
import BacterialTyper 
import config

###############   
def gettime (start_time):
    total_sec = time.time() - start_time
    m, s = divmod(int(total_sec), 60)
    h, m = divmod(m, 60)
    return h, m, s
###############   

###############   
def create_subfolder (name, path):
    ## create subfolder  ##	
	subfolder_path = path + "/" + name
	access_rights = 0o755
	
    # define the access rights
	try:
		os.mkdir(subfolder_path, access_rights)
	except OSError:  
	   	print ("\tDirectory %s already exists" % subfolder_path)
	else:  
		print ("\tSuccessfully created the directory %s " % subfolder_path)
	
	print ("")
	return subfolder_path
###############   
    
###############  
def create_folder (path):
    ## create subfolder  ##	
	access_rights = 0o755
	
    # define the access rights
	try:
		os.mkdir(path, access_rights)
	except OSError:  
	   	print ("\tDirectory %s already exists" %path)
	else:  
		print ("\tSuccessfully created the directory %s " %path)
	
	print ("")
	return path
###############  

###############	
def timestamp (start_time_partial):
	h,m,s = gettime(start_time_partial)
	print_sepLine("-", 25)
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print_sepLine("-", 25)
	return time.time()
############### 

############### 
def get_symbolic_link (sample_list, path_to_samples, directory):
	for samplex in sample_list:
		sample_path = path_to_samples + '/' + samplex
		cmd = 'ln -s %s %s' %(sample_path, directory)
		system_call(cmd)

	files2return = os.listdir(directory)
	return files2return
###############

###############
def system_call(cmd):	
	## call system
	## send command
	print (colored("[** System: %s **]" % cmd, 'green'))
	try:
		subprocess.check_output(cmd, shell = True)
		return ('OK')
	except subprocess.CalledProcessError as err:
		print (colored(err.output, 'red'))
		return ('FAIL')

###############
def wget_download(url, path):
	print ('+ Downloading: ', url)
	wget.download(url, path)
	print ('\n')
###############

###############
def subset_fasta(ident, fasta, out):
	output_FASTA = open(out, 'w')	
	for record in SeqIO.parse(fasta, "fasta"):
		all_id = record.description
		species_search = re.search(r"%s" % ident, all_id)
		if species_search:
			head = ">" + all_id + "\n"
			output_FASTA.write(head)
			output_FASTA.write(str(record.seq))
			output_FASTA.write("\n")

###############

###############
def extract(fileGiven):
	xtract(fileGiven, all=True)
###############

###############
def makeblastdb(DBname, fasta):
	## generate blastdb for genome
	makeblastDBexe = config.EXECUTABLES['makeblastdb']
	
	if (os.path.isfile(DBname + '.nhr')):
		print ("+ BLAST database is already generated...")
	else:
		cmd_makeblast = "%s -in %s -input_type fasta -dbtype %s -out %s" %(makeblastDBexe, fasta, 'nucl', DBname)
		code = system_call(cmd_makeblast)

		if (code == 'FAIL'):
			print (colored('****ERROR: Some error happened during the makeblastDB command', 'red'))
			print (cmd_makeblast)
			exit()
	
###############

###############	
def blastn(outFile, DBname, fasta, threads):
	# blastn plasmids vs contigs
	blastnexe = config.EXECUTABLES['blastn']
	cmd_blastn = "%s -db %s -query %s -out %s -evalue 1e-20 -outfmt \'6 std qlen slen\' -num_threads %s" %(blastnexe, DBname, fasta, outFile, threads )
	codeBlastn = system_call(cmd_blastn)
	
	if (codeBlastn == 'FAIL'):
		print (colored('****ERROR: Some error happened during the blastn command', 'red'))
		print (cmd_blastn)
		exit()
	
###############

###############
def progbar(curr, total, full_progbar):
	frac = (curr/total)
	filled_progbar = round(frac*full_progbar)
	print ('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')
	sys.stdout.flush()
###############

###############
def concat_fasta(dirFasta, Fasta):
	print ("+ Concatenating all information into one file...")	
	cmd = 'cat ' + dirFasta + '/*fna > ' + Fasta
	return(system_call(cmd))
###############

###############
def get_data(ID_file, SEP):	
	print ("+ Obtaining information from file: ", ID_file)
	data = pd.read_csv(ID_file, header=0, sep=SEP)
	print ("\n+ Data:")
	print (data)
	print ("\n\n")	
	return(data)
###############

###############
def get_number_lines(input_file):	
	with open(input_file) as foo:
		lines = len(foo.readlines())
	return (lines)
###############

###############
def pipeline_header():
	print_sepLine("#", 70)
	print ("BacterialTyper pipeline")
	print ("Jose F. Sanchez")
	print ("Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain")
	print_sepLine("#", 70)

###############

###############
def print_sepLine(char, num):
	string = char * num
	print (string)
###############


