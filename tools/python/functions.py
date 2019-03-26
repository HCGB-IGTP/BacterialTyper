#usr/bin/env python

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
	print ('--------------------------------')
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print ('--------------------------------')
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
	print ("[System: %s]" % cmd)
	try:
		subprocess.check_output(cmd, shell = True)
		return ('OK')
	except subprocess.CalledProcessError as err:
		print (err.output)
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

