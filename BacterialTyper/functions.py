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
from Bio import SeqIO
from termcolor import colored
from filehash import FileHash
import pandas as pd
import patoolib ## to extract
import shutil

## import configuration
import BacterialTyper 
import config

########################################################################
######## 					TIME								######## 					
########################################################################

###############   
def print_time ():
	now = datetime.now()
	date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
	print (date_time)

###############   
def gettime (start_time):
    total_sec = time.time() - start_time
    m, s = divmod(int(total_sec), 60)
    h, m = divmod(m, 60)
    return h, m, s

###############	
def timestamp (start_time_partial):
	h,m,s = gettime(start_time_partial)
	print_sepLine("-", 25)
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print_sepLine("-", 25)
	return time.time()
	
###############	
def print_time_stamp (out):
	timefile = open(out, 'w')    
	string2write = str(time.time())
	timefile.write(string2write)
	return()

###############	
def read_time_stamp (out):
	st_hd = open(out, 'r')
	st = st_hd.read()
	st_hd.close()
	stamp = datetime.fromtimestamp(float(st)).strftime('%Y-%m-%d %H:%M:%S')
	return(stamp)

###############	
def get_diff_time(stamp):
	time_today = time.time()
	elapsed = time_today - float(st)
	days_passed = int((elapsed/3600)/24)
	return(days_passed)


############################################################################
######## 					FILES/FOLDERS							######## 					
############################################################################

###############   
def create_subfolder (name, path):
    ## create subfolder  ##	
	subfolder_path = path + "/" + name
	access_rights = 0o755
	
    # define the access rights
	try:
		os.mkdir(subfolder_path, access_rights)
	except OSError:  
	   	#print ("\tDirectory %s already exists" % subfolder_path)
		return subfolder_path
	else:  
		print (colored("Successfully created the directory %s " % subfolder_path, 'yellow'))
	
	return subfolder_path
    
###############  
def create_folder (path):
    ## create subfolder  ##	
	access_rights = 0o755
	
    # define the access rights
	try:
		os.mkdir(path, access_rights)
	except OSError:  
	   	#print ("\tDirectory %s already exists" %path)
	   	return path
	else:  
		print (colored("Successfully created the directory %s " %path, 'yellow'))
	
	return path

############### 
def get_symbolic_link (sample_list, directory):
	for samplex in sample_list:
		cmd = 'ln -s %s %s' %(samplex, directory)
		system_call(cmd)

	files2return = os.listdir(directory)
	return files2return

#################
def get_fullpath_list(dir_given):
	return_path = []
	for root, dirs, files in os.walk(dir_given):
		for f in files:
			return_path.append(os.path.join(root,f))
	return return_path

#################
def printList2file(fileGiven, listGiven):
	file_hd = open(fileGiven, 'w')
	file_hd.write("\n".join(listGiven))
	file_hd.close()  

########################################################################
######## 					system call							######## 					
########################################################################

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
	print ('\t+ Downloading: ', url)
	wget.download(url, path)
	print ('\n')

###############
def check_md5sum(string, File):
	md5hasher = FileHash('md5')
	md5_file = md5hasher.hash_file(File)
	
	if (md5_file == string):
		#print (md5_file + '==' + string)
		return (True)
	else:
		return (False)

###############
def extract(fileGiven, out):
	cmd = patoolib.extract_archive(fileGiven, outdir=out, verbosity=0)
	## remove compress file
	print ("Remove compress file...")
	os.remove(fileGiven)
	print ("\n")

####################################################################
######## 					BLAST							######## 					
####################################################################

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
def blastn(outFile, DBname, fasta, threads):
	# blastn plasmids vs contigs
	blastnexe = config.EXECUTABLES['blastn']
	cmd_blastn = "%s -db %s -query %s -out %s -evalue 1e-20 -outfmt \'6 std qlen slen\' -num_threads %s" %(blastnexe, DBname, fasta, outFile, threads )
	codeBlastn = system_call(cmd_blastn)
	
	if (codeBlastn == 'FAIL'):
		print (colored('****ERROR: Some error happened during the blastn command', 'red'))
		print (cmd_blastn)
		exit()
	
########################################################################
######## 					FASTA files							######## 					
########################################################################

###############
def concat_fasta(dirFasta, Fasta):
	print ("+ Concatenating all information into one file...")	
	cmd = 'cat ' + dirFasta + '/*fna > ' + Fasta
	return(system_call(cmd))

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

############################################################################
######## 					AESTHETICS								######## 					
############################################################################
def pipeline_header():
	print ("\n")
	print_sepLine("#", 70)
	print('#', '{: ^66}'.format("BacterialTyper pipeline"), '#')
	print('#', '{: ^66}'.format("Jose F. Sanchez & Lauro Sumoy"), '#')
	print('#', '{: ^66}'.format("Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain"), '#')
	print_sepLine("#", 70)

###############
def print_sepLine(char, num):
	string = char * num
	print (string)

###############
def boxymcboxface(message):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
   	#print('-' * 79)
	print ('\n')
	print('|', '=' * 50, '|', sep='')
	print('|', '{: ^48}'.format(message), '|')
	print('|', '=' * 50, '|', sep='')
	print ('\n')
    #print('-' * 79)

###############
def progbar(curr, total, full_progbar):
	frac = (curr/total)
	filled_progbar = round(frac*full_progbar)
	print ('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')
	sys.stdout.flush()


########################################################################
######## 					Miscellaneous						######## 					
########################################################################

###############
def get_data(ID_file, SEP):	
	print ("+ Obtaining information from file: ", ID_file)
	data = pd.read_csv(ID_file, header=0, sep=SEP)
	print ("\n+ Data:")
	print (data)
	print ("\n\n")	
	return(data)

###############
def get_number_lines(input_file):	
	with open(input_file) as foo:
		lines = len(foo.readlines())
	return (lines)
  	    
#################
def optimize_threads(total, samples):
	cpu = int(int(total)/int(samples))
	
	if (cpu==0): ## 5 availabe cpus for 10 samples == 0 cpu
		cpu = 1
	
	return(cpu)

