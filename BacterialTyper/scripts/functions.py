#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez					##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Shared functions used along ``BacterialTyper`` pipeline.
With different purposes:

	- Print time stamps

	- Create system calls

	- Manage/list/arrange files/folders

	- Aesthetics

	- Manage fasta files

	- Other miscellaneous functions
"""
## useful imports
import time
import io
import os
import re
import subprocess
import sys
import wget
import json
from datetime import datetime
from Bio import SeqIO
from termcolor import colored
#from filehash import FileHash ## not conda supported
import pandas as pd
import patoolib ## to extract
import shutil

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.other_tools import tools

########################################################################
######## 					TIME								######## 					
########################################################################

###############   
def print_time ():
	"""Prints time stamp in human readable format: month/day/year, hour:minute:seconds."""
	now = datetime.now()
	date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
	print ('\t' + date_time + '\n')

###############   
def gettime (start_time):
	"""Obtains time stamp in human readable format: hour:minute:seconds from a time.time() format timestamp."""
	total_sec = time.time() - start_time
	m, s = divmod(int(total_sec), 60)
	h, m = divmod(m, 60)
	return h, m, s

###############	
def timestamp (start_time_partial):
	"""Prints a stamp of the time spent for a process in human readable format: hour:minute:seconds.
	Returns time in format time.time().
	"""
	h,m,s = gettime(start_time_partial)
	print_sepLine("-", 25, False)
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print_sepLine("-", 25, False)
	return time.time()

###############	
def print_time_stamp (out):
	"""Prints out timestamp in a file provided. Format: time.time()"""
	timefile = open(out, 'w')    
	string2write = str(time.time())
	timefile.write(string2write)
	return()

###############	
def read_time_stamp (out):
	"""Reads timestamp from a file provided. Format: time.time()"""
	st_hd = open(out, 'r')
	st = st_hd.read()
	st_hd.close()
	stamp = datetime.fromtimestamp(float(st)).strftime('%Y-%m-%d %H:%M:%S')
	return(stamp)

###############	
def get_diff_time(stamp):
	"""Obtains the time spent for a process in days given a stamp in time.time() format.
	Returns days passed since.
	"""

	time_today = time.time()
	elapsed = time_today - float(time_today)
	days_passed = int((elapsed/3600)/24)
	return(days_passed)

###############    
def create_human_timestamp():
	"""Generates human timestamp for the date of day in format (yearmonthday): e.g. 20191011"""
	now = datetime.now()
	timeprint = now.strftime("%Y%m%d")
	return timeprint
###############

############################################################################
######## 					FILES/FOLDERS							######## 					
############################################################################

###############
def is_non_zero_file(fpath):  
	# https://stackoverflow.com/a/15924160
	"""Returns TRUE/FALSE if file exists and non zero"""
	return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

###############
def outdir_project(outdir, project_mode, pd_samples, mode):
	"""
	"""
	# Group dataframe by sample name
	sample_frame = pd_samples.groupby(["new_name"])

	dict_outdir = {}	
	for name, cluster in sample_frame:
		if (project_mode):
			#print ("Create subdir for every sample: ", mode)
			sample_dir = create_subfolder('data', outdir)		

			## create sample
			sample_name_dir = create_subfolder(name, sample_dir)		

			## create subdir sub sample
			mode_name_dir = create_subfolder(mode, sample_name_dir)		
			dict_outdir[name] = mode_name_dir

		else:
			#print ("All samples share same folder")
			sample_name_dir = create_subfolder(name, outdir)		
			dict_outdir[name] = sample_name_dir

	return (dict_outdir)

###############
def outdir_subproject(outdir, pd_samples, mode):
	## we assume we want to create within a project dir a subdir
	# Group dataframe by sample name
	sample_frame = pd_samples.groupby(["name"])
	dict_outdir = {}	
	for name, cluster in sample_frame:
		mode_name_dir = create_subfolder(mode,os.path.join(outdir, name))		
		dict_outdir[name] = mode_name_dir

	return (dict_outdir)

###############
def create_subfolder (name, path):
	"""Create a subfolder named 'name' in directory 'path'. Returns path created."""
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
	"""Create a folder directory 'path'. Returns path created."""

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
	"""Creates symbolic links, using system call, for list of files given in directory provided"""
	for samplex in sample_list:
		cmd = 'ln -s %s %s' %(samplex, directory)
		system_call(cmd, returned=False)

	files2return = os.listdir(directory)
	return files2return

############### 
def get_symbolic_link_file (file2link, newfile):
	"""Creates symbolic link for a file into a new name file"""
	cmd = 'ln -s %s %s' %(file2link, newfile)
	system_call(cmd, returned=False)

#################
def get_fullpath_list(dir_given):
	"""Retrieve full absolute path for the files within a directory specified.

	:param dir_given: Directory to retrieve files
	:type dir_given: string

	:returns: List of absolute path files.
	"""
	return_path = []
	for root, dirs, files in os.walk(dir_given):
		for f in files:
			return_path.append(os.path.join(root,f))
	return return_path

#################
def printList2file(fileGiven, listGiven):
	"""Prints list given in the output file provided. One item per row."""
	file_hd = open(fileGiven, 'w')
	file_hd.write("\n".join(listGiven))
	file_hd.close()  

#################
def readList_fromFile(fileGiven):
	"""Reads list from the input file provided. One item per row."""
	# open file and read content into list
	lineList = [line.rstrip('\n') for line in open(fileGiven)]
	return (lineList)

######
def retrieve_matching_files(folder, string):
	"""Lists folder path provided and given a string to search, returns all files ending with the given string"""
	my_all_list = get_fullpath_list(folder)
	matching = [s for s in my_all_list if s.endswith(string)]
	return (matching)

#####
def file2dictionary(file2read, split_char):
	"""Read file and generate a dictionary"""
	d = {}
	with open(file2read) as f:
		for line in f:
			line = line.rstrip('\n')
			(key, val) = line.split(split_char)
			d[key] = val

	return(d)

def file2dataframe(file2read, names):
	## TODO: check if duplicated with get_data function
	"""Read csv file into pandas dataframe"""
	d = pd.read_csv(file2read, comment="#", names=names)
	return(d)


def chmod_rights(file, access_rights):
	"""Changes access permission of a file:
	
	Read additional information in:
	https://docs.python.org/3/library/os.html#os.chmod
	https://www.geeksforgeeks.org/python-os-chmod-method/
	"""
	access = os.chmod(file, access_rights, follow_symlinks=True)


########################################################################
######## 					system call							######## 					
########################################################################

###############
def system_call(cmd, returned=False, message=True):
	"""Generates system call using subprocess.check_output"""
	## call system
	## send command
	if (message):
		print (colored("[** System: %s **]" % cmd, 'magenta'))

	try:
		out = subprocess.check_output(cmd, shell = True)
		if (returned):
			return (out)
		return ('OK')
	except subprocess.CalledProcessError as err:
		if (returned):
			return (err.output)
		if (message):
			print (colored("** ERROR **", 'red'))
			print (colored(err.output, 'red'))
			print (colored("** ERROR **", 'red'))
		
		return ('FAIL')

###############
def wget_download(url, path):
	"""
	Downloads file in the given path.
	
	It uses python module wget to download given ``url`` in the ``path`` provided.
	
	:param url: File to download from http/ftp site.
	:param path: Absolute path to the folder where to save the downloaded file.
	:type url: string
	:type path: string  
	
	:returns: Print messages and generates download
	"""
	
	print ('\t+ Downloading: ', url)
	wget.download(url, path)
	print ('\n')

###############
## not conda supported
def check_md5sum(string, File):
	md5hasher = FileHash('md5')
	md5_file = md5hasher.hash_file(File)
	
	if (md5_file == string):
		#print (md5_file + '==' + string)
		return (True)
	else:
		return (False)

###############
def extract(fileGiven, out, remove=True):
	"""
	Extracts archived file
	
	This function extracts the file given in the ``out`` path provided.
	It uses the ``patoolib`` that is able to identify the type of file 
	and compression algorithm to use in each case.
	
	It also removes the compressed file using ``os`` module.
	
	:param fileGiven: Archived file to extract.
	:param out: Output name and absolute path for the extracted archived.
	:param remove: True/False for removing compressed file extracted
	
	:type fileGiven: string
	:type out: string
	:type remove: boolean
	
	"""
	## extract using patoolib
	patoolib.extract_archive(fileGiven, outdir=out, verbosity=0)
	
	if (remove):
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
	makeblastDBexe = set_config.get_exe('makeblastdb')
	
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
	blastnexe = set_config.get_exe('blastn')
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
	
	output_FASTA.close()

###############
def rename_fasta_seqs(fasta_file, name, new_fasta):
	"""Rename fasta sequences provided in file :file:`fasta_file` using id :file:`name`. Save results in file :file:`new_fasta` provided.
	
	Check for id character lenght do not exceed 37 characters as it might be a limitiation in further annotation and subsequent analysis. Read Prokka_ issue for further details: https://github.com/tseemann/prokka/issues/337.
	
	:param fasta_file: Absolute path to fasta file.
	:type fasta_file: string
	:param name: String to add every fasta sequence header.
	:type name: string
	:param new_fasta: Name for the new fasta file (Absolute path).
	:type new_fasta: string
	:return: Path to tabular delimited file containing conversion from all to new id for each sequence.
	:warnings: Returns FAIL if name is >37 characters.
	
	.. include:: ../../links.inc	 	
	"""

	output_FASTA = open(new_fasta, 'w')
	id_conversion = open(new_fasta + "_conversionID.txt", 'w')
	
	if len(name) > 37:
		print (colored("** ERROR **", 'red'))
		print (colored("BacterialTyper.functions.rename_fasta_seqs():: name id is > 37 characters.", 'red'))
		print (colored("** ERROR **", 'red'))
		return ('FAIL')
	
	counter_seqs = 0
	for record in SeqIO.parse(fasta_file, "fasta"):
		old_id = record.description
		counter_seqs += 1
		new_id = name + "_" + str(counter_seqs)
		head = ">" + new_id + "\n"
		output_FASTA.write(head)
		output_FASTA.write(str(record.seq))
		output_FASTA.write("\n")
		
		id_conversion.write(old_id + "\t" + new_id)
		id_conversion.write("\n")
	
	output_FASTA.close()
	id_conversion.close()
		
	return (new_fasta + "_conversionID.txt")
	

############################################################################
######## 					AESTHETICS								######## 					
############################################################################
def pipeline_header():
	"""
	Prints a common header for the pipeline including name, author, copyright and year.	    
	"""
	print ("\n")
	print_sepLine("#", 70, False)
	print('#', '{: ^66}'.format("BacterialTyper pipeline"), '#')
	print('#', '{: ^66}'.format("Jose F. Sanchez, Cristina Prat & Lauro Sumoy"), '#')
	print('#', '{: ^66}'.format("Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain"), '#')
	print_sepLine("#", 70, False)

###############
def print_sepLine(char, num, color):
	string = char * num
	if (color):
		print (colored(string, color))
	else:
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
def get_data(ID_file, SEP, options):	
	if options == 'index_col=0':
		data = pd.read_csv(ID_file, sep=SEP, index_col=0)
		return(data)
	else:
		data = pd.read_csv(ID_file, sep=SEP)
		return(data)
	## fix for another example if any
	
###############
def get_number_lines(input_file):	
	with open(input_file) as foo:
		lines = len(foo.readlines())
	foo.close()
	return (lines)

###############
def get_info_file(input_file):	
	with open(input_file) as f:
	    lines = f.read().splitlines() 
	f.close()
	return (lines)
  	    
#################
def optimize_threads(total, samples):
	cpu = int(int(total)/int(samples))
	
	if (cpu==0): ## 5 availabe cpus for 10 samples == 0 cpu
		cpu = 1
	
	return(cpu)

#################
def parse_sublist(lst, ind): 
	## extract elemnts of the list
	## Original Code: https://www.geeksforgeeks.org/python-get-first-element-of-each-sublist/
	return [item[ind] for item in lst]

##################
def decode(x):
	"""
	Python decode string method

	It converts bytes to string.

	:param x: String of text to decode
	:type x: string
	:returns: Text decoded

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)

		Give them credit accordingly.
	"""
	try:
		s = x.decode()
	except:
		return x

	return s

###################
def collect_info_run(out_folder, module, userinput, runInfo):
	"""
	Prints information of the job run in json format
	
	:param out_folder: Absolute path to folder to dump results
	:param module: Module or script name executed
	:param userinput: Contains information about the input parameters
	:param runInfo: Contains configuration detail information to dump
	
	:type out_folder: string 
	:type module: string 
	:type userinput: dict
	:type runInfo: dict
	
	:returns: Prints information in json format in the output folder provided.
	
	Exampple of infromation to include:
	
	userinput = {"filename":infiles, 
				"database":path_to_database,
				"file_format":file_format}
				
	runInfo = { "module":module, 
				"analysis":example,
				"date":date, 
				"time":time}
	
	Original idea extracted from https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/kmerfinder.py
	
	"""
	
	# Collect all data of run-analysis
	data = {service:{}}
	
	data[service]["user_input"] = userinput
	data[service]["run_info"] = runInfo
	
	# Save json output
	result_json_file = out_folder + "/" + module + ".json" 
	with open(result_json_file, "w") as outfile:  
	   json.dump(data, outfile)
	
def print_all_pandaDF(pd_df):
	pd.set_option('display.max_colwidth', None)
	pd.set_option('display.max_columns', None)
	print (pd_df)
