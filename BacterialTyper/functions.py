#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Common functions.
- time
- system call
- files/folders
- aesthetics
- fasta files
- miscellaneous
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

## import my modules
from BacterialTyper import config

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

	# Group dataframe by sample name
	sample_frame = pd_samples.groupby(["name"])
	
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
		mode_name_dir = create_subfolder(mode, outdir)		
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
	my_all_list = get_fullpath_list(folder)
	matching = [s for s in my_all_list if s.endswith(string)]
	return (matching)
	
#####
def file2dictionary(file2read):
	d = {}
	with open(file2read) as f:
		for line in f:
			(key, val) = line.split()
			d[key] = val

	return(d)
	
	
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
		print (colored("** ERROR **", 'red'))
		print (colored(err.output, 'red'))
		print (colored("** ERROR **", 'red'))
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
	makeblastDBexe = config.get_exe('makeblastdb')
	
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
	blastnexe = config.get_exe('blastn')
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
	"""
	Prints a common header for the pipeline including name, author, copyright and year.	    
	"""
	print ("\n")
	print_sepLine("#", 70, False)
	print('#', '{: ^66}'.format("BacterialTyper pipeline"), '#')
	print('#', '{: ^66}'.format("Jose F. Sanchez, Cristina Prat & Lauro Sumoy"), '#')
	print('#', '{: ^66}'.format("Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain"), '#')
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
def _access_check(fn, mode):
	## the original code belongs to shutil, slightly modified here
	# https://github.com/python/cpython/blob/master/Lib/shutil.py
	
	if os.path.isdir(fn):
		return False

	if os.path.exists(fn):
		if fn.endswith('.jar'):
			return True

		if os.access(fn, mode):
			return True
	
#################
def my_which(cmd, mode=os.F_OK | os.X_OK, path=None):
	## the original code belongs to shutil, slightly modified here
	# https://github.com/python/cpython/blob/master/Lib/shutil.py

	"""Given a command, mode, and a PATH string, return the path which
	conforms to the given mode on the PATH, or None if there is no such
	file.
	`mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
	of os.environ.get("PATH"), or can be overridden with a custom search
	path.
	"""
	# If we're given a path with a directory part, look it up directly rather
	# than referring to PATH directories. This includes checking relative to the
	# current directory, e.g. ./script
	if os.path.dirname(cmd):
		if _access_check(cmd, mode):
			return cmd
		return None

	use_bytes = isinstance(cmd, bytes)

	if path is None:
		path = os.environ.get("PATH", None)

	if path is None:
		try:
			path = os.confstr("CS_PATH")
		except (AttributeError, ValueError):
			# os.confstr() or CS_PATH is not available
			path = os.defpath
		# bpo-35755: Don't use os.defpath if the PATH environment variable is
		# set to an empty string

	# PATH='' doesn't match, whereas PATH=':' looks in the current directory
	if not path:
		return None

	if use_bytes:
		path = os.fsencode(path)
		path = path.split(os.fsencode(os.pathsep))
	else:
		path = os.fsdecode(path)
		path = path.split(os.pathsep)

	if sys.platform == "win32":
		# The current directory takes precedence on Windows.
		curdir = os.curdir

		if use_bytes:
			curdir = os.fsencode(curdir)
	
			if curdir not in path:
				path.insert(0, curdir)

			# PATHEXT is necessary to check on Windows.
			pathext = os.environ.get("PATHEXT", "").split(os.pathsep)

			if use_bytes:
				pathext = [os.fsencode(ext) for ext in pathext]
				# See if the given file matches any of the expected path extensions.
				# This will allow us to short circuit when given "python.exe".
				# If it does match, only test that one, otherwise we have to try
				# others.
			if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
				files = [cmd]
			else:
				files = [cmd + ext for ext in pathext]
	else:
		# On other platforms you don't have things like PATHEXT to tell you
		# what file suffixes are executable, so just pass on cmd as-is.
		files = [cmd]

	#print (files)
	#print (path)

	return_paths = [] ## modification
	seen = set()
	for dir in path:
		normdir = os.path.normcase(dir)
		#print ("Normdir: ", normdir)
		if not normdir in seen:
			seen.add(normdir)
			for thefile in files:
				name = os.path.join(dir, thefile)
				#print ("Name: ", name)
				if _access_check(name, mode):
					## return (name) ## previously, it would only return the first item
					return_paths.append(name) ## modification
	
	if (len(return_paths) >= 1):
		return return_paths
	else:
		return None

