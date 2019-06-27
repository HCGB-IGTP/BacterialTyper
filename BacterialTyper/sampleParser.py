#usr/bin/en python3
'''
This code prepares samples for further analysis.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from termcolor import colored
import shutil
import concurrent.futures

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

## todo check if 'link_files' folder exists within folder provided.

###############
def select_samples (list_samples, samples_prefix, pair=True, exclude=False, merge=False):
    
    #Get all files in the folder "path_to_samples"    
	sample_list = []
	for names in samples_prefix:
		for path_fastq in list_samples:	
			fastq = os.path.basename(path_fastq)
			samplename_search = re.search(r"(%s).*" % names, fastq)
			enter = ""
			if samplename_search:
				if (exclude): ## exclude==True
					enter = False
				else: ## exclude==True
					enter = True
			else:
				if (exclude): ## exclude==True
					enter = True
				else: ## exclude==True
					enter = False
					
			if enter:
				if fastq.endswith('.gz'):
					sample_list.append(path_fastq)
				elif fastq.endswith('fastq'):
					sample_list.append(path_fastq)
				else:
					print ("** ERROR: ", path_fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included')

	## discard duplicates if any
	non_duplicate_samples = list(set(sample_list))	
	discard_samples = []
	
	##	
	found = []

	if (pair):
		## initiate dataframe
		name_columns = ("samples", "R1", "R2")
	else:
		## initiate dataframe
		name_columns = ("samples", "read")

	## initiate dataframe
	df_samples = pd.DataFrame(columns=name_columns)

	for path_files in non_duplicate_samples:
		##	
		files = os.path.basename(path_files)
		trim_search = re.search(r".*trim.*", files)
		lane_search = re.search(r".*L\d+.*", files)

		## get name
		if (pair):
			if (trim_search):
				name_search = re.search(r"(.*)\_trim\_(R1|R2)(\.f.*q)(\..*){0,1}", files)
			else:
				## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
				if (lane_search):
					name_search = re.search(r"(.*)\_(L\d+)\_(R1|R2)\_(\d+)(\.f.*q)(\..*){0,1}", files)
				else:
					name_search = re.search(r"(.*)\_(R1|1|R2|2)(\.f.*q)(\..*){0,1}", files)
		else:
			name_search = re.search(r"(.*)(\.f.*q)(\..*){0,1}", files)
		
		if name_search:
			file_name = name_search.group(1)
			
			if not merge:
				if file_name in found:
					continue
			
			#print ("##")
			#print (files)
			
			if (pair):
				## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
				if (lane_search):
					lane_id = name_search.group(2)
					read_pair = name_search.group(3)
					lane_file = name_search.group(4)
					ext = name_search.group(5)
					gz = name_search.group(6)
				else:
					## could exist or not
					read_pair = name_search.group(2)
					ext = name_search.group(3)
					gz = name_search.group(4)
			else:
				ext = name_search.group(2)
				gz = name_search.group(3)
			
			dirN = os.path.dirname(path_files)
			
			#print (file_name)
			#print (read_pair)
			#print (ext)
			#print (gz)

			if (pair):
				## get second pair
				paired = ""
				R2_paired = ""
				
				if (trim_search):
					if (read_pair == 'R1'):
						paired = dirN + '/' + file_name + '_trim_R2' + ext
						R2_paired=True
					else:
						paired = dirN + '/' + file_name + '_trim_R1' + ext
						R2_paired=False
				else:
				
					#print (read_pair + ext + gz + '\n')
					## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
					if (lane_search):
						if (read_pair == 'R1'):
							paired = dirN + '/' + file_name + '_' + lane_id + '_R2_' + lane_file + ext
							R2_paired=True
						else:
							paired = dirN + '/' + file_name + '_' + lane_id + '_R1_' + lane_file + ext
							R2_paired=False
					else:
						if (read_pair == '1'):
							paired = dirN + '/' + file_name + '_2' + ext
							R2_paired=True
						elif (read_pair == 'R1'):
							paired = dirN + '/' + file_name + '_R2' + ext
							R2_paired=True
						elif (read_pair == 'R2'):
							paired = dirN + '/' + file_name + '_R1' + ext
							R2_paired=True
						else:
							paired = dirN + '/' + file_name + '_1' + ext
							R2_paired=False
		
			found.append(file_name) ## save retrieved samples
			
			## if gunzipped
			if gz:
				unzip = path_files
				if (pair):
					unzip = dirN + '/' + file_name + read_pair + ext
			
				if os.path.isfile(unzip):
					#discard_samples.append(path_files)
					if (pair):
						if os.path.isfile(paired):
							## both files are unzipped and available
							## save both files
							if R2_paired:
								df_samples.loc[len(df_samples)] = [file_name, unzip, paired]							
							else:
								df_samples.loc[len(df_samples)] = [file_name, paired, unzip]						
					else:
						df_samples.loc[len(df_samples)] = [file_name, path_files]							
				else:
					if (pair):
						## gunzipped
						## save both files
						if os.path.isfile(paired + gz):
							if R2_paired:
								df_samples.loc[len(df_samples)] = [file_name, path_files, paired + gz]
							else:
								df_samples.loc[len(df_samples)] = [file_name, paired + gz, path_files]	
					else:
						df_samples.loc[len(df_samples)] = [file_name, path_files]
			
			else:
				## not gunzipped files					
				if (pair):
					if os.path.isfile(paired):
						## save both files
						if R2_paired:
							df_samples.loc[len(df_samples)] = [file_name, path_files, paired]
						else:
							df_samples.loc[len(df_samples)] = [file_name, paired, path_files]							
				else:
					df_samples.loc[len(df_samples)] = [file_name, path_files]							

	## df_samples is a pandas dataframe containing info
	#print (df_samples)	
	number_samples = df_samples.index.size
	if (number_samples == 0):
		print (colored("\n**ERROR: No samples were retrieved. Check the input provided\n",'red'))
		exit()
	print (colored("\t" + str(number_samples) + " samples selected from the input provided...", 'yellow'))
	return (df_samples)


###############
def select_other_samples (list_samples, samples_prefix, mode, extensions, exclude=False):

	#print (list_samples, samples_prefix, mode, extensions)

    #Get all files in the folder "path_to_samples"    
	sample_list = []
	for names in samples_prefix:
		for path_file in list_samples:	
			f = os.path.basename(path_file)
			samplename_search = re.search(r"(%s).*" % names, f)
			enter = ""
			if samplename_search:
				if (exclude): ## exclude==True
					enter = False
				else: ## exclude==True
					enter = True
			else:
				if (exclude): ## exclude==True
					enter = True
				else: ## exclude==True
					enter = False
					
			if enter:
				if mode == 'annotation':
					if f.endswith(extensions):
						sample_list.append(path_file)
				else:
					if mode + '.' in f:
						if f.endswith(extensions):
							sample_list.append(path_file)
					
	## discard duplicates if any
	non_duplicate_samples = list(set(sample_list))	
	discard_samples = []

	## initiate dataframe
	name_columns = ("samples", "tag", "file")

	## initiate dataframe
	df_samples = pd.DataFrame(columns=name_columns)

	## iterate list
	for a_file in non_duplicate_samples:
		a_name = os.path.basename(a_file)

		if mode == 'annotation':
			file_name, ext = os.path.splitext(a_name)
			df_samples.loc[len(df_samples)] = [file_name, ext.split(".")[1], a_file]
	
		else:
			file_name = a_name.split("_" + mode)[0]
			df_samples.loc[len(df_samples)] = [file_name, mode, a_file]

	#print (non_duplicate_samples)
	number_samples = df_samples.index.size
	if (number_samples == 0):
		print (colored("\n**ERROR: No samples were retrieved. Check the input provided\n",'red'))
		exit()
	print (colored("\t" + str(number_samples) + " samples selected from the input provided...", 'yellow'))

	return (df_samples)
	

###############    
def gunzip_merge(outfile, list_files):
	
	list_files = list(list_files)
	list_files.sort()
	print ("\tMerging gz files into: ", outfile)
	print ("\tFiles: ", list_files)

	with open(outfile, 'wb') as wfp:
		for fn in list_files:
			with open(fn, 'rb') as rfp:
				shutil.copyfileobj(rfp, wfp)

	return()
	
###############    
def one_file_per_sample(dataFrame, outdir, threads):
	## merge sequencing files for sample, no matter of sector or lane generated.
	
	list_samples = set(dataFrame['samples'].tolist())
	print (colored("\t" + str(len(list_samples)) + " samples to be merged from the input provided...", 'yellow'))
	
	functions.create_folder(outdir)
	print ("+ Merging sequencing files for samples") 	

	pairs = ['R1', 'R2']
	sample_frame = dataFrame.groupby("samples")

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor: ## need to do 1 by one as there is a problem with the working directory
		for name, cluster in sample_frame: ## loop over samples
			print ("\n\tSample: ", name)
			## send for each sample
			commandsSent = { executor.submit(gunzip_merge, outdir + '/' + name + '_' + pair + '.fq.gz', set(cluster[pair].tolist())): pair for pair in pairs }
			
			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))

	return()
	
###############
def help_options():
	print ("\nUSAGE:\npython %s in_folder out_folder threads\n"  %os.path.abspath(sys.argv[0]))

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	in_folder = os.path.abspath(sys.argv[1])
	out_folder = os.path.abspath(sys.argv[2])
	threads = int(sys.argv[3])
	
	## get files path
	list_files = []
	
	for root, dirs, files in os.walk(in_folder):
		for f in files:
			list_files.append(os.path.join(root,f))
	
	names = ['.*']
	datareturn = select_samples(list_files, names, True, False, True)
	one_file_per_sample(datareturn, out_folder, threads)
	
######
if __name__== "__main__":
	main()

