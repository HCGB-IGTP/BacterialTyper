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

###############
def help_format():

	functions.boxymcboxface("Name format for samples")

	print ("Format for fastq files can be:\n")
	functions.print_sepLine("*",20,"red")
	print ("[1] Single end files")
	functions.print_sepLine("*",20,"red")
	print ("name.fastq.gz")
	print ("name.fastq")
	print ("name.fq")
	print ("\n")

	functions.print_sepLine("*",20,"red")
	print ("[2] Paired-end files")
	functions.print_sepLine("*",20,"red")
	print ("Read1 => name_1.fastq.gz\tname_R1.fastq.gz")
	print ("Read2 => name_2.fastq.gz\tname_R2.fastq.gz")
	print ("\n")

	functions.print_sepLine("*",55,"red")
	print ("[3] Containing Lane information (*L00x* and/or *00x*):")
	functions.print_sepLine("*",55,"red")
	print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
	print ("name_L00x_1.fastq.gz\tname_L00x_2.fastq.gz")
	print ("name_L00x_R1_00x.fastq.gz\tname_L00x_R1_00x.fastq.gz")
	print (" ** Merge lane files using module option merge module prep **")
	print ("\n")

	functions.print_sepLine("*",15,"red")
	print ("[4] Extensions:")
	functions.print_sepLine("*",15,"red")
	print ("name_L00x_R2.fastq\tname_L00x_R2.fq\nname_L00x_R2.fastq.gz\tname_L00x_R2.fq.gz")
	print ("\n")

###############
def get_fields(file_name_list, pair=True):

	## init dataframe
	name_columns = ("sample", "dirname", "name", "lane", "read_pair","lane_file","ext","gz")
	name_frame = pd.DataFrame(columns=name_columns)
	
	## loop through list
	for path_files in file_name_list:
		
		## get file name
		file_name = os.path.basename(path_files)
		dirN = os.path.dirname(path_files)

		##	
		trim_search = re.search(r".*trim.*", file_name)
		lane_search = re.search(r".*\_L\d+\_.*", file_name)

		## get name
		if (pair):
		
			## pair could be: R1|R2 or 1|2
			## lane should contain L00x			
	
			if (trim_search):
				name_search = re.search(r"(.*)\_trim\_(R1|1|R2|2)(\.f.*q)(\..*){0,1}", file_name)
			else:
				## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
				if (lane_search):
					name_search = re.search(r"(.*)\_(L\d+)\_(R1|1|R2|2)(.*)(\.f.*q)(\..*){0,1}", file_name)
				else:
					name_search = re.search(r"(.*)\_(R1|1|R2|2)(\.f.*q)(\..*){0,1}", file_name)
		else:
			name_search = re.search(r"(.*)(\.f.*q)(\..*){0,1}", file_name)
	
		### declare
		name= ""
		lane_id= ""
		read_pair= ""
		lane_file= ""
		ext= ""
		gz= ""
	
		if name_search:
			name = name_search.group(1)
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
	
			name_frame.loc [len(name_frame)] = (path_files, dirN, name, lane_id, read_pair, lane_file, ext, gz)
	
		else:
			print ("*** ATTENTION: Sample did not match the possible parsing options...")
			print (file_name)

	return (name_frame)

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
	
	## get fields
	name_frame_samples = get_fields(non_duplicate_samples, pair)	
	number_samples = name_frame_samples.index.size
	total_samples = set(name_frame_samples['name'].to_list())
	
	### get some stats
	if (number_samples == 0):
		print (colored("\n**ERROR: No samples were retrieved. Check the input provided\n",'red'))
		exit()
	print (colored("\t" + str(number_samples) + " files selected...", 'yellow'))
	print (colored("\t" + str(len(total_samples)) + " samples selected...", 'yellow'))
	if (pair):
		print (colored("\tPaired-end mode selected...", 'yellow'))
	else:
		print (colored("\tSingle end mode selected...", 'yellow'))
	
	## return info
	return (name_frame_samples)

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
	print ("\tMerging files into: ", outfile)
	#print ("\tFiles: ", list_files)

	with open(outfile, 'wb') as wfp:
		for fn in list_files:
			with open(fn, 'rb') as rfp:
				shutil.copyfileobj(rfp, wfp)

	return()
	
###############    
def one_file_per_sample(dataFrame, outdir, threads):
	## merge sequencing files for sample, no matter of sector or lane generated.
	
	list_samples = set(dataFrame['name'].tolist())
	print (colored("\t" + str(len(list_samples)) + " samples to be merged from the input provided...", 'yellow'))
	
	functions.create_folder(outdir)
	print ("+ Merging sequencing files for samples")

	##
	sample_frame = dataFrame.groupby(["name", "read_pair"])
	
	### get extension for files
	ext_list = dataFrame.ext.unique()
	gz_list = dataFrame.gz.unique()
	ext = ext_list[0] + gz_list[0] ## might generate a bug if several extension or some zip/unzip files provided

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor: ## need to do 1 by one as there is a problem with the working directory
		commandsSent = { executor.submit(gunzip_merge, outdir + '/' + name[0] + '_' + name[1] + ext, set(cluster["sample"].tolist())): name for name, cluster in sample_frame }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
				
							
	## return output name merged generated in dataframe
	name_columns = ("name", "dirname", "read_pair", "sample", "ext", "gz")
	name_frame = pd.DataFrame(columns=name_columns)

	## print to a file
	merge_details = outdir + '/merge_details.txt'
	merge_details_hd = open(merge_details, 'w')

	for name, cluster in sample_frame: ## loop over samples
		outfile = outdir + '/' + name[0] + '_' + name[1] + ext
		
		merge_details_hd.write("####################\n")		
		merge_details_hd.write("Sample: " + name[0] + '\n')
		merge_details_hd.write("Read: " + name[1] + '\n')
		merge_details_hd.write("Files:\n")
		merge_details_hd.write(",".join(cluster["sample"].tolist()))
		merge_details_hd.write('\n')
		merge_details_hd.write("####################\n")		
		
		name_frame.loc [len(name_frame)] = (name[0], outdir, name[1], outfile, ext_list[0], gz_list[0])

	merge_details_hd.close()
	return(name_frame)
	
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

