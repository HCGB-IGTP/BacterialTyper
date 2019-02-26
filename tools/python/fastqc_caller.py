#usr/bin/env python
'''
This code calls fastqc analysis and parses results generated.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
import fastqcparser
from sys import argv
from io import open

## import
thisDir = os.path.dirname(os.path.abspath(argv[0]))
sys.path.append(thisDir)
import functions
import config

######
class fastqcObject: 

    def __init__(self, sample=None, path=None, StatisticsFrame=None, StatusFrame=None):
        self.__sample=sample
        self.__path = path
        self.__StatusFrame = StatusFrame
        self.__StatisticsFrame = StatisticsFrame
        
    def print_header(self):
        if self.__sample:
            print("Working on sample: %s\nPath: %s " %(self.__sample, self.__path))
            
    def set_sample(self, sample):
        self.__sample = sample
        
    def get_sample(self):
        return self.__sample    

    def set_StatusFrame(self, by):
        self.__StatusFrame = by
        
    def get_StatusFrame(self):
        return self.__StatusFrame
    
    def set_StatisticsFrame(self, by):
        self.__StatisticsFrame = by
        
    def get_StatisticsFrame(self):
        return self.__StatisticsFrame
    
    def __repr__(self):
        return "fastqcObject('" + self.__path + "', Name: " +  str(self.__sample) + ")"

    def __str__(self):
        return "fastqcObject: " + self.__path + ", Name: " +  str(self.__sample)

######

def	help_options():
	print ("\nUSAGE: python %s folder file1 file2 name fastqc\n"  %os.path.abspath(argv[0]))

######
def call_fastqc(path, file1, file2, sample, fastqc_bin):	
	## call system for fastqc sample given
	name = functions.create_subfolder(sample, path)
	logFile = path + '/' + sample + '.log'
	
	if os.path.isfile(logFile):
		return ('OK')
	
	print ("+ Calling fastqc for samples...")	
	cmd_fastqc = '%s --extract -o %s %s %s > %s 2> %s' %(fastqc_bin, name, file1, file2, logFile, logFile)
	## send command	
	return (functions.system_call( cmd_fastqc ))
		
######	
def parse_fastqcFile(resultsfile):
	## parse fastqc_data.txt file 
	print ("+ Parsing results from fastqc analysis")
	# load file
	f = fastqcparser.FastQCParser(resultsfile)
	#gc_content = f.modules['Basic Statistics']['data'][-1][1]
	statistics = {
		'filename':f.filename,
		'encoding':f.encoding,
		'sequences':f.total_sequences,
		'filtered':f.filtered_sequences,
		'GC':f.modules['Basic Statistics']['data'][-1][1],
	}	
	status = {}	
	for values in f.modules.keys():
		status[values] = f.modules[values]['status']
		#print ('[ ' + values + ' ]: ' + f.modules[values]['status'])
	
	#print (f.filename), #print (f.file_type),	#print (f.encoding), #print (f.total_sequences), 
	#print (f.filtered_sequences), #print (f.sequence_length), #print ('%GC ' + str(gc_content))
	return (statistics, status)
	
######
def run_module_fastqc(path, file1, file2, sample):	
	## Arguments provided via ARGVs
	fastqc_bin = config.CONFIGURATION['fastqc']
	codeReturn = call_fastqc(path, file1, file2, sample, fastqc_bin)
	if codeReturn == 'FAIL':
		exit()
	path_to_sample = path + '/' + sample
	return path_to_sample

######
def get_files(path):
	files = os.listdir(path)
	fastqc_files = []
	## get files generated with summary information
	for f in files:
		if f.endswith('_fastqc'):
			tmp = path + '/' + f + '/fastqc_data.txt'
			if os.path.isfile(tmp):
				fastqc_files.append(tmp)
	return fastqc_files

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	path = os.path.abspath(argv[1])
	file1 = os.path.abspath(argv[2])
	file2 = os.path.abspath(argv[3])
	sample = argv[4]
	fastqc_bin = argv[5]

	##
	path_to_sample = run_module_fastqc(path, file1, file2, sample, fastqc_bin)
	fastqc_files = get_files(path_to_sample)
	
	for files in fastqc_files:
		parse_fastqcFile(files)
		
######

'''******************************************'''
if __name__== "__main__":
	main()








