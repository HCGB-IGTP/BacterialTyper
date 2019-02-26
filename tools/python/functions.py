#usr/bin/env python

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from datetime import datetime
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
	try:
		subprocess.check_output(cmd, shell = True)
		return ('OK')
	except subprocess.CalledProcessError as err:
		print (err.output)
		return ('FAIL')


