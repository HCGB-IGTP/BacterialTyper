#usr/bin/env python

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from datetime import datetime


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
