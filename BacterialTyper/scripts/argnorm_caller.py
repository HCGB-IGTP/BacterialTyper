#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:15:29 2023

@author: jsanchez
"""

from BacterialTyper.config import set_config
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.aesthetics_functions as HCGB_aes
from termcolor import colored

## argnorm is a package to standarize/normalize Antibiotic resistance gene (ARG) by mapping to the antibiotic resistance ontology (ARO) by CARD.
## https://github.com/BigDataBiology/argNorm
## installation available with pip: pip install https://github.com/BigDataBiology/argNorm/archive/refs/heads/main.zip
## and soon via pip/conda


###############
def help_argnorm():
    """
    Print help message information for argnorm
    """
    print (colored("\n\n***** argnorm help message *****\n\n", 'yellow'))
    print("""argnorm is a package to standarize/normalize Antibiotic resistance gene (ARG) by mapping to the antibiotic resistance ontology (ARO) by CARD""")

###############
def call_argnorm(tsv_file, out_file, softname="amrfinderplus", Debug=False):
    
    argnorm_v = set_config.check_package_version("argnorm", Debug=Debug)
    
    ## If it is available, argnorm is installed in path
    if (argnorm_v):
        argnorm_cmd = "argnorm %s -i %s -o %s" %(softname, tsv_file, out_file)
        code2return = HCGB_sys.system_call(argnorm_cmd)
        
        return(code2return)
    else:
        HCGB_aes.raise_and_exit("No argnorm package installed")



