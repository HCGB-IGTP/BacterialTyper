#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez                                            ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain        ##
##############################################################
"""
Generates a phylogenetic reconstruction
"""
## useful imports
import time
import io
import os
import re
import sys
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config