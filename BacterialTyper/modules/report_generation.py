#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Creates a report in HTML for each step of the process analyzed.
"""
## import useful modules
import os
import sys
import re
import time
from termcolor import colored

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config

## example report check: https://github.com/tseemann/nullarbor

## for staphyloccocus aureus: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0593-7
## for klebsiella: kleborate
## for micobacterium tuberculosis: MTBSeq

