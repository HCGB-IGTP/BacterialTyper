#!/usr/bin/env bash

export PYTHONPATH=""
python3 -m venv config/BacterialGenotyping_env
source config/BacterialGenotyping_env/bin/activate

## upgrade pip
pip install -I --isolated --upgrade pip

## install from requirements.txt file
pip install -r ./config/python_requirements.txt


##################################################
## Modules
##################################################
## pip install -I --isolated fastqcparser
## pip install -I --isolated configparser
## pip install -I --isolated numpy
## pip install -I --isolated scipy
## pip install -I --isolated cython
## pip install -I --isolated biopython
## pip install -I --isolated pandas
## pip install -I --isolated xlrd
## pip install -I --isolated xlwt
## pip install -I --isolated wget
## pip install -I --isolated xtract
## pip install -I --isolated ncbi-genome-download
## pip install -I --isolated 'matplotlib==2.2.4'
## pip install -I --isolated multiqc
## pip install -I --isolated termcolor 
## pip install -I --isolated ariba  ## dependencies: biopython, pysam, pymummer, dendropy, pyfastaq
## pip install -I --isolated openpyxl
## pip install -I --isolated xlsxwriter
## 
## pip freeze > requirements.txt
##################################################

