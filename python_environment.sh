#!/usr/bin/env bash
python3 -m venv BacterialGenotyping
source BacterialGenotyping/bin/activate
pip install --upgrade pip
pip install fastqcparser

pip install configparser
pip install numpy
pip install cython
pip install biopython
