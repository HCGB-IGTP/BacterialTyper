#!/usr/bin/env bash

export PYTHONPATH=""
python3 -m venv config/BacterialGenotyping_env
source config/BacterialGenotyping_env/bin/activate

## upgrade pip
pip install -I --isolated --upgrade pip

## install from requirements.txt file
pip install -r ./config/main/python_requirements.txt
