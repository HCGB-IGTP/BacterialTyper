#!/usr/bin/env bash

## create environment
export PYTHONPATH=""
python3 -m venv env/BacterialTyper

## activate env
source env/BacterialTyper/bin/activate

## upgrade pip
pip install -I --isolated --upgrade pip

## install python modules dependencies
pip install -I --isolated -r ./BacterialTyper/config/python_requirements.txt

export PYTHONPATH=$PYTHONPATH":"$PWD"/BacterialTyper"

export PATH=$PATH":"$PWD"/BacterialTyper.py"

## install additional perl and software dependencies
python BacterialTyper/config/install_dependencies.py