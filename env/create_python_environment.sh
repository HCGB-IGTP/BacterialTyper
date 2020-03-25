#!/usr/bin/env bash

echo ""
echo ""
echo "##### Create a virtual environment for BacterialTyper ####"
echo ""

export PYTHONPATH=""
echo "python3 -m venv env/BacterialTyper"
#python3 -m venv env/BacterialTyper
echo "Done..."
echo ""

echo "#### Activate virtual environment ####"
echo "source env/BacterialTyper/bin/activate"
source env/BacterialTyper/bin/activate
echo "Done..."
echo ""

## upgrade pip
echo "#### Upgrade pip package installer ####"
echo "pip install -I --isolated --upgrade pip"
pip install -I --isolated --upgrade pip
echo "Done...."
echo ""

# install BacterialTyper requirements

echo "#### Install pip package requirements ####"
echo "pip install -I --isolated --default-timeout=500 -r BacterialTyper/config/python_requirements.txt"
pip install -I --isolated --default-timeout=500 -r BacterialTyper/config/python_requirements.txt
echo "Done...."
echo ""

echo "Exit...."
echo ""
