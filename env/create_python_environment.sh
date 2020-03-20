#!/usr/bin/env bash

export PYTHONPATH=""
python3 -m venv env/BacterialTyper
source env/BacterialTyper/bin/activate

## upgrade pip
pip install -I --isolated --upgrade pip