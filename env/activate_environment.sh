#!/usr/bin/env bash

echo ""
echo ""
echo "##### Activate the virtual environment for BacterialTyper ####"
echo ""

echo "source env/BacterialTyper/bin/activate"
source env/BacterialTyper/bin/activate
echo "Done..."
echo ""

## export BacterialTyper
export PYTHONPATH="$PYTHONPATH:$PWD/BacterialTyper"

echo "Exit...."
echo ""
