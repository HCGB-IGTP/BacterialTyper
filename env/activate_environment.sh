#!/usr/bin/env bash

echo ""
echo ""
echo "##### Activate the virtual environment for BacterialTyper ####"
echo ""

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "source env/BacterialTyper/bin/activate"
source $DIR/BacterialTyper/bin/activate
echo "Done..."
echo ""

echo "Exit...."
echo ""
