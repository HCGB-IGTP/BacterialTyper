#!/usr/bin/env bash

cd BacterialTyper/docs/source/api/scripts
for i in `dir ../../../../BacterialTyper/`; do 
	name=(${i//.py/}); 
	file=$name".rst"; 
	echo ".. _"$name":" >> $file; 
	echo "" >> $file; 
	echo $name >> $file; 
	echo "========" >> $file; 
	echo "This script contains several functions. Here we show a graph" >> $file; 
	echo "representation of the different functions and relationships among them:" >> $file;

	echo ".. image:: ../../images/scripts_graph/"$name".png" >> $file;
		echo ":align: center" >> $file;

	echo ".. automodule:: BacterialTyper.scripts."$i >> $file; 
	echo "    :members:" >> $file; 
done;
