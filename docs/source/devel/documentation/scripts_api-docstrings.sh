#!/usr/bin/env bash

cd BacterialTyper/docs/source/api/scripts
for i in `dir ../../../../BacterialTyper/`; do 
	name=(${i//.py/}); 
	file=$name".rst"; 
	echo ".. _"$name":" >> $file; 
	echo "" >> $file; 
	echo $name >> $file; 
	echo "========" >> $file; 
	echo ".. automodule:: BacterialTyper.scripts."$i >> $file; 
	echo "    :members:" >> $file; 
done;
