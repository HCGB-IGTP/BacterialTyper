#!/usr/bin/env bash

cd BacterialTyper/docs/source/api/modules
for i in `dir ../../../../BacterialTyper/modules`; do 
	name=(${i//.py/}); 
	file=$name".rst"; 
	echo ".. _"$name":" >> $file; 
	echo "" >> $file; 
	echo $name >> $file; 
	echo "========" >> $file; 
	echo ".. automodule:: BacterialTyper.modules."$i >> $file; 
	echo "    :members:" >> $file; 
done;