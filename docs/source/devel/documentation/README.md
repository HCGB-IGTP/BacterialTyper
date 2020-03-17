# BacterialTyper API documentation

Documentation corresponding to python functions is included as docstrings within each module and submodule python scripts. Only, an rst file referring to each python file is necessary to access information. 

I would use shell to generate a loop to create the '.rst' file for each module and script:

modules:

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

scripts:

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
