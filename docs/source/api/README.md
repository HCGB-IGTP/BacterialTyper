# BacterialTyper API documentation

Documentation corresponding to python functions is included as docstrings within each module and submodule python script.

To generate the '.rst' file for each file module:

'''
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
'''

To generate the '.rst' file for each submodule file:

'''
cd BacterialTyper/docs/source/api/submodules
for i in `dir ../../../../BacterialTyper/`; do 
	name=(${i//.py/}); 
	file=$name".rst"; 
	echo ".. _"$name":" >> $file; 
	echo "" >> $file; 
	echo $name >> $file; 
	echo "========" >> $file; 
	echo ".. automodule:: BacterialTyper.submodules."$i >> $file; 
	echo "    :members:" >> $file; 
done;
'''
