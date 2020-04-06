for i in `dir ../../../../BacterialTyper/scripts/`; do name=(${i//.py/}); file=$name".rst"; echo ".. _"$name":" >> $file; echo "" >> $file; echo $name >> $file; echo "==========================================" >> $file; echo "This script contains several functions. Here we show a graph representation of the different functions and relationships among them:" >> $file; echo "" >> $file; echo ".. image:: ../../images/python_graph/"$name".png" >> $file; echo "    :align: center" >> $file; echo "" >> $file; echo ".. automodule:: BacterialTyper.scripts."$name >> $file; echo "    :members:" >> $file; echo "    :undoc-members:" >> $file; echo "" >> $file; echo ".. include:: ../../links.inc" >> $file;  done



