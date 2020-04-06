
cd ./BacterialTyper/docs/source/images/python_graph

for i in `dir ../../../../../BacterialTyper/scripts/`;
do
	echo "#########";
	echo "## Script: "$i;
	name=(${i//.py/});
	echo "##########";

	pyan --grouped -e --no-define --colored --dot ../../../../../BacterialTyper/scripts/$i > $name.dot; 
	dot -Tsvg -Gnewrank=true $name.dot > $name.svg;
done

