echo "#### Check BacterialTyper software dependencies ####"
echo "## Read dependencies from file: BacterialTyper/config/dependencies.csv"
echo "+ Get software name"
echo "+ Get minimun version required"
echo "+ Get regular expression to get version for each software"

##cat "BacterialTyper/config/dependencies.csv"

IFS=','
while read soft version_cmd get_version min_version soft_name
do
	echo "#################################"
	echo "## Check software: " $soft
	echo "#################################"
	echo "soft_name:" $soft_name
	echo "min_version:" $min_version
	echo ""
	

done < BacterialTyper/config/dependencies.csv



