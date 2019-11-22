## BacterialTyper: a bioinformatics pipeline for the integrative analysis of bacterial WGS 

We present BacterialTyper, a pipeline for the analysis of bacterial WGS data that integrates and facilitates the interpretation of results generated from multiple software analysis. It is capable of processing and identifying bacterial strains, identifying resistance and virulence genes, and generating data for outbreak analysis using WGS data from isolated microbial cultured colonies. The design of this bioinformatic tool allows comparing samples with an internal database (previously identified samples) and external databases.  

The pipeline is written in Python with a modular architecture and based on open-source software and databases engines. Multiple tasks are performed by each of several modules including: preparation of raw data; assembly and annotation; bacterial strain identification; mobile genetic elements identification (plasmids, putative pathogenicity islands or phage insertions regions); generation of a virulence and resistance profile; clustering based on sequence similarity; phylogenetic analysis; integration of metadata, etc. The tool allows to compare samples with previously identified samples (collected and internal database) but it also uses, and updates periodically, external databases from different sources. 


![Workflow](docs/source/images/workflow/all.jpg "BacterialTyper pipeline")
