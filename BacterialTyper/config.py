EXECUTABLES = {
	## Executables
	'blastn': '/software/debian-8/bin/blastn',
	'makeblastdb': '/software/debian-8/bin/makeblastdb',
	'spades': '/imppc/labs/lslab/jsanchez/software/SPAdes-3.13.0-Linux/bin/spades.py',
	'fastqc': '/software/debian-8/bin/fastqc',
	'busco_bin': '',
	'Phispy_folder': '/imppc/labs/lslab/jsanchez/git_repo/PhiSpy',
	'prokka': '/software/debian-8/bin/prokka',
	'trimmomatic':'/soft/bio/trimmomatic/trimmomatic.jar',
	'Rscript': '/usr/bin/Rscript',
	'java':'',
	'kma':'/imppc/labs/lslab/jsanchez/git_repo/KmerFinder/kma/kma',
	'plasmidID':'/imppc/labs/lslab/jsanchez/git_repo/plasmidID/plasmidID.sh'
}

PARAMETERS = {
	## Parameters
	'threads':2,
}

DATA = {
	## trimmomatic
	'trimmomatic_adapters':'/imppc/labs/lslab/share/data/references/Trimmomatic_adapters.fa',

	## database generated
	'database':'/home/labs/lslab/jsanchez/DATA/Saures_test/BacterialTyper_test/database',

	## BUSCO datasets
	'busco_bacteria': '/imppc/labs/lslab/jsanchez/software/BUSCO/datasets/bacteria_odb9',
	'busco_firmicutes': '/imppc/labs/lslab/jsanchez/software/BUSCO/datasets/firmicutes_odb9',

	## Virulence Factor Database
	#'VFDB_VFs':'',
	#'VFDB_Comparative_tables':'',
	#'VFDB_set_folder':'',

	## MLSTar
	#'MLSTar_profile_folder': '',
	#'MLSTar_sequence_folder': '',
}


