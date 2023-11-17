#conda create -y -n BT_dev; conda activate BT_dev
conda install -y -c conda-forge mamba
mamba install -y -c r -c bioconda -c r -c conda-forge perl-bio-searchio-hmmer  \
    perl-bioperl trimmomatic fastqc spades snippy prokka kma \
    agrvate staphopia-sccmec\
    kraken2 bowtie2 ncbi-amrfinderplus  \
    bracken krakentools mlst 
    

mamba install -y spyder-kernels=2.4

## install argnorm
## pip install ....githubzip

## busco=5.5.0; sepp=4.4.0
## missing busco, etc

## Install pip packages under development
echo "Install card-trick"
cd /imppc/labs/lslab/jsanchez/git_repo/HCGB_git/card_trick
sh ./devel/pypi/test_module.sh

echo "Install HCGB"
cd /imppc/labs/lslab/jsanchez/git_repo/HCGB_git/HCGB_python_functions
sh ./devel/pypi/test_module.sh

echo "Install BacterialTyper"
cd /imppc/labs/lslab/jsanchez/git_repo/HCGB_git/BacterialTyper 
sh ./devel/pypi/test_module.sh

