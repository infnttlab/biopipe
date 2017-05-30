
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod ugo+x Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
source .bashrc 
pip install -U snakemake
conda config --add channels bioconda


