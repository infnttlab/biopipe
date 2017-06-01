
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod ugo+x Miniconda3-latest-Linux-x86_64.sh 
bash ./Miniconda3-latest-Linux-x86_64.sh -b
source ~/.bashrc 
PATH=~/miniconda3/bin:$PATH
export PATH
conda update conda -y
pip install -U psutil
pip install -U snakemake
pip install -U PyYAML
conda config --add channels bioconda


