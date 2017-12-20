wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod ugo+x Miniconda3-latest-Linux-x86_64.sh 
bash ./Miniconda3-latest-Linux-x86_64.sh -b
source ~/.bashrc
export PATH=~/miniconda3/bin:$PATH 
conda update conda -y
pip install -U psutil
pip install -U snakemake
pip install -U PyYAML
pip install matplotlib
pip install pandas
pip pysam
pip install ipython
conda config --add channels bioconda
rm ~/biopipe/Miniconda3-latest-Linux-x86_64.sh
chmod ugo+x script.sh
echo " " >> ~/.bashrc
echo "PATH=~/miniconda3/bin:$PATH"  >> ~/.bashrc

# for different architectures

#mkdir ~/data_ref
#mkdir ~/data_ref/normal_sample
#mkdir ~/data_ref/tumor_sample

#scp bio8:~/data_ref/normal_sample  ~/data_ref/normal_sample
#scp bio8:~/data_ref/tumor_sample ~/data_ref/tumor_sample