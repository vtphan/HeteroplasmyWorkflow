# install dependencies
conda install -y -c anaconda biopython
conda install -y -c bokeh bokeh
conda install -y -c bokeh flexx

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install -y bwa
conda install -y samtools
conda install -y bzip2
