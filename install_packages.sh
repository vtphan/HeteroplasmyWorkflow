# install dependencies

# Install biopython, bokeh, flexx
conda install -y -c anaconda biopython
conda install -y -c bokeh bokeh
# conda install -y -c bokeh flexx
conda install -y -c conda-forge flexx

# Set the proper channels to install 
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

# Install BWa, samtools and bzip2
conda install -y bwa
conda install -y samtools
conda install -y bzip2
