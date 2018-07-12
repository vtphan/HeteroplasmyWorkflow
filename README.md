# icHET: Exploratory Visualization of Cytoplasmic Heteroplasmy

## Heteroplasmy Workflow

A workflow (designed to be run on high-performance clusters) for detecting heteroplasmies across multiple genomic samples.

#### REQUIREMENTS:
1. Python 3 packages:
	- Biopython
	- Bokeh
	- Flexx

2. Other tools:
	- Bwa (http://bio-bwa.sourceforge.net/)
	- Samtools (http://samtools.sourceforge.net/)

#### INSTALL REQUIRED PACKAGES VIA ANACONDA:

You can use Anaconda distribution for easier installation.

1. Install Anaconda:

	- Download the appropriate .sh file from https://www.anaconda.com/download/
	- In the directory with the .sh file, run the .sh file using the following commands:
		- Make executable if needed:  ```chmod 755 SampleFileName.sh```
		- Run installer script: ```./SAMPLEFILENAME.sh```

2. Install required packaged:
	- Biopython : ```conda install -c anaconda biopython```
	- Bokeh : ```conda install -c bokeh bokeh```
	- Flexx : ```conda install -c bokeh flexx```

3. Install BWA, Samtools:
	- Set the proper channels to install bioconda:

	```
	conda config --add channels defaults
	conda config --add channels conda-forge
	conda config --add channels bioconda
	```

	- Install BWA, Samtools and Bzip2
	```
	conda install bwa
	conda install samtools
	conda install bzip2
	```

#### CONFIGURATION: 
You need to specify paths to your data in a configuration file. See config.txt for example.

##### config.txt:

1. Required inputs:

- READS_DIR: path to reads directory.
- REF: path to reference genomes. This is the concatenated of all genomes (nuclear DNA, mitochondrial genome, chloroplast genome).
- LOG_FILE: path to log file.
- OUTPUT_DIR: path to output directory.
- cp_ref : path to chloroplast genome.
- cp_annotation : path to chloroplast annotation file.
- mt_ref: path to mitochondria genome.
- mt_annotation: path to mitochondri annotation file.
- mitochondria: mitochondria sequence IDs. This can be a list, separated by commas.
- chloroplast: chloroplast sequence IDs. This can be a list, separated by commas.

It is not neccessary to use single or double quote for these paths. The workflow will generate the OUTPUT_DIR if it doesn't exists.

If there are no input sequence IDs for mitochondria or chloroplast, the program terminates. 

2. Optional inputs:
- DIST: name of distance function used to compute conservation scores of heteroplasmic sites (hellinger or consine). Default = hellinger distance.
- alignment_quality: quality threshold for SAMtools to filter alignments. Default = 20.
- score_threshold: threshold for conservation scores of heteroplasmic sites to be shown in visualization. Default = 10.
- percentage_threshold: threshold for base percentage of heteroplasmic sites to be shown in visualization. Default = 0.05


##### Readids.txt:
Text file contains all input reads ID you want to run. Each line is reserved for only one ID. The output plots the samples by the ordering of read names in this file.


#### HOW TO RUN:

```python run_hpc.py config.txt readids.txt```

- config.txt: configuration file
- readids.txt: read IDs file.

#### VISUALIZATION:

The program will output both visualization for mitochondria and chloroplast if users gives paths to chloroplast and mitochondrial genomes, annotation files, as well as sequence IDs.

Outputs for mitochondria and chloroplast will be separated into OUTPUT_DIR/mitochondria and OUTPUT_DIR/chloroplast directories.


