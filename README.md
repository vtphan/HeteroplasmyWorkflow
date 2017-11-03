A workflow (designed to be run on high-performance clusters) for detecting heteroplasmies across multiple genomic samples.

REQUIREMENTS:
1. Python 3 packages:
- Biopython
- Bokeh
- Flexx

You can use Anaconda distribution (https://www.anaconda.com/download) for easier installation.
Install required packaged using Anaconda:
- Biopython : conda install -c anaconda biopython
- Bokeh : conda install -c bokeh bokeh
- Flexx : conda install -c bokeh flexx

2. Other tools:
- Bwa (http://bio-bwa.sourceforge.net/)
- Samtools (http://samtools.sourceforge.net/)


CONFIGURATION: 
You need to specify paths to your data in a configuration file. See config.txt for example.

config.txt:
- READS_DIR: path to reads directory.
- REF_DIR: path to reference genomes directory.
- ANNOTATION: path to Annotation file.
- LOG_FILE: path to log file.
- OUTPUT_DIR: path to output directory.
- DIST: name of distance function used to compute conservation scores of heteroplasmic sites (hellinger or consine).

It is not neccessary to use single or double quote for these paths.

Readids.txt:
Text file contains all input reads ID you want to run. Each line is reserved for only one ID.


HOW TO RUN:

python workflow.py config.txt readids.txt

- config.txt: configuration file
- readids.txt: read IDs file.


