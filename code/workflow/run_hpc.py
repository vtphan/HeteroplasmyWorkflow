import subprocess
import os
import sys
import datetime
import time
import random
from configparser import ConfigParser
from datetime import datetime
import s02_hpc_align, run_hpc_het
import multiprocessing

def check_exist(cmd, thing):
    try:
        subprocess.check_output('%s %s' % (cmd, thing), shell=True)
    except subprocess.CalledProcessError:
        print("Error: did not find %s in path." % thing)
        sys.exit(0)

def log_error(cmd, exec_output, exec_error):
    import datetime
    with open(LOG_FILE, 'a') as f:
        f.write('time: %s\ncmd: %s\noutput: %s\nexec error:%s\n' % (str(datetime.datetime.now()), cmd, exec_output, exec_error))
                
def log_final(no_error, argv):
    log_output = os.path.join(SCRIPT_DIR, 'log_align_analyze_sort.txt')
    with open(log_output, 'a') as f:
        f.write('%s %s %s %s\n' % (no_error, argv[0], argv[1], str(datetime.now())))

if len(sys.argv) != 3:
    print('Usage: python', sys.argv[0], 'config_file.txt','read_file.txt')
    sys.exit(0)

#--------------------------------------------------------------
# read defaults file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open('defaults.ini'))
default_dist = config.get('defaults', 'DIST')
default_score_threshold = config.get('defaults', 'score_threshold')
default_percentage_threshold = config.get('defaults', 'percentage_threshold')
default_alignment_quality = config.get('defaults', 'alignment_quality')
default_chloroplast = config.get('defaults', 'chloroplast')
default_mitochondria = config.get('defaults', 'mitochondria')

#--------------------------------------------------------------
# read version
#--------------------------------------------------------------
with open('VERSION', 'r') as f:
    line = f.readline()
    version = line.strip()

# get output_day for all output files
output_day = str(datetime.now()).split(" ")[0].replace(",","")

# make info for output filename
output_info = "_v"+version + "_" + output_day

#--------------------------------------------------------------
# read config file
#--------------------------------------------------------------
config.readfp(open(sys.argv[1]))
ref = config.get('config', 'REF')
READS_DIR = config.get('config', 'READS_DIR')
OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')
LOG_FILE = config.get('config', 'LOG_FILE')
cp_ref = config.get('config', 'cp_ref')
mt_ref = config.get('config', 'mt_ref')
cp_annotation = config.get('config', 'cp_annotation')
mt_annotation = config.get('config', 'mt_annotation')


try:
    chloroplast = config.get('config', 'chloroplast')
except:
    chloroplast = default_chloroplast

try:
    mitochondria = config.get('config', 'mitochondria')
except:
    mitochondria = default_mitochondria

if chloroplast == 'None' and mitochondria == 'None':
    print('No sequence ID input for chloroplast or mitochondrial genome.')
    exit()

# read optional parameters
try:
    dist = config.get('config', 'DIST')
except:
    dist = default_dist

try:
    alignment_quality = config.get('config', 'alignment_quality')
except:
    alignment_quality = default_alignment_quality

try:
    score_threshold = config.get('config', 'score_threshold')
except:
    score_threshold = default_score_threshold

try:
    percentage_threshold = config.get('config', 'percentage_threshold')
except:
    percentage_threshold = default_percentage_threshold

check_exist('ls', cp_ref)
check_exist('ls', mt_ref)
check_exist('ls', cp_annotation)
check_exist('ls', mt_annotation)

logf = open(LOG_FILE,'w')
logf.write("icHET\n")
logf.close()

#--------------------------------------------------------------
reads = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        if ',' in line:
            read = line.strip().split(",")[0]
        else:
            read = line.strip()
        reads.append(read)
        # reads = f.readlines()

n_reads = len(reads)

SCRIPT_DIR = os.getcwd()
print("icHET: Exploratory Visualization of Cytoplasmic Heteroplasmy")

output = 'None'

###########################################################
# check if OUTPUT_DIR exists
###########################################################
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
else:
    ans = input("\nOutput directory exists!!! Overwrite? (Y to continue, N to exit): ")
    if ans in ['n','N','No','no']:
        print("\nOutput exists! Please change the OUTPUT_DIR in config file and re-run the program.\n")
        exit()
    else:
        cmd = 'rm -rf '+OUTPUT_DIR
        try:
            output = subprocess.check_call(cmd, shell=True)
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())
        os.makedirs(OUTPUT_DIR)
        print("\nOverwrite OUTPUT_DIR.\n")

start_time = time.time()
###########################################################
# 01_bwa
###########################################################
check_exist('which', 'bwa')
check_exist('which', 'samtools')
# import datetime

if os.path.exists(ref + '.bwt'):
    print('\nIndex exists. Skip indexing by bwa.')
else:
    print('\nIndex', ref)
    cmd = 'bwa index %s' % ref
    with open(LOG_FILE, 'a') as f:
        f.write('%s\n%s\n' % (str(datetime.now()), cmd))

    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())

index_time = time.time()
print("Index time: ", index_time-start_time)

###########################################################
# 02_alignment
# 02_filter_by_samtools
###########################################################
if chloroplast != 'None':
    cp_out = os.path.join(OUTPUT_DIR, 'chloroplast')
    if not os.path.exists(cp_out):
        os.makedirs(cp_out)

if mitochondria != 'None':
    mt_out = os.path.join(OUTPUT_DIR, 'mitochondria')
    if not os.path.exists(mt_out):
        os.makedirs(mt_out)


print("Run hpc_align")
# hpc_align.main(sys.argv[1], sys.argv[2])

P = multiprocessing.Pool()
jobs = []
for r in reads:
    kw = {'config_file': sys.argv[1], 'read_ID': r}
    print(kw)
    jobs.append(P.apply_async(s02_hpc_align.process, (), kw))

P.close()
P.join()
   
align_filter_time = time.time()
print("Alignment and Filtering time:", align_filter_time-index_time)


# ###########################################################
# run partial workflow for chloroplast
# ###########################################################

if chloroplast != 'None':
    cp_out = os.path.join(OUTPUT_DIR,'chloroplast')
    params = {
        'ref': cp_ref,
        'annotation': cp_annotation,
        'dist': dist,
        'read_file': sys.argv[2],
        'out_html_name': 'chloroplast'+output_info+'.html',
        'random_id': "1",
        'read_dir': READS_DIR,
        'output_dir': cp_out,
        'log_file': LOG_FILE,
        'alignment_quality': alignment_quality,
        'score_threshold': score_threshold,
        'percentage_threshold': percentage_threshold
    }

    run_hpc_het.process(params)

else:
    print("No sequence ID input for chloroplast genome.")

# ###########################################################
# run partial workflow for mitochondria
# ###########################################################
if mitochondria != 'None':
    mt_out = os.path.join(OUTPUT_DIR,'mitochondria')

    params = {
        'ref': mt_ref,
        'annotation': mt_annotation,
        'dist': dist,
        'read_file': sys.argv[2],
        'out_html_name': 'mitochondria'+output_info+'.html',
        'random_id': "1",
        'read_dir': READS_DIR,
        'output_dir': mt_out,
        'log_file': LOG_FILE,
        'alignment_quality': alignment_quality,
        'score_threshold': score_threshold,
        'percentage_threshold': percentage_threshold
    }

    run_hpc_het.process(params)
else:
    print("No sequence ID input for mitochondrial genome.")

