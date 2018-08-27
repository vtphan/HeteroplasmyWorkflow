import subprocess
import os
import sys
import datetime
import time
import random
from configparser import ConfigParser
from datetime import datetime

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

#--------------------------------------------------------------

with open(sys.argv[2], 'r') as f:
    reads = f.readlines()

n_reads = len(reads)

SCRIPT_DIR = os.getcwd()
print("HETEROPLASMY")

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


random_id = random.randint(1,999999)
split_read_dir = os.path.join(OUTPUT_DIR,"out/")
if not os.path.exists(split_read_dir):
    os.makedirs(split_read_dir)

# split read id to separated files
for i in range(0, n_reads):
    read_file = os.path.join(split_read_dir, 'readids'+str(i)+'_'+str(random_id)+'.txt')
    with open(read_file,'w') as rf:
        rf.write(reads[i])

check_exist('ls', ref)

output = 'None'
no_error = True

# make bash file
job_name = 'HTPLASMY_JOB'+str(random_id)
fname = 'heteroplasmy_align'+str(random_id)+'.sh'
bash_file = os.path.join(SCRIPT_DIR, fname)
with open(bash_file, 'w') as bf:
    bf.write('#!/bin/sh \n')
    bf.write('#PBS -l nodes=1:default:ppn=1 \n')
    bf.write('#PBS -l walltime=72:00:00 \n')
    # bf.write('#PBS -A COMP')
    bf.write('#PBS -N '+job_name)                             
    bf.write('#PBS -t 0-'+str(n_reads-1)+' \n')
    bf.write('cd '+SCRIPT_DIR+' \n')
    bf.write('python 02_hpc_align.py '+sys.argv[1]+' '+split_read_dir+'readids${PBS_ARRAYID}_'+str(random_id)+'.txt \n')

check_exist('ls', bash_file)
print('\nAlign and filter reads...')
cmd = 'qsub '+bash_file


try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

check = True
while check:
    cmd = 'qstat | grep "JOB'+str(random_id)+'"'
    output = subprocess.check_output(cmd, shell=True)
    if output:
        parts = output.split()
        if b'C' in parts[-2]:
            check = False
    else:
        check = False

align_filter_time = time.time()
print("Alignment and Filtering time:", align_filter_time-index_time)

script = "run_hpc"
# ###########################################################
# run partial workflow for chloroplast
# ###########################################################
if chloroplast != 'None':
    partial_workflow = os.path.join(SCRIPT_DIR, 'run_hpc_het.py')
    check_exist('ls', partial_workflow)
    cp_out = os.path.join(OUTPUT_DIR,'chloroplast')
    params = [cp_ref, cp_annotation, dist, sys.argv[2], 'chloroplast'+output_info+'.html', str(random_id), READS_DIR, cp_out, LOG_FILE, alignment_quality, score_threshold, percentage_threshold, script]
    cmd = 'python run_hpc_het.py %s' %(" ".join(params))
    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())
else:
    print("No sequence ID input for chloroplast genome.")

# ###########################################################
# run partial workflow for mitochondria
# ###########################################################
if mitochondria != 'None':
    partial_workflow = os.path.join(SCRIPT_DIR, 'run_hpc_het.py')
    check_exist('ls', partial_workflow)
    mt_out = os.path.join(OUTPUT_DIR,'mitochondria')
    params = [mt_ref, mt_annotation, dist, sys.argv[2], 'mitochondria'+output_info+'.html', str(random_id), READS_DIR, mt_out, LOG_FILE, alignment_quality, score_threshold, percentage_threshold, script]
    cmd = 'python run_hpc_het.py %s' %(" ".join(params))
    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())
else:
    print("No sequence ID input for mitochondrial genome.")