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
        with open(LOG_FILE, 'a') as f:
                f.write('time: %s\ncmd: %s\noutput: %s\nexec error:%s\n' % (str(datetime.datetime.now()), cmd, exec_output, exec_error))
                
def log_final(no_error, argv):
    log_output = os.path.join(SCRIPT_DIR, 'log_align_analyze_sort.txt')
    with open(log_output, 'a') as f:
        f.write('%s %s %s %s\n' % (no_error, argv[0], argv[1], str(datetime.datetime.now())))

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
cutoff = config.get('config', 'cutoff')
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
        print("\nOutput exists! Please change the OUTPUT_DIR in config file and re-run the program.")
        exit()
    else:
        ans = input("Remove all existing files in "+OUTPUT_DIR+"? (Y to remove, N to re-use these files)")
        if ans in ['y','Y','Yes','yes']:
            cmd = 'rm -rf '+OUTPUT_DIR
            try:
                output = subprocess.check_call(cmd, shell=True)
            except:
                no_error = False
                log_error(cmd, output, sys.exc_info())
            os.makedirs(OUTPUT_DIR)
            print("\nOverwrite OUTPUT_DIR.")
        if ans in ['n','N','No','no']:
            print("The workflow will re-use the existing files in "+OUTPUT_DIR+".")

start_time = time.time()
###########################################################
# 01_bwa
###########################################################
check_exist('which', 'bwa')
check_exist('which', 'samtools')

if os.path.exists(ref + '.bwt'):
    print('Index exists. Skip indexing by bwa.')
else:
    print('Index', ref)
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
for line in reads:
    read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
    read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
    check_exist('ls', read1)
    check_exist('ls', read2)

    name = read1.split('/')[-1].split('_R1')[0]
    out_sam = os.path.join(OUTPUT_DIR, name + '.sam')
    out_filtered_sam = os.path.join(OUTPUT_DIR, name + '_f2_F0x900_q' + alignment_quality + '.sam')

    # 02_alignment
    if os.path.exists(out_sam):
        print('Alignment might have been done already.  Skip bwa.')
    else:
        cmd = 'bwa mem %s %s %s' % (ref, read1, read2)
        try:
            output = subprocess.check_call(cmd, shell=True, stdout=open(out_sam, 'w'))
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    alignment_time = time.time()
    print("Alignment time for ", line.strip(), ": ", alignment_time - start_time)

    # 02_filter_by_samtools
    # only keep primary alignments
    if os.path.exists(out_filtered_sam):
        print('Alignment might have been filtered already.  Skip samtools.')
    else:
        print("Filter bwa's output")
        cmd = 'samtools view -f 2 -F 0x900 -q %s %s' % (alignment_quality, out_sam)
        try:
            ouptut = subprocess.check_call(cmd, shell=True, stdout=open(out_filtered_sam, 'w'))
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    print('Filter alignments for chloroplast and mitochondrial genomes.')
    cmd = 'python filter_samfiles_cp_mt.py %s %s %s %s' %(out_filtered_sam, OUTPUT_DIR, chloroplast, mitochondria)
    try:
        output = subprocess.check_output(cmd, shell=True)
        # output.wait()
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())            

    filter_time = time.time()
    print("Filter time for ", line.strip(), ": ", filter_time-alignment_time)


    print ("Finished %s. " %(line))


script = "run_local"
random_id = "0"
# ###########################################################
# run partial workflow for chloroplast
# ###########################################################
if chloroplast != 'None':
    partial_workflow = os.path.join(SCRIPT_DIR, 'run_hpc_het.py')
    check_exist('ls', partial_workflow)
    cp_out = os.path.join(OUTPUT_DIR,'chloroplast')

    # make temp parameters file
    param_file = os.join.path(OUTPUT_DIR, "temp_params.txt")
    params = [cp_ref, cp_annotation, dist, sys.argv[2], 'chloroplast'+output_info+'.html', str(random_id), READS_DIR, cp_out, LOG_FILE, alignment_quality, score_threshold, percentage_threshold, script, str(cutoff)]
    f = open(param_file,'w')
    for item in params:
        f.write(item+'\n')
    f.close()

    # cmd = 'python run_hpc_het.py %s' %(" ".join(params))
    cmd = 'python run_hpc_het.py temp_params.txt'
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
    params = [mt_ref, mt_annotation, dist, sys.argv[2], 'mitochondria'+output_info+'.html', str(random_id), READS_DIR, mt_out, LOG_FILE, alignment_quality, score_threshold, percentage_threshold, script, cutoff]
    
    # make temp parameters file
    param_file = os.join.path(OUTPUT_DIR, "temp_params.txt")
    params = [cp_ref, cp_annotation, dist, sys.argv[2], 'chloroplast'+output_info+'.html', str(random_id), READS_DIR, cp_out, LOG_FILE, alignment_quality, score_threshold, percentage_threshold, script, str(cutoff)]
    f = open(param_file,'w')
    for item in params:
        f.write(item+'\n')
    f.close()

    # cmd = 'python run_hpc_het.py %s' %(" ".join(params))
    cmd = 'python run_hpc_het.py temp_params.txt'
    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())
else:
    print("No sequence ID input for mitochondrial genome.")

