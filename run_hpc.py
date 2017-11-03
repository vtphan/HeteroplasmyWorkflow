import subprocess
import os
import sys
import datetime
from configparser import ConfigParser


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

#------------------------------------------------------------------------

if len(sys.argv) != 3:
    print('Usage: python', sys.argv[0], 'config_file.txt','read_file.txt')
    sys.exit(0)

#--------------------------------------------------------------
# read config file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open(sys.argv[1]))

OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')
ref = config.get('config', 'REF_DIR')
LOG_FILE = config.get('config', 'LOG_FILE')
#--------------------------------------------------------------

with open(sys.argv[2], 'r') as f:
    reads = f.readlines()

n_reads = len(reads)

SCRIPT_DIR = os.getcwd()

# split read id to separated files
for i in range(0, n_reads):
    read_file = os.path.join(SCRIPT_DIR, 'readids'+str(i)+'.txt')
    with open(read_file,'w') as rf:
        rf.write(reads[i])

check_exist('ls', ref)

output = 'None'
no_error = True

# build index
prog = os.path.join(SCRIPT_DIR, 'build_index.py')
cmd = 'python %s %s ' %(prog, sys.argv[1])
print(cmd)
try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())


# make bash file
bash_file = os.path.join(SCRIPT_DIR, 'heteroplamy_submit.sh')
with open(bash_file, 'w') as bf:
    bf.write('#!/bin/sh \n')
    bf.write('#PBS -l nodes=1:default:ppn=1 \n')
    bf.write('#PBS -l walltime=72:00:00 \n')
    # bf.write('#PBS -A COMP')
    bf.write('#PBS -N HTPLASMY_JOB')
    bf.write('#PBS -t 0-'+str(n_reads-1)+' \n')
    bf.write('cd '+SCRIPT_DIR+' \n')
    bf.write('python hpc_1.py '+sys.argv[1]+' readids${PBS_ARRAYID}.txt \n')

# run 1st phase
check_exist('ls', bash_file)
print('Processing reads...')
cmd = 'qsub '+bash_file

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())


if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

check = True
while check:
    cmd = 'qstat | grep "HTPLASMY_JOB"'
    output = subprocess.check_output(cmd, shell=True)
    if output:
        parts = output.split()
        if b'C' in parts[-2]:
            check = False
    else:
        check = False
   

# print (finished_jobs)
print ('Finished computing heteroplasmy scores.\n')


# move all csv files to csv directory
csv_dir = os.path.join(OUTPUT_DIR, "csv")
if not os.path.exists(csv_dir):
    os.makedirs(csv_dir)

cmd = 'mv '+OUTPUT_DIR+'/*.csv '+csv_dir

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

2nd phase
hpc2 = os.path.join(SCRIPT_DIR, 'hpc_2.py')
check_exist('ls', hpc2)

cmd = 'python %s %s' %(hpc2, sys.argv[1])

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

