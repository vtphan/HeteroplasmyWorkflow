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

if len(sys.argv) != 3:
    print('Usage: python', sys.argv[0], 'config_file.txt','read_file.txt')
    sys.exit(0)

#--------------------------------------------------------------
# read config file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open('defaults.ini'))
default_dist = config.get('defaults', 'DIST')
config.readfp(open(sys.argv[1]))

ref = config.get('config', 'REF')
annotation = config.get('config', 'ANNOTATION')
try:
    dist = config.get('config', 'DIST')
except:
    dist = default_dist

READS_DIR = config.get('config', 'READS_DIR')
OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')
LOG_FILE = config.get('config', 'LOG_FILE')
#--------------------------------------------------------------

with open(sys.argv[2], 'r') as f:
    reads = f.readlines()

n_reads = len(reads)

SCRIPT_DIR = os.getcwd()
print("HETEROPLASMY")
###########################################################
# 01_bwa
###########################################################
check_exist('which', 'bwa')
check_exist('which', 'samtools')

output = 'None'
if os.path.exists(ref + '.bwt'):
    print('Index exists. Skip indexing by bwa.')
else:
    print('Index', ref)
    cmd = 'bwa index %s' % ref
    with open(LOG_FILE, 'a') as f:
        f.write('%s\n%s\n' % (str(datetime.datetime.now()), cmd))

    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())

###########################################################
# 02_alignment
# 02_filter_by_samtools
###########################################################
# split read id to separated files
for i in range(0, n_reads):
    read_file = os.path.join(SCRIPT_DIR, 'readids'+str(i)+'.txt')
    with open(read_file,'w') as rf:
        rf.write(reads[i])

check_exist('ls', ref)
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

output = 'None'
no_error = True

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
    bf.write('python hpc_align.py '+sys.argv[1]+' readids${PBS_ARRAYID}.txt \n')

check_exist('ls', bash_file)
print('\nAlign reads...')
cmd = 'qsub '+bash_file

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

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
   
###########################################################
# 03_compute_heteroplasmy likelihood
# 04_sort_sites
###########################################################
heteroplasmy_likelihood = os.path.join(SCRIPT_DIR, 'heteroplasmy_likelihood.py')
sort_candidates = os.path.join(SCRIPT_DIR, 'sort_candidates.py')
check_exist('ls', heteroplasmy_likelihood)
check_exist('ls', sort_candidates)
check_exist('ls', annotation)
read_file = open(sys.argv[2])

csv_dir = os.path.join(OUTPUT_DIR, "csv")
if not os.path.exists(csv_dir):
    os.makedirs(csv_dir)

for line in read_file:
    read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
    read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
    name = read1.split('/')[-1].split('_R1')[0]
    out_csv = os.path.join(csv_dir, name+'_f2_q20.csv')
    out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_q20.sam')
    no_error = True

    output = 'None'

    if os.path.exists(out_csv):
        print("Scores might have been calculated already. Skip calculating.")
    else:
        print("Calculate heteroplasmy scores")
        cmd = 'python %s %s %s %s' % (heteroplasmy_likelihood,ref,out_filtered_sam,annotation)
        print(cmd)
        try:
            output = subprocess.check_call(cmd, shell=True, stdout=open(out_csv,'w'))
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    # 04_sort_sites
    print("Sort scores")
    cmd = 'python %s %s' % (sort_candidates,out_csv)
    print(cmd)
    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())

# print (finished_jobs)
print ('Finished computing heteroplasmy scores.\n')

###########################################################
# 05_select_sites
###########################################################
print('Select heteroplasmy sites.')
select_sites = os.path.join(SCRIPT_DIR, 'select_sites.py')
check_exist('ls', select_sites)
# run select_sites.py
result_dir = os.path.join(OUTPUT_DIR,"Result")
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

cp_heteroplasmy = os.path.join(result_dir,"cp_heteroplasmy.csv")
cmd = 'python %s %s > %s' %(select_sites, csv_dir, cp_heteroplasmy)
print(cmd)

output = 'None'
try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

###########################################################
# 06_compute_site_conservation
###########################################################
# run location_conservation.py
print('\nCompute site conservation.')
location_conservation = os.path.join(SCRIPT_DIR, 'location_conservation.py')
check_exist('ls', location_conservation)

cp_conserved = os.path.join(result_dir, "cp_conserved_"+dist+".csv")

cmd = 'python %s %s %s > %s' % (location_conservation, cp_heteroplasmy, dist, cp_conserved)
print(cmd)

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

###########################################################
# 07_plot
###########################################################
# run plot_heteroplasmy.py
print('\nPlot heteroplasmies.')
plot_heteroplasmy = os.path.join(SCRIPT_DIR, 'plot_heteroplasmy.py')
check_exist('ls',plot_heteroplasmy)

genome_name = '"Daucus carota chloroplast genome"'
out_html = os.path.join(OUTPUT_DIR,"cp.html")
cmd = 'python %s %s %s %s %s %s' %(plot_heteroplasmy, genome_name, annotation, cp_heteroplasmy, cp_conserved, out_html)
print(cmd)

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

print("\nSuccess!\n")
print("Vizualization file : ", out_html)