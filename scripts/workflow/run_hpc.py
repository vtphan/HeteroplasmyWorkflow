import subprocess
import os
import sys
import datetime
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
# read config file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open('defaults.ini'))
default_dist = config.get('defaults', 'DIST')
default_score_threshold = config.get('defaults', 'score_threshold')
default_percentage_threshold = config.get('defaults', 'percentage_threshold')

# get version
with open('VERSION', 'r') as f:
    line = f.readline()
    version = line.strip()

# get output_day for all output files
output_day = str(datetime.now()).split(" ")[0].replace(",","")

# make info for output filename
output_info = "_v"+ver + "_" + output_day

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


try:
    score_threshold = config.get('config', 'score_threshold')
except:
    score_threshold = default_score_threshold


try:
    percentage_threshold = config.get('config', 'percentage_threshold')
except:
    percentage_threshold = default_percentage_threshold

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

###########################################################
# 01_bwa
###########################################################
check_exist('which', 'bwa')
check_exist('which', 'samtools')

if os.path.exists(ref + '.bwt'):
    print('\nIndex exists. Skip indexing by bwa.')
else:
    print('\nIndex', ref)
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
random_id = random.randint(1,999999)
split_read_dir = os.path.join(OUTPUT_DIR,"out/")
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
print('\nAlign reads...')
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
   
###########################################################
# 03_compute_heteroplasmy likelihood
# 04_sort_sites
###########################################################
heteroplasmy_likelihood = os.path.join(SCRIPT_DIR, '03_heteroplasmy_likelihood.py')
sort_candidates = os.path.join(SCRIPT_DIR, '04_sort_candidates.py')
check_exist('ls', heteroplasmy_likelihood)
check_exist('ls', sort_candidates)
check_exist('ls', annotation)
read_file = open(sys.argv[2])

csv_dir = os.path.join(OUTPUT_DIR, "csv")
if not os.path.exists(csv_dir):
    os.makedirs(csv_dir)

# make bash file
job_name = 'HTPLASMY_SCR'+str(random_id)
fname = 'heteroplasmy_score'+str(random_id)+'.sh'
bash_score = os.path.join(SCRIPT_DIR, fname)
with open(bash_score, 'w') as bf:
    bf.write('#!/bin/sh \n')
    bf.write('#PBS -l nodes=1:default:ppn=1 \n')
    bf.write('#PBS -l walltime=72:00:00 \n')
    # bf.write('#PBS -A COMP')
    bf.write('#PBS -N '+job_name+'\n')                             
    bf.write('cd '+SCRIPT_DIR+' \n')
    
    for line in read_file:
        read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
        read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
        name = read1.split('/')[-1].split('_R1')[0]
        out_csv = os.path.join(csv_dir, name+'_f2_q20.csv')
        out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_q20.sam')
        no_error = True

        output = 'None'

        # 03_compute heteroplasmy likelihood
        bf.write('echo "Calculate heteroplasmy scores"\n')
        cmd = 'python %s %s %s %s > %s' % (heteroplasmy_likelihood,ref,out_filtered_sam,annotation, out_csv)
        bf.write(cmd+'\n')
        
        # 04_sort_sites
        bf.write('echo "Sort scores"\n')
        cmd = 'python %s %s' % (sort_candidates,out_csv)
        bf.write(cmd+'\n')
     
check_exist('ls', bash_score)
print('\nCalculate and sort heteroplasmy scores...')
cmd = 'qsub '+bash_score

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

check = True
while check:
    cmd = 'qstat | grep "SCR'+str(random_id)+'"'
    output = subprocess.check_output(cmd, shell=True)
    if output:
        parts = output.split()
        if b'C' in parts[-2]:
            check = False
    else:
        check = False

# print (finished_jobs)
print ('Finished computing heteroplasmy scores.\n')

###########################################################
# clean up output files
###########################################################
hpc_out_dir = os.path.join(OUTPUT_DIR,"hpc_out/")
os.makedirs(hpc_out_dir)
cmd = 'mv '+SCRIPT_DIR+'/HTPLASMY_*'+str(random_id)+'* '+hpc_out_dir

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

bash_dir = os.path.join(OUTPUT_DIR,"bash_out/")
os.makedirs(bash_dir)

cmd = 'mv '+SCRIPT_DIR+'/*'+str(random_id)+'.sh '+bash_dir
try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())


###########################################################
# 05_select_sites
###########################################################
print('Select heteroplasmy sites.')
select_sites = os.path.join(SCRIPT_DIR, '05_select_sites.py')
check_exist('ls', select_sites)
# run select_sites.py
result_dir = os.path.join(OUTPUT_DIR,"Result")
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

cp_het_filename = "chloroplast_heteroplasmy"+output_info+".csv"
cp_heteroplasmy = os.path.join(result_dir,cp_het_filename)
cmd = 'python %s %s %s %s > %s' %(select_sites, csv_dir, score_threshold, percentage_threshold, cp_heteroplasmy)
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
location_conservation = os.path.join(SCRIPT_DIR, '06_location_conservation.py')
check_exist('ls', location_conservation)

cp_conserved_filename = "chloroplast_conserved_"+dist+output_info+".csv"
cp_conserved = os.path.join(result_dir, cp_conserved_filename)

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
plot_heteroplasmy = os.path.join(SCRIPT_DIR, '07_plot_heteroplasmy.py')
check_exist('ls',plot_heteroplasmy)

genome_name = '"Daucus carota chloroplast genome"'

outfilename = "chloroplast"+output_info+".html"

out_html = os.path.join(OUTPUT_DIR, outfilename)
cmd = 'python %s %s %s %s %s %s' %(plot_heteroplasmy, genome_name, annotation, cp_heteroplasmy, cp_conserved, out_html)
print(cmd)

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

print("\nSuccess!\n")
print("Vizualization file : ", out_html)

