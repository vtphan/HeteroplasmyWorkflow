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
                f.write('time: %s\ncmd: %s\noutput: %s\nexec error:%s\n' % (str(datetime.now()), cmd, exec_output, exec_error))
                
def log_final(no_error, argv):
    log_output = os.path.join(SCRIPT_DIR, 'log_align_analyze_sort.txt')
    with open(log_output, 'a') as f:
        f.write('%s %s %s %s\n' % (no_error, argv[0], argv[1], str(datetime.now())))

if len(sys.argv) != 13:
    print('Usage: python', sys.argv[0], 'ref', 'annotation', 'dist', 'read_file', 'output.html', 'random_id', 'READS_DIR', 'output_dir', 'log_file', 'alignment_quality', 'score_threshold', 'percentage_threshold')
    sys.exit(0)

ref = sys.argv[1]
annotation = sys.argv[2]
dist = sys.argv[3]
read_file = sys.argv[4]
out_html_name = sys.argv[5]
random_id = sys.argv[6]
READS_DIR = sys.argv[7]
OUTPUT_DIR = sys.argv[8]
LOG_FILE = sys.argv[9]
alignment_quality = sys.argv[10]
score_threshold = sys.argv[11]
percentage_threshold = sys.argv[12]


# read version
with open('VERSION','r') as f:
    line = f.readline()
    version = float(line.strip())

#--------------------------------------------------------------
SCRIPT_DIR = os.getcwd()
print("\nComputing scores")
# print("Version: "+str(version))

output = 'None'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

   
###########################################################
# 03_compute_heteroplasmy likelihood
# 04_sort_sites
###########################################################
heteroplasmy_likelihood = os.path.join(SCRIPT_DIR, '03_heteroplasmy_likelihood.py')
sort_candidates = os.path.join(SCRIPT_DIR, '04_sort_candidates.py')
check_exist('ls', heteroplasmy_likelihood)
check_exist('ls', sort_candidates)
check_exist('ls', annotation)

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
    
    with open(read_file,'r') as f:
        for line in f:
            read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
            read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
            name = read1.split('/')[-1].split('_R1')[0]
            # out_csv = os.path.join(csv_dir, name+'_f2_q'+alignment_quality+'.csv')
            out_csv = os.path.join(csv_dir, name+'_f2_F0x900_q'+alignment_quality+'.csv')
            # out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_q'+alignment_quality+'.sam')
            out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_F0x900_q'+alignment_quality+'.sam')
            no_error = True

            output = 'None'

            bf.write('echo "Calculate heteroplasmy scores "'+name+'\n')
            cmd = 'python %s %s %s %s > %s' % (heteroplasmy_likelihood,ref,out_filtered_sam,annotation, out_csv)
            bf.write(cmd+'\n')
            
            # 04_sort_sites
            bf.write('echo "Sort scores "'+name+'\n')
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
if not os.path.exists(hpc_out_dir):
    os.makedirs(hpc_out_dir)
cmd = 'mv '+SCRIPT_DIR+'/HTPLASMY_*'+str(random_id)+'* '+hpc_out_dir

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

bash_dir = os.path.join(OUTPUT_DIR,"bash_out/")
if not os.path.exists(bash_dir):
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

if 'chloroplast' in out_html_name:
    cp_heteroplasmy = os.path.join(result_dir,"chloroplast_heteroplasmy.csv")
if 'mitochondria' in out_html_name:
    cp_heteroplasmy = os.path.join(result_dir,"mitochondria_heteroplasmy.csv")

cmd = 'python %s %s %s %s %s > %s' %(select_sites, csv_dir, score_threshold, percentage_threshold, read_file, cp_heteroplasmy)
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

if 'chloroplast' in out_html_name:
    cp_conserved = os.path.join(result_dir, "chloroplast_conserved_"+dist+".csv")
if 'mitochondria' in out_html_name:
    cp_conserved = os.path.join(result_dir, "mitochondria_conserved_"+dist+".csv")

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

# genome_name = '"Daucus carota chloroplast genome"'
if 'chloroplast' in out_html_name:
    genome_name = '"Daucus carota chloroplast genome"'
if 'mitochondria' in out_html_name:
    genome_name = '"Daucus carota mitochondrial genome"'

out_html = os.path.join(OUTPUT_DIR, out_html_name)
cmd = 'python %s %s %s %s %s %s' %(plot_heteroplasmy, genome_name, annotation, cp_heteroplasmy, cp_conserved, out_html)
print(cmd)

try:
    output = subprocess.check_call(cmd, shell=True)
except:
    no_error = False
    log_error(cmd, output, sys.exc_info())

print("\nSuccess!\n")
print("Vizualization file : ", out_html)

