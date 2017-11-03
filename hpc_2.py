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

if len(sys.argv) != 2:
    print('Usage: python', sys.argv[0], 'config_file.txt')
    sys.exit(0)

#--------------------------------------------------------------
# read config file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open(sys.argv[1]))

LOG_FILE = config.get('config', 'LOG_FILE')
OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')
annotation = config.get('config', 'ANNOTATION')
#--------------------------------------------------------------

SCRIPT_DIR = os.getcwd()
csv_dir = os.path.join(OUTPUT_DIR, "csv")
###########################################################
# 05_select_sites
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
# run location_conservation.py
location_conservation = os.path.join(SCRIPT_DIR, 'location_conservation.py')
check_exist('ls', location_conservation)

dist = config.get('config', 'DIST')
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
# run plot_heteroplasmy.py
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

print("Success!")
print("Output file : ", out_html)
