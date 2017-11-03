import os
import subprocess
import sys
import datetime
from configparser import ConfigParser

#--------------------------------------------------------------
# read config file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open(sys.argv[1]))

LOG_FILE = config.get('config', 'LOG_FILE')

ref = config.get('config', 'REF_DIR')

#--------------------------------------------------------------

output = 'None'
# 01_bwa
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
  
