This folder stores all the scripts of the workflow.

#### Run workflow on a cluster
Use this command to run the workflow on the cluster:

`run_hpc.py config.txt readids.txt`

where:
- config.txt is a configuration file. See examples/config.txt
- readids.txt is a file which contains all sample IDs. See examples/readids.txt

#### Run workflow on a single server
Use this command to run the workflow on a single server:

`run_local.py config.txt readids.txt` 
