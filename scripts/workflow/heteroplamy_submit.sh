#!/bin/sh 
#PBS -l nodes=1:default:ppn=1 
#PBS -l walltime=72:00:00 
#PBS -N HTPLASMY_JOB#PBS -t 0-1 
cd /home/dpham2/carrot/scripts 
python hpc_align.py ../config.txt readids${PBS_ARRAYID}.txt 
