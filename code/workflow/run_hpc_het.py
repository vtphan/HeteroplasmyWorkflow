import subprocess
import os
import sys
import datetime
import random
from configparser import ConfigParser
from datetime import datetime
import s03_heteroplasmy_likelihood, s04_sort_candidates, s05_select_sites, s06_location_conservation
import multiprocessing

def check_exist(cmd, thing):
    try:
        subprocess.check_output('%s %s' % (cmd, thing), shell=True)
    except subprocess.CalledProcessError:
        print("Error: did not find %s in path." % thing)
        sys.exit(0)

def log_error(cmd, exec_output, exec_error, LOG_FILE):
        with open(LOG_FILE, 'a') as f:
                f.write('time: %s\ncmd: %s\noutput: %s\nexec error:%s\n' % (str(datetime.now()), cmd, exec_output, exec_error))
                
def log_final(no_error, argv):
    log_output = os.path.join(SCRIPT_DIR, 'log_align_analyze_sort.txt')
    with open(log_output, 'a') as f:
        f.write('%s %s %s %s\n' % (no_error, argv[0], argv[1], str(datetime.now())))

def process(params):
    ref = params['ref']
    annotation = params['annotation']
    dist = params['dist']
    read_file = params['read_file']
    out_html_name = params['out_html_name']
    random_id = params['random_id']
    READS_DIR = params['read_dir']
    OUTPUT_DIR = params['output_dir']
    LOG_FILE = params['log_file']
    alignment_quality = params['alignment_quality']
    score_threshold = params['score_threshold']
    percentage_threshold = params['percentage_threshold']

    # print(ref)
    # print(annotation)
    # print(dist)
    # print(read_file)
    # print(READS_DIR)
    # print(OUTPUT_DIR)
    # print(LOG_FILE)
    # print(alignment_quality)
    # print(score_threshold)
    # print(percentage_threshold)

    # read version
    with open('VERSION','r') as f:
        line = f.readline()
        version = float(line.strip())

    # #--------------------------------------------------------------
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
    check_exist('ls', annotation)

    csv_dir = os.path.join(OUTPUT_DIR, "csv")
    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    print("Compute heteroplasmy likelihood")

    P = multiprocessing.Pool()
    jobs = []
    with open(read_file, 'r') as f:
        for line in f:
            read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
            read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
            name = read1.split('/')[-1].split('_R1')[0]
            # name = line.strip()
            out_csv = os.path.join(csv_dir, name+'_f2_F0x900_q'+alignment_quality+'.csv')
            out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_F0x900_q'+alignment_quality+'.sam')
            no_error = True
            output = 'None'

            kw = {
                'ref': ref,
                'out_filtered_sam': out_filtered_sam,
                'annotation': annotation,
                'out_csv': out_csv,
            }

            jobs.append(P.apply_async(s03_heteroplasmy_likelihood.process, (), kw))

    P.close()
    P.join()

    # Sort score
    P = multiprocessing.Pool()
    jobs = []
    with open(read_file, 'r') as f:
        for line in f:
            read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
            read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
            name = read1.split('/')[-1].split('_R1')[0]
            # name = line.strip()
            out_csv = os.path.join(csv_dir, name+'_f2_F0x900_q'+alignment_quality+'.csv')
            
            kw2 = {
                'out_csv': out_csv
            }
            jobs.append(P.apply_async(s04_sort_candidates.process, (), kw2))

    P.close()
    P.join()

    print ('Finished computing heteroplasmy scores.\n')

    ###########################################################
    # 05_select_sites
    ###########################################################
    print('Select heteroplasmy sites.')
    # run select_sites.py
    result_dir = os.path.join(OUTPUT_DIR,"Result")
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    organellar_type = None
    if 'chloroplast' in out_html_name:
        organellar_type = 'chloroplast'
    if 'mitochondria' in out_html_name:
        organellar_type = 'mitochondria'

    select_sites_inputs = {
        'csv_dir' : csv_dir,
        'score_threshold': score_threshold,
        'percentage_threshold': percentage_threshold,
        'name_list' : None,
        'organellar_type': organellar_type,
        'result_dir': result_dir
    }

    het_file = s05_select_sites.process(select_sites_inputs)

    ###########################################################
    # 06_compute_site_conservation
    ###########################################################
    # run location_conservation.py
    print('\nCompute site conservation.')

    cp_conserved = None
    if organellar_type == 'chloroplast':
        cp_conserved = os.path.join(result_dir, "chloroplast_conserved_"+dist+".csv")
    if organellar_type == 'mitochondria':
        cp_conserved = os.path.join(result_dir, "mitochondria_conserved_"+dist+".csv")

    location_conservation_inputs = {
        'het_file': het_file,
        'func': dist,
        'output': cp_conserved
    }

    s06_location_conservation.main(location_conservation_inputs)

    ###########################################################
    # 07_plot
    ###########################################################
    # run plot_heteroplasmy.py
    print('\nPlot heteroplasmies.')
    plot_heteroplasmy = os.path.join(SCRIPT_DIR, 's07_plot_heteroplasmy.py')
    check_exist('ls',plot_heteroplasmy)

    # genome_name = '"Daucus carota chloroplast genome"'
    if organellar_type == 'chloroplast':
        genome_name = '"Daucus carota chloroplast genome"'
    if organellar_type == 'mitochondria':
        genome_name = '"Daucus carota mitochondrial genome"'

    out_html = os.path.join(OUTPUT_DIR, out_html_name)
    cmd = 'python %s %s %s %s %s %s' %(plot_heteroplasmy, genome_name, annotation, het_file, cp_conserved, out_html)
    print(cmd)
    print()

    try:
        output = subprocess.check_call(cmd, shell=True)
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info(), LOG_FILE)

    print("\nSuccess!\n")
    print("Vizualization file : ", out_html)

if __name__ == '__main__':
    if len(sys.argv) != 13:
        print('Usage: python', sys.argv[0], 'ref', 'annotation', 'dist', 'read_file', 'output.html', 'random_id', 'READS_DIR', 'output_dir', 'log_file', 'alignment_quality', 'score_threshold', 'percentage_threshold')
        sys.exit(0)

    params = {
        'ref': sys.argv[1],
        'annotation': sys.argv[2],
        'dist': sys.argv[3],
        'read_file': sys.argv[4],
        'out_html_name': sys.argv[5],
        'random_id': sys.argv[6],
        'READS_DIR': sys.argv[7],
        'OUTPUT_DIR': sys.argv[8],
        'LOG_FILE': sys.argv[9],
        'alignment_quality': sys.argv[10],
        'score_threshold': sys.argv[11],
        'percentage_threshold': sys.argv[12],
    }

    process(params)
