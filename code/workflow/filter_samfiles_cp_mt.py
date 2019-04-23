import sys
import os


def filter_reads_for_cp_mt(filtered_sam, SAMdir, chloroplast, mitochondria):
    genome_name = {}
    if mitochondria != 'None':
        genome_name['mitochondria'] = mitochondria.split(",")
    if chloroplast != 'None':
        genome_name['chloroplast'] = chloroplast.split(",")

    if len(genome_name) > 0:
        # get alignment of mt and cp

        outfilename = filtered_sam.split("/")[-1]

        if chloroplast != 'None':
            cp_out = os.path.join(SAMdir, 'chloroplast')
            if not os.path.exists(cp_out):
                os.makedirs(cp_out)
            cp_outfile = os.path.join(cp_out, outfilename)
            cpf = open(cp_outfile, 'w')

        if mitochondria != 'None':
            mt_out = os.path.join(SAMdir, 'mitochondria')
            if not os.path.exists(mt_out):
                os.makedirs(mt_out)
            mt_outfile = os.path.join(mt_out, outfilename)
            mtf = open(mt_outfile, 'w')

        with open(filtered_sam, "rU") as file:
            for line in file:
                if line[0] != '@':
                    items = line.split('\t')
                    ref_name = items[2]
                    sam_flag, ref_pos, cigar, read, phred = int(items[1]), int(items[3]) - 1, items[5], items[9], items[10]

                    if 'S' not in cigar:
                        name_parts = ref_name.split("|")
                        for n in name_parts:
                            if chloroplast != 'None':
                                if n in genome_name['chloroplast']:
                                    cpf.write(line)
                                    break

                            if mitochondria != 'None':
                                if n in genome_name['mitochondria']:
                                    mtf.write(line)
                                    break

        if chloroplast != 'None':
            cpf.close()
        if mitochondria != 'None':
            mtf.close()


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: python', sys.argv[0], 'filtered_sam_file', 'SAM_files_directory', 'chloroplast_ID', 'mitochondria_ID')
        sys.exit(0)

    filtered_sam = sys.argv[1]  # filtered SAM file
    SAMdir = sys.argv[2]  # dir contains all sam files
    chloroplast = sys.argv[3]
    mitochondria = sys.argv[4]

    filter_reads_for_cp_mt(filtered_sam, SAMdir, chloroplast, mitochondria)
