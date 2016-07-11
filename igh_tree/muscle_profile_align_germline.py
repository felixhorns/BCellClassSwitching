import sys
import os
import time
import uuid
import json
from itertools import tee, izip_longest
import subprocess

def main(argv):

    start_time = time.time()

    infile = argv[1]
    infile_field_dict = argv[2]
    outfile_lineages = argv[3]
    path_to_muscle_bin = argv[4]

    wdir = os.path.dirname(infile)
    wdir_original = os.path.dirname(infile_field_dict)

    with open(infile_field_dict, 'rb') as f:
        field_dict = json.load(f)

    alignments_dir = wdir_original + "/alignments/"
    alignments_dir_scratch = wdir + "/alignments/"

    for d in [alignments_dir, alignments_dir_scratch]:
        mkdir(d)

    outfile_lineages_header = ["uid", "alignment"]
    out_lineages = open(outfile_lineages, 'w')
    write_header(out_lineages, outfile_lineages_header)

    muscle(infile, field_dict, wdir, alignments_dir, path_to_muscle_bin, out_lineages)

    out_lineages.close()

    print "Done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

def muscle(infile, field_dict, wdir, alignments_dir, path_to_muscle_bin, out_lineages):

    with open(infile, 'rU') as f:
        for line in f:

            if line == None:
                print "Empty input file"
                break

            lineage_uid = line.rstrip()
            muscle_lineage(lineage_uid, field_dict, wdir, alignments_dir, path_to_muscle_bin, out_lineages)

    return None

def muscle_lineage(lineage_uid, field_dict, wdir, alignments_dir, path_to_muscle_bin, out_lineages):

    infile_seqs = alignments_dir + "/" + str(lineage_uid) + ".afa.gapless"
    infile_germline = alignments_dir + "/" + str(lineage_uid) + ".germline_VJ.fa"

    outfile = alignments_dir + "/" + str(lineage_uid) + ".afa.gapless.germline_VJ"

    result = muscle_profile_profile(infile_seqs, infile_germline, outfile, path_to_muscle_bin)

    if result == 0:
        vals = [lineage_uid, 1]
        out_lineages.write("\t".join(map(str, vals)) + "\n")
    else:
        print "Error: muscle returned exit status 0"
        print "lineage uid = ", lineage_uid

    return None

def muscle_profile_profile(infile1, infile2, outfile, path_to_muscle_bin):
    cmd = path_to_muscle_bin + " -profile -in1 " + infile1 + " -in2 " + infile2 + " -out " + outfile + " -quiet -maxiters 2 -diags"
    result = subprocess.call(cmd, shell=True)
    return result

def mkdir(path):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST:
                pass
            else: raise
    return None

def write_header(out, h): out.write("\t".join(h)+"\n")

if __name__ == "__main__":
    main(sys.argv)
