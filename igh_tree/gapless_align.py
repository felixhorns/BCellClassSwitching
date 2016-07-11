import sys
import os
import time
import json
from itertools import tee, izip_longest

def main(argv):

    start_time = time.time()

    infile = argv[1]
    infile_field_dict = argv[2]

    wdir = os.path.dirname(infile)
    wdir_original = os.path.dirname(infile_field_dict)

    with open(infile_field_dict, 'rb') as f:
        field_dict = json.load(f)

    alignments_dir = wdir_original + "/alignments/"
    alignments_dir_scratch = wdir + "/alignments/"

    for d in [alignments_dir, alignments_dir_scratch]:
        mkdir(d)

    align(infile, field_dict, wdir, alignments_dir_scratch)

    print "Done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

def align(infile, field_dict, wdir, alignments_dir):
    
    lineage = []

    with open(infile, 'rU') as f:
        for line, next_line in pairwise(f):

            if line == None:
                print "Empty input file"
                break

            lineage.append(line.rstrip())

            if next_line == None:
                # end of file                                     
                align_lineage(lineage, field_dict, wdir, alignments_dir)
                break

            this_lineage_uid = line.split("\t")[field_dict["lineage_uid"]]
            next_lineage_uid = next_line.split("\t")[field_dict["lineage_uid"]]

            if next_lineage_uid != this_lineage_uid:
                # end of lineage
                align_lineage(lineage, field_dict, wdir, alignments_dir)
                lineage = []

    return None

def align_lineage(lineage, field_dict, wdir, alignments_dir):
    """ Aligns sequences using a very simple gapless alignment approach.
    We assume that the boundaries of the D sequences are aligned (these are CDR3 start and end).
    Thus, we simply need to pad the V sequences and the J sequences to the same length,
    and concatenate them with the D sequence."""

    ### Do not process if lineage has only 1 sequence
    if len(lineage) == 1: return None

    ### Parse inputs
    lineage_uid = lineage[0].split("\t")[field_dict["lineage_uid"]]

    labels = []
    V_seqs = []
    D_seqs = []
    J_seqs = []

    for line in lineage:

        vals = line.split("\t")

        V_seq = vals[field_dict["V_seq"]]
        D_seq = vals[field_dict["D_seq"]]
        J_seq = vals[field_dict["J_seq"]]

        V_seqs.append(V_seq)
        D_seqs.append(D_seq)
        J_seqs.append(J_seq)

        label = "~".join(map(str, [vals[field_dict["sequence_uid"]], vals[field_dict["subisotype"]],
                                   vals[field_dict["abundance"]]]))
        labels.append(label)

    ### Find maximum V length and pad all V sequences to that length
    max_V_len = len(max(V_seqs, key=len))

    V_seqs_padded = []
    for s in V_seqs:
        s_padded = s.rjust(max_V_len, '-')
        V_seqs_padded.append(s_padded)

    ### Find maximum J length and pad all J sequences to that length
    max_J_len = len(max(J_seqs, key=len))

    J_seqs_padded = []
    for s in J_seqs:
        s_padded = s.ljust(max_J_len, '-')
        J_seqs_padded.append(s_padded)

    ### Make VDJ sequences
    VDJ_seqs = []

    for V_seq, D_seq, J_seq in zip(V_seqs_padded, D_seqs, J_seqs_padded):
        VDJ_seq = V_seq + D_seq + J_seq
        VDJ_seqs.append(VDJ_seq)

    ### Write sequences to output
    outfile = alignments_dir + "/" + str(lineage_uid) + ".afa.gapless"

    with open(outfile, 'w') as out:
        for label, s in zip(labels, VDJ_seqs):
            out.write(">" + label + "\n")
            out.write(s + "\n")

    return None

def mkdir(path):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST:
                pass
            else: raise
    return None

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip_longest(a, b)


if __name__ == "__main__":
    main(sys.argv)

