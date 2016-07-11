import sys
import os
import time
import uuid
import json
from itertools import tee, izip_longest
import subprocess
import random

def main(argv):

    start_time = time.time()

    infile = argv[1]
    infile_field_dict = argv[2]
    outfile_lineages = argv[3]

    path_to_muscle_bin = argv[4]
    path_to_V_germline_ref_file = argv[5]
    path_to_J_germline_ref_file = argv[6]

    wdir = os.path.dirname(infile)
    wdir_original = os.path.dirname(infile_field_dict)

    with open(infile_field_dict, 'rb') as f:
        field_dict = json.load(f)

    out_lineages = open(outfile_lineages, 'w')
    outfile_lineages_header = ["uid", "alignment"]
    write_header(out_lineages, outfile_lineages_header)

    alignments_dir = wdir_original + "/alignments/"
    alignments_dir_scratch = wdir + "/alignments/"

    for d in [alignments_dir, alignments_dir_scratch]:
        mkdir(d)

    align(infile, field_dict, wdir, alignments_dir_scratch, out_lineages, path_to_muscle_bin,
          path_to_V_germline_ref_file, path_to_J_germline_ref_file)

    out_lineages.close()

    print "Done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

def align(infile, field_dict, wdir, alignments_dir_scratch, out_lineages, path_to_muscle_bin,
          path_to_V_germline_ref_file, path_to_J_germline_ref_file):

    lineage = []

    with open(infile, 'rU') as f:
        for line, next_line in pairwise(f):

            if line == None:
                print "Empty input file"
                break

            lineage.append(line.rstrip())

            if next_line == None:
                # end of file                                     
                align_lineage(lineage, field_dict, wdir, alignments_dir_scratch, out_lineages, path_to_muscle_bin,
                              path_to_V_germline_ref_file, path_to_J_germline_ref_file)
                break

            this_lineage_uid = line.split("\t")[field_dict["lineage_uid"]]
            next_lineage_uid = next_line.split("\t")[field_dict["lineage_uid"]]

            if next_lineage_uid != this_lineage_uid:
                # end of lineage
                align_lineage(lineage, field_dict, wdir, alignments_dir_scratch, out_lineages, path_to_muscle_bin,
                              path_to_V_germline_ref_file, path_to_J_germline_ref_file)
                lineage = []

    return None

def align_lineage(lineage, field_dict, wdir, alignments_dir_scratch, out_lineages, path_to_muscle_bin,
                  path_to_V_germline_ref_file, path_to_J_germline_ref_file):

    ### Do not process if lineage has only 1 sequence
    if len(lineage) == 1: return None

    ### Subsample lineage if it is too large
    len_lineage_cutoff = 10000
    if len(lineage) > len_lineage_cutoff:
        lineage = subsample(lineage, N = len_lineage_cutoff)

    ### Parse inputs
    lineage_uid = lineage[0].split("\t")[field_dict["lineage_uid"]]

    V_germline = lineage[0].split("\t")[field_dict["V_germline"]]
    J_germline = lineage[0].split("\t")[field_dict["J_germline"]]

    labels = []
    V_seqs = []
    J_seqs = []
    VDJ_seqs = []

    for line in lineage:

        vals = line.split("\t")

        V_seq = vals[field_dict["V_seq"]]
        J_seq = vals[field_dict["J_seq"]]
        VDJ_seq = V_seq + vals[field_dict["D_seq"]] + J_seq

        V_seqs.append(V_seq)
        J_seqs.append(J_seq)
        VDJ_seqs.append(VDJ_seq)

        label = "~".join(map(str, [vals[field_dict["sequence_uid"]],
                                   clean_str(vals[field_dict["visits"]]),
                                   clean_str(vals[field_dict["abundances"]]),
                                   clean_str(vals[field_dict["isotypes"]])]))
        labels.append(label)

    ### Trim sequences to same length

    # min_length = len(min(VDJ_seqs, key=len))
    # VDJ_seqs_trimmed = [s[:min_length] for s in VDJ_seqs]

    VDJ_seqs_trimmed = VDJ_seqs # no trimming

    ### Prune duplicate sequences
    VDJ_seqs_trimmed_nodup = []
    labels_nodup = []

    seqs_seen = set()

    for seq, label in zip(VDJ_seqs_trimmed, labels):
        if seq not in seqs_seen:
            VDJ_seqs_trimmed_nodup.append(seq)
            labels_nodup.append(label)

    ### Write sequences and germline sequence to temporary fasta file for alignment
    temp_fasta_file = wdir + "/" + str(uuid.uuid4()) + ".fa.temp"
    write_to_fasta(VDJ_seqs_trimmed_nodup, labels_nodup, temp_fasta_file)

    ### Make multiple sequence alignment
    gapopen = -1.0 * float(len(VDJ_seqs_trimmed_nodup))
    temp_aln_file = alignments_dir_scratch + lineage_uid + ".afa.muscle.temp"
    muscle(temp_fasta_file, temp_aln_file, path_to_muscle_bin, gapopen=gapopen)

    ### Get germline sequence
    V_germline_seq = str(get_seq_by_id(path_to_V_germline_ref_file, V_germline)).upper()
    J_germline_seq = str(get_seq_by_id(path_to_J_germline_ref_file, J_germline)).upper()

    ### Truncate germlines to longest V and J sequences
    longest_V_seq = max(V_seqs, key=len)
    V_germline_seq_truncated = truncate_germline(V_germline_seq, longest_V_seq, wdir)

    longest_J_seq = max(J_seqs, key=len)
    J_germline_seq_truncated = truncate_germline(J_germline_seq, longest_J_seq, wdir)

    ### Make full germline sequence and write to temporary fasta file for alignment
    VJ_germline_seq = V_germline_seq_truncated + J_germline_seq_truncated
    write_to_fasta([VJ_germline_seq], ["germline"], temp_fasta_file)

    ### Profile-profile alignment to align germline to all sequences
    aln_file = alignments_dir_scratch + lineage_uid + ".afa.muscle"
    muscle_profile_profile(temp_aln_file, temp_fasta_file, aln_file, path_to_muscle_bin)

    ### Add gap to germline to make it the same length
#    length_gap = min_length - len(V_germline_seq_truncated) - len(J_germline_seq_truncated)
#    VJ_germline_seq = V_germline_seq_truncated + str('-' * length_gap) + J_germline_seq_truncated

    ### Set alignment flag in SQL database
    vals = [lineage_uid, 1]
    out_lineages.write("\t".join(map(str, vals)) + "\n")

    os.remove(temp_fasta_file)
    os.remove(temp_aln_file)

    return None

def muscle(infile_name, outfile_name, path_to_muscle_bin, gapopen=float(-1000.0)):
    cmd = path_to_muscle_bin + " -in " + infile_name + " -out " + outfile_name + " -quiet -maxiters 2 -diags -gapopen " + str(gapopen) + " -center 0.0 "
    subprocess.call(cmd, shell=True)
    return outfile_name

def muscle_profile_profile(infile1, infile2, outfile, path_to_muscle_bin):
    cmd = path_to_muscle_bin + " -profile -in1 " + infile1 + " -in2 " + infile2 + " -out " + outfile + " -quiet -maxiters 2 -diags"
    subprocess.call(cmd, shell=True)
    return outfile

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def blastn2seq(query, subject, wdir, **blast_params):
    query_record = SeqRecord(Seq(str(query)), id="query")
    subject_record = SeqRecord(Seq(str(subject)), id="subject")
    query_file_name = wdir + "/" + str(uuid.uuid4()) + ".fasta"
    subject_file_name = wdir + "/" + str(uuid.uuid4()) + ".fasta"
    out = wdir + "/" + str(uuid.uuid4()) + ".blastn_out.xml~"
    SeqIO.write(query_record, query_file_name, "fasta")
    SeqIO.write(subject_record, subject_file_name, "fasta")
    NcbiblastnCommandline(query=query_file_name, subject=subject_file_name, outfmt=5, out=out, **blast_params)()
    os.remove(query_file_name)
    os.remove(subject_file_name)
    return out

def truncate_germline(germline_seq, seq, wdir):

    germline_seq = str(germline_seq).upper()
    seq = str(seq).upper()

    blast_params = {'evalue': 1000,
                    'word_size': 7,
                    'gapopen': 10,
                    'gapextend': 2,
                    'penalty': -3,
                    'reward': 1}

    blast_output = blastn2seq(germline_seq, seq, wdir, **blast_params)
    blast_result_record = NCBIXML.read(open(blast_output, 'rU'))

    best_evalue = 1000
    best_hsp = None
    best_primer_i = -1

    evalue_warning_threshold = 0.001

    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:
            if float(hsp.expect) < best_evalue:
                best_evalue = float(hsp.expect)
                best_hsp = hsp

    if best_hsp == None:
        print "Error: no HSP"
        print
        quit()

    if best_evalue > evalue_warning_threshold:
        print "Warning: best E-value >", evalue_warning_threshold
        print "Best E value = ", str(best_evalue)
        print

    truncated_germline_seq = germline_seq[best_hsp.query_start-1:best_hsp.query_end]

    os.remove(blast_output)

    return truncated_germline_seq

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

def clean_str(s): return s.replace(",", "+")

def write_to_fasta(seqs, labels, outfile):
    with open(outfile, 'w') as out:
        for seq, label in zip(seqs, labels):
            out.write(">" + label + "\n")
            out.write(seq + "\n")
    return None

def write_header(out, h): out.write("\t".join(h)+"\n")

def get_seq_by_id(fasta_file_name, seq_id):
    seq = ''
    with open(fasta_file_name,'r') as f:
        for record in SeqIO.parse(f,'fasta'):
            if record.id == seq_id:
                seq = record.seq
    f.close()
    return seq

def subsample(lineage, N):

    lineage_subsampled = []

    p = float(N) / float(len(lineage)) # probability p to keep each sequence

    for line in lineage:
        if random.random() < p:
            lineage_subsampled.append(line)

    return lineage_subsampled

if __name__ == "__main__":
    main(sys.argv)
