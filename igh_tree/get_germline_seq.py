import sys
import os
import time
import json
import uuid
import errno
from itertools import tee, izip_longest
from Bio.Blast.Applications import NcbiblastnCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def main(argv):

    start_time = time.time()

    infile = argv[1]
    infile_field_dict = argv[2]
    path_to_V_germline_ref_file = argv[3]
    path_to_J_germline_ref_file = argv[4]

    wdir = os.path.dirname(infile)
    wdir_original = os.path.dirname(infile_field_dict)

    with open(infile_field_dict, 'rb') as f:
        field_dict = json.load(f)

    alignments_dir = wdir_original + "/alignments/"
    alignments_dir_scratch = wdir + "/alignments/"

    for d in [alignments_dir, alignments_dir_scratch]:
        mkdir(d)

    get_germline_seq(infile, field_dict, wdir, alignments_dir_scratch,
                     path_to_V_germline_ref_file, path_to_J_germline_ref_file)

    print "Done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

    return None

def get_germline_seq(infile, field_dict, wdir, alignments_dir,
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
                get_germline_seq_lineage(lineage, field_dict, wdir, alignments_dir,
                                         path_to_V_germline_ref_file, path_to_J_germline_ref_file)
                break

            this_lineage_uid = line.split("\t")[field_dict["lineage_uid"]]
            next_lineage_uid = next_line.split("\t")[field_dict["lineage_uid"]]

            if next_lineage_uid != this_lineage_uid:
                # end of lineage
                get_germline_seq_lineage(lineage, field_dict, wdir, alignments_dir,
                                         path_to_V_germline_ref_file, path_to_J_germline_ref_file)
                lineage = []

    return None

def get_germline_seq_lineage(lineage, field_dict, wdir, alignments_dir,
                             path_to_V_germline_ref_file, path_to_J_germline_ref_file):

    ### Do not process if lineage has only 1 sequence
    if len(lineage) == 1: return None

    ### Parse inputs
    lineage_uid = lineage[0].split("\t")[field_dict["lineage_uid"]]

    vals = lineage[0].split("\t")

    V_seq = vals[field_dict["V_seq"]]
    J_seq = vals[field_dict["J_seq"]]

    V_germline = vals[field_dict["V_germline"]]
    J_germline = vals[field_dict["J_germline"]]

    ### Get germline sequences
    V_germline_seq = str(get_seq_by_id(path_to_V_germline_ref_file, V_germline)).upper()
    J_germline_seq = str(get_seq_by_id(path_to_J_germline_ref_file, J_germline)).upper()

    ### Align V germline sequence to V sequence, truncate
    V_blast_output_file = wdir + "/" + str(uuid.uuid4()) + ".blastn_out.xml~"
    V_germline_seq_trimmed = truncate_germline(V_germline_seq, V_seq, V_blast_output_file, wdir, lineage_uid)

    ### Align J germline sequence to J sequence, truncate
    J_blast_output_file = wdir + "/" + str(uuid.uuid4()) + ".blastn_out.xml~"
    J_germline_seq_trimmed = truncate_germline(J_germline_seq, J_seq, J_blast_output_file, wdir, lineage_uid)
    
    ### Make full VJ sequence
    VJ_germline_seq = V_germline_seq_trimmed + J_germline_seq_trimmed

    ### Write sequence to output
    outfile = alignments_dir + "/" + str(lineage_uid) + ".germline_VJ.fa"

    with open(outfile, 'w') as out:
        out.write(">germline" + "\n")
        out.write(VJ_germline_seq + "\n")

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

def get_seq_by_id(fasta_file_name, seq_id):
    seq = ''
    with open(fasta_file_name,'r') as f:
        for record in SeqIO.parse(f,'fasta'):
            if record.id == seq_id:
                seq = record.seq
    f.close()
    return seq

def truncate_germline(germline_seq, seq, blast_output, wdir, lineage_uid):

    germline_seq = str(germline_seq).upper()
    seq = str(seq).upper()

    blast_params = {'evalue': 1000,
                    'word_size': 7,
                    'gapopen': 10,
                    'gapextend': 2,
                    'penalty': -3,
                    'reward': 1}

    blast_output = blastn2seq(germline_seq, seq, blast_output, wdir, **blast_params)
    blast_result_record = NCBIXML.read(open(blast_output, 'rU'))

    best_evalue = 100000
    best_hsp = None
    best_primer_i = -1

    evalue_warning_threshold = 0.001

    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:
            if float(hsp.expect) < best_evalue:
                best_evalue = float(hsp.expect)
                best_hsp = hsp

    if best_hsp == None:
        print "Warning: no HSP"
        print "Lineage uid = ", str(lineage_uid)
        print
        return germline_seq

    if best_evalue > evalue_warning_threshold:
        print "Warning: best E-value >", evalue_warning_threshold
        print "Best E value = ", str(best_evalue)
        print "Lineage uid = ", str(lineage_uid)
        print

    truncated_germline_seq = germline_seq[best_hsp.query_start-1:best_hsp.query_end]

    os.remove(blast_output)

    return truncated_germline_seq

def blastn2seq(query, subject, out="blastn_out.xml~", wdir=".", **blast_params):

    query_record = SeqRecord(Seq(str(query)), id="query")
    subject_record = SeqRecord(Seq(str(subject)), id="subject")

    u = str(uuid.uuid4())
    query_file = wdir + "/" + u + ".query.fasta"
    subject_file = wdir + "/" + u + ".subject.fasta"

    SeqIO.write(query_record, query_file, "fasta")
    SeqIO.write(subject_record, subject_file, "fasta")

    NcbiblastnCommandline(query=query_file, subject=subject_file, outfmt=5, out=out, **blast_params)()

    os.remove(query_file)
    os.remove(subject_file)

    return out

if __name__ == "__main__":
    main(sys.argv)

