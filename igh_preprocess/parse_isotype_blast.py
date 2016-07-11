"""
Parses output of blastn to call isotypes.

Inputs:
    - parsed_igblast file (previously split)
    - blastn output file  (previously split)

Outputs:
    - parsed_igblast_isotypes.{n}
    - losses_parse_isotype_blast.txt.{n}
    - seq_ids_no_isotype_call.txt.{n}

Example:                
python /datastore/rfhorns/Bcell_class_switching_prelim/scripts/igh_pipeline/parse_isotype_blast.py ./parsed_igblast ./isotype_blastn.out
"""

import sys
import os
import json
from Bio.Blast import NCBIXML
import numpy as np

import matplotlib
from socket import gethostname
hostname = gethostname()
matplotlib.use('Agg')
from matplotlib import pyplot as plt

path_to_isotype_to_color_dict = "/datastore/rfhorns/resources/immune/isotype_to_color_dict.json"

def main(argv):

    ##### Input
    in_parsed_igblast = argv[1]
    in_isotype_blastn_output = argv[2]
    
    # Split file number
    print in_parsed_igblast
    split_num = in_parsed_igblast.split('.')[-1]
    print split_num
    ##### Output
    output_dir = os.path.dirname(in_parsed_igblast)
    out = open("%s/parsed_igblast_isotypes.txt.%s" % (output_dir,split_num), 'w')

    ##### Parsing
    isotype_data = parse_isotype_blastn(in_isotype_blastn_output)

    isotypes, log_best_isotype_evalues, isotype_ambiguities, error_seq_ids = parse_parsed_igblast(in_parsed_igblast, isotype_data, out)

    report_errors(error_seq_ids, output_dir, split_num)

    return None

def parse_isotype_blastn(in_isotype_blastn_output):

    isotype_data = {}

    blast_records = NCBIXML.parse(open(in_isotype_blastn_output, 'rU'))

    for record in blast_records:
        seq_id, i = parse_blast_record(record)
        if i is not None: isotype_data[seq_id] = i

    return isotype_data

def parse_blast_record(record):

    seq_id = record.query

    # Get high scoring pairs
    hsps = [] 
    for alignment in record.alignments:
        subject_isotype = str(alignment.accession).split("_")[0]
        evalue = float(alignment.hsps[0].expect)
        score = float(alignment.hsps[0].score)
        hsps.append((subject_isotype, evalue, score))

    hsps_sorted_by_evalue = sorted(hsps, key=lambda tup: tup[1]) # Sort by E value
    hsps_sorted_by_evalue_unique = find_unique_keys(hsps_sorted_by_evalue) # Keep only best hit for each isotype

    isotype = "?"

    try:
        isotype = hsps_sorted_by_evalue_unique[0][0]
        best_isotype_evalue = hsps_sorted_by_evalue_unique[0][1]
        best_isotype_score = hsps_sorted_by_evalue_unique[0][2]
    except:
        isotype = "?"
        best_isotype_evalue = "NaN"
        best_isotype_score = "NaN"

    try:
        second_best_isotype = hsps_sorted_by_evalue_unique[1][0]
        second_best_isotype_evalue = hsps_sorted_by_evalue_unique[1][1]
        second_best_isotype_score = hsps_sorted_by_evalue_unique[1][2]
    except:
        second_best_isotype = "None"
        second_best_isotype_evalue = "NaN"
        second_best_isotype_score = "NaN"

    if best_isotype_score != "NaN" and second_best_isotype_score != "NaN":
        ambiguity = float(second_best_isotype_score) / float(best_isotype_score)
    else:
        ambiguity = "NaN"

    if isotype == "?":
        isotype_data = None
    else:
        isotype_data = [isotype, best_isotype_evalue, best_isotype_score,
                        second_best_isotype, second_best_isotype_evalue, second_best_isotype_score,
                        ambiguity]

    return seq_id, isotype_data

def find_unique_keys(L):
    unique_keys = set()
    L_unique = []
    for item in L:
        key = item[0]
        if key not in unique_keys:
            unique_keys.add(key)
            L_unique.append(item)
    return L_unique

def parse_parsed_igblast(in_parsed_igblast, isotype_data, out):

    # Keep track of errors
    error_seq_ids = {}
    error_seq_ids["no_isotype_call"] = []

    # Keep track of values for plotting distributions
    isotypes = []
    log_best_isotype_evalues = []
    isotype_ambiguities = []

    with open(in_parsed_igblast, 'rU') as f:

        for line in f:

            vals = line.rstrip().split("\t")
            seq_id = vals[0]

            try:

                i = isotype_data[seq_id]
                vals.insert(2, i[0]) # isotype
                vals.extend(i[1:])
                out.write("\t".join(map(str, vals))+"\n")

                if i[1] != "NaN" and i[-1] != "NaN":
                    isotypes.append(i[0])
                    log_best_isotype_evalues.append(np.log10(float(i[1])))
                    isotype_ambiguities.append(float(i[-1]))

            except:

                error_seq_ids["no_isotype_call"].append(seq_id)

    return isotypes, log_best_isotype_evalues, isotype_ambiguities, error_seq_ids

def report_errors(error_seq_ids, output_dir, split_num):

    with open("%s/losses_parse_isotype_blast.txt.%s" % (output_dir,split_num), 'w') as f:
        f.write("parse_isotype_blast.py"+"\t"+"no_isotype_call"+"\t"+str(len(error_seq_ids["no_isotype_call"]))+"\n")

    with open("%s/seq_ids_no_isotype_call.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in error_seq_ids["no_isotype_call"]:
            f.write(str(seq_id)+"\n")

    return None

def touch(f, times=None):
    with open(f, 'a'):
        os.utime(f, times)

if __name__ == "__main__":
    main(sys.argv)
