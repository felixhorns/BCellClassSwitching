"""
Parses output of IgBlast to custom format with sequences plus annotation.

Inputs:
    - directory containing IgBlast output file (igblast_out)
    - path to file containing isotype match dictionaries
    - cutoff on C sequence length min, e.g. 80
    - cutoff on C sequence length max, e.g. 120
    - cutoff on E value for V (log10 of value), e.g. -5.0
    - cutoff on E value for J (log10 of value), e.g. -5.0

Outputs:
    *** Note that these all have a split number extension 
    *** e.g. parsed_igblast.txt.{n}
    - parsed_igblast
    - losses_parse_igblast.txt
    - seq_ids_no_hits_found.txt
    - seq_ids_no_CDR3_start.txt
    - seq_ids_no_CDR3_end.txt
    - seq_ids_no_alignment_summary.txt
    - seq_ids_V_Evalue_greater_than_cutoff.txt
    - seq_ids_J_Evalue_greater_than_cutoff.txt
    - seq_ids_len_C_too_small_large.txt
    - seq_ids_no_V_call.txt
    - seq_ids_no_J_call.txt

Example:
python ~/Bcell_class_switching_prelim/scripts/parse_igblast.py . ~/Bcell_class_switching_prelim/resources/isotype_C_sequences.txt 80 120 -5.0 -5.0
"""

from __future__ import division
import sys
import os
import json
import regex
from distance import hamming
from Bio.Seq import Seq
import numpy as np

verbose = False

def main(argv):

    ##### Input
    in_sequences_abundances = argv[1]
    in_igblast_out = argv[2]
    C_primer_match_dict_file = argv[3]
    in_sequences_reads_per_molecule = argv[4]

    # Split file number
    split_num = in_sequences_abundances.split('.')[-1]
    if verbose: 
        print split_num
    # Filtering parameters
    len_C_cutoff_min = int(argv[5]) # 80
    len_C_cutoff_max = int(argv[6]) # 120
    log_V_Evalue_cutoff = float(argv[7]) # -5.0
    log_J_Evalue_cutoff = float(argv[8]) # -5.0

    ##### Output
    output_dir = os.path.dirname(in_sequences_abundances)
    out = open("%s/parsed_igblast.txt.%s" % (output_dir,split_num), 'w')
    out_C_seq_file = open("%s/C_seqs_without_primer.fasta.%s" % (output_dir,split_num), 'w')

    ##### Begin parsing
    C_primer_match_dict = load_dict(C_primer_match_dict_file) # Load matching dictionary for calling isotype
    num_records_parsed_successfully = 0

    global error_counts
    error_counts = {}
    error_counts["no_hits_found"] = 0
    error_counts["no_CDR3_start"] = 0
    error_counts["no_CDR3_end"] = 0
    error_counts["no_alignment_summary"] = 0
    error_counts["V_Evalue_greater_than_cutoff"] = 0
    error_counts["J_Evalue_greater_than_cutoff"] = 0
    error_counts["len_C_too_small_large"] = 0
    error_counts["no_V_call"] = 0
    error_counts["no_J_call"] = 0

    global seq_ids
    seq_ids = {}
    seq_ids["no_hits_found"] = []
    seq_ids["no_CDR3_start"] = []
    seq_ids["no_CDR3_end"] = []
    seq_ids["no_alignment_summary"] = []
    seq_ids["V_Evalue_greater_than_cutoff"] = []
    seq_ids["J_Evalue_greater_than_cutoff"] = []
    seq_ids["len_C_too_small_large"] = []
    seq_ids["no_V_call"] = []
    seq_ids["no_J_call"] = []

    global log_V_Evalues
    global log_D_Evalues
    global log_J_Evalues

    log_V_Evalues = []
    log_D_Evalues = []
    log_J_Evalues = []

    ### Load original query sequences
    sequences = {}

    file_length = get_file_length(in_sequences_abundances)

    line_num = 0
    with open(in_sequences_abundances, 'rU') as f:
        while line_num < file_length:
            seq_id = f.readline().rstrip().split(">")[1]
            seq = f.readline().rstrip()
            sequences[seq_id] = seq
            line_num += 2

    ### Load read per molecule counts
    reads_per_molecule = {}
    
    file_length = get_file_length(in_sequences_reads_per_molecule)

    line_num = 0
    with open(in_sequences_reads_per_molecule, 'rU') as f:
        while line_num < file_length:
            seq_id = f.readline().rstrip().split(">")[1]
            my_reads_per_molecule = f.readline().rstrip()
            reads_per_molecule[seq_id] = my_reads_per_molecule
            line_num += 2

    ### Parse igblast output
    with open(in_igblast_out, 'rU') as f:

    # concatenated igblast file no longer has headers or footers
    #    readlines(f, n = 15) # skip header

        while True:

            line = f.readline()
            if "Query=" in line:
                seq_id = line.rstrip().split(" ")[1]
                seq = sequences[seq_id]
                my_reads_per_molecule = reads_per_molecule[seq_id]
                num_records_parsed_successfully += parse_record(f, seq_id, seq, out, out_C_seq_file, C_primer_match_dict,
                                                                len_C_cutoff_min, len_C_cutoff_max,
                                                                log_V_Evalue_cutoff, log_J_Evalue_cutoff,
                                                                my_reads_per_molecule)
            else:
                if verbose:
                    print "Done!!"
                    print "Parsed successfully:", num_records_parsed_successfully
                    print "Errors:", error_counts
                break

    out.close()


    ### Write losses to log
    with open("%s/losses_parse_igblast.txt.%s" % (output_dir,split_num), 'w') as f:
        f.write("parse_igblast.py"+"\t"+"parsed_igblast"+"\t"+str(num_records_parsed_successfully)+"\n")
        f.write("parse_igblast.py"+"\t"+"no_hits_found"+"\t"+str(error_counts["no_hits_found"])+"\n")
        f.write("parse_igblast.py"+"\t"+"no_CDR3_start"+"\t"+str(error_counts["no_CDR3_start"])+"\n")
        f.write("parse_igblast.py"+"\t"+"no_CDR3_end"+"\t"+str(error_counts["no_CDR3_end"])+"\n")
        f.write("parse_igblast.py"+"\t"+"no_alignment_summary"+"\t"+str(error_counts["no_alignment_summary"])+"\n")
        f.write("parse_igblast.py"+"\t"+"V_Evalue_greater_than_cutoff"+"\t"+str(error_counts["V_Evalue_greater_than_cutoff"])+"\n")
        f.write("parse_igblast.py"+"\t"+"J_Evalue_greater_than_cutoff"+"\t"+str(error_counts["J_Evalue_greater_than_cutoff"])+"\n")
        f.write("parse_igblast.py"+"\t"+"len_C_too_small_large"+"\t"+str(error_counts["len_C_too_small_large"])+"\n")
        f.write("parse_igblast.py"+"\t"+"no_V_call"+"\t"+str(error_counts["no_V_call"])+"\n")
        f.write("parse_igblast.py"+"\t"+"no_J_call"+"\t"+str(error_counts["no_J_call"])+"\n")

    ### Write seq_ids of sequences lost to each error to files
    with open("%s/seq_ids_no_hits_found.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["no_hits_found"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_no_CDR3_start.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["no_CDR3_start"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_no_CDR3_end.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["no_CDR3_end"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_no_alignment_summary.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["no_alignment_summary"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_V_Evalue_greater_than_cutoff.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["V_Evalue_greater_than_cutoff"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_J_Evalue_greater_than_cutoff.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["J_Evalue_greater_than_cutoff"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_len_C_too_small_large.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["len_C_too_small_large"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_no_V_call.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["no_V_call"]:
            f.write(str(seq_id)+"\n")

    with open("%s/seq_ids_no_J_call.txt.%s" % (output_dir,split_num), 'w') as f:
        for seq_id in seq_ids["no_J_call"]:
            f.write(str(seq_id)+"\n")

def get_file_length(f):
    L = 0
    with open(f, 'rU') as infile:
        for line in infile: L += 1
    return L

def load_dict(f):
    d = {}
    with open(f, 'rU') as infile:
        for line in infile:
            try:
                d[line.split()[1].strip()] = line.split('\t')[0]
            except:
                continue
    return d

def readlines(f, n):
    for i in range(n): f.readline()

def readlines_to_end_of_record(f):
    while True:
        line = f.readline()
        if "Effective search space used" in line:
            readlines(f, n = 2)
            break
    return None

def parse_record(f, seq_id, seq, out, out_C_seq_file, C_primer_match_dict, len_C_cutoff_min, len_C_cutoff_max, log_V_Evalue_cutoff, log_J_Evalue_cutoff, my_reads_per_molecule):

    abundance = seq_id.split("_")[-1]

    seq = seq.upper()

    V = "?"
    V_Evalue = "NaN"
    D = "?"
    D_Evalue = "NaN"
    J = "?"
    J_Evalue = "NaN"
    stop_codon = ""
    productive = ""

    while True:
        line = f.readline()
        if "separated by a comma" in line: break
        if "IGHV" in line: V, V_Evalue = get_germline_gene_Evalue(line)
        if "IGHD" in line: D, D_Evalue = get_germline_gene_Evalue(line)
        if "IGHJ" in line: J, J_Evalue = get_germline_gene_Evalue(line)
        if "No hits found" in line:
            if verbose: print "Warning: No hits found. seq_id", seq_id, ". Skipping."
            error_counts["no_hits_found"] += 1
            seq_ids["no_hits_found"].append(seq_id)
            readlines_to_end_of_record(f)
            return 0

    line = f.readline()

    # Filter out sequences with Evalues worse than cutoff for V and J hits, or no call for V and J
    if V_Evalue != "NaN":
        if np.log10(float(V_Evalue)) > log_V_Evalue_cutoff:
            if verbose: print 'Warning: log10(V_Evalue) >', log_V_Evalue_cutoff, '. seq_id', seq_id, ". Skipping."
            error_counts["V_Evalue_greater_than_cutoff"] += 1
            seq_ids["V_Evalue_greater_than_cutoff"].append(seq_id)
            readlines_to_end_of_record(f)
            return 0

    if J_Evalue != "NaN":
        if np.log10(float(J_Evalue)) > log_J_Evalue_cutoff:
            if verbose: print 'Warning: log10(J_Evalue) >', log_J_Evalue_cutoff, '. seq_id', seq_id, ". Skipping."
            error_counts["J_Evalue_greater_than_cutoff"] += 1
            seq_ids["J_Evalue_greater_than_cutoff"].append(seq_id)
            readlines_to_end_of_record(f)
            return 0

    if V == "?":
        if verbose: print 'Warning: V == "?". seq_id', seq_id, ". Skipping."
        error_counts["no_V_call"] += 1
        seq_ids["no_V_call"].append(seq_id)
        readlines_to_end_of_record(f)
        return 0

    if J == "?":
        if verbose: print 'Warning: J == "?". seq_id', seq_id, ". Skipping."
        error_counts["no_J_call"] += 1
        seq_ids["no_J_call"].append(seq_id)
        readlines_to_end_of_record(f)
        return 0

    stop_codon, productive = get_stop_codon_productive(line)

    while True:
        line = f.readline()
        if "Alignment summary" in line: break
        if "Alignments" in line:
            if verbose: print "Warning: alignment summary did not exist. seq_id " + seq_id + ". Skipping."
            error_counts["no_alignment_summary"] += 1
            seq_ids["no_alignment_summary"].append(seq_id)
            readlines_to_end_of_record(f)
            return 0

    # Parse FWR/CDR boundaries
    boundaries, alignment_length = parse_boundaries(f)

    readlines(f, n = 5)

    # Parse alignment to get mutations
    query_start_position, query_reversed, query_end_position, mutation_positions, germline_bases, derived_bases, V_germline_identity = parse_alignment(f)

#    if query_reversed: seq = reverse_complement(seq) # TODO not sure how this could happen

    # Split sequence into regions
    leader_seq = seq[:query_start_position-1]
    VDJ_seq = seq[query_start_position-1:query_end_position-1]
    C_seq = seq[query_end_position-1:]

    # Call CDR3 start
    if boundaries[5] == -1:
        # CDR3 start has not been called, so try to call it
        nearest_boundary_to_CDR3_start = find_nearest_boundary_to_CDR3_start(boundaries) # start with the nearest boundary
        CDR3_start_with_respect_to_nearest_boundary = call_CDR3_start(VDJ_seq[nearest_boundary_to_CDR3_start:])
        if CDR3_start_with_respect_to_nearest_boundary != -1:
            boundaries[5] = CDR3_start_with_respect_to_nearest_boundary + nearest_boundary_to_CDR3_start

    if boundaries[5] == -1:
        if verbose: print "Warning: could not find CDR3_start. seq_id", seq_id, ". Skipping."
        error_counts["no_CDR3_start"] += 1
        seq_ids["no_CDR3_start"].append(seq_id)
        readlines_to_end_of_record(f)
        return 0

    # Call CDR3 end
    CDR3_end_with_respect_to_VDJ_seq = call_CDR3_end(VDJ_seq, boundaries[5]-(query_start_position-1))

    if CDR3_end_with_respect_to_VDJ_seq != -1:
        boundaries[6] = CDR3_end_with_respect_to_VDJ_seq + (query_start_position - 1)

    if boundaries[6] == -1:
        if verbose: print "Warning: could not find CDR3_end. seq_id", seq_id, ". Skipping."
        error_counts["no_CDR3_end"] += 1
        seq_ids["no_CDR3_end"].append(seq_id)
        readlines_to_end_of_record(f)
        return 0

    # Split sequence into V, D and J
    V_seq = seq[query_start_position-1:boundaries[5]-1]
    D_seq = seq[boundaries[5]-1:boundaries[6]-1]
    J_seq = seq[boundaries[6]-1:query_end_position-1]

    primer_isotype, C_seq_without_primer = call_C_primer(C_seq, C_primer_match_dict)

    AA = translate_AA(seq, CDR3_start=boundaries[5])

    mutation_density = len(mutation_positions) / len(VDJ_seq)

    len_V = len(V_seq)
    len_D = len(D_seq)
    len_J = len(J_seq)
    len_C = len(C_seq)

    region_lengths = calc_region_lengths(boundaries)

    boundaries_in_VDJ_seq_indexing = convert_boundaries_to_VDJ_seq_indexing(boundaries)

    # Check for failure to pass basic filters

    # Length of C
    if len_C < len_C_cutoff_min or len_C > len_C_cutoff_max:
        if verbose: print "Warning: length of C sequence <", len_C_cutoff_min, "or >", len_C_cutoff_max, ". seq_id", seq_id, ". Skipping."
        error_counts["len_C_too_small_large"] += 1
        seq_ids["len_C_too_small_large"].append(seq_id)
        readlines_to_end_of_record(f)
        return 0

    # Write record to parsed output file
    vals = [seq_id, abundance, V, D, J,
            V_Evalue, D_Evalue, J_Evalue,
            V_seq, D_seq, J_seq, C_seq,
            len_V, len_D, len_J, len_C,
            ",".join(map(str, boundaries_in_VDJ_seq_indexing)),
            ",".join(map(str, region_lengths)),
            stop_codon, productive, AA,
            ",".join(mutation_positions), ",".join(germline_bases), ",".join(derived_bases),
            mutation_density, V_germline_identity, leader_seq,
            my_reads_per_molecule, primer_isotype]

    out.write('\t'.join(map(str, vals)))
    out.write('\n')

    # Write C sequence to call isotype
    out_C_seq_file.write(">"+seq_id+"\n")
    out_C_seq_file.write(C_seq_without_primer+"\n")

    # Keep track of distributions
    if V_Evalue != "NaN": log_V_Evalues.append(np.log10(float(V_Evalue)))
    if D_Evalue != "NaN": log_D_Evalues.append(np.log10(float(D_Evalue)))
    if J_Evalue != "NaN": log_J_Evalues.append(np.log10(float(J_Evalue)))

    # Write output to standard out if verbose
    if verbose:
        for v in vals: print v
        print

    # Read to the end of the record
    readlines_to_end_of_record(f)

    return 1

def get_germline_gene_Evalue(line):
    vals = line.rstrip().split()
    gene = vals[0].split("lcl|")[1]
    E = vals[2]
    return gene, E

def get_stop_codon_productive(line):
    vals = line.rstrip().split("\t")

    if vals[4] == 'No':
        stop_codon = 0
    else:
        stop_codon = 1

    if vals[6] == 'No':
        productive = 0
    else:
        productive = 1

    return stop_codon, productive

def get_start(line):
    vals = line.rstrip().split("\t")
    return int(vals[1])

def get_end(line):
    vals = line.rstrip().split("\t")
    return int(vals[2])

def parse_boundaries(f):

    boundaries = [-1, -1, -1, -1, -1, -1, -1]

    while True:
        line = f.readline()
        if "FR1-IMGT" in line: boundaries[0] = get_start(line)
        if "CDR1-IMGT" in line: boundaries[1] = get_start(line)
        if "FR2-IMGT" in line: boundaries[2] = get_start(line)
        if "CDR2-IMGT" in line: boundaries[3] = get_start(line)
        if "FR3-IMGT" in line:
            # use FR3 end to get CDR3 start because sometimes CDR3 line is not present
            boundaries[4] = get_start(line)
            boundaries[5] = get_end(line) + 1
        if "CDR3-IMGT" in line: boundaries[5] = get_start(line)
        if "Total" in line:
            alignment_length = int(line.rstrip().split("\t")[3])
            break

    return boundaries, alignment_length

def parse_alignment(f):

    query_start_position = -1 # end of leader sequence, start of query sequence
    query_reversed = False
    query_end_position = -1 # end of alignment, start of constant region

    mutation_positions = []
    germline_bases = []
    derived_bases = []

    V_germline_identity = ""

    query = ""
    germline = ""
    germline_gene = ""

    germline_gene_in_previous_chunk = ""

    while True:

        line = f.readline().replace('\n', '')

        if line == "": continue

        if "Lambda" in line: break

        if "Query" in line:
            vals = line.split()
            query = vals[2]
            if query_start_position == -1: query_start_position = int(vals[1])
            if "reversed" in line: query_reversed = True
            query_end_position = int(vals[3])

        if "IGH" in line:

            # In rare cases (one observed so far), there is no start or last position.
            # We allow for the possibility that this occurs and assume start = last = 1.
            if len(line.split()) == 7:
                # Normal case
                germline = line.split()[5]
                germline_gene = line.split()[0] # V, D, or J
                start_position_in_germline = int(line.split()[4])
                last_position_in_germline = int(line.split()[6])
            elif len(line.split()) == 5:
                # No start or last position
                germline = line.split()[4]
                germline_gene = line.split()[0] # V, D, or J
                start_position_in_germline = 1
                last_position_in_germline = 1

            if germline_gene == "V": V_germline_identity = float(line.split()[1].replace("%", ""))/100

            if germline_gene == germline_gene_in_previous_chunk:
                entered_alignment = True
            else:
                entered_alignment = False

            m, g, d = call_mutations(query, germline, germline_gene, start_position_in_germline, last_position_in_germline, entered_alignment)

            mutation_positions.extend(m)
            germline_bases.extend(g)
            derived_bases.extend(d)

            germline_gene_in_previous_chunk = germline_gene

    return query_start_position, query_reversed, query_end_position, mutation_positions, germline_bases, derived_bases, V_germline_identity

def call_mutations(query, germline, germline_gene, start_position_in_germline, last_position_in_germline, entered_alignment=False):

    mutation_positions = []
    germline_bases = []
    derived_bases = []

    position_in_germline = start_position_in_germline - 1

    for q, g in zip(query, germline):

        if g == "-" and not entered_alignment: continue
        if g != "-" and not entered_alignment: entered_alignment = True # we have found the first position in the alignment

        if g in ["A", "T", "C", "G", "."]: position_in_germline += 1

        if q != g and g != ".":
            mutation_positions.append(germline_gene + str(position_in_germline))
            germline_bases.append(g)
            derived_bases.append(q)

        if position_in_germline == last_position_in_germline: break

    return mutation_positions, germline_bases, derived_bases

def call_CDR3_start(VDJ_seq):

    CDR3_start_anchor_sequence = 'TATTACTGT'
    minimum_match_distance = 2
    CDR3_start = -1

    for i in xrange(0, len(VDJ_seq) - len(CDR3_start_anchor_sequence) - 1):
        
        d = hamming(VDJ_seq[i:i+len(CDR3_start_anchor_sequence)], CDR3_start_anchor_sequence)
       
        if d <= minimum_match_distance:
            CDR3_start = i + len(CDR3_start_anchor_sequence) + 1
            minimum_match_distance = d


    return CDR3_start

def call_CDR3_end(VDJ_seq, CDR3_start):

    CDR3_end_anchor_sequence = 'CTGGGG'
    minimum_match_distance = 1
    CDR3_end = -1

    for i in xrange(CDR3_start, len(VDJ_seq) - len(CDR3_end_anchor_sequence) + 1):
        try:
            d = hamming(VDJ_seq[i:i+len(CDR3_end_anchor_sequence)], CDR3_end_anchor_sequence)
        except:
            print 'expected two strings of the same length for hamming'
            d=10
        if d <= minimum_match_distance:
            CDR3_end = i + 1
            minimum_match_distance = d

    return CDR3_end    
    
def call_C_primer(C_seq, C_primer_match_dict):

    primer_isotype = 'no_primer'
    primer_position = -1
    C_seq_without_primer = C_seq
    C_seq_length = len(C_seq)
    # series of ifs necessary for correct execution (as opposed to if / else)
    for primer_seq, primer_name in C_primer_match_dict.items():
        primer_position = C_seq[-25:].find(primer_seq)  # look only in the last 25nt for the primer seq
        if primer_position == -1:  # try to find a shortened version of the primer (for 12nt barcodes as opposed to 8)
            primer_position = C_seq[-25:].find(primer_seq[:-4])
        if primer_position == -1:  # still no exact primer match, moving to fuzzy primer search
            match = regex.findall('(%s){e<=2}'%primer_seq[:-4],C_seq[-25:],regex.BESTMATCH)  # fuzzy match allowing 2 errors
            if match != []:
                primer_position = C_seq[-25:].find(match[0])  # regex.findall returns a list
        if primer_position != -1:
            primer_isotype = primer_name
            # primer search is done on only the last 30nt of the constant
            # sequence, thus the match location must be reindexed relative to
            # the entire constant sequence
            adj_primer_pos = primer_position + C_seq_length - 25
            C_seq_without_primer = C_seq[:adj_primer_pos]
            break

    return primer_isotype, C_seq_without_primer

def translate_AA(seq, CDR3_start):
    if CDR3_start == -1: return "" # do not translate if we could not find CDR3 start anchor
    frame = CDR3_start % 3
    start_position = frame + 2
    num_trailing_bases = (len(seq) - start_position) % 3
    if num_trailing_bases > 0:
        seq_to_translate = seq[start_position:-num_trailing_bases]
    else:
        seq_to_translate = seq[start_position:]
    AA = Seq(seq_to_translate).translate()
    return AA

def calc_region_lengths(boundaries):
    lengths = [-1, -1, -1, -1, -1, -1]
    for i, l in enumerate(lengths):
        if boundaries[i+1] != -1 and boundaries[i] != -1:
            lengths[i] = boundaries[i+1] - boundaries[i]
    return lengths

def find_nearest_boundary_to_CDR3_start(boundaries):

    for i in range(0, 5)[::-1]:
        if boundaries[i] > -1:
            return boundaries[i] - 1

    return 0

def convert_boundaries_to_VDJ_seq_indexing(boundaries):

    # If the first value is not -1, 0, or 1, then correct boundary indexing with an offset
    if boundaries[0] > 1:

        offset = boundaries[0] - 1

        # Move all boundaries that are not -1 by offset
        boundaries_in_VDJ_seq_indexing = []

        for v in boundaries:

            if v != -1:
                v_in_new_index = v - offset
            else:
                v_in_new_index = -1
            
            boundaries_in_VDJ_seq_indexing.append(v_in_new_index)

        return boundaries_in_VDJ_seq_indexing

    else:

        return boundaries

if __name__ == "__main__":
    main(sys.argv)
