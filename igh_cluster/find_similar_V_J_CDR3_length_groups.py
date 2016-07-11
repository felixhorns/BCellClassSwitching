""" Finds groups of sequences with the same V, J, and CDR3 length.
"""

from __future__ import division
import sys
import os
import time
import zipfile

global verbose
verbose = False

field_dict = {"id":0, "abundance":1, "V":3, "J":5, "V_seq":9, "D_seq":10, "J_seq":11, "boundaries":17, "region_lengths":18, "quality_scores":-1}

def main(argv):

    # Inputs
    path_to_clustering_input_file = argv[1] # clustering_inputs.txt
    path_to_V_allele_collapse_dict_file = argv[2] # resources/Vsegments_20150201_IGHV70_primer_amplicons_identicals_collapsed.txt

    # Outputs
    output_dir = os.path.dirname(path_to_clustering_input_file)

    V_allele_collapse_dict = load_allele_collapse_dict(path_to_V_allele_collapse_dict_file)

    # Find groups of sequences sharing similar V, J, and CDR3 length
    start_time = time.time()

    print "Starting..."

    records, groups = find_groups(path_to_clustering_input_file, V_allele_collapse_dict)

    # Write records to look up rest of data associated with each sequence after clustering
    write_records(records, output_dir)

    # Write groups into individual files for clustering
    write_groups(groups, output_dir)

    time_elapsed = time.time() - start_time

    print "Wall clock time elapsed", time_elapsed
    print "Done!!"

def load_allele_collapse_dict(f):
    """ Loads a dictionary used to collapse different alleles of the 
    same V gene onto a single number. Each set of alleles that will
    be collapsed to the same gene are listed on a single line starting
    with ">" and separated by ",".
    """
    d = {}
    i = 0
    for line in open(f):
        i += 1
        if line[0] == '>':
            alleles = line.strip().replace('>','').split(',')
            for a in alleles:
                d[a] = i
    return d

def find_groups(path_to_clustering_input_file, V_allele_collapse_dict):
    """ Builds 2 dictionaries:
    (1) Records, containing every line of the input parsed igblast files,
    indexed by the unique id of the line (library_id plus sequence_id).
    This will be used to retrieve all the data associated with a sequence
    after we have clustered the sequences.
    (2) Unclustered sequences, containing every sequence in the input
    parsed igblast files, indexed by the key of its group. Each group
    is defined by the V and J germline calls, and the CDR3 length. Each
    element of the group has its sequence, uid (for retrieving its full
    data later), region start and end (for knowing which region of the
    sequence to compare to another), and abundance and quality scores
    (for filtering based on quality).
    """

    records = []
    groups = {}

    with open(path_to_clustering_input_file, 'rU') as f:

        for line in f:

            lib_id = line.split()[0]
            path_to_parsed_igblast_quals = line.split()[1]

            with open(path_to_parsed_igblast_quals+"/parsed_igblast_isotypes_quals") as f2:
                
                for record in f2:

                    vals = record.rstrip().split("\t")

                    seq_id = vals[field_dict["id"]]
                    abundance = int(vals[field_dict["abundance"]])
                    V = vals[field_dict["V"]]
                    J = vals[field_dict["J"]]
                    V_seq = vals[field_dict["V_seq"]]
                    D_seq = vals[field_dict["D_seq"]]
                    J_seq = vals[field_dict["J_seq"]]
                    boundaries = vals[field_dict["boundaries"]]
                    region_lengths = vals[field_dict["region_lengths"]]
                    quality_scores = vals[field_dict["quality_scores"]].split("~")

                    uid = lib_id + "-" + seq_id

                    vals.insert(0, lib_id)
                    vals.insert(0, uid)

                    record_with_uid = "\t".join(vals)
                    records.append(record_with_uid)
                    
                    # Key for clustering, which defines group
                    V_id = V_allele_collapse_dict[V]
                    CDR3_length = region_lengths.split(",")[-1]
                    key = str(V_id) + "_" + J + "_" + str(CDR3_length)

                    # Get region of sequence used for comparison
                    # Truncate quality score string to region used for comparison

                    # VDJ
                    VDJ_seq = V_seq + D_seq + J_seq
                    qual = get_first_truncated_quality_score(quality_scores, 0, len(VDJ_seq))

                    # CDR3
                    CDR3_start = int(boundaries.split(",")[5]) + 6
                    CDR3_end = CDR3_start - 6 + int(CDR3_length)
                    # CDR3_end = int(boundaries.split(",")[6]) - 1

                    CDR3_seq = VDJ_seq[CDR3_start:CDR3_end]
                    CDR3_qual = get_first_truncated_quality_score(quality_scores, CDR3_start, CDR3_end)

                    # Templated region
                    templated_seq = VDJ_seq[:CDR3_start] + VDJ_seq[CDR3_end:]
                    templated_qual = get_first_truncated_quality_score(quality_scores, 0, CDR3_start) + get_first_truncated_quality_score(quality_scores, CDR3_end, len(VDJ_seq))

                    # Templated region before and after CDR3 separately
                    templated_seq_before_CDR3 = VDJ_seq[:CDR3_start]
                    templated_qual_before_CDR3 = get_first_truncated_quality_score(quality_scores, 0, CDR3_start)
                    templated_seq_after_CDR3 = VDJ_seq[CDR3_end:]
                    templated_qual_after_CDR3 = get_first_truncated_quality_score(quality_scores, CDR3_end, len(VDJ_seq))

                    uid_abundance_seq_qual = [uid, abundance, VDJ_seq, qual, CDR3_seq, CDR3_qual, templated_seq, templated_qual,
                                              templated_seq_before_CDR3, templated_qual_before_CDR3, templated_seq_after_CDR3, templated_qual_after_CDR3]

                    try:
                        groups[key].append(uid_abundance_seq_qual)
                    except:
                        groups[key] = [uid_abundance_seq_qual]

    return records, groups

def truncate_quality_scores(quality_scores, start, end):
    truncated = []
    for q in quality_scores:
        truncated.append(q[start:end])
    return truncated

def get_first_truncated_quality_score(quality_scores, start, end):
    q = quality_scores[0]
    return q[start:end]

def write_groups(groups, output_dir):

    zf = zipfile.ZipFile(output_dir + "/groups.zip.temp", 'w', allowZip64=True)

    for key, group in groups.items():

        name = key + ".txt"
        name = name.replace("*", "-")

        vals = []

        for seq in group:
            vals.append("~".join(map(str, seq)))

        s = "}".join(vals)

        try:
            zf.writestr(name, s)
        except:
            pass

    zf.close()

    return None

def write_records(records, output_dir):

    output = output_dir + "/seq_records.txt.temp"

    with open(output, 'w') as out:
        for record in records:
            out.write(record+"\n")

    return None

if __name__ == "__main__":
    main(sys.argv)


