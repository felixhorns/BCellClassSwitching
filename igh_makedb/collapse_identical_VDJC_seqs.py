import sys
import os

global verbose
verbose = True

field_dict = {"lineage_uid":0, "visit":2, "abundance":4, "isotype":5, "V_seq":12, "D_seq":13, "J_seq":14, "C_seq":15}

def main(argv):

    infile = argv[1]
    outfile = argv[2]

    out = open(outfile, 'w')

    lineage = []

    with open(infile, 'rU') as f:

        for line, next_line in pairwise(f):

            lineage.append(line)

            if next_line == None:
                collapse_identical_sequences_in_lineage(lineage, out)
                break

            lineage_uid = line.split("\t")[field_dict["lineage_uid"]]
            next_lineage_uid = next_line.split("\t")[field_dict["lineage_uid"]]

            if next_lineage_uid != lineage_uid:
                collapse_identical_sequences_in_lineage(lineage, out)
                lineage = []

    out.close()

    print "Done!!"
                
def collapse_identical_sequences_in_lineage(lineage, out):

    unique_records = {}
    visits = {}
    abundances = {}
    isotypes = {}

    # Group records by unique sequences
    for record in lineage:

        vals = record.rstrip().split("\t")

        visit = vals[field_dict["visit"]]
        abundance = vals[field_dict["abundance"]]
        isotype = vals[field_dict["isotype"]]
        V_seq = vals[field_dict["V_seq"]]
        D_seq = vals[field_dict["D_seq"]]
        J_seq = vals[field_dict["J_seq"]]
        C_seq = vals[field_dict["C_seq"]]

        VDJC_seq = V_seq + D_seq + J_seq + C_seq

        try:
            unique_records[VDJC_seq].append(record)
            visits[VDJC_seq].append(visit)
            abundances[VDJC_seq].append(abundance)
            isotypes[VDJC_seq].append(isotype)
        except:
            unique_records[VDJC_seq] = [record]
            visits[VDJC_seq] = [visit]
            abundances[VDJC_seq] = [abundance]
            isotypes[VDJC_seq] = [isotype]

    # Write record for each unique sequence with summary
    for unique_seq in unique_records.keys():

        my_visits = visits[unique_seq]
        my_records = unique_records[unique_seq]
        my_abundances = abundances[unique_seq]
        my_isotypes = isotypes[unique_seq]

        # Sort by visits
        s = sorted(zip(my_visits, my_records, my_abundances, my_isotypes), key=lambda tuple: tuple[0])
        my_visits_sorted = []
        my_records_sorted = []
        my_abundances_sorted = []
        my_isotypes_sorted = []
        for v, r, a, i in s:
            my_visits_sorted.append(v)
            my_records_sorted.append(r)
            my_abundances_sorted.append(a)
            my_isotypes_sorted.append(i)

        # Make summaries
        visits_out = ",".join(my_visits_sorted)
        abundances_out = ",".join(map(str, my_abundances_sorted))
        isotypes_out = ",".join(map(str, my_isotypes))

        # Write records
        is_first_seen = 1

        for record in my_records_sorted:

            vals = record.rstrip().split("\t")
            vals.append(visits_out)
            vals.append(abundances_out)
            vals.append(isotypes_out)
            vals.append(is_first_seen)
            out.write("\t".join(map(str, vals)) + "\n")            

            is_first_seen = 0

    return None

from itertools import tee, izip_longest
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip_longest(a, b)

if __name__ == "__main__":
    main(sys.argv)
