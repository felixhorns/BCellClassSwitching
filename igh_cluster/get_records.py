import sys
import os

infile_clusters = sys.argv[1] # clusters.txt.temp
infile_records = sys.argv[2] # seq_records.txt.temp

output_dir = os.path.dirname(infile_clusters)
out = open(output_dir + "/sequences_lineages.txt", 'w')

# Load record dictionary
records = {}
with open(infile_records, 'rU') as f:
    for line in f:
        vals = line.rstrip().split("\t")
        records[vals[0]] = "\t".join(vals[1:])

# Get record for each sequence
with open(infile_clusters, 'rU') as f:

    lineage_uid = 1

    for line in f:

        uids = line.rstrip().split(",")

        lineage_num_seqs = len(uids)

        for uid in uids:

            record = records[uid].rstrip().split("\t")

            record.insert(0, lineage_uid)
            record.insert(1, lineage_num_seqs)

            out.write("\t".join(map(str, record)) + "\n")

        lineage_uid += 1

out.close()
