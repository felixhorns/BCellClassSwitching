import sys
import os

infile_seedfile = sys.argv[1] # seedfile_igh_cluster_original_dirs_subsampled.txt
infile_params = sys.argv[2] # all_cluster_params.txt

outfile = os.path.dirname(infile_seedfile) + "/seedfile_igh_cluster.txt"
out = open(outfile, 'w')

# Read parameter sets
parameter_sets = []
with open(infile_params, 'rU') as f:
    for line in f:
        parameter_sets.append(line.rstrip())

# Generate seedfile
with open(infile_seedfile, 'rU') as f:
    for line in f:
        for i, param_set in enumerate(parameter_sets):
            vals = line.rstrip().split("\t")
            new_output_dir = vals[2].rstrip("/") + "/clustering_params" + str(i + 1) + "/"
            vals[2] = new_output_dir
            out.write("\t".join(vals) + "\t" + param_set + "\n")

out.close()
