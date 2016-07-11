import sys
import os

infile = sys.argv[1]

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

with open(infile) as f:
    for line in f:

        vals = line.rstrip().split('\t')

        handles = vals[0].split(",")
        paths = vals[1].split(",")
        output_dir = vals[2].strip()
        
        clustering_method = vals[-5]
        clustering_CDR3_similarity_cutoff = vals[-4]
        clustering_VDJ_similarity_cutoff = vals[-3]
        clustering_templated_similarity_cutoff = vals[-2]
        clustering_Q_cutoff = vals[-1]

        # Create output dir
        if not os.path.exists(output_dir): os.makedirs(output_dir)

        # Write list of input files into dir
        with open(output_dir+"/clustering_inputs.txt", 'w') as out:
            for h, p in zip(handles, paths):
                out.write(h.strip()+"\t"+p.strip()+"\n")

        # Write parameters into dir
        with open(output_dir+"/clustering_params.txt", 'w') as out:
            out.write(clustering_method + "\n")
            out.write(clustering_CDR3_similarity_cutoff + "\n")
            out.write(clustering_VDJ_similarity_cutoff + "\n")
            out.write(clustering_templated_similarity_cutoff + "\n")
            out.write(clustering_Q_cutoff + "\n")

        # Touch everything in dir (so that existing files will not be recreated)
        files = [output_dir + "/" + f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))]
        for f in files:
            touch(f)
