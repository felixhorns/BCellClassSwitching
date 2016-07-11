import sys
import os
import errno
import time
from itertools import tee, izip_longest
import numpy as np
import copy
import json
from Bio import SeqIO, AlignIO
from igh_mst import find_mst_from_alignment

global verbose
verbose = False

def main(argv):

    start_time = time.time()

    ### Parse arguments
    infile = argv[1]
    infile_field_dict = argv[2]
    path_to_edmonds_bin = argv[3]
    constrain_by_isotype = int(argv[4])

    i_worker = int(infile.split(".")[-1])

    wdir = os.path.dirname(infile)
    wdir_original = os.path.dirname(infile_field_dict)

    alignments_dir = wdir_original + "/alignments/"
    alignments_dir_scratch = wdir + "/alignments/"

    with open(infile_field_dict, 'rb') as f:
        field_dict = json.load(f)

    ### Set output file names
    if constrain_by_isotype == 1:

        outfile_sequences = wdir + "/mst.sequences.temp." + str(i_worker)
        outfile_lineages = wdir + "/mst.lineages.temp." + str(i_worker)
        
        outfile_sequences_header = ["uid", "mst_ancestor", "mst_ancestor_is_germline", "mst_ancestor_sequence_uid",
                                    "mst_ancestor_isotype", "mst_dist_to_ancestor",
                                    "mst_steps_to_root", "mst_dist_to_root", "mst_steps_to_root_normalized",
                                    "mst_dist_to_root_normalized", "mst_num_descendants", "mst_centrality"]
        outfile_lineages_header = ["uid", "abundance", "mst", "mst_max_centrality",
                                   "mst_num_descendants_from_germline"]
        msts_dir = wdir_original + "/msts/"
        msts_dir_scratch = wdir + "/msts/"

    elif constrain_by_isotype == 0:

        outfile_sequences = wdir + "/mst_no_isotype_constraint.sequences.temp." + str(i_worker)
        outfile_lineages = wdir + "/mst_no_isotype_constraint.lineages.temp." + str(i_worker)

        outfile_sequences_header = ["uid", "mst_no_isotype_constraint_ancestor", "mst_no_isotype_constraint_ancestor_is_germline",
                                    "mst_no_isotype_constraint_ancestor_sequence_uid",
                                    "mst_no_isotype_constraint_ancestor_isotype",
                                    "mst_no_isotype_constraint_dist_to_ancestor",
                                    "mst_no_isotype_constraint_steps_to_root", "mst_no_isotype_constraint_dist_to_root",
                                    "mst_no_isotype_constraint_steps_to_root_normalized",
                                    "mst_no_isotype_constraint_dist_to_root_normalized",
                                    "mst_no_isotype_constraint_num_descendants",
                                    "mst_no_isotype_constraint_centrality"]
        outfile_lineages_header = ["uid", "abundance", "mst_no_isotype_constraint",
                                   "mst_no_isotype_constraint_max_centrality",
                                   "mst_no_isotype_constraint_num_descendants_from_germline"]
        msts_dir = wdir_original + "/msts_no_isotype_constraint/"
        msts_dir_scratch = wdir + "/msts_no_isotype_constraint/"

    else:

        print "Error: constrain by isotype parameter must be either 0 or 1."
        quit()

    for d in [msts_dir, msts_dir_scratch]:
        mkdir(d)

    ### Initialize output files
    out_sequences = open(outfile_sequences, 'w')
    out_lineages = open(outfile_lineages, 'w')

    write_header(out_sequences, outfile_sequences_header)
    write_header(out_lineages, outfile_lineages_header)

    msts(infile, field_dict, wdir, alignments_dir, path_to_edmonds_bin,
         out_sequences, out_lineages, msts_dir_scratch, constrain_by_isotype)

    out_sequences.close()
    out_lineages.close()

    print "Done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

def msts(infile, field_dict, wdir, alignments_dir, path_to_edmonds_bin,
         out_sequences, out_lineages, msts_dir_scratch, constrain_by_isotype):

    with open(infile, 'rU') as f:
        for line in f:

            if line == None:
                print "Empty input file"
                break

            lineage_uid = line.rstrip()
            mst_lineage(lineage_uid, field_dict, wdir, alignments_dir, path_to_edmonds_bin,
                        out_sequences, out_lineages, msts_dir_scratch, constrain_by_isotype)

    return None

def mst_lineage(lineage_uid, field_dict, wdir, alignments_dir, path_to_edmonds_bin,
                out_sequences, out_lineages, msts_dir_scratch, constrain_by_isotype):

    ### Load alignment and parse into sequences, isotypes, abundances, and names
    alnfile = alignments_dir + "/" + str(lineage_uid) + ".afa.gapless.germline_VJ"
    aln = AlignIO.read(alnfile, 'fasta')

    seqs = []
    isotypes = []
    abundances = []
    uids = []

    root = -1

    for i, record in enumerate(aln):

        if record.id != "germline":
            uid = record.id.split("~")[0]
            isotype = record.id.split("~")[1]
            abundance = record.id.split("~")[2]
        else:
            uid = "NA"
            isotype = "germline"
            abundance = "NA"
            root = i

        seqs.append(str(record.seq))
        isotypes.append(isotype)
        abundances.append(abundance)
        uids.append(uid)

    if root < 0:
        print "Error: germline sequence not found"
        quit()

    ### Find MST
    mst, root, d_constrained, c, d, seqs = find_mst_from_alignment(seqs, isotypes, uids, root, wdir, path_to_edmonds_bin,
                                                                   constrain_by_isotype=constrain_by_isotype)
    mst_reversed = _reverse(mst)


    ### Calculate features
    # Distances
    N = len(isotypes) # number of sequences in lineage, including germline if rooted on germline
    steps_to_root, distances_to_root = calc_steps_and_distances_to_root(mst_reversed, N, root)
    steps_to_root_normalized = normalize_by_max(steps_to_root)
    distances_to_root_normalized = normalize_by_max(distances_to_root)

    # Number of descendants and centrality
    num_descendants = count_descendants(mst, N, isotypes)
    centralities = [float(i)/(N - 1) for i in num_descendants]
    max_centrality = max(centralities)

    # Count number of descendants from germline
    num_descendants_from_germline = count_descendants_from_germline(mst, isotypes, root)
    
    if verbose:
        print "*** Inputs ***"
        print "Number of sequences", N
        print seqs
        print isotypes
        print abundances
        print
        print "*** MST ***"
        print mst
        print root
        print isotypes
        print d
        print c
        print d_constrained
        print seqs
        print
        print "*** Features ***"
        print steps_to_root
        print steps_to_root_normalized
        print distances_to_root
        print distances_to_root_normalized
        print num_descendants
        print centralities
        print max_centrality
        print num_descendants_from_germline
        print
        print "*** Outputs ***"
    
    ### Write outputs

    ### MST
    labels = get_labels(uids, isotypes, abundances)
    mst_labeled = label_mst(mst, labels)

    path_to_mst_scratch = msts_dir_scratch + "/" + str(lineage_uid) + ".mst"
    with open(path_to_mst_scratch, 'w') as f:
        json.dump(mst_labeled, f)

    if verbose:
        print mst_labeled

    ### Sequences
    null_str = "NULL"

    for i, u in enumerate(uids):

        if u == "NA":
            # germline, root node
            continue

        # find ancestor
        ancestor, dist_to_anc = mst_reversed[i].popitem()
        mst_reversed[i][ancestor] = dist_to_anc # put item back

        # check whether ancestor is germline
        if ancestor == root:
            # ancestor is germline
            ancestor_sequence_uid = null_str
            ancestor_subisotype = null_str
            ancestor_is_germline = 1
        else:
            ancestor_sequence_uid = uids[ancestor]
            ancestor_subisotype = isotypes[ancestor]
            ancestor_is_germline = 0

        # Write to file        
        vals = [u, str(ancestor), ancestor_is_germline, ancestor_sequence_uid, ancestor_subisotype,
                dist_to_anc, steps_to_root[i], distances_to_root[i],
                steps_to_root_normalized[i], distances_to_root_normalized[i],
                num_descendants[i], centralities[i]]
        
        out_sequences.write("\t".join(map(str, vals)) + "\n")

        if verbose:
            print vals

    ### Lineages

    # Find total abundance of lineage
    total_abundance = 0
    for a in abundances:
        if a != "NA":
            total_abundance += int(a)

    vals = [lineage_uid, total_abundance, 1, max_centrality, num_descendants_from_germline]
    out_lineages.write("\t".join(map(str, vals)) + "\n")

    return None

def write_header(out, h): out.write("\t".join(h)+"\n")

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

def _reverse(graph):
    """ Constructs the reverse representation of the graph
    (destination -> source).
    """
    r = {}
    for src in graph:
        for (dst,c) in graph[src].items():
            if dst in r:
                r[dst][src] = c
            else:
                r[dst] = { src : c }
    return r

def calc_steps_and_distances_to_root(mst_reversed, N, root):
    """ Calculates number of steps and total distance to the root.
    Calls a recursive helper function, which walks the tree until
    it finds the root.  As it walks, it adds 1 to the step count and
    the distance to the immediate ancestor to the distance for each
    recursive call.
    """
    steps = []
    distances = []

    for i in range(0,N):
        s, d = _calc_steps_and_distances_to_root(mst_reversed, i, root)
        steps.append(s)
        distances.append(d)

    return steps, distances

def _calc_steps_and_distances_to_root(mst_reversed, i, root):
    """ Helper function for calc_steps_and_distances_to_root.
    Walks the tree and sums the number of steps and the total
    distance to the root vertex.
    """
    if i == root:
        return 0, 0
    else:
        anc, dist_to_anc = mst_reversed[i].popitem()
        mst_reversed[i][anc] = dist_to_anc # put item back
        s, d = _calc_steps_and_distances_to_root(mst_reversed, anc, root)
        return s+1, d+dist_to_anc

def normalize_by_max(x):
    """ Normalizes a vector by the maximum value.
    """
    m = max(x)
    if m == 0:
        return [1.0 for i in x]
    else:
        return [float(i)/m for i in x]

def count_descendants(mst, N, isotypes):
    counts = [0] * N
    for ancestor in mst.keys():
        counts[ancestor] = len(mst[ancestor].keys())
    return counts

def count_descendants_from_germline(mst, isotypes, germline_index):
    return len(mst[germline_index].keys())

def gini_coeff(x):
    # requires all values in x to be zero or positive numbers,
    # otherwise results are undefined
    n = len(x)
    s = x.sum()
    r = np.argsort(np.argsort(-x)) # calculates zero-based ranks
    return 1 - (2.0 * (r*x).sum() + s)/(n*s)

def calc_delta_feature(mst_reversed, N, root, isotypes, features):

    delta_feature = []

    # When rooting tree on germline, we must offset the indexing
    # into the feature list by 1.  This is because the vertexes
    # of the MST corresponding to the sequences are numbered
    # 1, ..., N, and the germline is 0.
    if "germline" in isotypes:
        index_offset = -1
    else:
        index_offset = 0

    for i in range(0,N):

        if i == root:
            delta_feature.append(missing_data_str)
        else:
            anc, dist_to_anc = mst_reversed[i].popitem()
            mst_reversed[i][anc] = dist_to_anc # put item back

            d = 0.
            if isotypes[anc] == "germline": # If ancestor is germline, then delta_feature should not exist
                d = missing_data_str
            else:
                try:
                    d = float(features[i + index_offset]) - float(features[anc + index_offset])
                except ValueError:
                    d = missing_data_str

            delta_feature.append(d)

    return delta_feature

def count_isotypes(lineage):
    counts = {"IgM":0, "IgD":0, "IgG":0, "IgA":0, "IgE":0}
    for line in lineage:
        subisotype = line.rstrip().split(sep)[field_dict["subisotype"]]
        for key in counts:
            if key in subisotype:
                counts[key] += 1
    return counts

def get_labels(uids, isotypes, abundances):

    labels = {}
    
    for i, (uid, isotype, abundance) in enumerate(zip(uids, isotypes, abundances)):
        labels[i] = str(i) + "," + str(uid) + "," + isotype + "," + str(abundance)

    return labels
    
def label_mst(mst, labels):
        
    mst_labeled = {}

    for ancestor in mst.keys():

        mst_labeled[labels[ancestor]] = {}

        for child, distance in mst[ancestor].items():

            mst_labeled[labels[ancestor]][labels[child]] = distance

    return mst_labeled

if __name__ == "__main__":
    main(sys.argv)
