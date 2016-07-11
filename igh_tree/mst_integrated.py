import sys
import os
import errno
import time
from itertools import tee, izip_longest
import numpy as np
import copy
import json
from Bio import SeqIO
from igh_mst import find_mst

global verbose
verbose = False

def main(argv):

    ### Parse arguments
    infile_query = argv[1]
    infile_query_field_dict = argv[2] # used to get original working dir (not scratch)
    i_worker = int(argv[3].split(".")[-1])
    path_to_V_germline_ref_file = argv[4]
    path_to_J_germline_ref_file = argv[5]
    path_to_muscle_bin = argv[6]
    path_to_edmonds_bin = argv[7]
    constrain_by_isotype = int(argv[8])

    wdir = os.path.dirname(infile_query)
    wdir_original = os.path.dirname(infile_query_field_dict)

    with open(infile_query_field_dict, 'rb') as f:
        field_dict = json.load(f)

    ### Set output file names
    if constrain_by_isotype == 1:

        outfile_sequences = wdir + "/mst.sequences.temp." + str(i_worker)
        outfile_lineages = wdir + "/mst.lineages.temp." + str(i_worker)
        
        outfile_sequences_header = ["uid", "mst_ancestor", "mst_ancestor_is_germline", "mst_ancestor_sequence_uid",
                                    "mst_ancestor_isotype", "mst_ancestor_uid_seq_str", "mst_dist_to_ancestor",
                                    "mst_steps_to_root", "mst_dist_to_root", "mst_steps_to_root_normalized",
                                    "mst_dist_to_root_normalized", "mst_num_descendants", "mst_centrality"]
        outfile_lineages_header = ["uid", "abundance", "gini_abundance", "mst", "mst_max_centrality",
                                   "mst_num_descendants_from_germline"]
        msts_dir = wdir_original + "/msts/"
        msts_dir_scratch = wdir + "/msts/"

    elif constrain_by_isotype == 0:

        outfile_sequences = wdir + "/mst_no_isotype_constraint.sequences.temp." + str(i_worker)
        outfile_lineages = wdir + "/mst_no_isotype_constraint.lineages.temp." + str(i_worker)

        outfile_sequences_header = ["uid", "mst_no_isotype_constraint_ancestor", "mst_no_isotype_constraint_ancestor_is_germline",
                                    "mst_no_isotype_constraint_ancestor_sequence_uid",
                                    "mst_no_isotype_constraint_ancestor_isotype", "mst_no_isotype_constraint_ancestor_uid_seq_str",
                                    "mst_no_isotype_constraint_dist_to_ancestor",
                                    "mst_no_isotype_constraint_steps_to_root", "mst_no_isotype_constraint_dist_to_root",
                                    "mst_no_isotype_constraint_steps_to_root_normalized",
                                    "mst_no_isotype_constraint_dist_to_root_normalized",
                                    "mst_no_isotype_constraint_num_descendants",
                                    "mst_no_isotype_constraint_centrality"]
        outfile_lineages_header = ["uid", "abundance", "gini_abundance", "mst_no_isotype_constraint",
                                   "mst_no_isotype_constraint_max_centrality",
                                   "mst_no_isotype_constraint_num_descendants_from_germline"]
        msts_dir = wdir_original + "/msts_no_isotype_constraint/"
        msts_dir_scratch = wdir + "/msts_no_isotype_constraint/"

    else:

        print "Error: constrain by isotype parameter must be either 0 or 1."
        quit()

    mkdir(msts_dir)
    mkdir(msts_dir_scratch)

    ### Find MSTs
    start_time_mst = time.time()

    out_sequences = open(outfile_sequences, 'w')
    out_lineages = open(outfile_lineages, 'w')

    write_header(out_sequences, outfile_sequences_header)
    write_header(out_lineages, outfile_lineages_header)

    make_msts(infile_query, field_dict, wdir, path_to_muscle_bin, path_to_edmonds_bin,
              out_sequences, out_lineages, msts_dir_scratch,
              path_to_V_germline_ref_file, path_to_J_germline_ref_file, constrain_by_isotype)

    out_sequences.close()
    out_lineages.close()

    print "Done!!"
    print "mst elapsed time (wall clock):", time.time() - start_time_mst

def make_msts(infile, field_dict, wdir, path_to_muscle_bin, path_to_edmonds_bin,
              out_sequences, out_lineages, msts_dir_scratch,
              path_to_V_germline_ref_file, path_to_J_germline_ref_file, constrain_by_isotype):

    lineage = []

    with open(infile, 'rU') as f:
        for line, next_line in pairwise(f):

            if line == None:
                print "Empty input file"
                break

            lineage.append(line.rstrip())

            if next_line == None:
                # end of file
                make_mst_for_lineage(lineage, field_dict, wdir, path_to_muscle_bin, path_to_edmonds_bin,
                                     out_sequences, out_lineages, msts_dir_scratch,
                                     path_to_V_germline_ref_file, path_to_J_germline_ref_file, constrain_by_isotype)
                break

            this_lineage_uid = line.split("\t")[field_dict["lineage_uid"]]
            next_lineage_uid = next_line.split("\t")[field_dict["lineage_uid"]]

            if next_lineage_uid != this_lineage_uid:
                # end of lineage
                make_mst_for_lineage(lineage, field_dict, wdir, path_to_muscle_bin, path_to_edmonds_bin,
                                     out_sequences, out_lineages, msts_dir_scratch,
                                     path_to_V_germline_ref_file, path_to_J_germline_ref_file, constrain_by_isotype)
                lineage = []

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

def make_mst_for_lineage(lineage, field_dict, wdir, path_to_muscle_bin, path_to_edmonds_bin,
                         out_sequences, out_lineages, msts_dir_scratch,
                         path_to_V_germline_ref_file, path_to_J_germline_ref_file, constrain_by_isotype):

    ### Do not process if lineage has only 1 sequence
    if len(lineage) == 1: return None

    ### Parse inputs
    lineage_uid = lineage[0].split("\t")[field_dict["lineage_uid"]]
    V_germline = lineage[0].split("\t")[field_dict["V_germline"]]
    D_germline = lineage[0].split("\t")[field_dict["D_germline"]]
    J_germline = lineage[0].split("\t")[field_dict["J_germline"]]

    sequence_uids = []
    sequence_string_uids = []
    VDJ_seqs = []
    subisotypes = []
    abundances = []

    for line in lineage:

        vals = line.split("\t")

        VDJ_seq = vals[field_dict["V_seq"]] + vals[field_dict["D_seq"]] + vals[field_dict["J_seq"]]

        sequence_uids.append(vals[field_dict["sequence_uid"]])
        sequence_string_uids.append(vals[field_dict["sequence_string_uid"]])
        VDJ_seqs.append(VDJ_seq)
        subisotypes.append(vals[field_dict["subisotype"]])
        abundances.append(int(vals[field_dict["abundance"]]))


    ### Get germline sequence
    V_germline_seq = str(get_seq_by_id(path_to_V_germline_ref_file, V_germline)).upper()
    J_germline_seq = str(get_seq_by_id(path_to_J_germline_ref_file, J_germline)).upper()
    VJ_germline_seq = V_germline_seq + J_germline_seq

    ### Find MST
    mst, root, isotypes_with_root, d_constrained, c, d, seqs = find_mst(VDJ_seqs, subisotypes, VJ_germline_seq,
                                                                        wdir, path_to_muscle_bin, path_to_edmonds_bin,
                                                                        constrain_by_isotype=constrain_by_isotype)
    mst_reversed = _reverse(mst)

    ### Calculate features
    # Distances
    N = len(isotypes_with_root) # number of sequences in lineage, including germline if rooted on germline
    steps_to_root, distances_to_root = calc_steps_and_distances_to_root(mst_reversed, N, root)
    steps_to_root_normalized = normalize_by_max(steps_to_root)
    distances_to_root_normalized = normalize_by_max(distances_to_root)

    # Number of descendants and centrality
    num_descendants = count_descendants(mst, N, isotypes_with_root)
    centralities = [float(i)/len(lineage) for i in num_descendants]
    max_centrality = max(centralities)

    # Count number of descendants from germline
    num_descendants_from_germline = count_descendants_from_germline(mst, isotypes_with_root)
    
    # Gini coefficients
    gini_abundance = gini_coeff(np.array(abundances))

    ### Unused features
    # Change in mutations and sigma between sequences along tree
    # delta_Focused_Sigma_CDR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, Focused_Sigma_CDRs)
    # delta_Focused_Sigma_FWR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, Focused_Sigma_FWRs)

    # delta_Observed_CDR_R = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, Observed_CDR_Rs)
    # delta_Observed_CDR_S = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, Observed_CDR_Ss)
    # delta_Observed_FWR_R = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, Observed_FWR_Rs)
    # delta_Observed_FWR_S = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, Observed_FWR_Ss)
    
    # delta_CDR1Focused_Sigma_CDR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, CDR1Focused_Sigma_CDRs)
    # delta_CDR2Focused_Sigma_CDR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, CDR2Focused_Sigma_CDRs)
    # delta_CDR3Focused_Sigma_CDR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, CDR3Focused_Sigma_CDRs)
    # delta_FWR1Focused_Sigma_FWR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, FWR1Focused_Sigma_FWRs)
    # delta_FWR2Focused_Sigma_FWR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, FWR2Focused_Sigma_FWRs)
    # delta_FWR3Focused_Sigma_FWR = calc_delta_feature(mst_reversed, N, root, isotypes_with_root, FWR3Focused_Sigma_FWRs)

    # Class switching rates
    # class_switch_rates, class_switch_rates_to_alternative_isotypes, arrival_rates, num_descendants_by_isotype = calc_class_switch_rates(mst_reversed, N, root, isotypes, isotypes_with_root)

    # # Calculate mutation rates
    # muts, muts_3mer_templated, muts_3mer_untemplated, muts_5mer_templated, muts_5mer_untemplated = count_muts_lineage(mst_reversed, N, aln, AAs, V_seqs, J_seqs, codon_position_to_n_or_s, temp_fasta_file_name)
    # (contexts_3mer_templated, contexts_3mer_untemplated,
    #     contexts_5mer_templated, contexts_5mer_untemplated,
    #     contexts_3mer_FWR, contexts_3mer_CDR,
    #     contexts_5mer_FWR, contexts_5mer_CDR) = count_contexts_lineage(mst_reversed, N, aln, AAs, V_seqs, J_seqs, regions_list, codon_position_to_n_or_s, temp_fasta_file_name)

    if verbose:
        print "*** Inputs ***"
        print "Number of sequences", str(len(sequence_uids))
        print sequence_uids
        print sequence_string_uids
        print VDJ_seqs
        print subisotypes
        print abundances
        print VJ_germline_seq
        print
        print "*** MST ***"
        print mst
        print root
        print isotypes_with_root
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
        print gini_abundance
        print
        print "*** Outputs ***"
    
    ### Write to outputs

    # MST
    labels = get_labels(sequence_uids, isotypes_with_root, abundances)
    mst_labeled = label_mst(mst, labels)

    path_to_mst_scratch = msts_dir_scratch + "/" + str(lineage_uid) + ".mst"
    with open(path_to_mst_scratch, 'w') as f:
        json.dump(mst_labeled, f)

    if verbose:
        print mst_labeled

    null_str = "NULL"

    # Sequences

    # Set offsets for indexing
    # When rooting tree on germline, we must offset the indexing
    # into the isotypes list by 1.  This is because the vertexes
    # of the MST corresponding to the sequences are numbered
    # 1, ..., N, and the germline is 0.
    if "germline" in isotypes_with_root:
        index_offset = -1
        root = -1 # Root is germline
    else:
        index_offset = 0

    for i, u in enumerate(sequence_uids):
        
        # Identify ancestor
        if i + 1 == root:
            ancestor = "root"
            dist_to_anc = null_str
        else:
            ancestor, dist_to_anc = mst_reversed[i + 1].popitem()
            mst_reversed[i + 1][ancestor] = dist_to_anc # put item back

        # Find features of ancestor
        if isotypes_with_root[ancestor] == 'germline':

            ancestor_sequence_uid = null_str
            ancestor_subisotype = null_str
            ancestor_VDJ_seq = null_str
            ancestor_sequence_string_uid = null_str
            ancestor_is_germline = 1

        else:
        
            ancestor_sequence_uid = sequence_uids[ancestor + index_offset]
            ancestor_subisotype = subisotypes[ancestor + index_offset]
            ancestor_VDJ_seq = VDJ_seqs[ancestor + index_offset]
            ancestor_sequence_string_uid = sequence_string_uids[ancestor + index_offset]
            ancestor_is_germline = 0

        # Write to file        
        vals = [u, str(ancestor), ancestor_is_germline, ancestor_sequence_uid, ancestor_subisotype,
                ancestor_sequence_string_uid, dist_to_anc, steps_to_root[i + 1], distances_to_root[i + 1],
                steps_to_root_normalized[i + 1], distances_to_root_normalized[i + 1],
                num_descendants[i + 1], centralities[i + 1]]
        
        out_sequences.write("\t".join(map(str, vals)) + "\n")

        if verbose:
            print vals

    # Lineages
    vals = [lineage_uid, sum(abundances), gini_abundance, 1, max_centrality, num_descendants_from_germline]
    out_lineages.write("\t".join(map(str, vals)) + "\n")

    if verbose:
        print vals
        print
        print

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

def count_descendants(mst, N, isotypes_with_root):
    counts = [0] * N
    for ancestor in mst.keys():
        counts[ancestor] = len(mst[ancestor].keys())
    return counts

def count_descendants_from_germline(mst, isotypes_with_root):
    if "germline" in isotypes_with_root:
        return len(mst[0].keys())
    else:
        return missing_data_str

def gini_coeff(x):
    # requires all values in x to be zero or positive numbers,
    # otherwise results are undefined
    n = len(x)
    s = x.sum()
    r = np.argsort(np.argsort(-x)) # calculates zero-based ranks
    return 1 - (2.0 * (r*x).sum() + s)/(n*s)

def calc_delta_feature(mst_reversed, N, root, isotypes_with_root, features):

    delta_feature = []

    # When rooting tree on germline, we must offset the indexing
    # into the feature list by 1.  This is because the vertexes
    # of the MST corresponding to the sequences are numbered
    # 1, ..., N, and the germline is 0.
    if "germline" in isotypes_with_root:
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
            if isotypes_with_root[anc] == "germline": # If ancestor is germline, then delta_feature should not exist
                d = missing_data_str
            else:
                try:
                    d = float(features[i + index_offset]) - float(features[anc + index_offset])
                except ValueError:
                    d = missing_data_str

            delta_feature.append(d)

    return delta_feature

def calc_class_switch_rates(mst_reversed, N, root, isotypes, isotypes_with_root):

    mst_reversed = copy.deepcopy(mst_reversed) # copy graph, so we do not edit object

    # When rooting tree on germline, we must offset the indexing
    # into the isotypes list by 1.  This is because the vertexes
    # of the MST corresponding to the sequences are numbered
    # 1, ..., N, and the germline is 0.
    if "germline" in isotypes_with_root:
        index_offset = -1
    else:
        index_offset = 0

    # Count number of ancestor -> descendant edges for each pair of isotypes
    num_descendants_by_isotype = np.zeros((len(isotype_to_int),len(isotype_to_int)))
    for i in range(0,N):
        if i == root: continue
        ancestor, _ = mst_reversed[i].popitem()
        if isotypes_with_root[ancestor] == 'germline': continue
        num_descendants_by_isotype[isotype_to_int[isotypes[ancestor + index_offset]]][isotype_to_int[isotypes[i + index_offset]]] += 1

    # Normalize to obtain rates

    # Rate of class switching to any alternative isotype (units are per child)
    # Sum of off diagonals (which equals 1 - diagonal) divided by row sum
    r = 1. - np.divide(num_descendants_by_isotype.diagonal(), num_descendants_by_isotype.sum(axis=1))

    # Rate of class switching to each alternative isotype, given a switch event (units are per event)
    # Off diagonals divided by row sum of off diagonals
    row_sums = num_descendants_by_isotype.sum(axis=1)
    r_alternative_isotypes = num_descendants_by_isotype / (row_sums - num_descendants_by_isotype.diagonal())[:,np.newaxis]

    # Rate of arrival from each isotype, given descendant isotype
    col_sums = num_descendants_by_isotype.sum(axis=0)
    r_arrival = num_descendants_by_isotype / col_sums[np.newaxis,:]

    return r, r_alternative_isotypes, r_arrival, num_descendants_by_isotype

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
    
    uids_copy = list(uids)
    uids_copy.insert(0, "NA")

    abundances_copy = list(abundances)
    abundances_copy.insert(0, "NA")

    for i, (uid, isotype, abundance) in enumerate(zip(uids_copy, isotypes, abundances_copy)):
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
