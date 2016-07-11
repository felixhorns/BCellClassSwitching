""" Finds clusters of similar sequences.
"""

from __future__ import division
import sys
import os
import time
import zipfile
import numpy as np

global verbose
verbose = False

def main(argv):

    # Inputs
    infile_groups = argv[1] # groups.zip
    infile_group_worker_assignments = argv[2] # group_worker_assignments.txt
    i_worker = int(argv[3].split('.')[-1]) # clusters.txt.temp.{n}
    infile_clustering_params = argv[4] # clustering_params.txt

    # Outputs
    output_dir = os.path.dirname(infile_groups)
    outfile_CDR3_distances = output_dir + "/CDR3_distances.txt.temp." + str(i_worker)
    outfile_templated_distances = output_dir + "/templated_distances.txt.temp." + str(i_worker)
    outfile_min_CDR3_distances = output_dir + "/min_CDR3_distances.txt.temp." + str(i_worker)
    outfile_min_templated_distances = output_dir + "/min_templated_distances.txt.temp." + str(i_worker)

    outfile_CDR3_similarities = output_dir + "/CDR3_similarities.txt.temp." + str(i_worker)
    outfile_templated_similarities = output_dir + "/templated_similarities.txt.temp." + str(i_worker)
    outfile_min_CDR3_similarities = output_dir + "/min_CDR3_similarities.txt.temp." + str(i_worker)
    outfile_min_templated_similarities = output_dir + "/min_templated_similarities.txt.temp." + str(i_worker)


    # Get input filenames
    with open(infile_group_worker_assignments, 'rU') as f:
        for i, line in enumerate(f):
            if i == i_worker:
                filenames = line.rstrip().split("\t")

    # Get clustering parameters
    with open(infile_clustering_params, 'rU') as f:
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        Q_cutoff = int(f.readline().rstrip()) # mismatches with Q scores lower than this value are not counted (Q = -log10(probability of error)) (e.g. 5)

    # Calculate pairwise distances within groups
    start_time = time.time()

    (pairwise_CDR3_distances, pairwise_templated_distances, min_CDR3_distances, min_templated_distances,
     pairwise_CDR3_similarities, pairwise_templated_similarities, min_CDR3_similarities, min_templated_similarities) = calc_distances(infile_groups, filenames, Q_cutoff)

    write_list(pairwise_CDR3_distances, outfile_CDR3_distances)
    write_list(pairwise_templated_distances, outfile_templated_distances)
    write_list(min_CDR3_distances, outfile_min_CDR3_distances)
    write_list(min_templated_distances, outfile_min_templated_distances)

    write_list(pairwise_CDR3_similarities, outfile_CDR3_similarities)
    write_list(pairwise_templated_similarities, outfile_templated_similarities)
    write_list(min_CDR3_similarities, outfile_min_CDR3_similarities)
    write_list(min_templated_similarities, outfile_min_templated_similarities)

    time_elapsed = time.time() - start_time

    print "Wall clock time elapsed", time_elapsed
    print "Done!!"

def calc_distances(infile, filenames, Q_cutoff):

    pairwise_CDR3_distances = []
    pairwise_templated_distances = []
    min_CDR3_distances = []
    min_templated_distances = []

    pairwise_CDR3_similarities = []
    pairwise_templated_similarities = []
    min_CDR3_similarities = []
    min_templated_similarities = []

    zf = zipfile.ZipFile(infile)

    for filename in filenames:

        try:
            f = zf.read(filename)
        except KeyError:
            print 'ERROR: Did not find %s in zip file' % filename
        
        lines = f.split("}")
        group = [l.split("~") for l in lines]

        # add group name
        pairwise_CDR3_distances.append(filename)
        pairwise_templated_distances.append(filename)
        min_CDR3_distances.append(filename)
        min_templated_distances.append(filename)

        pairwise_CDR3_similarities.append(filename)
        pairwise_templated_similarities.append(filename)
        min_CDR3_similarities.append(filename)
        min_templated_similarities.append(filename)

        d1, d2, d3, d4, s1, s2, s3, s4 = calc_distances_on_group(group, Q_cutoff)

        pairwise_CDR3_distances.extend(d1)
        pairwise_templated_distances.extend(d2)
        min_CDR3_distances.extend(d3)
        min_templated_distances.extend(d4)

        pairwise_CDR3_similarities.extend(s1)
        pairwise_templated_similarities.extend(s2)
        min_CDR3_similarities.extend(s3)
        min_templated_similarities.extend(s4)

        result = (pairwise_CDR3_distances, pairwise_templated_distances, min_CDR3_distances, min_templated_distances,
                  pairwise_CDR3_similarities, pairwise_templated_similarities, min_CDR3_similarities, min_templated_similarities)

    return result

def calc_distances_on_group(group, Q_cutoff):

    N = len(group)

    if N == 1: return [], [], [], [], [], [], [], []

    # initialize lists of distances and similarities
    pairwise_CDR3_distances = []
    pairwise_templated_distances = []
    pairwise_CDR3_similarities = []
    pairwise_templated_similarities = []

    # initialize distance and similarity matrices
    CDR3_distances = np.zeros((N, N), np.int16)
    CDR3_distances.fill(32000)

    templated_distances = np.zeros((N, N), np.int16)
    templated_distances.fill(32000)

    CDR3_similarities = np.zeros((N, N), np.float16)
    CDR3_similarities.fill(-np.inf)

    templated_similarities = np.zeros((N, N), np.float16)
    templated_similarities.fill(-np.inf)

    # calculate distances and similarities
    for i, x in enumerate(group):
        for j, y in enumerate(group):

            if i >= j: continue

            CDR3_distance, CDR3_similarity = similarity2seqs_CDR3(x, y, Q_cutoff)
            templated_distance, templated_similarity = similarity2seqs_templated(x, y, Q_cutoff)

            CDR3_distances[i,j] = CDR3_distance
            CDR3_distances[j,i] = CDR3_distance

            templated_distances[i,j] = templated_distance
            templated_distances[j,i] = templated_distance

            CDR3_similarities[i,j] = CDR3_similarity
            CDR3_similarities[j,i] = CDR3_similarity

            templated_similarities[i,j] = templated_similarity
            templated_similarities[j,i] = templated_similarity

            pairwise_CDR3_distances.append(CDR3_distance)
            pairwise_templated_distances.append(templated_distance)

            pairwise_CDR3_similarities.append(CDR3_similarity)
            pairwise_templated_similarities.append(templated_similarity)

    min_CDR3_distances = list(np.amin(CDR3_distances, axis=0))
    min_templated_distances = list(np.amin(templated_distances, axis=0))

    max_CDR3_similarities = list(np.amax(CDR3_similarities, axis=0))
    max_templated_similarities = list(np.amax(templated_similarities, axis=0))

    result = (pairwise_CDR3_distances, pairwise_templated_distances, min_CDR3_distances, min_templated_distances,
              pairwise_CDR3_similarities, pairwise_templated_similarities, max_CDR3_similarities, max_templated_similarities)

    return result

def similarity2seqs_CDR3(x, y, Q_cutoff):

    _, x_abundance, _, _, x_seq, x_qual, _, _, _, _, _, _ = tuple(x)
    _, y_abundance, _, _, y_seq, y_qual, _, _, _, _, _, _ = tuple(y)

    d, s = similarity2seqs_same_length(x_abundance, x_seq, x_qual,
                                       y_abundance, y_seq, y_qual,
                                       Q_cutoff)

    return d, s

def similarity2seqs_templated(x, y, Q_cutoff):

    # before CDR3
    _, x_abundance, _, _, _, _, _, _, x_seq, x_qual, _, _ = tuple(x)
    _, y_abundance, _, _, _, _, _, _, y_seq, y_qual, _, _ = tuple(y)

    max_len_before_CDR3 = max([len(x_seq), len(y_seq)])
    x_seq.rjust(max_len_before_CDR3, "-")
    y_seq.rjust(max_len_before_CDR3, "-")

    d1, s1 = similarity2seqs_same_length(x_abundance, x_seq, x_qual,
                                         y_abundance, y_seq, y_qual,
                                         Q_cutoff)

    # after CDR3
    _, x_abundance, _, _, _, _, _, _, _, _, x_seq, x_qual = tuple(x)
    _, y_abundance, _, _, _, _, _, _, _, _, y_seq, y_qual = tuple(y)

    max_len_after_CDR3 = max([len(x_seq), len(y_seq)])
    x_seq.ljust(max_len_after_CDR3, "-")
    y_seq.ljust(max_len_after_CDR3, "-")

    d2, s2 = similarity2seqs_same_length(x_abundance, x_seq, x_qual,
                                         y_abundance, y_seq, y_qual,
                                         Q_cutoff)


    # weight by length
    s = s1 * (float(max_len_before_CDR3)/float(max_len_before_CDR3 + max_len_after_CDR3)) + s2 * (float(max_len_after_CDR3)/float(max_len_before_CDR3 + max_len_after_CDR3))

    d = d1 + d2

    return d, s

def similarity2seqs_ungapped_alignment(cluster_abundance, cluster_seq, cluster_qual,
                                       candidate_abundance, candidate_seq, candidate_qual,
                                       Q_cutoff):

    similarities = []

    # Set offsets
    ungapped_alignment_offset_window_size = 5

    offsets = range(-ungapped_alignment_offset_window_size, ungapped_alignment_offset_window_size + 1)

    for o in [0, -4, 4]:
        if o in offsets:
            offsets.remove(o)

    offsets = [0, 4, -4] + offsets # test these offsets first because they are most likely (with 8 and 12 nt barcodes)

    # Do alignment for each offset
    for offset in offsets:

        if offset >= 0:

            seq1 = cluster_seq
            qual1 = cluster_qual
            abundance1 = cluster_abundance

            seq2 = candidate_seq
            qual2 = candidate_qual
            abundance2 = candidate_abundance

        else:

            offset = -offset

            seq1 = candidate_seq
            qual1 = candidate_qual
            abundance1 = candidate_abundance

            seq2 = cluster_seq
            qual2 = cluster_qual
            abundance2 = cluster_abundance

        # Slide sequences and truncate to same length
        s1 = seq1[0:min(len(seq1), len(seq2)-offset)]
        q1 = qual1[0:min(len(seq1), len(seq2)-offset)]

        s2 = seq2[offset:min(len(seq1), len(seq2)-offset)+offset]
        q2 = qual2[offset:min(len(seq1), len(seq2)-offset)+offset]

        similarity = similarity2seqs_same_length(abundance1, s1, q1, abundance2, s2, q2, Q_cutoff)

        similarities.append(similarity)

    return max(similarities)

def similarity2seqs_same_length(cluster_abundance, cluster_seq, cluster_qual,
                                candidate_abundance, candidate_seq, candidate_qual,
                                Q_cutoff):

    if len(cluster_seq) == 0: return 0, 1.0

    mismatches = 0

    for s_cluster, s_candidate, q_cluster, q_candidate in zip(cluster_seq, candidate_seq, cluster_qual, candidate_qual):

        if s_candidate != s_cluster and s_candidate != "-" and s_cluster != "-":

            if candidate_abundance > 1 and cluster_abundance > 1:

                # if abundance of the candidate > 1, then we trust that it's a real sequence, and thus a real mismatch
                mismatches += 1

            elif candidate_abundance == 1 and cluster_abundance > 1:

                Q_candidate = ord(q_candidate) - 33
                if Q_candidate > Q_cutoff: mismatches += 1

            elif candidate_abundance > 1 and cluster_abundance == 1:

                Q_cluster = ord(q_cluster) - 33
                if Q_cluster > Q_cutoff: mismatches += 1

            else:
                
                Q_candidate = ord(q_candidate) - 33
                Q_cluster = ord(q_cluster) - 33
                if Q_cluster > Q_cutoff and Q_candidate > Q_cutoff: mismatches += 1

        else:

            # match
            # if the matching base is low quality, then we count it as a mismatch anyway

            if candidate_abundance == 1 and cluster_abundance > 1:

                Q_candidate = ord(q_candidate) - 33
                if Q_candidate <= Q_cutoff: mismatches += 1

            elif candidate_abundance > 1 and cluster_abundance == 1:

                Q_cluster = ord(q_cluster) - 33
                if Q_cluster <= Q_cutoff: mismatches += 1

            else:

                Q_candidate = ord(q_candidate) - 33
                Q_cluster = ord(q_cluster) - 33
                if Q_cluster <= Q_cutoff or Q_candidate <= Q_cutoff: mismatches += 1

    d = float(mismatches) / float(len(cluster_seq))
    similarity = 1.0 - d

    return mismatches, similarity

def write_list(L, outfile):

    out = open(outfile, 'w')

    for val in L:
        out.write(str(val) + "\n")

    out.close()

    return None

if __name__ == "__main__":
    main(sys.argv)
