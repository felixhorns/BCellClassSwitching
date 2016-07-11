import sys
import os
import numpy as np
import scipy.sparse
from scipy.misc import comb
from itertools import combinations, chain
import math
from scipy.spatial.distance import squareform
import sys
sys.path.append('/datastore/dcroote/resources/scipy-2015-cython-tutorial/exercises/02-typing/01-hamming-distance')
import hamming_cython_solution


import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
from matplotlib import pyplot as plt

plt.style.use('ggplot')

def main(argv):
    infile = argv[1]
    outfile = argv[2]
    output_dir = os.path.dirname(outfile)
    D_cutoff = int(argv[3])
    q_threshold = int(argv[4])

    # normally this step of the pipeline is cripplingly slow due to the pairwise dist calc that 
    # requires comparing every base and every base call quality score
    # A faster version that does not take into account the quality score and includes a cython optimized
    # hamming dist calculation can be used instead
    # Profiling: 71s to 2s
    try:
        fast_pwdist = int(argv[5])
    except:
        fast_pwdist = 0

    out = open(outfile,'w')

    held_lines = []  # list for holding all reads in a barcode
    # (see add_line_to_dict function for explicit definition) 
    
    # for histogram plotting
    pairwise_distances = {}  # key: # of seqs in a barcode, val: list of pairwise distances

    total_lines = get_number_of_lines(infile)

    if total_lines == 0:
        out.close()
        return 0

    with open(infile,'rU') as f:
        
        # read the first line to get current barcode and read, then loop
        line = f.readline().split()
        held_lines = add_line_to_list(held_lines,line)
        current_bar = line[2]+line[3]

        for line_num,line in enumerate(f):

            line=line.split()
            bar = line[2]+line[3]

            if bar == current_bar:
                # if the read belongs to the current barcode, add then continue
                held_lines = add_line_to_list(held_lines,line)

                # case for last line matching previous barcode
                if line_num == total_lines-2:  # -2 from reading in first line outside of loop
                    # compare then write out
                    clustered_bar,pairwise_distances = cluster_within_barcode(held_lines,pairwise_distances,D_cutoff,q_threshold,fast_pwdist)
                    for cluster in clustered_bar:
                        write_new_bar(out,cluster)

            else:

                if len(held_lines) == 1:
                    # current barcode has only 1 read so append 1 to groupID then write out
                    for read in held_lines:
                        read['groupID'] = '%s_%s' % (read['groupID'],1)
                        write_new_bar(out,read)
                
                else:
                    # execute pairwise comparison
                    clustered_bar,pairwise_distances = cluster_within_barcode(held_lines,pairwise_distances,D_cutoff,q_threshold,fast_pwdist)

                    # write lines
                    for cluster in clustered_bar:
                        write_new_bar(out,cluster)

                # move on to next barcode
                held_lines = []
                current_bar = bar
                held_lines = add_line_to_list(held_lines,line)

                # case for last line has a unique barcode
                if line_num == total_lines-2:
                    #write out
                    for read in held_lines:
                        read['groupID'] = '%s_%s' % (read['groupID'],1)
                        write_new_bar(out,read)

    
    out.close()
        
    # Make plots of pairwise distances
    # Commented because this is very slow and causes errors (width or height > 32768 crashes renderer)
#    fig, ax = plot_pairwise_distances(pairwise_distances)
#    fig.savefig("%s/pairwise_distances.png" % output_dir)

    print "Done!!"

def get_number_of_lines(filein):
    with open(filein,'rU') as f:
        c=0
        for x in f:
            c+=1
    return c


def add_line_to_list(lines,line):
    '''adds the current line to the running list of lines in a barcode.
    each element in the list is a dict containing all read information
    '''
    lines.append({'bar1': line[2], 'bar2':line[3], 'seq':line[4], '0':line[0],
                           'groupID':line[1], 'qual':line[5], 'uid':line[6], 'loc':line[7]})
    return lines
    

def cluster_within_barcode(lines,pairwise_distances,D_cutoff,q_threshold,fast_pwdist):
    ''' takes a set of reads belonging to a given barcode and subdivides them into unique barcodes 
        if substantial differences exist between between sequences, as determined by a mismatch 
        difference of *** bases
    '''

    if fast_pwdist == 1:
        D = pairwise_dist_fast(lines)
    else:
        D = pairwise_dist_calc(lines,q_threshold)

    #iu1 = np.triu_indices(len(lines), 1)
    #D_unique = list(D[iu1])

    # no plotting so dictionary generation useless
    #if len(lines) in pairwise_distances:
    #    pairwise_distances[len(lines)].extend(D_unique)
    #else:
    #    pairwise_distances[len(lines)] = D_unique

    # Prune graph by distance cutoff
    C = D <= D_cutoff

    # Find connected components
    csrm = scipy.sparse.csr_matrix(C)
    n_components, labels = scipy.sparse.csgraph.connected_components(scipy.sparse.csr_matrix(C), directed=False)

    seqs_per_new_bar = {}
    for x in labels:
        if x not in seqs_per_new_bar:
            seqs_per_new_bar[x] = 1
        else:
            seqs_per_new_bar[x] += 1

    for cnt,read in enumerate(lines):
        if n_components > 1:
            # if sequences cluster into more than one cluster
            read['groupID'] = '%s_%s_%s' % (read['groupID'],labels[cnt],seqs_per_new_bar[labels[cnt]])
        else:
            # reads cluster together so no subdivision of groupID, but still
            # append number of reads to groupID
            read['groupID'] = '%s_%s' % (read['groupID'],len(lines))

    # sorting reads in the barcode so they're grouped together by determine_consensus
    lines = sorted(lines, key=lambda x: x['groupID'])

    return lines, pairwise_distances


def pairwise_dist_calc(lines,qthresh):
    ''' returns numpy distance matrix from pairwise comparison of reads in dictionary dict_lines
        compares only bases within reads above a quality threshold to ensure that only reads that are 
        truly different are divided into separate barcodes
    '''

    dist_matrix,rows,cols = initialize_dist_matrix(lines)

    for cnt,read in enumerate(combinations(lines,2)):

        seq1 = read[0]['seq']
        seq2 = read[1]['seq']
        qual1 = read[0]['qual']
        qual2 = read[1]['qual']
        
        mismatches=0
        L=0
        for s1, q1, s2, q2 in zip(seq1, qual1, seq2, qual2):
            q1_int = ord(q1) - 33
            q2_int = ord(q2) - 33
            if q1_int > qthresh and q2_int > qthresh:
                L += 1
                if s1 != s2: mismatches += 1

        dist_matrix[rows[cnt],cols[cnt]] = mismatches
        dist_matrix[cols[cnt],rows[cnt]] = mismatches  # symmetric
        
    return dist_matrix


def pairwise_dist_fast(lines):

    num_seqs = len(lines)
    num_comparisons = num_seqs*(num_seqs-1)/2
    dists = np.zeros(num_comparisons)

    for cnt, read in enumerate(combinations(lines,2)):
        dists[cnt] = hamming_cython_solution.hamming_loop(read[0]['seq'],read[1]['seq'])

    return squareform(dists)

        
def initialize_dist_matrix(lines):
    ''' initializes an NxN numpy array where N is the number of sequences and additionally 
        returns the row and column indices for indexing mismatch values
    '''
    num_seqs = len(lines)
    m = np.zeros([num_seqs,num_seqs])
    r,c = np.triu_indices_from(m,1)  #1 is for offset (dont want diagonal)
    return m,r,c


def write_new_bar(outfile,line):
    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (line['0'],line['groupID'], 
                                                  line['bar1'], line['bar2'],
                                                  line['seq'], line['qual'],
                                                  line['uid'], line['loc']))


def plot_pairwise_distances(pairwise_distances):

    num_plots = len(pairwise_distances.keys())
    fig, ax = plt.subplots(num_plots, 1, figsize=(4,num_plots*4))

    for i, (N, distances) in enumerate(pairwise_distances.items()):
        # calculate the number of barcodes for that have the given number of reads
        nbarcodes = len(distances) / (math.factorial(N) / (2*math.factorial(N-2)))
        median_barcode = np.median(distances)
        #shift the bins by +1 to show the 0 bin
        distances = [x+1 for x in distances]

        ax[i].hist(distances, bins=range(int(max(distances))+1))
        ax[i].set_xscale("log")
        ax[i].set_xlabel("Distance")
        ax[i].set_ylabel("Number of pairs")
        ax[i].set_title(str(N) + " reads: " + str(nbarcodes) + " barcodes\nMedian distance: " + str(median_barcode))
    plt.tight_layout()

    return fig, ax


if __name__ == "__main__":
    main(sys.argv)
