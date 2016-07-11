""" Finds clusters of similar sequences.
"""

from __future__ import division
import sys
import os
import time
import zipfile
import numpy as np
from itertools import combinations,product
from scipy.spatial.distance import pdist,squareform,cdist
import scipy.sparse

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
    outfile = output_dir + "/clusters.txt.temp." + str(i_worker)

    # Get input filenames
    filenames = []
    with open(infile_group_worker_assignments, 'rU') as f:
        for i, line in enumerate(f):
            if i == i_worker:
                filenames = line.rstrip().split("\t")

    if filenames == []:
        # If filenames is empty, then this worker has no groups to process
        write_clusters([], outfile)
        quit()

    # Get clustering parameters
    with open(infile_clustering_params, 'rU') as f:
        method = f.readline().rstrip() # CDR3, VDJ, hybrid
        CDR3_similarity_cutoff = float(f.readline().rstrip()) # minimum similarity in CDR3 to be added to the cluster (e.g. 0.9)
        VDJ_similarity_cutoff = float(f.readline().rstrip()) # minimum similarity in VDJ to be added to the cluster (e.g. 0.95)
        templated_similarity_cutoff = float(f.readline().rstrip()) # minimum similarity in templated to be added to the cluster (e.g. 0.95)
        Q_cutoff = int(f.readline().rstrip()) # mismatches with Q scores lower than this value are not counted (Q = -log10(probability of error)) (e.g. 5)

    if method not in ["CDR3", "VDJ", "hybrid"]:
        print 'Error: method must be one of: "CDR3" "VDJ" "hybrid"'
        quit()

    # Clustering
    start_time = time.time()

    clusters = do_clustering(infile_groups, filenames, method, CDR3_similarity_cutoff, VDJ_similarity_cutoff,
                             templated_similarity_cutoff, Q_cutoff)
    write_clusters(clusters, outfile)

    time_elapsed = time.time() - start_time

    print "Wall clock time elapsed", time_elapsed
    print "Done!!"

def do_clustering(infile, filenames, method, CDR3_similarity_cutoff, VDJ_similarity_cutoff,
                  templated_similarity_cutoff, Q_cutoff):

    clusters = []

    zf = zipfile.ZipFile(infile)

    for filename in filenames:

        try:
            f = zf.read(filename)
        except KeyError:
            print 'ERROR: Did not find %s in zip file' % filename
        
        lines = f.split("}")
        group = [l.split("~") for l in lines]

        c = vector_dist(group, method, CDR3_similarity_cutoff,
                                   VDJ_similarity_cutoff, templated_similarity_cutoff, Q_cutoff)
        clusters.extend(c)

    return clusters


def vector_dist(unclustered_seqs, method, CDR3_similarity_cutoff,VDJ_similarity_cutoff, templated_similarity_cutoff, Q_cutoff):

    seqs_data = np.array(unclustered_seqs)

    max_len_before_CDR3 = len(max(seqs_data[:,8],key=len))
    max_len_after_CDR3  = len(max(seqs_data[:,10],key=len))
    max_CDR3_len        = len(max(seqs_data[:,4], key=len))

    before_CDR3 = numeric_array(seqs_data[:,8],seqs_data[:,9],max_len_before_CDR3,Q_cutoff,direction='reverse')
    after_CDR3  = numeric_array(seqs_data[:,10],seqs_data[:,11],max_len_after_CDR3,Q_cutoff)
    CDR3        = numeric_array(seqs_data[:,4],seqs_data[:,5],max_CDR3_len,Q_cutoff)
    
    seq_names = seqs_data[:,0]
    del seqs_data
    
    if method == 'hybrid' and CDR3_similarity_cutoff == templated_similarity_cutoff:
        
        full_seq   = np.hstack((before_CDR3,CDR3,after_CDR3))
        del before_CDR3,CDR3,after_CDR3
    
        if seq_names.shape[0] > 10000:
            print seq_names.shape[0]
            # too many sequences to operate on in memory: operate row-wise
            set_list=[]
            fullsplit = np.array_split(full_seq,100)
            split_lens = [len(x) for x in fullsplit]
            
            del full_seq
            
            for i in range(len(fullsplit)):
                bool_within = pairwise_within_group(fullsplit[i],i,CDR3_similarity_cutoff,split_lens)
                set_list = set_list + bool_within
                
                if i != len(fullsplit)-1:
                    # distance between groups unless this is the last group
                    between_set = pairwise_between_groups(fullsplit,i,split_lens,CDR3_similarity_cutoff)
                    set_list = set_list + between_set

                # collapse sets every round if num seqs > 10000
                if seq_names.shape[0] > 10000:
                    set_list = merge_sets(set_list)
                    
            # otherwise collapse at the end
            if seq_names.shape[0] <= 10000:
                set_list = merge_sets(set_list)

            # replace indices with seq labels
            for cnt,aset in enumerate(set_list):
                aset=list(aset)
                for c,val in enumerate(aset):
                    aset[c] = seq_names[val]
                
                set_list[cnt]=set(aset)

            return set_list

        else:
            # number of sequences is below threshold- operate on in 1 pass
            dist_all = pairwise_within_group(full_seq, -1, CDR3_similarity_cutoff)
            
            d = squareform(dist_all)

            # connected components from a sparse matrix
            n_components, labels = scipy.sparse.csgraph.connected_components(scipy.sparse.csr_matrix(d), directed=False)
            

            clusters_out = [set() for k in xrange(n_components)]
            # add seqs to their respective clusters
            for cnt,n in enumerate(labels):
                clusters_out[n].add(seq_names[cnt])
    
            return clusters_out
    


def numeric_array(seq_array,qual_array,size,Q_cutoff,direction='forward'):
    '''
    Converts sequence array to numeric array for faster pairwise comparison. Variable lengths,
    "N"s, and low quality bases are masked using np.inf
    '''
    
    num_seqs = seq_array.shape[0]

    convert_dict = { 'A':1, 'T':2, 'C':3, 'G':4, 'N':np.inf}
    
    # initialize array with np.inf and replace each base
    # if above quality score threshold
    mod_seq = np.empty((num_seqs,size))
    mod_seq.fill(np.inf)
    
    for row in xrange(num_seqs):
        for cnt,(s,q) in enumerate(zip(seq_array[row],qual_array[row])):
            if ord(q) - 33 < Q_cutoff:
                pass
            else:
                if direction=='reverse':
                    # reverse direction handles case where the sequences 
                    # in before_CDR3 are aligned by their ends, therefore
                    # flip the seqs and pad the tail of the array if the length
                    # of the seq is shorter than the max seq length
                    mod_seq[row,len(seq_array[row])-cnt-1] = convert_dict[s]
                else:
                    mod_seq[row,cnt] = convert_dict[s]   
    return mod_seq

def pairwise_within_group(full_seq,i,CDR3_similarity_cutoff,split_lens=None):
    d_all      = pdist(full_seq,'hamming')
    inf_all    = pdist(np.isinf(full_seq),'hamming')
    finite_all = np.isfinite(full_seq).sum(axis=1)

    # combinations is faster than pdist w/ a lambda
    pdf_all=np.empty((len(finite_all)/2*(len(finite_all)-1)))
    for c,fin in enumerate(combinations(finite_all,2)):
        pdf_all[c] = min(fin)

    norm_dist_all  = ((d_all-inf_all)/pdf_all*full_seq.shape[1]) < (1-CDR3_similarity_cutoff)
    if i==-1:
        # case for single group
        return norm_dist_all
    else:
        # case for split computation
        sets = sets_within(norm_dist_all,sum(split_lens[:i]))
        return sets 
    
def pairwise_between_groups(fullsplit,i_own,split_lens,CDR3_similarity_cutoff):
    i_others = len(split_lens)
    if i_own != i_others-1:
        # if this is not the last group, stack groups being compared against
        fso = np.vstack(fullsplit[i_own+1:i_others])
    else:
        fso = fullsplit[i_others-1]
        
    bool_dis = cdist(fullsplit[i_own],fso,'hamming').flatten()
    bool_inf = cdist(np.isinf(fullsplit[i_own]),np.isinf(fso),'hamming').flatten()
    finite_own = np.isfinite(fullsplit[i_own]).sum(axis=1)
    finite_others = np.isfinite(fso).sum(axis=1)
    
    pdf_all=np.empty(len(bool_dis))
    for c,fin in enumerate(product(finite_own,finite_others)):
        pdf_all[c] = min(fin)

    norm_dist_all  = (bool_dis-bool_inf)/pdf_all*fullsplit[0].shape[1]
    bool_all = norm_dist_all < (1-CDR3_similarity_cutoff)
    
    # given boolean array, find sequences belonging in a cluster, row-wise
    bool_all = bool_all.reshape(fullsplit[i_own].shape[0],fso.shape[0])
    sets=[]
    
    col_offset = sum(split_lens[:i_own])+fullsplit[i_own].shape[0]
        
    for cnt,row in enumerate(bool_all):
        row_offset = cnt+sum(split_lens[:i_own])

        sets_from_group = set(np.add(np.nonzero(row)[0],col_offset))
        sets_from_group.add(row_offset)
        sets.append(sets_from_group)

    return sets

def merge_sets(data):
    """Check each number in our arrays only once, merging when we find
    a number we have seen before.
    """

    bins = range(len(data))
    nums = dict()

    for r, row in enumerate(data):
        for num in row:
            if num not in nums:
                # New number: tag it with a pointer to this row's bin
                nums[num] = r
                continue
            else:
                dest = locatebin(bins, nums[num])
                if dest == r:
                    continue # already in the same bin

                if dest > r:
                    dest, r = r, dest   # always merge into the smallest bin

                data[dest].update(data[r]) 
                data[r] = None
                # Update indices to reflect the move
                bins[r] = dest
                r = dest 

    # Filter out the empty bins
    have = [ m for m in data if m ]
    return have

def locatebin(bins, n):
    """
    Find the bin where list n has ended up: Follow bin references until
    we find a bin that has not moved.
    """
    while bins[n] != n:
        n = bins[n]
    return n

def sets_within(bool_array,global_offset):
    '''
    Given a boolean pairwise distance array, find connected components and return a list of cluster sets. 
    Global offset used to correct sequence numbering for the given group.
    '''
    sq = squareform(bool_array)
    n_components, labels = scipy.sparse.csgraph.connected_components(scipy.sparse.csr_matrix(sq), directed=False)
    clusters_out = [set() for k in xrange(n_components)]
    # add seqs to their respective clusters
    for cnt,n in enumerate(labels):
        clusters_out[n].add(cnt+global_offset)
        
    return clusters_out

def write_clusters(clusters, outfile):

    out = open(outfile, 'w')

    for cluster in clusters:
        out.write(",".join(cluster) + "\n")

    out.close()

    return None

if __name__ == "__main__":
    main(sys.argv)
