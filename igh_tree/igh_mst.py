"""
Stuff
"""

import os
import numpy as np
import subprocess
import uuid
import time
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline


def find_mst_from_alignment(seqs, isotypes, uids, root, wdir, path_to_edmonds_bin, constrain_by_isotype=True):
    """ Finds minimum spanning tree starting with sequence alignment. """

    ### Initialize temporary files
    u = str(uuid.uuid4())
    temp_graph_file_name = wdir+'/mst.' + u + ".graph.temp"
    temp_mst_output_file_name = wdir+'/mst.' + u + ".mst.temp"

    ### Calculate distance matrix
    d = dist_dna(seqs)

    ### Constrain possible ancestral relationships using isotypes
    isotype_to_order_map = _load_isotype_to_order_map(constrain_by_isotype)
    c = _get_constraint_matrix(isotypes, isotype_to_order_map)
    d_constrained = np.multiply(c, d)
    d_constrained = _remove_negative_zero_edges(c, d_constrained)

    ### Find MST using Edmond's algorithm
    _write_graph_to_csv(d_constrained, root, temp_graph_file_name)
    mst = edmonds(temp_graph_file_name, temp_mst_output_file_name, path_to_edmonds_bin)

    ### Clean up
    os.remove(temp_graph_file_name)
    os.remove(temp_mst_output_file_name)

    return mst, root, d_constrained, c, d, seqs

def find_mst(seqs, isotypes, germline_seq, wdir, path_to_muscle_bin, path_to_edmonds_bin, constrain_by_isotype=True, root_on_germline=True):
    """ Finds minimum spanning tree on aligned sequences constrained by isotypes. """

    ### Initialize temporary files
    u = str(uuid.uuid4())
    temp_fasta_file_name = wdir+'/mst.' + u + ".fa.temp"
    temp_aln_file_name = wdir+'/mst.' + u + ".afa.temp"
    temp_graph_file_name = wdir+'/mst.' + u + ".graph.temp"
    temp_mst_output_file_name = wdir+'/mst.' + u + ".mst.temp"

    # Write sequences to temporary fasta file for alignment
    temp_fasta_file = open(temp_fasta_file_name, 'w')

    for i, (seq, isotype) in enumerate(zip(seqs, isotypes)):
        temp_fasta_file.write('>' + str(i+1).zfill(6) + ',' + str(isotype) + '\n')
        temp_fasta_file.write(seq + '\n')

    temp_fasta_file.write('>' + str(0).zfill(6) + ',' + 'germline' + '\n')
    temp_fasta_file.write(germline_seq + '\n')

    temp_fasta_file.close()

    ### Align sequences
    align_muscle(temp_fasta_file_name, temp_aln_file_name, path_to_muscle_bin)
    sort_fasta_file_by_id(temp_aln_file_name)
    aln = AlignIO.read(temp_aln_file_name, 'fasta')
    aligned_seqs = aln_to_seq_list(aln)

    ### Calculate distance matrix
    d = dist_dna(aligned_seqs, count_gaps=True)

    ### Choose root
    if root_on_germline:
        # root on germline
        isotypes_with_root = list(isotypes)
        isotypes_with_root.insert(0, 'germline')
        root = 0
    else:
        # root on sequence nearest to germline
        root = _find_root_mst(d, isotypes, isotype_to_order)
        d = np.delete(d, 0, 0) # delete distances to germ line
        d = np.delete(d, 0, 1)
        isotypes_with_root = list(isotypes)
        isotypes_with_root[root] = 'root'

    ### Constrain possible ancestral relationships using isotypes
    isotype_to_order_map = _load_isotype_to_order_map(constrain_by_isotype)
    c = _get_constraint_matrix(isotypes_with_root, isotype_to_order_map)
    d_constrained = np.multiply(c, d)
    d_constrained = _remove_negative_zero_edges(c, d_constrained)

    ### Find MST using Edmond's algorithm
    _write_graph_to_csv(d_constrained, root, temp_graph_file_name)
    mst = edmonds(temp_graph_file_name, temp_mst_output_file_name, path_to_edmonds_bin)

    ### Clean up
#    os.remove(temp_fasta_file_name)
#    os.remove(temp_aln_file_name)
#    os.remove(temp_graph_file_name)
#    os.remove(temp_mst_output_file_name)

    return mst, root, isotypes_with_root, d_constrained, c, d, seqs

def align_muscle(infile_name, outfile_name, path_to_muscle_bin, gap_open_score=float(-1000.0)):
    cmd = path_to_muscle_bin + " -in " + infile_name + " -out " + outfile_name + " -gapopen " + str(gap_open_score) + " -diags -maxiters 2 -quiet -verbose"
    start_time = time.time()
    #print start_time, "Aligning sequences using MUSCLE..."
    #print cmd
    subprocess.call(cmd, shell=True)
    #print "Done!!"
    #print "Elapsed time (wall clock):", time.time() - start_time
    return outfile_name

def sort_fasta_file_by_id(fasta_file_name):
    sortedList = []
    with open(fasta_file_name,'r') as f:
        l = SeqIO.parse(f,'fasta')
        sortedList = [f for f in sorted(l, key=lambda x: x.id)]
    with open(fasta_file_name,'w') as out:
        SeqIO.write(sortedList, out, 'fasta')
    return None

def aln_to_seq_list(aln):
    return [str(a.seq) for a in aln]

def dist_dna_2seqs(x, y, count_gaps=False):
    x = x.upper()
    y = y.upper()
    if count_gaps:
        d = sum(ele_x != ele_y and not (ele_x == 'N' or ele_y == 'N') for ele_x, ele_y in zip(x,y))
    else:
        d = sum(ele_x != ele_y and not (ele_x == '-' or ele_y == '-' or ele_x == 'N' or ele_y == 'N') for ele_x, ele_y in zip(x,y))
    return d

def dist_dna(seqs, count_gaps=False):
    #print "Calculating pairwise distances between aligned sequences..."
    N = len(seqs)
    d = np.zeros((N,N))
    for i, ele_1 in enumerate(seqs):
        for j, ele_2 in enumerate(seqs):
            if j >= i:
                break
            dist = dist_dna_2seqs(ele_1, ele_2, count_gaps)
            d[i,j] = dist
            d[j,i] = dist
    return d

def _find_root_mst(d, isotypes, isotype_to_order):
    """ Finds the best root for a graph, given a
    distance matrix d and the isotypes of the sequences.
    The root found is the one with the minimum distance
    to the germ line sequence, which is assumed to be
    given by d[0,i] for sequence i, among the sequences
    with the minimum order isotype.  For example,
    among sequences with isotypes [IgM, IgA, IgM], and
    distances [7, 2, 5], the root is 2.
    """
    orders = [isotype_to_order[i] for i in isotypes]
    min_order = min(orders)
    candidate_roots = [int(i == min_order) for i in orders]
    candidate_root_distances_to_germline = [dist if k > 0 else -1 for k, dist in zip(candidate_roots, d[0][1:])]
    root = np.where(candidate_root_distances_to_germline == min([i for i in candidate_root_distances_to_germline if i >= 0.]))[0][0]
    return root

def _load_isotype_to_order_map(constrain_by_isotype):

    if constrain_by_isotype:

        m = {'germline':{'germline':1, 'root':1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'root':{'germline':-1, 'root':1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgM':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgD':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgA':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':-1, 'IgG1':-1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgE':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':-1, 'IgA':1, 'IgE':1, 'IgG3':-1, 'IgG1':-1, 'IgA1':-1, 'IgG2':-1, 'IgG4':-1, 'IgA2':1},
             'IgG3':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG1':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':-1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgA1':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':-1, 'IgG1':-1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG2':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':-1, 'IgG1':-1, 'IgA1':-1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG4':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':-1, 'IgG1':-1, 'IgA1':-1, 'IgG2':-1, 'IgG4':1, 'IgA2':1},
             'IgA2':{'germline':-1, 'root':-1, 'IgM':-1, 'IgD':-1, 'IgG':-1, 'IgA':1, 'IgE':-1, 'IgG3':-1, 'IgG1':-1, 'IgA1':-1, 'IgG2':-1, 'IgG4':-1, 'IgA2':1}}

    else:

        m = {'germline':{'germline':1, 'root':1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'root':{'germline':-1, 'root':1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgM':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgD':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgA':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgE':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG3':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG1':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgA1':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG2':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgG4':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1},
             'IgA2':{'germline':-1, 'root':-1, 'IgM':1, 'IgD':1, 'IgG':1, 'IgA':1, 'IgE':1, 'IgG3':1, 'IgG1':1, 'IgA1':1, 'IgG2':1, 'IgG4':1, 'IgA2':1}}
    
    return m

def _get_constraint_matrix(isotypes, isotype_to_order_map):
    """ Gets constraint matrix for a list
    of isotypes.  The constraint matrix states
    whether sequence i can be an ancestor of sequence j.
    If i can be an ancestor of j, then c[i,j] = 1,
    otherwise c[i,j] = -1.
    """

    N = len(isotypes)
    c = np.zeros((N,N))

    for i, ele_i in enumerate(isotypes):
        for j, ele_j in enumerate(isotypes):

            if i == j: continue

            c[i,j] = isotype_to_order_map[ele_i][ele_j]
            c[j,i] = isotype_to_order_map[ele_j][ele_i]

    return c

def _remove_negative_zero_edges(c, d):
    # Edge should be removed due to constraint, but weight is zero.
    # So it will not be set to negative by c * d.
    # We correct this by setting such cases equal to -1.
    assert c.shape == d.shape
    for i in range(0, c.shape[0]):
        for j in range(0, c.shape[1]):
            if c[i][j] == -1 and d[i][j] == 0:
                d[i][j] = -1
    return d

def _write_graph_to_csv(d, root, graph_file_name):
    """ Writes a graph to plain text format, which can be read
    by the Edmond's algorithm program.  Takes a distance matrix d
    containing the edge weights, with the weight on the edge from
    (source, destination) is d[source, destination].  The root of
    the graph is an integer in the range [0,N-1], where N is the
    number of vertices in the graph. """

    graph_file = open(graph_file_name,'w')
    graph_file.write(str(d.shape[0])+'\n') # First line is the number of vertices in the graph N
    graph_file.write(str(root))     # Second line is the root for the directed MST

    # Each remaining line is a single edge in the format
    # source, destination, weight
    for i in range(0, d.shape[0]):
        for j in range(0, d.shape[1]):
            if i != j and d[i,j] >= 0:
                graph_file.write('\n'+str(i)+"\t"+str(j)+"\t"+str(d[i,j]))

    graph_file.close()
    return graph_file_name

def edmonds(graph_file_name, result_file_name, edmonds_bin_path):
    """ Finds minimum spanning tree (MST) of graph specified in the
    plain text file called graph_file_name.  MST is written to
    result_file_name, then is parsed into a dictionary structure,
    which is returned.  Calls executable at edmonds_bin_path to find
    the MST.
    """
    result_file_name = _edmonds(graph_file_name, result_file_name, edmonds_bin_path) # Run Edmond's algorithm to find mst

    # Read result file and create graph
    mst = {}
    for line in open(result_file_name, 'r'):
        i, j, d = line.rstrip().split('\t')
        i, j, d = int(i), int(j), float(d)

        if i in mst:
            mst[i][j] = d
        else:
            mst[i] = {j:d}

    return mst

def _edmonds(graph_file_name, result_file_name, edmonds_bin_path):
    """ Helper function that calls Edmond's algorithm executable.
    """
    cmd = " ".join([edmonds_bin_path, graph_file_name, result_file_name])
    subprocess.call(cmd, shell=True)
    return result_file_name

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

