# Visualization of sequencing reads.
# Samples reads, makes alignment, builds phylogeny.
#
# Input: fastq file containing reads
# Output: figure of alignment with phylogeny with similar reads clustered

import sys
import os
import random
from ete2 import PhyloTree, TreeStyle

def file_len(fileName):
    i = 0
    for line in open(fileName): i+=1
    return i

from Bio.Align.Applications import MuscleCommandline
def align_muscle(infile_name, outfile_name, gapopen=-1000.0):
    muscle_program='/datastore/dcroote/resources/muscle3.8.31'
    cline = MuscleCommandline(muscle_program, input=infile_name, out=outfile_name, gapopen=gapopen)
    cline()
    return outfile_name

from cogent import LoadSeqs, DNA
from cogent.app.fasttree import build_tree_from_alignment as fasttree_build_tree
def build_tree_FT(aligned_fasta_file_name):
    # Builds approximate maximum likelihood tree using aligned sequences
    aln = LoadSeqs(aligned_fasta_file_name, format='fasta')
    tree = fasttree_build_tree(aln, DNA)
    return tree, aln

def get_label(s):
    """ Splits Illumina read name into a sensible label. """
    flow_cell_lane_id_x_y = s.split(" ")[0].split(":")[-4:]
    # index_seq = s.split(" ")[1].split(":")[-1]
    # label_items = flow_cell_lane_id_x_y + [index_seq]
    label_items = flow_cell_lane_id_x_y
    label = "_".join(label_items)
    return label

def truncate_first_n_bases(s, n): return s[n:]

verbose = True

# infile = "/Users/lime/Data/research/quake/Bcell_class_switching_prelim/fastq/141009_M00361_0225_000000000-ABC22/018-091_V1_70__IL5463-1/test.fastq"
# infile = "/Users/lime/Data/research/quake/Bcell_class_switching_prelim/fastq/141009_M00361_0225_000000000-ABC22/018-091_V1_70__IL5463-1/018-091_V1_70__IL5463-1_ATCACG-TAGTGC_L001_R1_001.fastq"
# infile = "/Users/lime/Data/research/quake/Bcell_class_switching_prelim/fastq/141009_M00361_0225_000000000-ABC22/018-091_V1_70__IL5463-1/018-091_V1_70__IL5463-1_ATCACG-TAGTGC_L001_R1_001.fastq"

infile = sys.argv[1]

N = 100 # number of reads to sample

outfile = infile + ".read_alignment.png"

adapters = {"D701D712adapterL":"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "D701D712adapterR":"ATCTCGTATGCCGTCTTCTGCTTG",
            "D501D508adapterL":"AATGATACGGCGACCACCGAGATCTACAC",
            "D501D508adapterR":"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"}

adapter_reverse_complements = {"D701D712adapterLRC":"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC",
                               "D701D712adapterRRC":"CAAGCAGAAGACGGCATACGAGAT",
                               "D501D508adapterLRC":"GTGTAGATCTCGGTGGTCGCCGTATCATT",
                               "D501D508adapterRRC":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}

adapters_read1 = {"D701D712adapterL":"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                  "D701D712adapterR":"ATCTCGTATGCCGTCTTCTGCTTG"}

adapters_read2 = {"D501D508adapterLRC":"GTGTAGATCTCGGTGGTCGCCGTATCATT",
                  "D501D508adapterRRC":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}

known_sequences = adapters_read1

### Load read sequences
if verbose: print "loading read sequences..."

seqs = []
names = []

L = file_len(infile)
i = 0

with open(infile, 'rU') as f:
    while i < L:
        name = f.readline().rstrip()
        seq = f.readline().rstrip()
        f.readline()
        quality = f.readline().rstrip()

        seqs.append(truncate_first_n_bases(seq, 8))
        names.append(name)

        i += 4

### Sample reads
if verbose: print "sampling reads..."

P = float(N) / float(len(seqs)) # probability of keeping read

temp_file_name = os.path.splitext(infile)[0] + ".subsampled.fa"

sampled_seqs = []
sampled_names = []
j = 0

with open(temp_file_name, 'w') as f:
    for seq, name in zip(seqs, names):
        if random.random() < P:
            sampled_seqs.append(seq)
            sampled_names.append(name)
            f.write(">"+str(j)+"\n") # running number
            # f.write(">"+get_label(name)+"\n")
            f.write(seq+"\n")
        j += 1

### Add known sequences (adapters)
with open(temp_file_name, 'a') as f:
    for name, seq in known_sequences.items():
        f.write(">"+str(name)+"\n")
        f.write(seq+"\n")

### Align
if verbose: print "aligning..."
aln_file_name = os.path.splitext(temp_file_name)[0] + ".afa"
align_muscle(temp_file_name, aln_file_name, gapopen=-1000.0)

### Build tree
if verbose: print "building tree..."
tree, aln = build_tree_FT(aln_file_name)
### Show in pretty format
pretty_tree = PhyloTree(str(tree), alignment=aln_file_name, alg_format="fasta")

pretty_tree.ladderize()

ts = TreeStyle()

pretty_tree.render(outfile, tree_style=ts)

### Clean up your mess
os.remove(temp_file_name)
os.remove(aln_file_name)

### TODO
# highlight adapter rows
# root on adapter?

# tweak alignment parameters so it's better
# understand muscle alignment score. how can we rationally change parameters to improve?

# consider NJ. pairwise distances. how would you overlay alignment?
# would need pairwise alignment and distance for N^2 ~ 10,000 sequences
# what we have is simple and achieves the desired outcome

