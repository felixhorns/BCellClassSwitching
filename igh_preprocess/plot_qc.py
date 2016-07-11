""" Plots QC figures for immune repertoire sequencing IgH pipeline.
"""

import sys
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
from matplotlib import pyplot as plt
from matplotlib import cm

# add to sys path for importing
# markup.py to produce HTML output
sys.path.append('/datastore/dcroote/allergy/gitbranches/master')
import markup

##### Input
infile_indexed_reads_sorted = sys.argv[1]
infile_parsed_igblast  = sys.argv[2]

# Combined_R1.fastq
# Combined_R1.trimmed.fastq
# cutadapt.1.out
# cutadapt.2.out
# consensus_complete.txt
# out.extendedFrags.fastq
# out.notCombined_1.fastq
# out.notCombined_2.fastq
# sequences_abundances.fasta
# losses_parse_igblast.txt
# losses_parse_isotype_blast.txt
# parsed_igblast_isotypes_quals

##### Output
output_dir = os.path.dirname(infile_indexed_reads_sorted) + "/"
# output_statistics.txt
# output_statistics_plotting.txt
# reads_per_molecule_hist.png
# abundance_hist.png
# isotype_counts.png
# V_counts.png
# J_counts.png
# output_statistics
# isotype_counts.txt
# CDR3_lengths.png
# losses_pipeline.png
# index.html 
####

field_dict = {"id":0, "abundance":1, "isotype":2, "V":3, "D":4, "J":5, "region_lengths":18}

# Number of molecules, reads per molecule
reads_per_molecule = {}
isotype_counts_indexed_reads = {}

with open(infile_indexed_reads_sorted, 'rU') as f:
    for line in f:
        values = line.rstrip().split('\t')
        molecular_index = values[1]
        isotype = values[-1]

        try:
            reads_per_molecule[molecular_index] += 1
        except:
            reads_per_molecule[molecular_index] = 1

        try:
            isotype_counts_indexed_reads[isotype] += 1
        except:
            isotype_counts_indexed_reads[isotype] = 1

num_indexed_reads = sum(reads_per_molecule.values())
num_molecules = len(reads_per_molecule.keys())

# Sequence data
abundances = []
V_counts = {}
J_counts = {}
isotype_counts = {}
isotype_counts_abund = {}
CDR3_lengths = []

with open(infile_parsed_igblast, 'rU') as f:
    for line in f:
        values = line.rstrip().split('\t')
        abundance = int(values[field_dict["abundance"]])
        V = values[field_dict["V"]]
        J = values[field_dict["J"]]
        isotype = values[field_dict["isotype"]]
        CDR3_length = int(values[field_dict["region_lengths"]].split(",")[-1])

        abundances.append(abundance)

        try:
            V_counts[V] += 1
        except:
            V_counts[V] = 1

        try:
            J_counts[J] += 1
        except:
            J_counts[J] = 1

        try:
            isotype_counts[isotype] += 1
            isotype_counts_abund[isotype] += abundance
        except:
            isotype_counts[isotype] = 1
            isotype_counts_abund[isotype] = abundance

        CDR3_lengths.append(CDR3_length)

def aggregate_split_losses(infile):
    '''Concatenating losses files e.g. losses_parse_igblast.txt.{n} results in 
    tabular loss data that must be summed:
    parse_igblast.py        parsed_igblast  6083
    parse_igblast.py        no_hits_found   16
    parse_igblast.py        no_CDR3_start   35
        etc ...
    # now it repeats for another split file
    parse_igblast.py        parsed_igblast  6064
    parse_igblast.py        no_hits_found   19 
    parse_igblast.py        no_CDR3_start   25
        etc...    ''' 
    losses_dict = {}
    with open(infile) as f:
        for line in f:
            line = line.strip().split()
            if line[1] not in losses_dict:
                losses_dict[line[1]] = int(line[2])
            else:
                losses_dict[line[1]] += int(line[2])
    return losses_dict
    
### Compile number of reads/molecules/sequences through pipeline and losses at each step
def compile_counts(output_dir):

    # Each count is a 3-ple (rule, description, number)
    counts = []
    losses_plot = []

    # Count reads/molecules/sequences at each step
    # input
    c = count_fastq(output_dir+"/Combined_R1.fastq")
    t = ("input", "Reads", c)
    counts.append(t)
    losses_plot.append(t)

    # trim adapters
    c = count_fastq(output_dir+"/Combined_R1.trimmed.fastq")
    t = ("trim_adapters", "Reads", c)
    counts.append(t)
    losses_plot.append(t)

    trim_adapt_total_cnt=0
    lines = grep("Trimmed:", output_dir+"/cutadapt.1.out")
    for line in lines:
        c = line.split("Trimmed: ")[1].split(" times")[0]
        seq = line.split("Sequence: ")[1].split(";")[0]
        t = ("trim_adapters", "Read 1, adapter " + seq + " trimmed", c)
        counts.append(t)
        trim_adapt_total_cnt+=int(c)

    lines = grep("Trimmed:", output_dir+"/cutadapt.2.out")
    for line in lines:
        c = line.split("Trimmed: ")[1].split(" times")[0]
        seq = line.split("Sequence: ")[1].split(";")[0]
        t = ("trim_adapters", "Read 2, adapter " + seq + " trimmed", c)
        counts.append(t)
        trim_adapt_total_cnt+=int(c)

    t = ("trim_adapters", "Adapter trimming", trim_adapt_total_cnt)
    losses_plot.append(t)

    # index reads
    c = num_indexed_reads
    t = ("index_reads", "Reads", c)
    counts.append(t)

    c = num_molecules
    t = ("index_reads", "Molecules", c)
    counts.append(t)
    losses_plot.append(t)

    c = num_indexed_reads - num_molecules
    t = ("index_reads", "Repetitive barcodes", c)
    losses_plot.append(t)

    # determine consensus
    c = count_fastq(output_dir+"/consensus_complete.txt")
    t = ("determine_consensus", "Molecules", c)
    counts.append(t)
    losses_plot.append(t)

    c = count_fastq(output_dir+"/out.extendedFrags.fastq")
    t = ("determine_consensus", "number of molecules, out.extendedFrags.fastq", c)
    counts.append(t)

    c = count_fastq(output_dir+"/out.notCombined_1.fastq")
    t = ("determine_consensus", "Reads, out.notCombined_1.fastq", c)
    t1 = ("determine_consensus", "Not merged", c)
    counts.append(t)
    losses_plot.append(t1)

    c = count_fastq(output_dir+"/out.notCombined_2.fastq")
    t = ("determine_consensus", "number of reads, out.notCombined_2.fastq", c)
    counts.append(t)

    # combine identical
    cmol_combi = sum_abundances(output_dir+"sequences_abundances.FASTA")
    t = ("combine_identical", "Molecules", cmol_combi)
    counts.append(t)

    cseq_combi = count_fasta(output_dir+"sequences_abundances.FASTA")
    t = ("combine_identical", "Sequences", cseq_combi)
    counts.append(t)
    losses_plot.append(t)

    c=cmol_combi - cseq_combi
    t = ("combine_identical", "Abundant sequences", c)
    losses_plot.append(t)

    # parse_igblast
    # losses dictionary organized by step (key), loss (value)
    loss_parse_ig = aggregate_split_losses(output_dir+"losses_parse_igblast.txt")
    parse_seq_losses=0
    for key in loss_parse_ig:
        if key == "parsed_igblast":
            # value for this key is not losses, but instead total input sequences
            c = loss_parse_ig[key] 
            t = ("parse_igblast", "Sequences", c)
            counts.append(t)
            losses_plot.append(t)
        else:
            c = loss_parse_ig[key]
            t = ("parse_igblast", "Sequences, " + key, c)
            counts.append(t)
            parse_seq_losses+=int(c)

    t = ("parse_igblast", "Igblast losses", parse_seq_losses)
    losses_plot.append(t)

    # parse_isotype_blast
    loss_parse_iso = aggregate_split_losses(output_dir+"losses_parse_isotype_blast.txt")
    for key in loss_parse_iso:
        c = loss_parse_iso[key]
        t = ("parse_isotype_blast", "number of sequences, " + key, c)
        counts.append(t)

    # parse_quality_scores
    c = count_lines(output_dir+"parsed_igblast_isotypes_quals")
    t = ("parse_quality_scores", "Sequences", c)
    counts.append(t)
    losses_plot.append(t)

    # Write statistics to file
    with open(output_dir+"/output_statistics", 'w') as out_statistics:
        for i in counts:
            out_statistics.write("\t".join(map(str, i))+"\n")

    with open(output_dir+"/output_statistics_plotting.txt", "w") as out_stats_plotting:
        for i in losses_plot:
            out_stats_plotting.write("\t".join(map(str, i))+"\n")

    return None

def count_lines(infile):
    count = 0
    with open(infile, 'rU') as f:
        for line in f:
            count += 1
    return count

def count_fastq(infile):
    count = 0
    with open(infile, 'rU') as f:
        for line in f:
            count += 1
    return int(count/4.0)

def count_fasta(infile):
    count = 0
    with open(infile, 'rU') as f:
        for line in f:
            count += 1
    return int(count/2.0)

def grep(s, infile):
    lines = []
    with open(infile, 'rU') as f:
        for line in f:
            if re.search(s, line):
                lines.append(line)
    return lines

def sum_abundances(infile):
    total = 0
    with open(infile, 'rU') as f:
        for line in f:
            if ">" in line:
                abundance = int(line.rstrip().split("_")[-1])
                total += abundance
    return total

compile_counts(output_dir)

### Write statistics to files
with open(output_dir+"/isotype_counts.txt", 'w') as out:
    isotypes = ['IgM', 'IgD', 'IgG3', 'IgG1', 'IgA1', 'IgG2', 'IgG4', 'IgE', 'IgA2']
    out.write('isotype,uniq_seq,total_molec\n')
    for isotype in isotypes:
        try:
            out.write('%s,%s,%s\n' % (isotype,isotype_counts[isotype],isotype_counts_abund[isotype]))
            #out.write(str(isotype_counts[isotype])+'\n')
        except:
            out.write('%s,0,0\n' % isotype)
            #out.write('0\n')

### Make plots
def plot_count_dict(ax, d, yscale_log=False, tick_params_labelsize=12):

    points = [(k, d[k]) for k in d.keys()]
    sorted_points = sorted(points, key=lambda tup:tup[0])
    labels = [tup[0] for tup in sorted_points]
    x = range(1,len(labels)+1)
    y = [tup[1] for tup in sorted_points]

    ax.bar(x, y)

    x_labels = [x_i + 0.5 for x_i in x]
    plt.xticks(x_labels, labels, rotation=90)
    if tick_params_labelsize != 12:
        plt.tick_params(axis='x', which='major', labelsize=tick_params_labelsize)
        plt.tick_params(axis='x', which='minor', labelsize=tick_params_labelsize)
    ax.set_ylabel("Number of sequences")
    ax.set_xlim(min(x)-1, max(x)+1)
    ax.set_ylim(bottom=0.1)
    if yscale_log: ax.set_yscale('log')
    plt.tight_layout()
    return None

def plot_histogram(ax, l, bins=30, xlabel="", ylabel="Density", xscale_log=False, yscale_log=False):
    ax.hist(l, bins=bins, histtype='bar')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_ylim(ax.get_ylim()[0]*0.1, ax.get_ylim()[1]) # does something glitchy to display
    if xscale_log: ax.set_xscale('log')
    if yscale_log: ax.set_yscale('log')
    plt.tight_layout()
    return None

figsize=(6,6)

fig, ax = plt.subplots(1, 1, figsize=figsize)
reads_mols_bins = range(11)+range(11,101,10)+range(101,1001,100)+range(1001,10000,1000)
plot_histogram(ax, reads_per_molecule.values(), xlabel="Reads per molecule", ylabel="Molecules", yscale_log=True,bins=reads_mols_bins,xscale_log=True)
ax.set_ylim(bottom=0.1)  # lower log floor to 0.1 to properly depict histogram bars with height = 10^0 (1)
fig.savefig(output_dir+"/reads_per_molecule_hist.png", dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(1, 1, figsize=figsize)
log_abundances = [np.log10(a) for a in abundances]
plot_histogram(ax, log_abundances, xlabel="Log10 Abundance (molecules per sequence)", ylabel="Sequences", yscale_log=True)
ax.set_ylim(bottom=0.1)  # lower log floor to 0.1 to properly depict histogra
fig.savefig(output_dir+"/abundance_hist.png", dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(1, 1, figsize=figsize)
plot_count_dict(ax, isotype_counts, yscale_log=True)
fig.savefig(output_dir+"/isotype_counts.png", dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(1, 1, figsize=(10,4))
plot_count_dict(ax, V_counts, yscale_log=True, tick_params_labelsize=4)
fig.savefig(output_dir+"/V_counts.png", dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(1, 1, figsize=(6,4))
plot_count_dict(ax, J_counts, yscale_log=True)
fig.savefig(output_dir+"/J_counts.png", dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(1, 1, figsize=figsize)
bins = max(CDR3_lengths) - min(CDR3_lengths)
plot_histogram(ax, CDR3_lengths, bins=bins, xlabel="Length of CDR3 (bp)", ylabel="Sequences", yscale_log=False)
fig.savefig(output_dir+"/CDR3_lengths.png", dpi=300, bbox_inches='tight')

df = pd.io.parsers.read_csv(output_dir+'/output_statistics_plotting.txt',sep='\t',header=None)
dft=pd.DataFrame()

for x in list(set(df[0])):  # dataframe manipulation to turn rows into columns
    b=df[df[0]==x].T         # to enable stacked bar plotting
    b.drop(0,axis=0,inplace=True)
    b.columns=b.iloc[0,:]
    b.drop(1,axis=0,inplace=True)
    b.index=[str(x)]
    dft = pd.concat([dft,b])

# Setting some nice colors
num_cols = len(dft.columns)
cs=cm.Set1(np.arange(num_cols)/10.)
# Plot the pipeline steps in the correct order

newi=['input','trim_adapters','index_reads','determine_consensus','combine_identical','parse_igblast','parse_quality_scores']

fig,ax = plt.subplots(1,1,figsize=(8,5))
dft=dft/1e6  # numbers in millions of reads
colorder=['Reads','Molecules','Sequences','Adapter trimming','Repetitive barcodes','Not merged','Abundant sequences','Igblast losses']
dft = dft[colorder]
ax = dft.ix[list(newi)].plot(kind='bar',stacked=True,ax=ax,color=cs,alpha=0.7)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel('Million reads')
fig.savefig(output_dir+"/losses_pipeline.png", dpi=300, bbox_inches='tight')

##### HTML output 
title = 'QC Plots'
header = '%s' % infile_parsed_igblast
page = markup.page()
page.init(title=title, header=header)
page.br()
paragraphs = ['Pipeline losses','CDR3 Lengths','Isotype counts','J Counts','V counts','Abundance Histogram','Reads per Molecule Histogram']
page.p(paragraphs[0])
page.img(width=800,height=500,src='losses_pipeline.png')
page.p(paragraphs[1])
page.img(width=500,height=500,src='CDR3_lengths.png')
page.p(paragraphs[2])
page.img(width=500,height=500,src='isotype_counts.png')
page.p(paragraphs[3])
page.img(width=600,height=400,src='J_counts.png')
page.p(paragraphs[4])
page.img(width=900,height=360,src='V_counts.png')
page.p(paragraphs[5])
page.img(width=500,height=500,src='abundance_hist.png')
page.p(paragraphs[6])
page.img(width=500,height=500,src='reads_per_molecule_hist.png')

# Write
with open('%s/index.html' % output_dir,'w') as qcout:
    qcout.write('%s' % page)


