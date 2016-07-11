import sys
import os

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
from matplotlib.colors import LogNorm
from matplotlib.ticker import AutoMinorLocator
from pylab import *

def main(argv):

    infile_clustering_inputs = argv[1] # clustering_inputs.txt
    infile_CDR3_distances = argv[2] # CDR3_distances.txt
    infile_templated_distances = argv[3] # templated_distances.txt

    # Load identity to germline

    field_dict = {"abundance":1, "isotype": 2, "mut_density":25, "V_germline_identity":26}

    abundances = []
    isotypes = []
    mut_density = []
    V_germline_identity = []

    with open(infile_clustering_inputs, 'rU') as f:
        for line in f:

            input_dir = line.rstrip().split("\t")[1]
            infile = input_dir + "/parsed_igblast_isotypes_quals"

            with open(infile, 'rU') as f2:
                for line2 in f2:

                    vals = line2.rstrip().split("\t")
                    
                    abundances.append(int(vals[field_dict["abundance"]]))
                    isotypes.append(vals[field_dict["isotype"]])
                    mut_density.append(float(vals[field_dict["mut_density"]]))
                    V_germline_identity.append(float(vals[field_dict["V_germline_identity"]]))

    # Load identity

    CDR3_identity = []
    with open(infile_CDR3_distances, 'rU') as f:
        for line in f:
            CDR3_identity.append(float(line.rstrip()))

    templated_identity = []
    with open(infile_templated_distances, 'rU') as f:
        for line in f:
            templated_identity.append(float(line.rstrip()))

    # Make plots

    minorLocator   = AutoMinorLocator(2)

    # Histogram of identity to V germline
    fig, ax = subplots()
    hist(V_germline_identity, bins=50)
    ax.set_yscale("log")
    xlim([0, 1])
    xlabel("Identity to V germline gene")
    ylabel("Sequences")
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.grid(b=True, which='minor', color='grey', linestyle='--')
    ax.xaxis.grid(b=True, which='major', color='grey', linestyle='--')
    ax.yaxis.grid(b=True, which='major', color='grey', linestyle='--')
    savefig("V_germline_identity.pdf", bbox_inches='tight')

    # Histogram of mutation density (includes D if called and J)
    fig, ax = subplots()
    hist(mut_density, bins=50)
    ax.set_yscale("log")
    xlim([0, 1])
    xlabel("Mutation density (V, J, D if called)")
    ylabel("Sequences")
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.grid(b=True, which='minor', color='grey', linestyle='--')
    ax.xaxis.grid(b=True, which='major', color='grey', linestyle='--')
    ax.yaxis.grid(b=True, which='major', color='grey', linestyle='--')
    savefig("mutation_density.pdf", bbox_inches='tight')

    # 2D histogram of CDR3, templated identity
    fig, ax = subplots()
    hist2d(CDR3_identity, templated_identity, bins=50, norm=LogNorm())
    colorbar()
    xlim([0, 1])
    ylim([0, 1])
    xlabel("CDR3 identity")
    ylabel("Templated identity")
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    grid(b=True, which='minor', color='grey', linestyle='--')
    grid(b=True, which='major', color='grey', linestyle='--')
    savefig("CDR3_templated_identity.pdf", bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv)
