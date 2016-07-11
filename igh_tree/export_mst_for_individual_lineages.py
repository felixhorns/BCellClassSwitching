import sys
import os
import time
import json
import numpy as np

global verbose
verbose = True

isotype_to_color_file = "/datastore/rfhorns/resources/immune/isotype_to_color_dict.json"
with open(isotype_to_color_file, 'rU') as f:
    isotype_to_color = json.load(f)

def main(argv):

    ### Parse arguments
    infile_query = argv[1]
    infile_query_field_dict = argv[2] # used to get original working dir (not scratch)

    wdir = os.path.dirname(infile_query)
    wdir_original = os.path.dirname(infile_query_field_dict)

    with open(infile_query_field_dict, 'rb') as f:
        field_dict = json.load(f)

    ### Set output file names
    msts_dir = wdir_original + "/msts/"
    msts_dir_scratch = wdir + "/msts/"

    start_time_mst = time.time()

    export_msts(infile_query, field_dict, msts_dir_scratch)

    print "Done!!"
    print "export_mst_for_lineages elapsed time (wall clock):", time.time() - start_time_mst

def export_msts(infile, field_dict, output_dir):

    with open(infile, 'rU') as f:
        for line in f:
            lineage_uid = line.rstrip().split("\t")[field_dict["lineage_uid"]]
            export_mst(lineage_uid, output_dir)

    return None

def export_mst(lineage_uid, output_dir):

    mst_file = output_dir + "/" + lineage_uid + ".mst"

    with open(mst_file, 'rU') as f:
        mst = json.load(f)

    output_edges = mst_file + ".edges.cytoscape"
    out_edges = open(output_edges, 'w')

    out_header_vals = ["source", "type", "target", "distance"]
    out_edges.write("\t".join(map(str, out_header_vals)) + "\n")

    nodes = set()

    for parent, children in mst.items():
        for child, dist in children.items():

            source = parent.split(",")[0]
            target = child.split(",")[0]

            source_isotype = parent.split(",")[2]
            if source_isotype == "germline":
                t = "g"
            else:
                t = "x"

            vals = [source, t, target, int(dist)]
            out_edges.write("\t".join(map(str, vals)) + "\n")

            nodes.add(parent)
            nodes.add(child)

    out_edges.close()

    output_nodes = mst_file + ".nodes.cytoscape"
    write_nodes(nodes, output_nodes, lineage_uid)

    return None

def write_nodes(nodes, output_nodes, lineage_uid):

    out_nodes = open(output_nodes, 'w')

    out_header_vals = ["seq_id", "uid", "isotype", "abundance", "color", "size", "lineage_uid"]
    out_nodes.write("\t".join(map(str, out_header_vals)) + "\n")

    for node in nodes:

        vals = node.split(",")
        
        seq_id = vals[0]
        uid = vals[1]
        isotype = vals[2]
        abundance = vals[3]

        color = isotype_to_color[isotype]
        size = abundance_to_size(abundance)

        out_vals = [seq_id, uid, isotype, abundance, color, size, lineage_uid]
        out_nodes.write("\t".join(map(str, out_vals)) + "\n")

    out_nodes.close()

    return None

def abundance_to_size(a):
    try:
        return 5 * np.log10(float(a)) + 10
    except:
        return 5

if __name__ == "__main__":
    main(sys.argv)
