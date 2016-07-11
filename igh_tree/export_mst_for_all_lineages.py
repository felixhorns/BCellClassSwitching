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

    start_time_mst = time.time()

    export_msts(infile_query, field_dict, msts_dir, wdir_original)

    print "Done!!"
    print "export_mst_for_lineages elapsed time (wall clock):", time.time() - start_time_mst

def export_msts(infile, field_dict, msts_dir, output_dir):

    # Open and write headers into output files
    output_nodes = output_dir + "/mst.all_lineages.nodes.cytoscape"
    out_nodes = open(output_nodes, 'w')
    out_header_vals = ["seq_id", "uid", "isotype", "abundance", "color", "size", "lineage_uid"]
    out_nodes.write("\t".join(map(str, out_header_vals)) + "\n")

    output_edges = output_dir + "/mst.all_lineages.edges.cytoscape"
    out_edges = open(output_edges, 'w')
    out_header_vals = ["source", "type", "target", "distance"]
    out_edges.write("\t".join(map(str, out_header_vals)) + "\n")

    id_offset = 0

    with open(infile, 'rU') as f:
        for line in f:
            lineage_uid = line.rstrip().split("\t")[field_dict["lineage_uid"]]
            id_offset = export_mst(lineage_uid, msts_dir, out_nodes, out_edges, id_offset)

    out_nodes.close()
    out_edges.close()

    return None

def export_mst(lineage_uid, msts_dir, out_nodes, out_edges, id_offset):

    mst_file = msts_dir + "/" + lineage_uid + ".mst"

    with open(mst_file, 'rU') as f:
        mst = json.load(f)

    nodes = set()

    for parent, children in mst.items():
        for child, dist in children.items():

            source = int(parent.split(",")[0]) + id_offset
            target = int(child.split(",")[0]) + id_offset

            source_isotype = parent.split(",")[2]
            if source_isotype == "germline":
                t = "g"
            else:
                t = "x"

            target_isotype = child.split(",")[2]

            if not (source_isotype == target_isotype and int(dist) == 0):
                vals = [source, t, target, int(dist)]
                out_edges.write("\t".join(map(str, vals)) + "\n")

            nodes.add(parent)
            nodes.add(child)

    next_id_offset = write_nodes(nodes, out_nodes, id_offset, lineage_uid)

    return next_id_offset

def write_nodes(nodes, out_nodes, id_offset, lineage_uid):

    seq_ids = []

    for node in nodes:

        vals = node.split(",")
        
        seq_id = int(vals[0]) + id_offset
        uid = vals[1]
        isotype = vals[2]
        abundance = vals[3]

        color = isotype_to_color[isotype]
        size = abundance_to_size(abundance)

        out_vals = [seq_id, uid, isotype, abundance, color, size, lineage_uid]
        out_nodes.write("\t".join(map(str, out_vals)) + "\n")

        seq_ids.append(seq_id)

    next_id_offset = max(seq_ids) + 1

    return next_id_offset

def abundance_to_size(a):
    try:
        return 5 * np.log10(float(a)) + 10
    except:
        return 5

if __name__ == "__main__":
    main(sys.argv)
