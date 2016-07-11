import sys
import os
import json
from itertools import tee, izip_longest
from string import digits
import time
import math

from uid import *

def main(argv):
    
    # Parameters
    infile = argv[1] # sequences_lineages file to read
    infile_original_dir = argv[2] # sequences lineages file original location to get directory where data is located
    path_to_patient_name_to_uid_map = argv[3] # /resources/patient_name_to_uid_map.json
    path_to_patient_data = argv[4] # /resources/patient_data.csv
    path_to_year_visit_str_to_uid_map = argv[5] # /resources/year_visit_str_to_uid_map.json
    path_to_subsampling_p_to_uid_map = argv[6] # /resources/subsampling_p_to_uid_map.json
    path_to_clustering_params_to_uid_map = argv[7] # /resources/clustering_params_to_uid_map.json
    path_to_sequences_lineages_field_dict = argv[8] # resources/sequences_lineages_field_dict.json

    start_time = time.time()
    
    # Get working directory
    wdir = os.path.dirname(infile) # directory containing input files for working
    wdir_original = os.path.dirname(infile_original_dir) # directory where data is located permanently

    # Patient
    patient_name = get_patient_name(wdir)
    patient_uid = get_patient_uid(patient_name, path_to_patient_name_to_uid_map)
    make_patients_sql_input(patient_name, patient_uid, path_to_patient_data, wdir)

    # Subsampling and clustering
    p, replicate, clustering_params = make_libs_clusterings_sql_input(patient_uid,
                                                                      path_to_year_visit_str_to_uid_map,
                                                                      path_to_subsampling_p_to_uid_map,
                                                                      path_to_clustering_params_to_uid_map, wdir)

    # Sequences and lineages
    make_sequences_sequence_strings_quality_scores_lineages_sql_input(infile, patient_uid, patient_name,
                                                                      path_to_year_visit_str_to_uid_map,
                                                                      p, path_to_subsampling_p_to_uid_map, replicate,
                                                                      clustering_params,
                                                                      path_to_clustering_params_to_uid_map,
                                                                      path_to_sequences_lineages_field_dict,
                                                                      wdir, wdir_original)

    # UIDs
    clustering_uid = get_clustering_uid(patient_uid, p, path_to_subsampling_p_to_uid_map, replicate, 
                                        clustering_params, path_to_clustering_params_to_uid_map)
    write_uids(patient_uid, clustering_uid, wdir)

    print "Done!!"
    print "Elapsed time (wall clock):", time.time() - start_time

def load_dict(infile):
    with open(infile, 'rU') as f:
        d = json.load(f)
    return d

def get_patient_name(wdir):
    with open(wdir+"/lib_info_combined.txt") as f:
        f.readline()
        name = f.readline().rstrip()
    return name

def write_uids(patient_uid, clustering_uid, wdir):
    with open(wdir+'/uids.txt', 'w') as out:
        out.write(str(patient_uid)+"\n")
        out.write(str(clustering_uid))
    return None

# Functions to make the SQL input files
def make_patients_sql_input(patient_name, patient_uid, path_to_patient_data, wdir):

    with open(path_to_patient_data, 'rU') as f:
        for line in f:
            vals = line.rstrip().split("\t")
            if vals[0] == patient_name:
                out_vals = vals
                out_vals.insert(0, patient_uid)
                with open(wdir+"/patients.txt.temp", 'w') as out:
                    out.write("\t".join(map(str, out_vals))+"\n")
                return None

    print "Error: patient name", patient_name, "not found in", path_to_patient_data
    quit()

def make_libs_clusterings_sql_input(patient_uid, path_to_year_visit_str_to_uid_map, path_to_subsampling_p_to_uid_map, path_to_clustering_params_to_uid_map, wdir):

    out_libs = open(wdir+"/libs.txt.temp", 'w')

    with open(wdir+"/lib_info_combined.txt", 'rU') as f:

        all_libs = []

        vals = []

        for line in f:
            
            if line in ["\n", "\r\n"]:

                # end of library
                del vals[0] # remove first instance of visit
                subsampling_p = vals[4]
                replicate = vals[-1]
                year_visit_uid = get_year_visit_uid(vals[1], path_to_year_visit_str_to_uid_map)
                uid = get_lib_uid(patient_uid, year_visit_uid, subsampling_p, path_to_subsampling_p_to_uid_map, replicate)

                vals.insert(0, uid)
                vals.insert(1, patient_uid)
                vals.insert(3, year_visit_uid)
                out_libs.write("\t".join(map(str, vals)) + "\n")
                all_libs.append(vals)
                vals = []
            else:
                vals.append(line.rstrip())

    out_libs.close()

    original_num_seqs = sum([ int(vals[8]) for vals in all_libs ])
    target_num_seqs = sum([ float(vals[9]) for vals in all_libs ])
    num_seqs = sum([ int(vals[10]) for vals in all_libs ])

    subsample_p = all_libs[0][7] # assume same p for all visits
    subsample_replicate = all_libs[0][-1] # assume same replicate for all visits

    with open(wdir+"/clustering_params.txt") as f:
        clustering_method = f.readline().rstrip()
        clustering_CDR3_similarity_cutoff = float(f.readline().rstrip())
        clustering_VDJ_similarity_cutoff = float(f.readline().rstrip())
        clustering_templated_similarity_cutoff = float(f.readline().rstrip())
        clustering_Q_cutoff = float(f.readline().rstrip())
        clustering_V_allele_collapse_dict_file = f.readline().rstrip()

    clustering_params = [clustering_method, clustering_CDR3_similarity_cutoff, clustering_VDJ_similarity_cutoff,
                         clustering_templated_similarity_cutoff, clustering_Q_cutoff]

    clustering_uid = get_clustering_uid(patient_uid, subsampling_p, path_to_subsampling_p_to_uid_map, replicate,
                                        clustering_params, path_to_clustering_params_to_uid_map)

    vals = [clustering_uid, patient_uid, subsample_p, subsample_replicate,
            original_num_seqs, target_num_seqs, num_seqs,
            clustering_method, clustering_CDR3_similarity_cutoff, clustering_VDJ_similarity_cutoff,
            clustering_templated_similarity_cutoff, clustering_Q_cutoff, clustering_V_allele_collapse_dict_file]

    with open(wdir+"/clusterings.txt.temp", 'w') as out_clusterings:
        out_clusterings.write("\t".join(map(str, vals)) + "\n")

    return subsample_p, subsample_replicate, clustering_params

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2,s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip_longest(a, b)

def split_year_visit(s):
    year = 0
    visit = 0
    if 'Y' in s:
        year = s.split('Y')[1][0] # assumes that year is single digit
    if 'V' in s:
        visit = s.split('V')[1] # assumes that visit is everything after 'V'
        visit = visit.replace("_Full", "") # remove _Full if present    
    return int(year), int(visit)

def get_isotype(s): return s.translate(None, digits)

def clean_float(f):
    try:
        if math.isnan(float(f)):
            return "NULL"
        else:
            return float(f)
    except:
        return "NULL"

def split_clean(s): return tuple([x if x != "-1" else "NULL" for x in s.split(",")])

def make_sequences_sequence_strings_quality_scores_lineages_sql_input(infile, patient_uid, patient_str,
                                                                      path_to_year_visit_str_to_uid_map,
                                                                      p, path_to_subsampling_p_to_uid_map, replicate,
                                                                      clustering_params,
                                                                      path_to_clustering_params_to_uid_map,
                                                                      path_to_sequences_lineages_field_dict,
                                                                      wdir, wdir_original):

    out_sequences = open(wdir+"/sequences.txt.temp", 'w')
    out_sequence_strings = open(wdir+"/sequence_strings.txt.temp", 'w')
    out_quality_scores = open(wdir+"/quality_scores.txt.temp", 'w')
    out_lineages = open(wdir+"/lineages.txt.temp", 'w')

    fields = load_dict(path_to_sequences_lineages_field_dict)

    lineage_num_unique_seqs = 0
    lineage_abundance = 0

    with open(infile, 'rU') as f:
        for line, next_line in pairwise(f):

            vals = line.rstrip().split("\t")

            seq_id = vals[fields["seq_id_abundance"]].split("_")[0]

            year_visit_uid = get_year_visit_uid(vals[fields["year_visit_str"]], path_to_year_visit_str_to_uid_map)
            year, visit = split_year_visit(vals[fields["year_visit_str"]])

            lib_uid = get_lib_uid(patient_uid, year_visit_uid, p, path_to_subsampling_p_to_uid_map, replicate)

            clustering_uid = get_clustering_uid(patient_uid, p, path_to_subsampling_p_to_uid_map, replicate,
                                                clustering_params,
                                                path_to_clustering_params_to_uid_map)

            uid = get_sequence_uid(patient_uid, year_visit_uid, p, path_to_subsampling_p_to_uid_map, replicate, clustering_params, path_to_clustering_params_to_uid_map, seq_id)
            sequence_string_uid = get_sequence_string_uid(patient_uid, year_visit_uid, seq_id)

            lineage_id = int(line.rstrip().split()[fields["lineage_id"]])
            lineage_uid = get_lineage_uid(patient_uid, p, path_to_subsampling_p_to_uid_map, replicate, clustering_params, path_to_clustering_params_to_uid_map, lineage_id)

            uid_work_group_4 = int(uid) % 4
            uid_work_group_16 = int(uid) % 16
            uid_work_group_64 = int(uid) % 64
            uid_work_group_128 = int(uid) % 128
            uid_work_group_192 = int(uid) % 192
            uid_work_group_384 = int(uid) % 384

            lineage_uid_work_group_4 = int(lineage_uid) % 4
            lineage_uid_work_group_16 = int(lineage_uid) % 16
            lineage_uid_work_group_64 = int(lineage_uid) % 64
            lineage_uid_work_group_128 = int(lineage_uid) % 128
            lineage_uid_work_group_192 = int(lineage_uid) % 192
            lineage_uid_work_group_384 = int(lineage_uid) % 384

            fwr1_start, cdr1_start, fwr2_start, cdr2_start, fwr3_start, cdr3_start, cdr3_end = split_clean(vals[fields["cdr_fwr_boundaries"]])
            fwr1_length, cdr1_length, fwr2_length, cdr2_length, fwr3_length, cdr3_length = split_clean(vals[fields["cdr_fwr_lengths"]])

            fwr1_AA_start, cdr1_AA_start, fwr2_AA_start, cdr2_AA_start, fwr3_AA_start, cdr3_AA_start, cdr3_AA_end = tuple(["NULL"] * 7)

            untemplated_start, untemplated_end = tuple(["NULL"] * 2)

            isotype = get_isotype(vals[fields["subisotype"]])
            subisotype = vals[fields["subisotype"]]
            second_best_subisotype = vals[fields["second_best_isotype"]]

            frame = "NULL"
            
            if vals[fields["reads_per_molecule"]] != "1":
                more_than_one_read = "1"
            else:
                more_than_one_read = "0"

            # Add line to sequences, sequence_strings, quality_scores
            out_vals = [uid, uid_work_group_4, uid_work_group_16, uid_work_group_64,
                        uid_work_group_128, uid_work_group_192, uid_work_group_384,
                        seq_id, sequence_string_uid, sequence_string_uid, lineage_uid,
                        lineage_uid_work_group_4, lineage_uid_work_group_16, lineage_uid_work_group_64,
                        lineage_uid_work_group_128, lineage_uid_work_group_192, lineage_uid_work_group_384,
                        patient_uid, patient_str, lib_uid, clustering_uid,
                        year_visit_uid, vals[fields["year_visit_str"]], year, visit,
                        vals[fields["V"]], vals[fields["D"]], vals[fields["J"]],
                        clean_float(vals[fields["V_Evalue"]]),
                        clean_float(vals[fields["D_Evalue"]]),
                        clean_float(vals[fields["J_Evalue"]]),
                        isotype, subisotype, vals[fields["isotype_primer"]],
                        clean_float(vals[fields["best_isotype_Evalue"]]),
                        clean_float(vals[fields["best_isotype_score"]]),
                        second_best_subisotype,
                        clean_float(vals[fields["second_best_isotype_Evalue"]]),
                        clean_float(vals[fields["second_best_isotype_score"]]),
                        clean_float(vals[fields["isotype_ambiguity"]]),
                        vals[fields["V_len"]], vals[fields["D_len"]], vals[fields["J_len"]], vals[fields["C_len"]],
                        fwr1_start, cdr1_start, fwr2_start, cdr2_start, fwr3_start, cdr3_start, cdr3_end,
                        fwr1_length, cdr1_length, fwr2_length, cdr2_length, fwr3_length, cdr3_length,
                        fwr1_AA_start, cdr1_AA_start, fwr2_AA_start, cdr2_AA_start, fwr3_AA_start, cdr3_AA_start, cdr3_AA_end,
                        untemplated_start, untemplated_end,
                        frame, vals[fields["stop_codon"]], vals[fields["productive"]],
                        vals[fields["mut_germline_positions"]], vals[fields["mut_germline_before"]],
                        vals[fields["mut_germline_after"]], vals[fields["mut_germline_density"]],
                        vals[fields["V_germline_identity"]],
                        vals[fields["abundance"]],
                        vals[fields["reads_per_molecule"]],
                        more_than_one_read,
                        vals[fields["visits"]], vals[fields["abundances"]], vals[fields["isotypes"]],
                        vals[fields["is_first_visit_seen"]]]

            # Additional column containing the original seq id (useful for Sanger reads with non-standardized names)
            if "original_seq_id" in fields.keys():
                out_vals.append(vals[fields["original_seq_id"]])

            out_sequences.write("\t".join(map(str, out_vals)) + "\n")

            out_vals = [sequence_string_uid, seq_id, patient_uid, lib_uid,
                        vals[fields["V_seq"]], vals[fields["D_seq"]], vals[fields["J_seq"]],
                        vals[fields["C_seq"]], vals[fields["leader_seq"]], vals[fields["AA"]]]
            out_sequence_strings.write("\t".join(map(str, out_vals)) + "\n")

            out_vals = [sequence_string_uid, seq_id, patient_uid, lib_uid,
                        vals[fields["quals"]]]
            out_quality_scores.write("\t".join(map(str, out_vals)) + "\n")

            # Keep track of number of sequences and abundance of lineage
            lineage_num_unique_seqs += 1
            lineage_abundance += int(vals[fields["abundance"]])

            # Decide if we should add the lineage
            if next_line == None:
                # end of file, add lineage
                clustering_uid = get_clustering_uid(patient_uid, p, path_to_subsampling_p_to_uid_map, replicate,
                                                    clustering_params, path_to_clustering_params_to_uid_map)

                out_vals = [lineage_uid, lineage_uid_work_group_4, lineage_uid_work_group_16, lineage_uid_work_group_64,
                            lineage_uid_work_group_128, lineage_uid_work_group_192, lineage_uid_work_group_384,
                            patient_uid, clustering_uid, p, replicate, wdir_original,
                            lineage_num_unique_seqs, lineage_abundance]

                out_lineages.write("\t".join(map(str, out_vals)) + "\n")

                break

            next_lineage_id = int(next_line.rstrip().split()[fields["lineage_id"]])

            if next_lineage_id != lineage_id:
                # end of lineage, add lineage
                clustering_uid = get_clustering_uid(patient_uid, p, path_to_subsampling_p_to_uid_map, replicate,
                                                    clustering_params, path_to_clustering_params_to_uid_map)

                out_vals = [lineage_uid, lineage_uid_work_group_4, lineage_uid_work_group_16, lineage_uid_work_group_64,
                            lineage_uid_work_group_128, lineage_uid_work_group_192, lineage_uid_work_group_384,
                            patient_uid, clustering_uid, p, replicate, wdir_original,
                            lineage_num_unique_seqs, lineage_abundance]

                out_lineages.write("\t".join(map(str, out_vals)) + "\n")
                lineage_num_unique_seqs = 0
                lineage_abundance = 0

    out_sequences.close()
    out_sequence_strings.close()
    out_lineages.close()

    return None

if __name__ == "__main__":
    main(sys.argv)
