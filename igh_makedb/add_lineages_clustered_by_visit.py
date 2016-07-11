from make_init_db_input import *


def add_lineages_clustered_by_visit_main(infile,
                                         path_to_patient_name_to_uid_map,
                                         path_to_patient_data,
                                         path_to_year_visit_str_to_uid_map,
                                         path_to_subsampling_p_to_uid_map,
                                         path_to_clustering_params_to_uid_map,
                                         path_to_sequences_lineages_field_dict,
                                         outfile):
    """
    Clustering sequences by visit is alternative to clustering sequences for all patient visits together. Given
    that the sequences are the same, the existing database structure can be updated with a single
    lineage_visit_uid column.
    :return Write uid / lineage_visit_uid data to a text file
    """

    # Get working directory
    wdir = os.path.dirname(infile)  # directory containing input files for working

    # Outfile
    out = open(outfile, 'w')

    # Patient
    patient_name = get_patient_name(wdir)
    patient_uid = get_patient_uid(patient_name, path_to_patient_name_to_uid_map)
    make_patients_sql_input(patient_name, patient_uid, path_to_patient_data, wdir)

    # Subsampling and clustering
    p, replicate, clustering_params = make_libs_clusterings_sql_input(patient_uid,
                                                                      path_to_year_visit_str_to_uid_map,
                                                                      path_to_subsampling_p_to_uid_map,
                                                                      path_to_clustering_params_to_uid_map, wdir)

    fields = load_dict(path_to_sequences_lineages_field_dict)

    m2 = load_line_dict_by_lines(path_to_clustering_params_to_uid_map)
    clustering_uid = m2[str(clustering_params)]

    # write header
    out.write('uid\tlineage_visit_uid\n')

    with open(infile, 'rU') as f:
        for line in f:

            vals = line.rstrip().split("\t")

            seq_id = vals[fields["seq_id_abundance"]].split("_")[0]

            year_visit_uid = get_year_visit_uid(vals[fields["year_visit_str"]], path_to_year_visit_str_to_uid_map)

            uid = get_sequence_uid_faster(patient_uid, year_visit_uid, clustering_uid, seq_id)

            lineage_id = int(line.rstrip().split()[fields["lineage_id"]])
            lineage_visit_uid = get_lineage_visit_uid_faster(patient_uid, clustering_uid, lineage_id, year_visit_uid)

            # write out
            out.write('%s\t%s\n' % (uid, lineage_visit_uid))

if __name__ == "__main__":
    add_lineages_clustered_by_visit_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                                         sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
