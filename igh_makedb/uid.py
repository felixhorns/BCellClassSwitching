import json

def load_dict(infile):
    with open(infile, 'rU') as f:
        d = json.load(f)
    return d

    if add_entry_for_one:
        # hacky way to add a mapping for "1.0" -> 1 (for p = 1.0, that is no subsampling)
        d[float(1.0)] = 1
        line_num = 2


def load_float_dict_by_lines(infile, add_entry_for_one=False):
    d = {}
    line_num = 1

    # hacky way to add "1.0" -> 1 (for p = 1.0 subsampling, that is no subsampling)
    if add_entry_for_one:
        d[float(1.0)] = 1
        line_num = 2

    with open(infile, 'rU') as f:
        for line in f:
            val = line.rstrip()
            d[float(val)] = line_num
            line_num += 1
    return d


def load_line_dict_by_lines(infile, add_entry_for_one=False):
    d = {}
    line_num = 1

    with open(infile, 'rU') as f:
        for line in f:
            vals = line.rstrip().split("\t")
            for i, v in enumerate(vals):
                try:
                    vals[i] = float(v)
                except:
                    continue
            d[str(vals)] = line_num
            line_num += 1
    return d


def get_patient_uid(name, path_to_patient_name_to_uid_map):
    m = load_dict(path_to_patient_name_to_uid_map)
    uid = m[name]
    return uid


def get_year_visit_uid(year_visit_str, path_to_year_visit_str_to_uid_map):
    m = load_dict(path_to_year_visit_str_to_uid_map)
    uid = m[year_visit_str]
    return uid


def get_lib_uid(patient_uid, year_visit_uid, subsampling_p, path_to_subsampling_p_to_uid_map, replicate):
    m = load_float_dict_by_lines(path_to_subsampling_p_to_uid_map, add_entry_for_one=True)
    p_uid = m[float(subsampling_p)]
    uid = str(patient_uid).zfill(3) + str(year_visit_uid).zfill(2) + str(p_uid).zfill(2) + str(replicate).zfill(2)
    return uid


def get_clustering_uid(patient_uid, subsampling_p, path_to_subsampling_p_to_uid_map, replicate,
                       clustering_params, path_to_clustering_params_to_uid_map):

    m = load_float_dict_by_lines(path_to_subsampling_p_to_uid_map, add_entry_for_one=True)
    p_uid = m[float(subsampling_p)]

    m2 = load_line_dict_by_lines(path_to_clustering_params_to_uid_map)
    clustering_uid = m2[str(clustering_params)]

    uid = str(patient_uid).zfill(3) + str(p_uid).zfill(2) + str(replicate).zfill(2) + str(clustering_uid).zfill(2)

    return uid


def get_sequence_uid(patient_uid, year_visit_uid, subsampling_p, path_to_subsampling_p_to_uid_map, replicate,
                     clustering_params, path_to_clustering_params_to_uid_map, seq_id):

    m = load_float_dict_by_lines(path_to_subsampling_p_to_uid_map, add_entry_for_one=True)
    p_uid = m[float(subsampling_p)]

    m2 = load_line_dict_by_lines(path_to_clustering_params_to_uid_map)
    clustering_uid = m2[str(clustering_params)]

    uid = str(patient_uid).zfill(3) + str(year_visit_uid).zfill(2) + str(p_uid).zfill(2) + str(replicate).zfill(2) + str(clustering_uid).zfill(2) + str(seq_id).zfill(8)

    return uid


def get_lineage_uid(patient_uid, subsampling_p, path_to_subsampling_p_to_uid_map, replicate,
                    clustering_params, path_to_clustering_params_to_uid_map, lineage_id):

    m = load_float_dict_by_lines(path_to_subsampling_p_to_uid_map, add_entry_for_one=True)
    p_uid = m[float(subsampling_p)]

    m2 = load_line_dict_by_lines(path_to_clustering_params_to_uid_map)
    clustering_uid = m2[str(clustering_params)]

    uid = str(patient_uid).zfill(3) + str(p_uid).zfill(2) + str(replicate).zfill(2) + str(clustering_uid).zfill(2) + str(lineage_id).zfill(10)

    return uid


def get_lineage_visit_uid(patient_uid, clustering_params, path_to_clustering_params_to_uid_map, lineage_id, year_visit_uid):
    """
    For use when sequences have been clustered on a per-visit basis rather than for all visits from a patient
    """

    m2 = load_line_dict_by_lines(path_to_clustering_params_to_uid_map)
    clustering_uid = m2[str(clustering_params)]

    uid = str(patient_uid).zfill(3) + str(year_visit_uid).zfill(2) + \
        str(clustering_uid).zfill(2) + str(lineage_id).zfill(10)

    return uid

def get_lineage_visit_uid_faster(patient_uid, clustering_uid, lineage_id, year_visit_uid):
    """
    For use when sequences have been clustered on a per-visit basis rather than for all visits from a patient
    """

    uid = str(patient_uid).zfill(2) + str(year_visit_uid).zfill(2) + \
        str(clustering_uid).zfill(2) + str(lineage_id).zfill(10)

    return uid


def get_sequence_uid_faster(patient_uid, year_visit_uid, clustering_uid, seq_id):
    """
    Also for clustering on a per-visit basis
    """
    uid = str(patient_uid).zfill(2) + str(year_visit_uid).zfill(2) + str(clustering_uid).zfill(2) + str(seq_id).zfill(8)

    return uid


def get_sequence_string_uid(patient_uid, year_visit_uid, seq_id):
    uid = str(patient_uid).zfill(3) + str(year_visit_uid).zfill(2) + str(seq_id).zfill(14)
    return uid
