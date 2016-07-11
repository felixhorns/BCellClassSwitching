import sys
import os

path_to_cleaned = sys.argv[1]
path_to_lib_info_dir_to_patient_visit_map = sys.argv[2]

wdir = os.path.normpath(path_to_cleaned).split("/cleaned")[0]

# Load map
m = {}
with open(path_to_lib_info_dir_to_patient_visit_map, 'rU') as f:
    for line in f:
        vals = line.rstrip().split("\t")
        m[vals[0]] = vals[1:]

# Look up library info
vals = m[wdir]
patient_name = vals[0]
visit_name = vals[1]

# Count number of sequences
num_seq = 0
with open(wdir+"/parsed_igblast_isotypes_quals", 'rU') as f:
    for line in f:
        num_seq += 1

# Write lib info file
with open(wdir+"/lib_info.txt", 'w') as out:
    out.write(patient_name + "\n")
    out.write(visit_name + "\n")
    out.write(wdir+"/parsed_igblast_isotypes_quals" + "\n")
    out.write("NA" + "\n")
    out.write("1.0" + "\n")
    out.write(str(num_seq) + "\n")
    out.write(str(num_seq) + "\n")
    out.write(str(num_seq) + "\n")
    out.write("1" + "\n")
