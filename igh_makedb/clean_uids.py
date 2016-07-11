""" Replaces non-numeric uids with numeric uids """
import sys
import json

infile = sys.argv[1]
infile_field_dict = sys.argv[2]
outfile = sys.argv[3]

# Load field dictionary
with open(infile_field_dict, 'rU') as f:
    field_dict = json.load(f)

out = open(outfile, 'w')

# Function that determines whether a string can be parsed as an int
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# Next uids
next_visit_uids = {}

# Function that returns next uid for given visit
def get_next_uid(visit, next_visit_uids):
    try:
        next_uid = next_visit_uids[visit]
        next_visit_uids[visit] += 1
    except:
        next_uid = 0
        next_visit_uids[visit] = 1
    return next_uid

# Parse file line by line
with open(infile, 'rU') as f:
    for line in f:
        
        vals = line.rstrip().split("\t")
        
        uid = vals[field_dict["seq_id_abundance"]].split("_")[0]

        if is_int(uid):

            vals.append(uid) # original uid is same as current uid

        else:

            visit = vals[field_dict["year_visit_str"]]
            clean_uid = get_next_uid(visit, next_visit_uids)
            vals[field_dict["seq_id_abundance"]] = clean_uid
            vals.append(uid)

            print "cleaned", uid, "to", clean_uid

        out.write("\t".join(map(str, vals)) + "\n")
