import sys
import os

input_parsed_igblast_isotypes = sys.argv[1]
input_sequences_quals = sys.argv[2]

output_dir = os.path.dirname(input_parsed_igblast_isotypes)
split_num  = input_parsed_igblast_isotypes.split('.')[-1]

out = open('%s/parsed_igblast_isotypes_quals.txt.%s' % (output_dir,split_num), 'w')

length_sequences_quals = 0
for line in open(input_sequences_quals, 'rU'): length_sequences_quals += 1

quals = {}

with open(input_sequences_quals, 'rU') as f:

    line_num = 0
    while line_num < length_sequences_quals:

        line = f.readline()

        if ">" in line:
            uid = line.strip().split(">")[1]
            num_quals = int(uid.split("_")[-1])
        else:
            print "Error in parsing quality scores: '>' not found in expected header line", line_num
            quit()

        line_num += 1

        quals_of_molecules = []
        for i in range(num_quals):
            qual = f.readline().strip()
            quals_of_molecules.append(qual)
            line_num += 1

            quals[uid] = quals_of_molecules

with open(input_parsed_igblast_isotypes, 'rU') as f:
    for line in f:
        vals = line.rstrip().split("\t")
        uid = vals[0]

        quals_of_molecules = quals[uid]

        out.write("\t".join(vals)+"\t"+"~".join(quals_of_molecules)+"\n")

out.close()
