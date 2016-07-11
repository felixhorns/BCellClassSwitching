import sys
import random

infile_R1 = sys.argv[1]
infile_R2 = sys.argv[2]

P = float(sys.argv[3]) # Subsampling probability

# Output files
outfile_R1 = infile_R1 + ".subsampled"
outfile_R2 = infile_R2 + ".subsampled"

out_R1 = open(outfile_R1, 'w')
out_R2 = open(outfile_R2, 'w')

# Count lines in input files
L = 0
with open(infile_R1, 'rU') as f:
    for line in f:
        L += 1

# Inputs
in_R1 = open(infile_R1, 'rU')
in_R2 = open(infile_R2, 'rU')

line_num = 0

while line_num < L:
    
    id1 = in_R1.readline()
    seq1 = in_R1.readline()
    plus1 = in_R1.readline()
    qual1 = in_R1.readline()

    id2 = in_R2.readline()
    seq2 = in_R2.readline()
    plus2 = in_R2.readline()
    qual2 = in_R2.readline()

    line_num += 4

    if random.random() < P:
        
        out_R1.write(id1)
        out_R1.write(seq1)
        out_R1.write(plus1)
        out_R1.write(qual1)

        out_R2.write(id2)
        out_R2.write(seq2)
        out_R2.write(plus2)
        out_R2.write(qual2)

out_R1.close()
out_R2.close()
