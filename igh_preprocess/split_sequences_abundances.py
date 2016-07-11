import sys
import os
import subprocess

def main(argv):

    infile = sys.argv[1]
    num_pieces = int(sys.argv[2])

    outfile_basename = infile
    output_dir = os.path.dirname(infile)

    L = get_length(infile)
    target_length_of_piece = float(L) / float(num_pieces)

    n = 0

    outfile = outfile_basename + "." + str(n)
    out = open(outfile, 'w')
    uids_split_on = open(output_dir + '/sequences_abundances_uids_split_on', 'w')

    i = 0
    j = 0
    with open(infile, 'rU') as f:
        while i < L:

            line1 = f.readline()
            line2 = f.readline()

            out.write(line1)
            out.write(line2)

            i += 2
            j += 2

            if j > target_length_of_piece:
                uids_split_on.write(line1.split('>')[1])
                out.close()
                n += 1 # next outfile
                j = 0
                outfile = outfile_basename + "." + str(n)
                out = open(outfile, 'w')
    
    uids_split_on.close()
    out.close()

def get_length(infile):
    i = 0
    with open(infile, 'rU') as f:
        for line in f:
            i += 1
    return i

if __name__ == "__main__":
    main(sys.argv)
