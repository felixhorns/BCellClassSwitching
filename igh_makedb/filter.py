import sys
import os

global verbose
verbose = True

field_dict = {"abundance": 4, "productive": 23, "reads_per_molecule": 31}

def main(argv):

    infile = argv[1]
    outfile = argv[2]
    outfile_nonproductive = argv[3]
    outfile_low_conf = argv[4]
    filter_nonproductive = int(argv[5])
    filter_low_conf = int(argv[6])

    out = open(outfile, 'w')
    out_nonproductive = open(outfile_nonproductive, 'w')
    out_low_conf = open(outfile_low_conf, 'w')

    with open(infile, 'rU') as f:

        for line in f:

            vals = line.rstrip().split("\t")

            productive = vals[field_dict["productive"]]
            abundance = int(vals[field_dict["abundance"]])
            reads_per_molecule = vals[field_dict["reads_per_molecule"]]

            if filter_nonproductive == 1:
                if productive == "0":
                    out_nonproductive.write(line)
                    continue

            if filter_low_conf == 1:
                if abundance == 1 and reads_per_molecule == "1":
                    out_low_conf.write(line)
                    continue

            out.write(line)

    out.close()
    out_nonproductive.close()
    out_low_conf.close()

    print "Done!!"

if __name__ == "__main__":
    main(sys.argv)
