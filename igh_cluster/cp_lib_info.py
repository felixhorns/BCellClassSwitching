import sys
import os

def main(argv):

    infile = argv[1] # do_clustering.in

    output_dir = os.path.dirname(infile)
    out = open(output_dir+"/lib_info_combined.txt", 'w')

    with open(infile, 'rU') as f:
        
        for line in f:
            
            lib_id = line.split()[0]
            path_to_parsed_igblast_quals = line.split()[1]

            path_to_lib_info = path_to_parsed_igblast_quals + "/lib_info.txt"

            if os.path.isfile(path_to_lib_info):

                with open(path_to_lib_info, 'rU') as f2:
                    out.write(lib_id+"\n")
                    for line in f2: out.write(line)
                    out.write("\n")

            else:

                print "Error:", path_to_lib_info, "does not exist"
                quit()

    out.close()

if __name__ == "__main__":
    main(sys.argv)
