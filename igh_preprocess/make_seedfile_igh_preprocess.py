import sys
import os

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

path = sys.argv[1]
seedfile_name = sys.argv[2]

with open(seedfile_name, 'w') as out:
    for lib_path in sorted(listdir_fullpath(path)):
        for data_path in sorted(listdir_fullpath(lib_path)):
            out.write(data_path+'\n')
