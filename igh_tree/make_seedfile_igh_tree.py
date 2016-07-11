import sys
import os

path = sys.argv[1]
target_file = sys.argv[2]
seedfile_name = sys.argv[3]

with open(seedfile_name, 'w') as out:
    for d in os.walk(path):
        if target_file in d[2]:
            out.write(d[0] + "\n")
