'''
Cleans split files
# e.g. parsed_igblast.{n}

    - first argument is number of splits
    - all subsequent arguments are absolute file paths for files 
        whose temp files should be deleted
'''

import sys
import subprocess

num_splits = int(sys.argv[1])

for arg in sys.argv[2:]:
    temp_files = ['%s.%s' % (arg,x) for x in range(num_splits)]
    temp_list = ' '.join('%s' % x for x in temp_files)
    cmd = 'rm %s' % temp_list
    process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = process.communicate()