"""
Parallelization of determining consensus
Bash splitting is not sufficient due to the potential for splitting to occur between molecules with the same sequence
that need to be used together in determining the consensus

Inputs:
- indexed_reads_sorted.txt
- desired parallelization (number of files to split into)

Outputs:
- indexed_reads_sorted.txt.n, where n is the zero-indexed split number (0 to
n-1)

"""

from __future__ import division
import sys
import os
from itertools import tee, izip_longest

##### Inputs
infile_indexed_reads_sorted = sys.argv[1]
num_splits = int(sys.argv[2])
output_dir = os.path.dirname(infile_indexed_reads_sorted)

##### Outputs
outfile_basename = infile_indexed_reads_sorted
outfiles = [open('%s.%s' % (outfile_basename,n),'w') for n in range(num_splits)]

def get_file_length(f):
    L = 0
    with open(f) as infile:
        for line in infile: L += 1
    return L

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip_longest(a, b)

file_len = get_file_length(infile_indexed_reads_sorted)
rough_split_len = round(int(file_len) / num_splits,0)

global_cnt=2
cnt=0
n=0  # used for tracking which file to write to
with open(infile_indexed_reads_sorted, 'rU') as f:

    for line, next_line in pairwise(f):
        mol_1_id = line.split('\t')[1]
        mol_2_id = next_line.split('\t')[1]

        if cnt < rough_split_len:
            outfiles[n].write(line)
            cnt+=1

        elif cnt >= rough_split_len:
            if mol_1_id != mol_2_id:
                # move to the next file b/c not splitting molecules
                outfiles[n].write(line)
                outfiles[n].close()
                n+=1  # move on to writing to the next file
                cnt=0  # reset the count

            elif mol_1_id == mol_2_id:
                outfiles[n].write(line)
                # n stays the same b/c don't want to split molecules
                cnt+=1

        if global_cnt == file_len:
            outfiles[n].write(next_line)
            break

        global_cnt+=1
        
outfiles[-1].close()
