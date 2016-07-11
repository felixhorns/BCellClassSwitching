import sys
import os
import zipfile
from itertools import izip
from collections import deque

infile = sys.argv[1] # "/groups.zip.temp"
num_workers = int(sys.argv[2])

outfile = os.path.dirname(infile) + "/group_worker_assignments.txt.temp"
out = open(outfile, 'w')

# Read file names and sizes
zf = zipfile.ZipFile(infile)

filenames = []
sizes = []

for f in zf.infolist():
    filenames.append(f.filename)
    sizes.append(int(f.file_size))

zf.close()

# Calculate loads
loads = [s**2 for s in sizes] # Load is N^2
target_load = sum(loads) / num_workers

# Assign loads to workers
sorted_lists = sorted(izip(filenames, loads), key=lambda x: x[1])
sorted_filenames, sorted_loads = [[x[i] for x in sorted_lists] for i in range(len(sorted_lists[0]))]

worker_queue = deque(range(0, num_workers))

my_inputs = {}
my_load = {}

while len(sorted_filenames) > 0:

    w = worker_queue.popleft()

    try:
        my_inputs[w].append(sorted_filenames.pop())
        my_load[w] += sorted_loads.pop()
    except:
        my_inputs[w] = [sorted_filenames.pop()]
        my_load[w] = sorted_loads.pop()

    if (my_load[w] >= target_load) and (len(worker_queue) > 0):
        pass
    else:
        worker_queue.append(w)

# Write input files for each worker
for _, inputs in my_inputs.items():
    out.write("\t".join(inputs) + "\n")

out.close()
