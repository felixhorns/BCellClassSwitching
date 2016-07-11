#########################################################################################
# authors: Chris Vollmers, Felix Horns, Derek Croote
#########################################################################################
# index_reads.py
# Splits reads into molecular barcode ("index") and remainder of read.
# Assigns a unique number to each barcode.
#########################################################################################

from __future__ import division
import time
import sys
import os

##### Inputs
infile1=sys.argv[1] # {dir}/Combined_R1.fastq
infile2=sys.argv[2] # {dir}/Combined_R2.fastq
print 'Combined input files: %s and %s' % (infile1,infile2)

L = int(sys.argv[3]) # length of molecular barcode (number of random bases at beginning of each read), e.g. 12

##### Outputs
output_dir = os.path.dirname(infile1) + "/"

out = open(output_dir+'indexed_reads.txt','w')
out_barcodes = open(output_dir+'all_barcodes','w')
#out_read1 = open(output_dir+'R1_prelim','w')
#out_read2 = open(output_dir+'R2_prelim','w')
######


def file_len(fileName):
    i = 0
    for line in open(fileName): i += 1
    return i

length = file_len(infile1)

if length == 0:
    print "Error: input file Combined_R1.fastq is empty"
    quit()

if length < 5000:
    print "Error: input file Combined_R1.fastq has less than 5000 lines. It will not be processed."
    quit()

indexes = {}
index_coverages = {}
next_uid = 1

in1=open(infile1,'rU')
in2=open(infile2,'rU')

start = time.time()

line_num = 0

while line_num < length:

    name1 = in1.readline()
    seq1 = in1.readline()
    plus_line1 = in1.readline()
    qual1 = in1.readline()

    name2 = in2.readline()
    seq2 = in2.readline()
    plus_line2 = in2.readline()
    qual2 = in2.readline()

    name = name1.strip().split(' ')[0]
    index1 = seq1[0:L]
    index2 = seq2[0:L]
    seq = seq1[L:].strip()+'~~~~'+seq2[L:].strip()
    qual = qual1[L:].strip()+'~~~~'+qual2[L:].strip()

    index = index1 + index2
    try:
        uid = indexes[index]
        index_coverages[index] += 1
    except:
        uid = next_uid
        indexes[index] = next_uid
        index_coverages[index] = 1
        next_uid += 1

    out.write("\t".join([str(0), str(uid), index1, index2, seq, qual, str(line_num), name]))
    out.write("\n")

    # individual read outputs for debugging
    #out_read1.write('@'+str(uid)+'_'+str(line_num)+'\n'+seq1+plus_line1+qual1)
    #out_read2.write('@'+str(uid)+'_'+str(line_num)+'\n'+seq2+plus_line2+qual2)

    line_num += 4

for index, cnts in index_coverages.iteritems():
    out_barcodes.write(str(indexes[index])+'\t'+str(cnts)+'\t'+index[0:L]+index[L:]+'\n')

