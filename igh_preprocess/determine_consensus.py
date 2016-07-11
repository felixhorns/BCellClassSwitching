from __future__ import division
import re
import sys
import os
import time
import math
import random


##### Input
infile=     		sys.argv[1]  # indexed_reads_sorted.txt.{n}
Single_Read_Length= 	int(sys.argv[2])
Frag_length=	    	int(sys.argv[3])
Frag_length_std_dev=	int(sys.argv[4])

##### Output
output_dir = os.path.dirname(infile)
num_split=  infile.split('.')[-1]
outfile=output_dir+'/consensus'  # format: {dir}/consensus_{file}

consensus1='%s_R1.txt.%s' % (outfile,num_split)  # also used for FLASH
consensus2='%s_R2.txt.%s' % (outfile,num_split)  # also used for FLASH

out1=           open(consensus1,'w')
out2=           open(consensus2,'w')
out_complete=   open('%s_complete.txt.%s' % (outfile,num_split),'w')
#qual_complete=  open('%s_qual_complete.txt.%s' % (outfile,num_split),'w')
#qual_combined=  open('%s_qual_combined.txt.%s' % (outfile,num_split),'w')
#out_all=        open('%s_all.txt.%s' % (outfile,num_split),'w')
# Output then input: output_dir+'out.extendedFrags.fastq

#####

def reverse_complement(sequence):
  Seq=''
  complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq


def determine_consensus_unequal_length(read):
    fast1=list()
    fast2=list()

    for sequence in read:
        seq1=sequence[0].split('~~~~')[0]
        qual1=sequence[1].split('~~~~')[0]

        seq2=sequence[0].split('~~~~')[1]
        qual2=sequence[1].split('~~~~')[1]


        fast1.append([seq1,qual1,sequence[2],sequence[3],sequence[4],'forward'])
        fast2.append([seq2,qual2,sequence[2],sequence[3],sequence[4],'forward'])


    final_1,qual_final_1,matched_reads_1, error_1,count1,Matrix1,Matrix_comb1,Length1=determine_consensus_equal_length(fast1,Single_Read_Length)
    final_2,qual_final_2,matched_reads_2, error_2,count2,Matrix2,Matrix_comb2,Length2=determine_consensus_equal_length(fast2,Single_Read_Length)

    return final_1,qual_final_1,final_2,qual_final_2


def determine_consensus_known_length(read,SeqPrep):

    fast1=list()
    fast2=list()

    for sequence in read:

        seq1=sequence[0].split('~~~~')[0]
        qual1=sequence[1].split('~~~~')[0]

        seq2=sequence[0].split('~~~~')[1]
        qual2=sequence[1].split('~~~~')[1]

        R2=reverse_complement(seq2)

        fast1.append([seq1+' '*(SeqPrep-len(seq1)),qual1+' '*(SeqPrep-len(seq1)),sequence[2],sequence[3],sequence[4],'forward'])
        fast1.append([' '*(SeqPrep-len(seq2))+str(R2),' '*(SeqPrep-len(seq2))+qual2[::-1],sequence[2],sequence[3],sequence[4],'reverse'])
        # for debugging R1 / R2 alignment:
        #out_all.write(seq1+'\n')
        #out_all.write('!'*(SeqPrep-len(seq2))+str(R2)+'\n')


    final_1,qual_final_1,matched_reads_1,error_1,count1,Matrix,Matrix_comb,Length=determine_consensus_equal_length(fast1,SeqPrep)
    return final_1,qual_final_1,Matrix,Matrix_comb,Length

def determine_consensus_equal_length(reads,Length):
    Consensus_Seq=''
    Consensus_Qual=''
    total_dis=0

    Matrix={}
    Matrix_comb={}

    for x in xrange(0,Length,1):

        Matrix[x]={'A':[],'T':[],'C':[],'G':[],'N':[]}
        Matrix_comb[x]={}
        Base_Count={'A':0,'T':0,'C':0,'G':0,'N':0}

        # Prior probability distribution, expressing complete ignorance:
        Proba={'A':0.25,'T':0.25,'C':0.25,'G':0.25}

        for read in reads:
            if x<len(read[0]):
            # Update probability distribution using Bayes theorem:
                if read[0][x]!=' ':
                    measured_base=read[0][x]
                    measured_q=ord(read[1][x])-33
                    Matrix[x][measured_base].append(measured_q)
                    Base_Count[measured_base]+=1
                if read[0][x] not in [' ','N']:
                    measured_base=read[0][x]
                    measured_q=ord(read[1][x])-33
                    pcond={}
                    for base in ['A','T','C','G']:
                        if measured_base==base:
                            pcond[base]=1-(10**(-measured_q/10.0))
                        else:
                            pcond[base]=(1.0/3)*(10**(-measured_q/10.0))

                    Proba_temp_A=Proba['A']
                    Proba_temp_T=Proba['T']
                    Proba_temp_C=Proba['C']
                    Proba_temp_G=Proba['G']
                    denominator=pcond['A']*Proba_temp_A+pcond['C']*Proba_temp_C+pcond['T']*Proba_temp_T+pcond['G']*Proba_temp_G
                    Proba['A']=pcond['A']*Proba_temp_A/denominator
                    Proba['T']=pcond['T']*Proba_temp_T/denominator
                    Proba['C']=pcond['C']*Proba_temp_C/denominator
                    Proba['G']=pcond['G']*Proba_temp_G/denominator

        # Record a combined quality score for each base possibility:
        for base in ['A','T','C','G']:
            if 1-Proba[base]!=0:
                Matrix_comb[x][base]=int(round(min(75,-10*math.log10(1-Proba[base])),0))
            else:
                Matrix_comb[x][base]=75

        # Choose base with largest probability as consensus, choosing at random if there is a tie:
        Winners=[(key,val) for key,val in Proba.iteritems() if val == max(Proba.values())]
        if len(Winners)==1:
            Consensus_Seq+=Winners[0][0]
            if 1-Winners[0][1]!=0:
                quality_score=-10*math.log10(1-Winners[0][1])
            else:
                quality_score=75
            Consensus_Qual+=chr(int(round(min(75,quality_score)+33,0)))
            total_dis=sum([val for key,val in Base_Count.iteritems() if key != Winners[0][0]])
        else:
            choice=random.randint(0,len(Winners)-1)
            Consensus_Seq+=Winners[choice][0]
            if 1-Winners[choice][1]!=0:
                quality_score=-10*math.log10(1-Winners[choice][1])
            else:
                quality_score=75
            Consensus_Qual+=chr(int(round(min(75,quality_score)+33,0)))
            total_dis=sum([val for key,val in Base_Count.iteritems() if key != Winners[choice][0]])

    count=len(reads)

    return Consensus_Seq,Consensus_Qual,reads,total_dis,count,Matrix,Matrix_comb,Length

def prelim_analysis(index,read1,path1,path2,path1_complete):
    SeqPrep=[]
    diversity=set()
    count=0
    for read in read1:

        diversity.add(read[0])
        count+=1
        bar1=read[2]
        bar2=read[3]
        SeqPrep=int(read[4])


    final_1,qual_final_1,final_2,qual_final_2=determine_consensus_unequal_length(read1)

    out1.write('@'+str(index)+'\n'+str(final_1)+'\n'+'+'+'\n'+str(qual_final_1)+'\n')
    out2.write('@'+str(index)+'\n'+str(final_2)+'\n'+'+'+'\n'+str(qual_final_2)+'\n')
    return path1, path2, path1_complete

def start_analysis(index,read1,path1,path2,path1_complete):

    SeqPrep=[]
    diversity=set()
    count=0
    for read in read1:

        diversity.add(read[0])
        count+=1
        bar1=read[2]
        bar2=read[3]
        SeqPrep=int(read[4])
    # print index,SeqPrep
    final_1,qual_final_1,Matrix,Matrix_comb,Length=determine_consensus_known_length(read1,SeqPrep)

    out_complete.write('@'+str(index)+'\n'+str(final_1)+'\n'+'+'+'\n'+str(qual_final_1)+'\n')

    # for debugging base calling and quality weighting:
    #qual_complete.write('@'+str(index)+'\t'+str(SeqPrep)+'\t')
    #qual_combined.write('@'+str(index)+'\t'+str(SeqPrep)+'\t')
    #for x in xrange(0,Length,1):
    #    for Base in ['A','T','C','G','N']:
    #         qual_complete.write(Base+':')
    #         Base_Scores=Matrix[x]
    #         for Base_Score in Base_Scores[Base]:
    #             qual_complete.write(str(Base_Score)+',')
    #         qual_complete.write('/')
    #    qual_complete.write('\t')
    #    for Base in ['A','T','C','G']:
    #         qual_combined.write(Base+':')
    #         qual_combined.write(str(Matrix_comb[x][Base])+',')
    #         qual_combined.write('/')
    #    qual_combined.write('\t')
    #qual_complete.write('\n')
    #qual_combined.write('\n')


random.seed(1)
path1=0
path1_complete=0
path2=0

molecule=''
counter=0


read1=[]
for line in open(infile):
    read=line
    a=read.strip().split('\t')
    index=a[1]
    if molecule=='':
        molecule=index
    bar1=a[2]
    bar2=a[3]
    base1=a[4]
    qual1=a[5]
    position1=-1

    if index==molecule:
        read1.append([base1,qual1,bar1,bar2,position1])
    else:
        prelim_analysis(molecule,read1,path1,path2,path1_complete)
        read1=[]
        read1.append([base1,qual1,bar1,bar2,position1])
        molecule=index

prelim_analysis(index,read1,path1,path2,path1_complete)

out1.close()
out2.close()

os.system('/datastore/dcroote/resources/FLASH-1.2.11/flash --suffix=%s -r %s -f %s -s %s \
          -d %s %s %s' % (num_split,Single_Read_Length,Frag_length,Frag_length_std_dev,output_dir,consensus1,consensus2))

Length=0
for line in open('%s/out.extendedFrags.fastq.%s' % (output_dir,num_split)):
    Length+=1

in1=open('%s/out.extendedFrags.fastq.%s' % (output_dir,num_split),'r')
counter=0
length_dict={}
while counter < Length:
    a=in1.readline().strip()
    b=in1.readline().strip()
    c=in1.readline().strip()
    d=in1.readline().strip()
    index=a.split('@')[1]
    read_length=len(b)
    length_dict[index]=read_length
    counter+=4

for line in open(infile):
    read=line
    a=read.strip().split('\t')
    index=a[1]
    if molecule=='':
        molecule=index
    bar1=a[2]
    bar2=a[3]
    base1=a[4]
    qual1=a[5]

    try:
        position1=length_dict[index]
    except:
        position1=-1

    if index==molecule:
        read1.append([base1,qual1,bar1,bar2,position1])
    else:
        if read1[0][4]!=-1:

            start_analysis(molecule,read1,path1,path2,path1_complete)
        read1=[]
        read1.append([base1,qual1,bar1,bar2,position1])
        molecule=index

if read1[0][4]!=-1:
    start_analysis(index,read1,path1,path2,path1_complete)
