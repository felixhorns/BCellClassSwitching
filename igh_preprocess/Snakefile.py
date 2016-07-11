import os
import re
import subprocess

##### Paths
RESOURCES=				'/local10G/rfhorns/resources'
ANACONDA=				RESOURCES+'/anaconda'

SCRIPTS=				'/b1_2/home/rfhorns/scripts/igh_preprocess'
configfile:				os.environ.get("CONFIGFILE")

INDEX_READS=				SCRIPTS+'/index_reads.py'
SPLIT_BARCODES=                         SCRIPTS+'/split_barcodes.py'
SPLIT_FOR_DETERMINE_CONSENSUS=		SCRIPTS+'/split_for_determine.py'
DETERMINE_CONSENSUS=			SCRIPTS+'/determine_consensus.py'
COMBINE_IDENTICALS=			SCRIPTS+'/combine_identicals.py'
SPLIT_SEQUENCES_ABUNDANCES=		SCRIPTS+'/split_sequences_abundances.py'
SPLIT_SEQUENCES_QUALS=			SCRIPTS+'/split_sequences_quals.py'
PARSE_IGBLAST=				SCRIPTS+'/parse_igblast.py'
PARSE_ISOTYPE_BLAST=			SCRIPTS+'/parse_isotype_blast.py'
PARSE_QUALITY_SCORES=			SCRIPTS+'/parse_sequences_quals.py'
QC=					SCRIPTS+'/plot_qc.py'
CLEAN=					SCRIPTS+"/clean_split.py"
MAKE_LIB_INFO=                          SCRIPTS+"/make_lib_info.py"
PHYLO=					SCRIPTS+'/make_read_alignment_phylogram.py'

IGBLAST_DIR=				'/b1_2/home/rfhorns/resources/igblast'
IGBLAST_BIN=				IGBLAST_DIR+'/ncbi-igblast-1.3.0/bin/igblastn'
IGBLAST_GERMLINE_DB_V=			IGBLAST_DIR+'/database/Vsegments_20150201.fasta'
IGBLAST_GERMLINE_DB_D=			IGBLAST_DIR+'/database/Dsegments_20150201.fasta'
IGBLAST_GERMLINE_DB_J=			IGBLAST_DIR+'/database/Jsegments_20150201.fasta'

BLASTN_BIN=				"/b1_2/home/rfhorns/resources/ncbi-blast-2.2.29+/bin/blastn"
BLASTN_DB=				'/b1_2/home/rfhorns/resources/ncbi-blast-2.2.29+/db/'

##### Parameters
workdir:				config["workdir"] + "/log"
SEEDFILE=				config["seedfile"]

num_pieces=				config["num_workers"]
scratch_storage=			config["scratch_storage"]

# Adapter trimming
ILLUMINA_D701_D712_ADAPTER_L=		config["ILLUMINA_D701_D712_ADAPTER_L"]
ILLUMINA_D701_D712_ADAPTER_R=       	config["ILLUMINA_D701_D712_ADAPTER_R"]
ILLUMINA_D501_D508_ADAPTER_L=        	config["ILLUMINA_D501_D508_ADAPTER_L"]
ILLUMINA_D501_D508_ADAPTER_R=           config["ILLUMINA_D501_D508_ADAPTER_R"]

# Index reads
barcode_length=				config["barcode_length"]

# Split barcodes
num_mismatches_cutoff=                  config["num_mismatches_cutoff"]
quality_score_cutoff=                   config["quality_score_cutoff"]

try:
    split_bar_no_qual=                  config["split_bar_no_qual"]
except:
    split_bar_no_qual=0

# Determine consensus:
read_length=				config["read_length"]
effective_read_length=			read_length - barcode_length
frag_length=				config["frag_length"]
effective_frag_length=			frag_length - barcode_length - barcode_length
frag_length_std_dev=			config["frag_length_std_dev"]

# Combine identicals
try:
    variable_barcode_len=                  config["variable_barcode_len"]
except:
    variable_barcode_len=0

# Parse igblast:
len_C_cutoff_min=			config["len_C_cutoff_min"]
len_C_cutoff_max=			config["len_C_cutoff_max"]
log_V_Evalue_cutoff=			config["log_V_Evalue_cutoff"]
log_J_Evalue_cutoff=			config["log_J_Evalue_cutoff"]
C_PRIMER_MATCH_DICT=			config["C_primer_match_dictionary"]

# Isotype blast:
isotype_blast_db=			config["isotype_blast_db"]
isotype_blast_evalue=			config["isotype_blast_evalue"]
isotype_blast_word_size=		config["isotype_blast_word_size"]
isotype_blast_gapopen=			config["isotype_blast_gapopen"]
isotype_blast_gapextend=		config["isotype_blast_gapextend"]
isotype_blast_penalty=			config["isotype_blast_penalty"]
isotype_blast_reward=			config["isotype_blast_reward"]

# Library information
LIB_INFO_FILE=                          config["lib_info_file"]

##### Function for loading seeds
def load_seeds(infile):
    seeds = []
    with open(infile, 'rU') as f:
        for line in f:
            seeds.append(line.rstrip())
    return seeds

if SEEDFILE is not None:
    SEEDS = load_seeds(SEEDFILE)

##### Functions for transferring files to/from cluster
def name_on_scratch(s, scratch):
    return scratch+"/"+os.path.basename(s)

def names_on_scratch(names, scratch):
    return [name_on_scratch(n, scratch) for n in names]

def cp_to_scratch(inputs, scratch):
    for i in inputs:
      cmd = "rsync -aW " + i + " " + name_on_scratch(i, scratch)
      subprocess.call(cmd, shell=True)
    return None

def cp_from_scratch(outputs, scratch):
    for o in outputs:
        cmd = "rsync -aW " + name_on_scratch(o, scratch) + " " + o
        subprocess.call(cmd, shell=True)
    return None

##### Rules

ruleorder: combine_identicals > sort_indexed_reads > make_lib_info > auto_merge

rule all:
  input: expand("{dir}/output_statistics", dir=SEEDS),
         expand("{dir}/lib_info.txt", dir=SEEDS)
  params: name='all', partition='general', mem='1024'

def get_all_files(d):
    return [d+"/"+f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def unzip_fastq(f):
    cmd = "gunzip " + f
    p = subprocess.Popen(cmd, shell=True)
    return None

def get_raw_fastqs(wildcards):
    # unzip all fastqs
    all_files = get_all_files(wildcards.dir)
    fastq_gzs = [f for f in all_files if ".fastq.gz" in f]
    for f in fastq_gzs: unzip_fastq(f)
    # find the raw fastqs
    all_files = get_all_files(wildcards.dir)
    raw_fastqs = []
    for f in all_files:
        if ("L00" in f) and ("R"+wildcards.R in f) and (os.path.splitext(f)[1] == ".fastq"):
            raw_fastqs.append(f)
    return sorted(raw_fastqs)

rule cat_raw_fastqs:
  input: get_raw_fastqs
  output: '{dir}/Combined_R{R,\d+}.fastq'
  threads: 1
  params: name='cat_raw_fastqs', partition=config["params_partition_cat_raw_fastqs"], mem=config["params_mem_cat_raw_fastqs"]
  shell: 'cat {input} > {output}'

rule trim_adapters:
  input: "{dir}/Combined_R1.fastq",
  	 "{dir}/Combined_R2.fastq"
  output: "{dir}/Combined_R1.trimmed.fastq",
  	  "{dir}/Combined_R2.trimmed.fastq",
	  "{dir}/cutadapt.1.info.out",
	  "{dir}/cutadapt.2.info.out",
	  "{dir}/cutadapt.1.out",
	  "{dir}/cutadapt.2.out",
	  temp("{dir}/tmp.1.fastq"),
	  temp("{dir}/tmp.2.fastq")
  threads: 1
  params: name="trim_adapters", partition=config["params_partition_trim_adapters"], mem=config["params_mem_trim_adapters"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           cutadapt --anywhere {ILLUMINA_D701_D712_ADAPTER_L} --anywhere {ILLUMINA_D701_D712_ADAPTER_R} --discard-trimmed --overlap 14 --info-file {output_on_scratch[2]} -o {output_on_scratch[6]} -p {output_on_scratch[7]} {input_on_scratch[0]} {input_on_scratch[1]} > {output_on_scratch[4]}")
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           cutadapt --anywhere {ILLUMINA_D501_D508_ADAPTER_L} --anywhere {ILLUMINA_D501_D508_ADAPTER_R} --discard-trimmed --overlap 14 --info-file {output_on_scratch[3]} -o {output_on_scratch[0]} -p {output_on_scratch[1]} {output_on_scratch[6]} {output_on_scratch[7]} > {output_on_scratch[5]}")
    cp_from_scratch(output, scratch)

rule index_reads:
  input: '{dir}/Combined_R1.trimmed.fastq', '{dir}/Combined_R2.trimmed.fastq'
  output: temp('{dir}/indexed_reads.txt'),
          '{dir}/all_barcodes',
  threads: 1
  params: name='index_reads', partition=config["params_partition_index_reads"], mem=config["params_mem_index_reads"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {INDEX_READS} {input_on_scratch[0]} {input_on_scratch[1]} {barcode_length}")
    cp_from_scratch(output, scratch)

rule sort_indexed_reads:
  input:  '{dir}/indexed_reads.txt'
  output: '{dir}/indexed_reads_sorted.TXT'
  threads: 1
  params: name='sort_indexed_reads', partition=config["params_partition_sort_indexed_reads"], mem=config["params_mem_sort_indexed_reads"]
  run:
    scratch = os.environ[scratch_storage]
    shell("sort -k2,2n -k7,7n -T {scratch} {input} > {output}")
    # -T refers to temp directory (/tmp/ runs out of space on large libraries)

rule split_for_determine_consensus:
  input: '{dir}/indexed_reads_sorted.TXT'
  output: ['{dir}/indexed_reads_sorted.TXT.%s' % x for x in range(num_pieces)]
  threads: 1
  params: name='split_for_determine_consensus', partition=config["params_partition_split_for_determine_consensus"], mem=config["params_mem_split_for_determine_consensus"]
  shell: "source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {SPLIT_FOR_DETERMINE_CONSENSUS} {input} {num_pieces}"
  
rule split_barcodes:
  input: '{dir}/indexed_reads_sorted.TXT.{n}'
  output: temp('{dir}/indexed_reads_sorted_barcodesplit.txt.{n}')
  threads: 1
  params: name='split_barcodes', partition=config["params_partition_split_barcodes"], mem=config["params_mem_split_barcodes"]
  run: 
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
          python {SPLIT_BARCODES} {input_on_scratch} {output_on_scratch} {num_mismatches_cutoff} {quality_score_cutoff} {split_bar_no_qual}")
    cp_from_scratch(output, scratch)

rule determine_consensus:
  input: '{dir}/indexed_reads_sorted_barcodesplit.txt.{n}'
  output: temp('{dir}/consensus_R1.txt.{n}'),
          temp('{dir}/consensus_R2.txt.{n}'),
          '{dir}/consensus_complete.txt.{n}',
          '{dir}/out.extendedFrags.fastq.{n}',
          '{dir}/out.notCombined_1.fastq.{n}',
          '{dir}/out.notCombined_2.fastq.{n}'
  threads: 1
  params: name='determine_consensus', partition=config["params_partition_determine_consensus"], mem=config["params_mem_determine_consensus"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {DETERMINE_CONSENSUS} {input_on_scratch} {effective_read_length} {effective_frag_length} {frag_length_std_dev}")
    cp_from_scratch(output, scratch)

rule combine_identicals:
  input: '{dir}/consensus_complete.txt'
  output:   '{dir}/sequences_molecule_ids.fasta',
            '{dir}/sequences_abundances.FASTA',
            '{dir}/sequences_reads_per_molecule.TXT',
            '{dir}/sequences_quals.TXT'
  threads: 1
  params: name='combine_identicals', partition=config["params_partition_combine_identicals"], mem=config["params_mem_combine_identicals"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {COMBINE_IDENTICALS} {input_on_scratch} {variable_barcode_len}")
    cp_from_scratch(output, scratch)

rule split_for_igblast:
  input: '{dir}/sequences_abundances.FASTA' 
  output: ['{dir}/sequences_abundances.FASTA.%s' % x for x in range(num_pieces)],
          '{dir}/sequences_abundances_uids_split_on'
  threads: 1
  params: name='split_for_igblast', partition=config["params_partition_split_for_igblast"], mem=config["params_mem_split_for_igblast"]
  run:
    scratch = os.environ[scratch_storage]
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {SPLIT_SEQUENCES_ABUNDANCES} {input} {num_pieces}")

rule split_sequences_quals:
  input: '{dir}/sequences_quals.TXT',
         '{dir}/sequences_abundances_uids_split_on'
  output: ['{dir}/sequences_quals.TXT.%s' % x for x in range(num_pieces)]
  threads: 1
  params: name='split_sequences_quals', partition=config["params_partition_split_sequences_quals"], mem=config["params_mem_split_sequences_quals"]
  shell: "source {ANACONDA}/bin/activate {ANACONDA} &&\
          python {SPLIT_SEQUENCES_QUALS} {num_pieces} {input[0]} {input[1]}"

rule igblast:
  input: '{dir}/sequences_abundances.FASTA.{n}'
  output: '{dir}/igblast_out.txt.{n}'
  threads: 1
  params: name='igblast', partition=config["params_partition_igblast"], mem=config["params_mem_igblast"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("cd {IGBLAST_DIR} && {IGBLAST_BIN} -query {input_on_scratch} -out {output_on_scratch} -outfmt 3 -num_alignments_V=1 -num_alignments_D=1 -num_alignments_J=1 -evalue=1e-20 -germline_db_V {IGBLAST_GERMLINE_DB_V} -germline_db_D {IGBLAST_GERMLINE_DB_D} -germline_db_J {IGBLAST_GERMLINE_DB_J} -domain_system imgt")
    cp_from_scratch(output, scratch)
    shell('sed "1,16d" {output} | head -n -19 > {wildcards.dir}/temp_igblast_out.{wildcards.n} &&\
           rm {output} && mv {wildcards.dir}/temp_igblast_out.{wildcards.n} {output}')
    
rule parse_igblast:
  input: '{dir}/sequences_abundances.FASTA.{n}',
         '{dir}/igblast_out.txt.{n}',
         '{dir}/sequences_reads_per_molecule.TXT'
  output: '{dir}/parsed_igblast.txt.{n}',
          '{dir}/C_seqs_without_primer.fasta.{n}',
          '{dir}/losses_parse_igblast.txt.{n}',
          '{dir}/seq_ids_no_hits_found.txt.{n}',
          '{dir}/seq_ids_no_CDR3_start.txt.{n}',
          '{dir}/seq_ids_no_CDR3_end.txt.{n}',
          '{dir}/seq_ids_no_alignment_summary.txt.{n}',
          '{dir}/seq_ids_V_Evalue_greater_than_cutoff.txt.{n}',
          '{dir}/seq_ids_J_Evalue_greater_than_cutoff.txt.{n}',
          '{dir}/seq_ids_len_C_too_small_large.txt.{n}',
          '{dir}/seq_ids_no_V_call.txt.{n}',
          '{dir}/seq_ids_no_J_call.txt.{n}'
  threads: 1
  params: name='parse_igblast', partition=config["params_partition_parse_igblast"], mem=config["params_mem_parse_igblast"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {PARSE_IGBLAST} {input_on_scratch[0]} {input_on_scratch[1]} {C_PRIMER_MATCH_DICT} {input_on_scratch[2]} \
          {len_C_cutoff_min} {len_C_cutoff_max} {log_V_Evalue_cutoff} {log_J_Evalue_cutoff}")
    cp_from_scratch(output, scratch)

rule isotype_blast:
  input: '{dir}/C_seqs_without_primer.fasta.{n}'
  output: '{dir}/isotype_blastn_out.xml.{n}'
  threads: 1
  params: name='isotype_blast', partition=config["params_partition_isotype_blast"], mem=config["params_mem_isotype_blast"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("export BLASTDB={BLASTN_DB}:$BLASTDB &&\
           {BLASTN_BIN} -db {isotype_blast_db} -query {input_on_scratch} -outfmt 5 -out {output_on_scratch} \
           -evalue {isotype_blast_evalue} -word_size {isotype_blast_word_size} \
           -gapopen {isotype_blast_gapopen} -gapextend {isotype_blast_gapextend} \
           -penalty {isotype_blast_penalty} -reward {isotype_blast_reward} -num_threads 1")
    cp_from_scratch(output, scratch)

rule parse_isotype_blast:
  input: '{dir}/parsed_igblast.txt.{n}', '{dir}/isotype_blastn_out.xml.{n}'
  output: '{dir}/parsed_igblast_isotypes.txt.{n}',
          '{dir}/losses_parse_isotype_blast.txt.{n}'
  threads: 1
  params: name='parse_isotype_blast', partition=config["params_partition_parse_isotype_blast"], mem=config["params_mem_parse_isotype_blast"]
  run:
    scratch = os.environ["LOCAL_SATA"]
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {PARSE_ISOTYPE_BLAST} {input_on_scratch[0]} {input_on_scratch[1]}")
    cp_from_scratch(output, scratch)

rule parse_quality_scores:
  input: '{dir}/parsed_igblast_isotypes.txt.{n}' , '{dir}/sequences_quals.TXT.{n}'
  output: '{dir}/parsed_igblast_isotypes_quals.txt.{n}'
  threads: 1
  params: name='parse_quality_scores', partition=config["params_partition_parse_quality_scores"], mem=config["params_mem_parse_quality_scores"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {PARSE_QUALITY_SCORES} {input_on_scratch[0]} {input_on_scratch[1]}")
    cp_from_scratch(output, scratch)

rule auto_merge:
  input: ["{dir}/{file}.{ext,[a-z]+}.%s" % x for x in range(num_pieces)]
  output: '{dir}/{file}.{ext,[a-z]+}'
  threads: 1
  params: name='auto_merge', partition=config["params_partition_auto_merge"], mem=config["params_mem_auto_merge"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("cat {input_on_scratch} > {output_on_scratch}")
    cp_from_scratch(output, scratch)

rule clean:
  input:  # parsing and blasts
          '{dir}/parsed_igblast_isotypes_quals.txt',
          '{dir}/parsed_igblast_isotypes.txt',
          '{dir}/parsed_igblast.txt',
          '{dir}/igblast_out.txt',
          '{dir}/isotype_blastn_out.xml',
          '{dir}/C_seqs_without_primer.fasta',
          '{dir}/losses_parse_igblast.txt',
          '{dir}/losses_parse_isotype_blast.txt',
          '{dir}/sequences_abundances.FASTA',
          '{dir}/sequences_quals.TXT',
          # determine_consensus
          '{dir}/consensus_complete.txt',
          '{dir}/out.notCombined_1.fastq',
          '{dir}/out.notCombined_2.fastq',
          '{dir}/out.extendedFrags.fastq',
          # parse_igblast
          '{dir}/seq_ids_no_hits_found.txt',
          '{dir}/seq_ids_no_CDR3_start.txt',
          '{dir}/seq_ids_no_CDR3_end.txt',
          '{dir}/seq_ids_no_alignment_summary.txt',
          '{dir}/seq_ids_V_Evalue_greater_than_cutoff.txt',
          '{dir}/seq_ids_J_Evalue_greater_than_cutoff.txt',
          '{dir}/seq_ids_len_C_too_small_large.txt',
          '{dir}/seq_ids_no_V_call.txt',
          '{dir}/seq_ids_no_J_call.txt'
  output: '{dir}/cleaned'
  threads: 1
  params: name='clean', partition=config["params_partition_clean"], mem=config["params_mem_clean"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} &&\
          python {CLEAN} {num_pieces} {input} &&\
          touch {output[0]}'

rule rename_parse_igblast_isotypes_quals:
  input: '{dir}/parsed_igblast_isotypes_quals.txt',
         '{dir}/cleaned'
  output: '{dir}/parsed_igblast_isotypes_quals'
  threads: 1
  params: name='rename_parse_igblast_isotypes_quals', partition='general', mem='2048'
  shell: 'mv {input[0]} {output} && touch {output}'
  # legacy rule because downstream analyses use this file

rule legacy_barcode_split:
  input: '{dir}/indexed_reads_sorted_barcodesplit.txt'
  output: '{dir}/indexed_reads_sorted_barcodesplit.TXT'
  threads: 1
  params: name="legacy_barcode_split", partition="general", mem="5300"
  shell: "rsync {input} {output}"

rule qc:
  input: '{dir}/indexed_reads_sorted_barcodesplit.TXT', #'{dir}/indexed_reads_sorted_barcodesplit.TXT',
         '{dir}/parsed_igblast_isotypes_quals',
         '{dir}/Combined_R1.fastq',
         '{dir}/Combined_R1.trimmed.fastq',
         '{dir}/cutadapt.1.out',
         '{dir}/cutadapt.2.out',
         '{dir}/consensus_complete.txt',
         '{dir}/out.extendedFrags.fastq',
         '{dir}/out.notCombined_1.fastq',
         '{dir}/out.notCombined_2.fastq',
         '{dir}/cleaned'
  output: '{dir}/reads_per_molecule_hist.png',
          '{dir}/abundance_hist.png',
          '{dir}/isotype_counts.png',
          '{dir}/V_counts.png',
          '{dir}/J_counts.png',
          '{dir}/output_statistics',
          '{dir}/isotype_counts.txt',
          '{dir}/CDR3_lengths.png',
          '{dir}/losses_pipeline.png',
          '{dir}/index.html'
  threads: 1
  params: name='qc', partition=config["params_partition_qc"], mem=config["params_mem_qc"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} &&\
          python {QC} {input[0]} {input[1]}'

rule make_lib_info:
  input: '{dir}/cleaned'
  output: '{dir}/lib_info.txt'
  threads: 1
  params: name='make_lib_info', partition=config["params_partition_make_lib_info"], mem=config["params_mem_make_lib_info"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} &&\
          python {MAKE_LIB_INFO} {input[0]} {LIB_INFO_FILE}'
