import os
import re
import subprocess

##### Paths
RESOURCES=	                                '/local10G/dcroote/resources'
ANACONDA=					RESOURCES+'/anaconda'

SCRIPTS=					'/b1_2/home/rfhorns/scripts/igh_tree'
configfile:					os.environ.get("CONFIGFILE")

UPDATE_DB=					SCRIPTS+'/update_db.py'
UPDATE_OR_ADD_FEATURES_TO_DB=			SCRIPTS+'/update_or_add_features_to_db.py'

MST_QUERY=					SCRIPTS+'/mst_query.py'
MST=						SCRIPTS+'/mst.py'

GAPLESS_ALIGN=                                  SCRIPTS+'/gapless_align.py'
GAPLESS_ALIGN_QUERY=                            SCRIPTS+'/gapless_align_query.py'
GET_GERMLINE_SEQ=                               SCRIPTS+'/get_germline_seq.py'
GET_GERMLINE_SEQ_QUERY=                         SCRIPTS+'/get_germline_seq_query.py'
MUSCLE_PROFILE_ALIGN_GERMLINE=                  SCRIPTS+'/muscle_profile_align_germline.py'
MUSCLE_PROFILE_ALIGN_GERMLINE_QUERY=            SCRIPTS+'/muscle_profile_align_germline_query.py'

EXPORT_MST_ALL_LINEAGES_QUERY=                  SCRIPTS+'/export_mst_for_all_lineages_query.py'
EXPORT_MST_ALL_LINEAGES=                        SCRIPTS+'/export_mst_for_all_lineages.py'
EXPORT_MST_NO_ISOTYPE_CONSTRAINT_ALL_LINEAGES_QUERY=SCRIPTS+'/export_mst_no_isotype_constraint_for_all_lineages_query.py'
EXPORT_MST_NO_ISOTYPE_CONSTRAINT_ALL_LINEAGES=      SCRIPTS+'/export_mst_no_isotype_constraint_for_all_lineages.py'

EXPORT_MST_INDIV_LINEAGES_QUERY=                SCRIPTS+'/export_mst_for_individual_lineages_query.py'
EXPORT_MST_INDIV_LINEAGES=                      SCRIPTS+'/export_mst_for_individual_lineages.py'

BLASTN_BIN=                                     "/b1_2/home/rfhorns/resources/ncbi-blast-2.2.29+/bin/blastn"
BLASTN_DIR=                                     "/b1_2/home/rfhorns/resources/ncbi-blast-2.2.29+/bin/"

##### Parameters
workdir:                                        config["workdir"] + "/log"
SEEDFILE=                                       config["seedfile"]

num_workers=					config["num_workers"]
scratch_storage=                                config["scratch_storage"]

path_to_V_germline_ref_file=			config["path_to_V_germline_ref_file"]
path_to_J_germline_ref_file=			config["path_to_J_germline_ref_file"]

path_to_muscle_bin=				config["path_to_muscle_bin"]
path_to_clustalo_bin=                           config["path_to_clustalo_bin"]
path_to_edmonds_bin=				config["path_to_edmonds_bin"]

db=						config["db"]
mysql_db_host=					config["mysql_db_host"]
mysql_db_user=					config["mysql_db_user"]
mysql_db_password=				config["mysql_db_password"]

alignment_filter_stringency=			config["alignment_filter_stringency"]

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

def cp_dir_from_scratch(scratch, dir_name, output):
    dir_name_on_scratch = scratch + "/" + dir_name.rstrip("/")
    output_dir = os.path.dirname(output)
    cmd = "rsync -aW " + dir_name_on_scratch + " " + output_dir + "/"
    subprocess.call(cmd, shell=True)
    return None

def replace_str(names, x, y): return [n.replace(x, y) for n in names]

def move(x, y):
    for x_i, y_i in zip(x, y):
    	cmd = "rsync -aW " + x_i + " " + y_i
	subprocess.call(cmd, shell=True)
    return None

##### Function for loading seeds
def load_seeds(infile):
    seeds = []
    with open(infile, 'rU') as f:
    	 for line in f:
	 	     seeds.append(line.rstrip())
    return seeds

if SEEDFILE is not None:
   SEEDS = load_seeds(SEEDFILE)

##### Rules

ruleorder: mst_worker > gapless_align_worker > get_germline_seq_worker > muscle_profile_align_germline_worker > auto_merge_header_safe

rule all:
  input: expand("{dir}/mst", dir=SEEDS),
  params: name='all', partition='general', mem='1024'

rule all_mst:
  input: expand("{dir}/mst", dir=SEEDS),
  params: name='all_mst', partition='general', mem='1024'

rule all_mst_no_isotype_constraint:
  input: expand("{dir}/mst_no_isotype_constraint", dir=SEEDS),
  params: name='all_mst_no_isotype_constraint', partition='general', mem='1024'

rule all_export_mst_for_all_lineages:
  input: expand("{dir}/mst.all_lineages.nodes.cytoscape", dir=SEEDS),
  params: name='all_export_mst_for_all_lineages:', partition='general', mem='1024'

rule all_export_mst_no_isotype_constraint_for_all_lineages:
  input: expand("{dir}/mst_no_isotype_constraint.all_lineages.nodes.cytoscape", dir=SEEDS),
  params: name='all_export_mst_no_isotype_constraint_for_all_lineages:', partition='general', mem='1024'

rule all_msa:
  input: expand("{dir}/msa", dir=SEEDS),
  params: name='all_msa', partition='general', mem='1024'

rule all_gapless:
  input: expand("{dir}/muscle_profile_align_germline", dir=SEEDS)
  params: name='all_gapless', partition='general', mem='1024'

rule auto_merge_header_safe:
  input: ["{dir}/{file}.{ext,[a-z]+}.%s" % x for x in range(num_workers)]
  output: '{dir}/{file}.{ext,[a-z]+}'
  threads: 1
  params: name='auto_merge_header_safe', partition=config["params_partition_auto_merge_header_safe"], mem=config["params_mem_auto_merge_header_safe"]
  run:
    shell("echo input: {input} &&\
           echo output: {output}")
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("head -n 1 {input_on_scratch[0]} > {output_on_scratch[0]}")
    for i in range(0, num_workers):
        shell('sed "1,1d" {input_on_scratch[%s]} >> {output_on_scratch[0]}' % i)
        continue
    cp_from_scratch(output, scratch)
    shell("rm {input}")

rule export_mst_for_all_lineages_query:
  input: "{dir}/uids.txt",
  	 "{dir}/init_db"
  output: temp("{dir}/export_mst_for_all_lineages.query.temp"),
  	  temp("{dir}/export_mst_for_all_lineages.query.field_dict.json.temp")
  threads: 1
  resources: mysql=1
  params: name='export_mst_for_all_lineages_query', partition="general", mem="4096"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {EXPORT_MST_ALL_LINEAGES_QUERY} {input[0]} {output[0]} {output[1]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db}")

rule export_mst_for_all_lineages:
  input: "{dir}/export_mst_for_all_lineages.query.temp",
         "{dir}/export_mst_for_all_lineages.query.field_dict.json.temp"
  output: "{dir}/mst.all_lineages.nodes.cytoscape",
  	  "{dir}/mst.all_lineages.edges.cytoscape"
  threads: 1
  params: name='export_mst_for_all_lineages', partition="general", mem="4096"
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {EXPORT_MST_ALL_LINEAGES} {input_on_scratch[0]} {input[1]}")
    cp_from_scratch(output, scratch)

rule export_mst_no_isotype_constraint_for_all_lineages_query:
  input: "{dir}/uids.txt",
  	 "{dir}/init_db"
  output: temp("{dir}/export_mst_no_isotype_constraint_for_all_lineages.query.temp"),
  	  temp("{dir}/export_mst_no_isotype_constraint_for_all_lineages.query.field_dict.json.temp")
  threads: 1
  resources: mysql=1
  params: name='export_mst_no_isotype_constraint_for_all_lineages_query', partition="general", mem="4096"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {EXPORT_MST_NO_ISOTYPE_CONSTRAINT_ALL_LINEAGES_QUERY} {input[0]} {output[0]} {output[1]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db}")

rule export_mst_no_isotype_constraint_for_all_lineages:
  input: "{dir}/export_mst_no_isotype_constraint_for_all_lineages.query.temp",
         "{dir}/export_mst_no_isotype_constraint_for_all_lineages.query.field_dict.json.temp"
  output: "{dir}/mst_no_isotype_constraint.all_lineages.nodes.cytoscape",
  	  "{dir}/mst_no_isotype_constraint.all_lineages.edges.cytoscape"
  threads: 1
  params: name='export_mst_no_isotype_constraint_for_all_lineages', partition="general", mem="4096"
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {EXPORT_MST_NO_ISOTYPE_CONSTRAINT_ALL_LINEAGES} {input_on_scratch[0]} {input[1]}")
    cp_from_scratch(output, scratch)



rule msa_query:
  input: "{dir}/uids.txt",
  	 "{dir}/init_db"
  output: temp("{dir}/msa.query.temp.{n}"),
  	  temp("{dir}/msa.query.field_dict.json.temp.{n}")
  threads: 1
  resources: mysql=1
  params: name='msa_query', partition=config["params_partition_msa_query"], mem=config["params_mem_msa_query"]
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {MSA_QUERY} {input[0]} {num_workers} {output[0]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db}")

rule msa_worker:
  input: "{dir}/msa.query.temp.{n}",
  	 "{dir}/msa.query.field_dict.json.temp.{n}"
  output: "{dir}/msa.lineages.temp.{n}"
  threads: 1
  params: name='msa_worker', partition=config["params_partition_msa_worker"], mem=config["params_mem_msa_worker"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)    
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           export PATH={BLASTN_DIR}:$PATH &&\
           python {MSA} {input_on_scratch[0]} {input[1]} {output_on_scratch[0]} {path_to_muscle_bin} \
           {path_to_V_germline_ref_file} {path_to_J_germline_ref_file}")
    cp_from_scratch(output, scratch)
    cp_dir_from_scratch(scratch, "alignments", output[0])
    shell("rm {input}")

rule msa:
  input: "{dir}/msa.lineages.temp"
  output: "{dir}/msa"
  threads: 1
  resources: mysql=20
  params: name='msa', partition="general", mem="1024"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {UPDATE_DB} {input[0]} lineages {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} &&\
           touch {output[0]}")

rule tree_query:
  input: "{dir}/uids.txt",
  	 "{dir}/muscle_profile_align_germline"
  output: temp("{dir}/tree.query.temp.{n}"),
  	  temp("{dir}/tree.query.field_dict.json.temp.{n}")
  threads: 1
  resources: mysql=1
  params: name='tree_query', partition=config["params_partition_tree_query"], mem=config["params_mem_tree_query"]
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {TREE_QUERY} {input[0]} {num_workers} {output[0]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db}")

rule gapless_align_query:
  input: "{dir}/uids.txt"
  output: temp("{dir}/gapless_align.query.temp.{n}"),
  	  temp("{dir}/gapless_align.query.field_dict.json.temp.{n}")
  threads: 1
  resources: mysql=1
  params: name='gapless_align_query', partition="general", mem="2048"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {GAPLESS_ALIGN_QUERY} {input[0]} {num_workers} {output[0]} {output[1]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} {alignment_filter_stringency}")

rule gapless_align_worker:
  input: "{dir}/gapless_align.query.temp.{n}",
  	 "{dir}/gapless_align.query.field_dict.json.temp.{n}"
  output: "{dir}/gapless_align.{n, [0-9]+}"
  threads: 1
  params: name='gapless_align_worker', partition=config["params_partition_gapless_align_worker"], mem=config["params_mem_gapless_align_worker"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)    
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {GAPLESS_ALIGN} {input_on_scratch[0]} {input[1]} &&\
           touch {output_on_scratch[0]}")
    cp_from_scratch(output, scratch)
    cp_dir_from_scratch(scratch, "alignments", output[0])

rule gapless_align:
  input: ["{dir}/gapless_align.%s" % x for x in range(num_workers)]
  output: "{dir}/gapless_align"
  threads: 1
  params: name='gapless_align', partition="general", mem="1024"
  shell: "touch {output[0]} &&\
          rm {input}"

rule get_germline_seq_query:
  input: "{dir}/uids.txt"
  output: temp("{dir}/get_germline_seq.query.temp.{n}"),
  	  temp("{dir}/get_germline_seq.query.field_dict.json.temp.{n}")
  threads: 1
  resources: mysql=1
  params: name='get_germline_seq_query', partition="general", mem="2048"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {GET_GERMLINE_SEQ_QUERY} {input[0]} {num_workers} {output[0]} {output[1]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} {alignment_filter_stringency}")

rule get_germline_seq_worker:
  input: "{dir}/get_germline_seq.query.temp.{n}",
  	 "{dir}/get_germline_seq.query.field_dict.json.temp.{n}"
  output: "{dir}/get_germline_seq.{n, [0-9]+}"
  threads: 1
  params: name='get_germline_seq_worker', partition=config["params_partition_get_germline_seq_worker"], mem=config["params_mem_get_germline_seq_worker"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)    
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           export PATH={BLASTN_DIR}:$PATH &&\
           python {GET_GERMLINE_SEQ} {input_on_scratch[0]} {input[1]} \
           {path_to_V_germline_ref_file} {path_to_J_germline_ref_file} &&\
           touch {output_on_scratch[0]}")
    cp_from_scratch(output, scratch)
    cp_dir_from_scratch(scratch, "alignments", output[0])

rule get_germline_seq:
  input: ["{dir}/get_germline_seq.%s" % x for x in range(num_workers)]
  output: "{dir}/get_germline_seq"
  threads: 1
  params: name='get_germline_seq', partition="general", mem="1024"
  shell: "touch {output[0]} &&\
          rm {input}"

rule muscle_profile_align_germline_query:
  input: "{dir}/uids.txt",
         "{dir}/gapless_align",
         "{dir}/get_germline_seq"
  output: temp("{dir}/muscle_profile_align_germline.query.temp.{n}"),
  	  temp("{dir}/muscle_profile_align_germline.query.field_dict.json.temp.{n}")
  threads: 1
  resources: mysql=1
  params: name='muscle_profile_align_germline_query', partition="general", mem="2048"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {MUSCLE_PROFILE_ALIGN_GERMLINE_QUERY} {input[0]} {num_workers} {output[0]} {output[1]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db}")

rule muscle_profile_align_germline_worker:
  input: "{dir}/muscle_profile_align_germline.query.temp.{n}",
  	 "{dir}/muscle_profile_align_germline.query.field_dict.json.temp.{n}"
  output: "{dir}/muscle_profile_align_germline.lineages.temp.{n}"
  threads: 1
  params: name='muscle_profile_align_germline_worker', partition=config["params_partition_muscle_profile_align_germline_worker"], mem=config["params_mem_muscle_profile_align_germline_worker"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)    
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {MUSCLE_PROFILE_ALIGN_GERMLINE} {input_on_scratch[0]} {input[1]} {output_on_scratch[0]} \
           {path_to_muscle_bin}")
    cp_from_scratch(output, scratch)
    cp_dir_from_scratch(scratch, "alignments", output[0])

rule muscle_profile_align_germline:
  input: "{dir}/muscle_profile_align_germline.lineages.temp"
  output: "{dir}/muscle_profile_align_germline"
  threads: 1
  resources: mysql=20
  params: name='muscle_profile_align_germline', partition="general", mem="1024"
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {UPDATE_DB} {input[0]} lineages {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} &&\
	   touch {output[0]}")

rule mst_query:
  input: "{dir}/uids.txt",
  	 "{dir}/muscle_profile_align_germline"
  output: temp("{dir}/mst.query.temp.{n}"),
  	  temp("{dir}/mst.query.field_dict.json.temp.{n}")
  threads: 1
  resources: mysql=1
  params: name='mst_query', partition=config["params_partition_mst_query"], mem=config["params_mem_mst_query"]
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           echo input: {input} && \
	   echo output: {output} && \
           python {MST_QUERY} {input[0]} {num_workers} {output[0]} {output[1]} \
           {mysql_db_host} {mysql_db_user} {mysql_db_password} {db}")

rule mst_worker:
  input: "{dir}/mst.query.temp.{n}",
  	 "{dir}/mst.query.field_dict.json.temp.{n}"
  output: "{dir}/mst.sequences.temp.{n}",
  	  "{dir}/mst.lineages.temp.{n}"
  threads: 1
  params: name='mst_worker', partition=config["params_partition_mst_worker"], mem=config["params_mem_mst_worker"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)    
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {MST} {input_on_scratch[0]} {input[1]} {path_to_edmonds_bin} 1")
    cp_from_scratch(output, scratch)
    cp_dir_from_scratch(scratch, "msts", output[0])

rule mst:
  input: "{dir}/mst.sequences.temp",
         "{dir}/mst.lineages.temp"
  output: "{dir}/mst"
  threads: 1
  resources: mysql=20
  params: name='mst', partition=config["params_partition_mst"], mem=config["params_mem_mst"]
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {UPDATE_DB} {input[0]} sequences {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} &&\
           python {UPDATE_DB} {input[1]} lineages {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} &&\
	   touch {output[0]}")

rule mst_no_isotype_constraint_worker:
  input: "{dir}/mst.query.temp.{n}",
  	 "{dir}/mst.query.field_dict.json.temp.{n}"
  output: "{dir}/mst_no_isotype_constraint.sequences.temp.{n}",
  	  "{dir}/mst_no_isotype_constraint.lineages.temp.{n}"
  threads: 1
  params: name='mst_no_isotype_constraint_worker', partition=config["params_partition_mst_worker"], mem=config["params_mem_mst_worker"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)    
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {MST} {input_on_scratch[0]} {input[1]} {path_to_edmonds_bin} 0")
    cp_from_scratch(output, scratch)
    cp_dir_from_scratch(scratch, "msts_no_isotype_constraint", output[0])

rule mst_no_isotype_constraint:
  input: "{dir}/mst_no_isotype_constraint.sequences.temp",
         "{dir}/mst_no_isotype_constraint.lineages.temp"
  output: "{dir}/mst_no_isotype_constraint"
  threads: 1
  resources: mysql=20
  params: name='mst_no_isotype_constraint', partition=config["params_partition_mst"], mem=config["params_mem_mst"]
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {UPDATE_DB} {input[0]} sequences {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} &&\
           python {UPDATE_DB} {input[1]} lineages {mysql_db_host} {mysql_db_user} {mysql_db_password} {db} &&\
	   touch {output[0]}")
