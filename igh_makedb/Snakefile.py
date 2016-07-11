import os
import re
import subprocess

##### Paths
RESOURCES=	                                '/local10G/rfhorns/resources'
ANACONDA=					RESOURCES+'/anaconda'

SCRIPTS=					'/datastore/rfhorns/scripts/igh_makedb'
configfile:					os.environ.get("CONFIGFILE")

FILTER=                                         SCRIPTS+'/filter.py'
COLLAPSE_IDENTICAL_SEQUENCES=			config["collapse_identical_sequences_script"]
CLEAN_UIDS=           			        SCRIPTS+'/clean_uids.py'
MAKE_INIT_DB_INPUT=           			SCRIPTS+'/make_init_db_input.py'
INIT_DB=					SCRIPTS+'/init_db.sh'

##### Parameters
workdir:                                        config["workdir"] + "/log"
SEEDFILE=                                       config["seedfile"]

num_workers=					64
scratch_storage=                                config["scratch_storage"]

patient_name_to_uid_map=			config["patient_name_to_uid_map"]
patient_data=					config["patient_data"]
year_visit_str_to_uid_map=			config["year_visit_str_to_uid_map"]
subsampling_p_to_uid_map=			config["subsampling_p_to_uid_map"]
clustering_params_to_uid_map=			config["clustering_params_to_uid_map"]
sequences_lineages_field_dict=			config["sequences_lineages_field_dict"]

filter_nonproductive=                           config["filter_nonproductive"]
filter_low_conf=                                config["filter_low_conf"]

db=						config["db"]
mysql_db_host=					config["mysql_db_host"]
mysql_db_user=					config["mysql_db_user"]
mysql_db_password=				config["mysql_db_password"]

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

rule all:
  input: expand("{dir}/init_db", dir=SEEDS),
  params: name='all', partition='general', mem='1024'

rule do_filter:
  input: "{dir}/sequences_lineages.txt"
  output: "{dir}/sequences_lineages.filtered.txt",
          "{dir}/sequences_lineages.nonproductive.txt",
          "{dir}/sequences_lineages.lowconf.txt"
  threads: 1
  params: name='filter', partition="general", mem="2048"
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} &&\
          python {FILTER} {input[0]} {output[0]} {output[1]} {output[2]}\
          {filter_nonproductive} {filter_low_conf}'

rule collapse_identical_sequences:
  input: "{dir}/sequences_lineages.filtered.txt"
  output: "{dir}/sequences_lineages.filtered.identicals_collapsed.txt"
  threads: 1
  params: name='collapse_identical_sequences', partition=config["params_partition_collapse_identical_sequences"], mem=config["params_mem_collapse_identical_sequences"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} &&\
  	  python {COLLAPSE_IDENTICAL_SEQUENCES} {input[0]} {output[0]}'

rule clean_uids:
  input: "{dir}/sequences_lineages.filtered.identicals_collapsed.txt"
  output: "{dir}/sequences_lineages.filtered.identicals_collapsed.uids_cleaned.txt"
  threads: 1
  params: name='clean_uids', partition=config["params_partition_clean_uids"], mem=config["params_mem_clean_uids"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} &&\
  	  python {CLEAN_UIDS} {input[0]} {sequences_lineages_field_dict} {output[0]}'

rule init_db:
  input: "{dir}/sequences_lineages.filtered.identicals_collapsed.uids_cleaned.txt",
  	 "{dir}/lib_info_combined.txt",
  	 "{dir}/clustering_params.txt",
  output: "{dir}/init_db",
          temp("{dir}/temp.sql"),
  	  "{dir}/uids.txt",
  	  temp("{dir}/patients.txt.temp"),
  	  temp("{dir}/clusterings.txt.temp"),
  	  temp("{dir}/libs.txt.temp"),
  	  temp("{dir}/sequences.txt.temp"),
  	  temp("{dir}/sequence_strings.txt.temp"),
  	  temp("{dir}/quality_scores.txt.temp"),
  	  temp("{dir}/lineages.txt.temp")
  threads: 1
  resources: mysql=1
  params: name='init_db', partition=config["params_partition_init_db"], mem=config["params_mem_init_db"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           python {MAKE_INIT_DB_INPUT} {input_on_scratch[0]} {input[0]} \
	   {patient_name_to_uid_map} {patient_data} {year_visit_str_to_uid_map} \
	   {subsampling_p_to_uid_map} {clustering_params_to_uid_map} {sequences_lineages_field_dict}")
    cp_from_scratch(output[2:], scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} &&\
           source {INIT_DB} {db} {input[0]} &&\
           mysql  --host={mysql_db_host} -u {mysql_db_user} --password={mysql_db_password} < {output[1]} &&\
           touch {output[0]}")
