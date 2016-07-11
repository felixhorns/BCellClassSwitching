import os
import subprocess

##### Paths
RESOURCES=      				'/local10G/rfhorns/resources'
ANACONDA=       				RESOURCES+'/anaconda'

SCRIPTS=                        		'/b1_2/home/rfhorns/scripts/igh_cluster'
configfile:                                     os.environ.get("CONFIGFILE")

CP_LIB_INFO=					SCRIPTS+'/cp_lib_info.py'
FIND_GROUPS=		 			SCRIPTS+'/find_similar_V_J_CDR3_length_groups.py'
BALANCE_LOAD=					SCRIPTS+'/balance_load.py'
DO_CLUSTERING=                     		SCRIPTS+'/cluster_lineages.py'
CALC_DISTANCES=					SCRIPTS+'/calc_pairwise_distances.py'
PLOT_DISTANCES=					SCRIPTS+'/plot_distances.py'
GET_RECORDS=					SCRIPTS+'/get_records.py'

##### Parameters
workdir:                                        config["workdir"] + "/log"
SEEDFILE=                                       config["seedfile"]

num_workers=					config["num_workers"]
scratch_storage=                                config["scratch_storage"]

clustering_V_allele_collapse_dict_file=         config["clustering_V_allele_collapse_dict_file"]

##### Function for loading seeds
def load_output_dirs(infile):
    output_dirs = []
    with open(infile, 'rU') as f:
        for line in f:
            output_dir = line.rstrip().split('\t')[2].strip()
            output_dirs.append(output_dir)
    return output_dirs

OUTPUT_DIRS = load_output_dirs(SEEDFILE)

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
rule all:
  input: expand("{dir}/sequences_lineages.txt", dir=OUTPUT_DIRS),
  	 expand("{dir}/lib_info_combined.txt", dir=OUTPUT_DIRS)
#	 expand("{dir}/CDR3_templated_identity.pdf", dir=OUTPUT_DIRS) # uncomment to plot identity to germline (>2 hrs per worker)
  params: name='all', partition='debug', mem='1024'

rule plot_all:
  input: expand("{dir}/CDR3_templated_identity.pdf", dir=OUTPUT_DIRS)
  params: name='plot_all', partition='general', mem='1024'

rule calc_distances_all:
  input: expand("{dir}/calc_distances", dir=OUTPUT_DIRS)
  params: name='calc_distances_all', partition='general', mem='1024'

rule make_lib_clustering_info:
  input: '{dir}/clustering_inputs.txt',
  	 '{dir}/clustering_params.txt'
  output: '{dir}/lib_info_combined.txt',
  threads: 1
  params: name="make_lib_clustering_info", partition=config["params_partition_make_lib_clustering_info"], mem=config["params_mem_make_lib_clustering_info"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} && \
  	  echo {clustering_V_allele_collapse_dict_file} >> {input[1]} && \
  	  python {CP_LIB_INFO} {input[0]}'

rule find_groups:
  input: '{dir}/clustering_inputs.txt'
  output: '{dir}/groups.zip.temp',
  	  '{dir}/seq_records.txt.temp'
  threads: 1
  params: name='find_groups', partition=config["params_partition_find_groups"], mem=config["params_mem_find_groups"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} && \
  	  python {FIND_GROUPS} {input[0]} {clustering_V_allele_collapse_dict_file}'

rule balance_load:
  input: '{dir}/groups.zip.temp'
  output: '{dir}/group_worker_assignments.txt.temp'
  threads: 1
  params: name='balance_load', partition=config["params_partition_balance_load"], mem=config["params_mem_balance_load"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} && \
  	  python {BALANCE_LOAD} {input[0]} {num_workers}'

rule do_clustering:
  input: '{dir}/groups.zip.temp',
  	 '{dir}/group_worker_assignments.txt.temp',
	 '{dir}/clustering_params.txt'
  output: temp('{dir}/clusters.txt.temp.{n}')
  threads: 1
  params: name='do_clustering', partition=config["params_partition_do_clustering"], mem=config["params_mem_do_clustering"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} && \
           python {DO_CLUSTERING} {input_on_scratch[0]} {input_on_scratch[1]} \
	   	  {output_on_scratch[0]} {input_on_scratch[2]}")
    cp_from_scratch(output, scratch)

rule cat_clusters:
  input: ['{dir}/clusters.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/clusters.txt'
  threads: 1
  params: name='cat_clusters', partition=config["params_partition_cat_clusters"], mem=config["params_mem_cat_clusters"]
  shell: "cat {input} > {output}"
#  	  rm {input}"

rule parse_clusters:
  input: '{dir}/clusters.txt',
  	 '{dir}/seq_records.txt.temp'
#	 '{dir}/groups.zip.temp',
#	 '{dir}/group_worker_assignments.txt.temp'
  output: '{dir}/sequences_lineages.txt'
  threads: 1
  params: name='parse_clusters', partition=config["params_partition_parse_clusters"], mem=config["params_mem_parse_clusters"]
  shell: 'source {ANACONDA}/bin/activate {ANACONDA} && \
  	  python {GET_RECORDS} {input[0]} {input[1]}'
#	  rm {input[1]} {input[2]} {input[3]}'

rule calc_distances:
  input: '{dir}/groups.zip.temp',
  	 '{dir}/group_worker_assignments.txt.temp',
	 '{dir}/clustering_params.txt'
  output: temp('{dir}/CDR3_distances.txt.temp.{n}'),
  	  temp('{dir}/templated_distances.txt.temp.{n}'),
          temp('{dir}/min_CDR3_distances.txt.temp.{n}'),
  	  temp('{dir}/min_templated_distances.txt.temp.{n}'),
          temp('{dir}/CDR3_similarities.txt.temp.{n}'),
  	  temp('{dir}/templated_similarities.txt.temp.{n}'),
          temp('{dir}/min_CDR3_similarities.txt.temp.{n}'),
  	  temp('{dir}/min_templated_similarities.txt.temp.{n}')
  threads: 1
  params: name='calc_distances', partition=config["params_partition_calc_distances"], mem=config["params_mem_calc_distances"]
  run:
    scratch = os.environ[scratch_storage]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    shell("source {ANACONDA}/bin/activate {ANACONDA} && \
           python {CALC_DISTANCES} {input_on_scratch[0]} {input_on_scratch[1]} \
	   	  {output_on_scratch[0]} {input_on_scratch[2]}")
    cp_from_scratch(output, scratch)

rule cat_CDR3_distances:
  input: ['{dir}/CDR3_distances.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/CDR3_distances.txt'
  threads: 1
  params: name='cat_CDR3_distances', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_templated_distances:
  input: ['{dir}/templated_distances.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/templated_distances.txt'
  threads: 1
  params: name='cat_templated_distances', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_min_CDR3_distances:
  input: ['{dir}/min_CDR3_distances.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/min_CDR3_distances.txt'
  threads: 1
  params: name='cat_min_CDR3_distances', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_min_templated_distances:
  input: ['{dir}/min_templated_distances.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/min_templated_distances.txt'
  threads: 1
  params: name='cat_min_templated_distances', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_CDR3_similarities:
  input: ['{dir}/CDR3_similarities.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/CDR3_similarities.txt'
  threads: 1
  params: name='cat_CDR3_similarities', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_templated_similarities:
  input: ['{dir}/templated_similarities.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/templated_similarities.txt'
  threads: 1
  params: name='cat_templated_similarities', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_min_CDR3_similarities:
  input: ['{dir}/min_CDR3_similarities.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/min_CDR3_similarities.txt'
  threads: 1
  params: name='cat_min_CDR3_similarities', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_min_templated_similarities:
  input: ['{dir}/min_templated_similarities.txt.temp.%s' % i for i in range(num_workers)]
  output: '{dir}/min_templated_similarities.txt'
  threads: 1
  params: name='cat_min_templated_similarities', partition='general', mem='1024'
  shell: "cat {input} > {output}"

rule cat_distances:
  input: '{dir}/clustering_inputs.txt',
  	 '{dir}/CDR3_distances.txt',
  	 '{dir}/templated_distances.txt',
  	 '{dir}/min_CDR3_distances.txt',
  	 '{dir}/min_templated_distances.txt',
  	 '{dir}/CDR3_similarities.txt',
  	 '{dir}/templated_similarities.txt',
  	 '{dir}/min_CDR3_similarities.txt',
  	 '{dir}/min_templated_similarities.txt'
  output: '{dir}/calc_distances'
  params: name='calc_distances', partition='general', mem='1024'
  shell: "touch {output[0]}"

rule plot_distances:
  input: '{dir}/clustering_inputs.txt',
  	 '{dir}/CDR3_distances.txt',
  	 '{dir}/templated_distances.txt',
  	 '{dir}/min_CDR3_distances.txt',
  	 '{dir}/min_templated_distances.txt',
  	 '{dir}/CDR3_similarities.txt',
  	 '{dir}/templated_similarities.txt',
  	 '{dir}/min_CDR3_similarities.txt',
  	 '{dir}/min_templated_similarities.txt'
  output: '{dir}/CDR3_templated_identity.pdf',
  	  '{dir}/V_germline_identity.pdf',
	  '{dir}/mutation_density.pdf'
  params: name='all', partition='general', mem='4096'
  shell: "source {ANACONDA}/bin/activate {ANACONDA} && \
  	  python {PLOT_DISTANCES} {input[0]} {input[1]} {input[2]}"
