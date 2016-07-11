#!/bin/bash

source /local10G/dcroote/resources/anaconda/envs/py3.4/bin/activate /local10G/dcroote/resources/anaconda/envs/py3.4/

MY_HOME=/local10G/rfhorns/Bcell_class_switching # set me to home directory for analysis
export CONFIGFILE=$MY_HOME/pipelines/igh_cluster/config.json # set to config file   

SNAKEFILE=/b1_2/home/rfhorns/scripts/igh_cluster/Snakefile.py

# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$MY_HOME/log/igh_cluster.$DATETIME.log

# Make seedfile
SEEDFILE=$MY_HOME/pipelines/igh_cluster/seedfile_igh_cluster.txt
SEEDFILE_BEFORE_CLUSTER_PARAMS=$MY_HOME/pipelines/igh_cluster/seedfile_igh_cluster_original_dirs.txt
CLUSTER_PARAMS=$MY_HOME/pipelines/igh_cluster/all_cluster_params.txt
MAKE_SEEDFILE=/b1_2/home/rfhorns/scripts/igh_cluster/make_seedfile_cluster.py
# python $MAKE_SEEDFILE $SEEDFILE_BEFORE_CLUSTER_PARAMS $CLUSTER_PARAMS

# Split seedfile and put directives for clustering into each output directory
SPLIT_SEEDFILE=/b1_2/home/rfhorns/scripts/igh_cluster/split_clustering_seedfile.py
# python $SPLIT_SEEDFILE $SEEDFILE

# Dry run snakemake
# snakemake all --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k -n

# Run snakemake
nohup snakemake all --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k > $LOGFILE &

echo Log is
echo $LOGFILE
echo
