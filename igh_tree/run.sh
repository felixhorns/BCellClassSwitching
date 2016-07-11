#!/bin/bash

source /local10G/rfhorns/resources/anaconda/envs/py3.4/bin/activate /local10G/rfhorns/resources/anaconda/envs/py3.4/

MY_HOME=/local10G/rfhorns/Bcell_class_switching_V1_v3 # set me to home directory for analysis
export CONFIGFILE=$MY_HOME/pipelines/igh_tree/config.json # set to config file   

SNAKEFILE=/b1_2/home/rfhorns/scripts/igh_tree/Snakefile.py

# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$MY_HOME/log/igh_tree.$DATETIME.log

# Make seedfile
SEEDFILE=$MY_HOME/pipelines/igh_tree/seedfile_igh_tree.txt
MAKE_SEEDFILE=/b1_2/home/rfhorns/scripts/igh_tree/make_seedfile_igh_tree.py
# python $MAKE_SEEDFILE $MY_HOME/lineages init_db $SEEDFILE

# Dry run snakemake
# snakemake all_mst -R mst_query --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k -n

# Run snakemake
nohup snakemake all_mst -R mst_query --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k > $LOGFILE &

# Run snakemake to make Cytoscape file for visualizing MSTs
# nohup snakemake all_export_mst_no_isotype_constraint_for_all_lineages --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k > $LOGFILE &

# nohup snakemake all_export_mst_for_all_lineages --forceall --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k > $LOGFILE &

echo Log is
echo $LOGFILE
echo
