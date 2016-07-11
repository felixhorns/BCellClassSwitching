#!/bin/bash

source /local10G/rfhorns/resources/anaconda/envs/py3.4/bin/activate /local10G/rfhorns/resources/anaconda/envs/py3.4/

MY_HOME=/local10G/rfhorns/Bcell_class_switching_V1_v3 # set me to home directory for analysis
export CONFIGFILE=$MY_HOME/pipelines/igh_makedb/config.json # set to config file   

SNAKEFILE=/datastore/rfhorns/scripts/igh_makedb/Snakefile.py

# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$MY_HOME/log/igh_makedb.$DATETIME.log

# Make seedfile
SEEDFILE=$MY_HOME/pipelines/igh_makedb/seedfile_igh_makedb.txt
MAKE_SEEDFILE=/datastore/rfhorns/scripts/igh_makedb/make_seedfile_igh_makedb.py
# python $MAKE_SEEDFILE $MY_HOME/lineages sequences_lineages.txt $SEEDFILE

# Dry run snakemake
# snakemake all --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files -j 100 -w 100 -k -n

# Run snakemake
nohup snakemake all --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files -j 100 -w 100 -k > $LOGFILE &

echo Log is
echo $LOGFILE
echo
