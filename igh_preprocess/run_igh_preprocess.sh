#!/bin/bash

source /local10G/dcroote/resources/anaconda/envs/py3.4/bin/activate /local10G/dcroote/resources/anaconda/envs/py3.4/

MY_HOME=/local10G/rfhorns/Bcell_class_switching # set me to home directory for analysis
export CONFIGFILE=$MY_HOME/pipelines/igh_preprocess/config.json # set to config file

SNAKEFILE=/b1_2/home/rfhorns/scripts/igh_preprocess/Snakefile.py

# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$MY_HOME/log/igh_preprocess.$DATETIME.log

# Make seedfile
SEEDFILE=$MY_HOME/pipelines/igh_preprocess/seedfile_igh_preprocess.txt
MAKE_SEEDFILE_SCRIPT=/b1_2/home/rfhorns/scripts/igh_preprocess/make_seedfile_igh_preprocess.py
# python $MAKE_SEEDFILE_SCRIPT $MY_HOME/fastq/ $SEEDFILE &> $MY_HOME/log/make_seedfile_igh_pipeline.log

# Dry run snakemake
# snakemake all --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k -n

# Run snakemake
nohup snakemake all --forceall --snakefile $SNAKEFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 100 -w 100 -k > $LOGFILE &

echo Log is
echo $LOGFILE
echo
