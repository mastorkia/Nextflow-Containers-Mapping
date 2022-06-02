#!/bin/bash
#
#$
#$ -cwd
#$ -j y
#$ -N illumina_nextflow
#$ -S /bin/bash
#$ -l h_vmem=80G

## backup dir format ##
temp_dir=$(date '+%d_%m_%y_')$RANDOM

module load /usr/bin/singularity
/share/data/software/nextflow/nextflow run Nextflow_v6.nf -with-trace -w temp_${temp_dir}

#done
