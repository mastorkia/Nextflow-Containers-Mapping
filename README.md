# Nextflow-Containers-Mapping
Small pipeline for QC, trimming and mapping using single-end fastq files. Use of containarized images for that.

### Pipeline steps
* If files are coming from NextSeq500 concatenate if not leave it as it is.
* Adaptor trimming using cutadapt
* Count of adaptors trimmed
* Quality Control step using FastQC and MultiQC
* Trimm by quality of the bases using Prinseq
* Download accession from NCBI and concatante
* Index genome with bowtie2
* Mapp reads and convert them from sam to bam  

### Run Nextflow pipeline.

Nextflow pipiline will be run by ```nextflow.sh``` file:

```bash
cat nextflow.sh 
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
/share/data/software/nextflow/nextflow run Nextflow.nf -with-trace -w temp_${temp_dir}

#done

```

Usage of ```nextflow.sh``` is as follows :

```bash
 qsub nextflow.sh 
```

In our home directory we will need the main pipeline ```Nextflow.nf``` and the config file ```nextflow.config``` (keep this name)and ```file_list``` and ```scripts``` folders in which the main utilities will be stored.

Since this pipeline consists on contanaraized images we have created images of different softwares using Docker and then transform them to singularity for HPC purposes. These singularity images are located in:

```bash
 /share/data/software/singularity/
```
