#!/bin/bash

# source and actiavte envs
source $HOME/.bash_profile

# output dir for the analysis
OUTPUT_DIRECTORY=~/LyrenosaPangenome/Snakemake_v6

# snakefile directory
SNAKEFILE_DIRECTORY=~/Software/Snakemake_Pipes/Genome-Assembly-Snakemake




# move into output directoru
mkdir -p $OUTPUT_DIRECTORY
cd $OUTPUT_DIRECTORY

# export variables used by the pipeline

# give the internal rules the output directory snakemake is using
export SNAKEMAKE_OUTPUTDIR=$OUTPUT_DIRECTORY

# give the internal rules the name for your anlysis, this is used as a prefix for compiled results
export ANALYSIS_NAME="test"

# activate env with snakemake, the yml should be in the envs directory
conda activate snakemake-tools

# create config and samples file
> config.ini
> samples.csv

# start building config ini
cat <<EOF >> config.ini

[PREPARING_HIFI_READS]
# variables which determine if we want to filter hifi reads first or not
filter_hifi_reads = False
#hifiadapterfiltDB_path = /path/to/hifiadapterfiltDB

[HIFIASM]
#use_genomescope_homcov = False

[QUALITY_CONTROL]
compleasm_download_path = /path/to/lineage/dir

[SCAFFOLDING]
#telomere_motif = TTTAGGG

[SLURM_ARGS]
# set mem (inMb), tmpdir, runtime (in minutes), allocation for jobs submitted
cpus_per_task = 16
mem_Mb = 32000
runtime = 120
slurm_partition = your_slurm_partition
slurm_account = your_slurm_account
tmpdir = /scratch

EOF

# start building sample csv file
cat <<EOF >> samples.csv
sample,HiFi,HiC1,HiC2,ploidy,chrom,haplomeSize,compleasm_lineage_odb
POPULATION1,/path/to/hifi/reads/{.bam or .fastq or.fastq.gz},,,2,8,210000000,brassicales_odb12
POPULATION2,/path/to/hifi/reads/{.bam or .fastq or .fastq.gz},/path/to/HiC1/reads/{.fastq.gz or .fq.gz},/path/to/HiC1/reads/{.fastq.gz or .fq.gz},2,8,210000000,brassicales_odb12
EOF

cd $SNAKEFILE_DIRECTORY

snakemake -p --executor slurm --directory $OUTPUT_DIRECTORY --jobs 10 --cores 96 --latency-wait 100 --restart-times 4 --sdm conda --conda-frontend conda


#snakemake --directory $BRAKER_DIR --dag | dot -Tsvg > dag.svg

# create regular dag with all samples
#snakemake --directory $BRAKER_DIR --dag | dot -Tpdf > dags/dag_regular.pdf

# create rule graph DAG
#snakemake --directory $BRAKER_DIR --rulegraph | dot -Tpdf > dag_rulegraph.pdf

# do a dry run to see how the pipe will run
#snakemake --directory $BRAKER_DIR -n -p

conda deactivate