# Genome-Assembly-Snakemake

## What is this workflow?
The aim of this snakemake workflow is to assemble genomes of any ploidy, and where HiC data is also available it will scaffold the assemblies. The final output for samples with HiC is a pretext map, and for samples without HiC it is the fasta files from hifiasm.

## Input Data
For this workflow we will need to fill in a config.ini file, and samples.csv file. Examples of both files can be found in the configs folder. 

### config ini
The config contains different parameters that are essential to the pipeline running properly. Anything that is commented out is not essential, but can be uncommented and given a value if needed. This workflow does not require you to give paths to individual fastq files, or paired end reads. 

```
[DEFAULT]
# samples file
samples_file= /path/to/samples/file


[PREPARING_HIFI_READS]
# variables which determine if we want to filter hifi reads first or not
filter_hifi_reads = False
hifiadapterfiltDB_path = /path/to/hifiadapterfiltDB

[HIFIASM]
hifi_reads = /path/to/bam/files/dir
hic1_reads = /path/to/HiC1_reads/dir # you can have more than one file in the directory for any given sample and it will be merged
hic2_reads = /path/to/HiC2_reads/dir
homcov = 1
use_genomescope_homcov = False

[QUALITY_CONTROL]
compleasm_download_path = /path/to/lineage/dir
compleasm_lineage = brassicales
compleasm_odb = odb12

[SCAFFOLDING]
telomere_motif = TTTAGGG

[SLURM_ARGS]
# set mem (inMb), tmpdir, runtime (in minutes), allocation for jobs submitted
cpus_per_task = 16
mem_Mb = 32000
runtime = 120
slurm_partition = pelle
slurm_account = slu2026-1-5
tmpdir = /scratch
```

`hifi_reads` : A path to a single directory containing hifi reads in .bam, .fastq, or .fastq.gz should be given.

`hic1_reads` : A path to a directory containing the population name as a subdirectory, and within that subdirectory, the hic1 reads. You can either merge your HiC data before hand and place it in this subdirectory, or have data from several runs in the subdirectory and the workflow will merge them together. But **DO NOT** keep both the individual files and the merged file together in the subdirectory, otherwise the individual files and the merged file will be merged again.

`hic2_reads` : A path to a directory containing the population name as a subdirectory, and within that subdirectory, the hic2 reads. You can either merge your HiC data before hand and place it in this subdirectory, or have data from several runs in the subdirectory and the workflow will merge them together. But **DO NOT** keep both the individual files and the merged file together in the subdirectory, otherwise the individual files and the merged file will be merged again.

`samples file` : Path to the samples.csv file. If not given then it will be assumed that it is in the output directory
