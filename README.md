# Genome-Assembly-Snakemake

## What is this workflow?
The aim of this snakemake workflow is to assemble genomes of any ploidy, and where HiC data is also available it will scaffold the assemblies. The final output for samples with HiC is a pretext map, and for samples without HiC it is the fasta files from hifiasm.

## Input Data
For this workflow we will need to fill in a config.ini file, and samples.csv file. Examples of both files can be found in the configs folder. 

### config ini
The config contains different parameters that are essential to the pipeline running properly. Anything that is commented out is not essential, but can be uncommented and given a value if needed. This workflow does not require you to give paths to individual fastq files, or paired end reads. 

```
[PREPARING_HIFI_READS]
# variables which determine if we want to filter hifi reads first or not
#filter_hifi_reads = False
#hifiadapterfiltDB_path = /path/to/hifiadapterfiltDB

[HIFIASM]
# use_genomescope_homcov = False

[QUALITY_CONTROL]
compleasm_download_path = /path/to/lineage/dir

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

### samples csv
The csv file contains columns that allow for specific operations to be performed on different populations in the analysis. Columns that should not be left empty are highlighted in bold. An example is provided below, and a description of the expected values in each column.

```
sample,HiFi,HiC1,HiC2,ploidy,chrom,haplomeSize,compleasm_lineage_odb
PopA,/absolute/path/to/bam/file,/absolute/path/to/HiC1/file,/absolute/path/to/HiC2/file,2,8,210000000,brassicales_odb12
```

`sample` : The name of a population/sample. This name does not have to be present in the file name of the hifi or hic reads

`HiFi` : The absolute path to a file with hifi reads. Valid file extensions are **.fastq**, .**fastq.gz**, **.fq**, **.fq.gz**, and **.bam**

`HiC1` : There are types of information that can be provided here. Either a path to a directory containing HiC reads, or a space separated list of file(s) with HiC reads can be given. If a space separated list of files is given, the order of the files for HiC1 and HiC2 does not matter. Where more than one file is provided, the data will be merged. Valid extensions for fastq files are **.fastq.gz** or **.fq.gz**

`HiC2` : There are types of information that can be provided here. Either a path to a directory containing HiC reads, or a space separated list of file(s) with HiC reads can be given. If a space separated list of files is given, the order of the files for HiC1 and HiC2 does not matter. Where more than one file is provided, the data will be merged. Valid extensions for fastq files are **.fastq.gz** or **.fq.gz**

`ploidy` : The ploidy of the population

`chrom` : The number of chromosomes per haplotype

`haplomeSize` : The expected haplome size (in bp)

`compleasm_lineage_db` : The lineage and odb to use for compleasm
