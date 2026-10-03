import configparser
import sys
import os
import pandas as pd

# set config ini path
# path is relative and set relative to working directory
# you can set environemntal variable called CONFIG_INI_PATH if the file 'config.ini' is not in your working directory
configIniPath = os.environ.get('CONFIG_INI_PATH','config.ini')

# load in config file, 
configIni = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))

# read in config gile
configIni.read(configIniPath)

# build snakemake config
### DEFAULTS ###
config['ploidy'] = configIni['DEFAULT'].getint('ploidy',fallback=2)

### PREPARING HIFI READS
config['filter_hifi_reads'] = configIni['PREPARING_HIFI_READS'].getboolean('filter_hifi_reads',fallback=False)
config['hifiadapterfiltDB_path'] = configIni['PREPARING_HIFI_READS'].get('hifiadapterfiltDB_path',fallback='')

### HIFIASM
config['hifi_reads'] = configIni['HIFIASM'].get('hifi_reads',fallback='')
config['hic1_reads'] = configIni['HIFIASM'].get('hic1_reads',fallback='')
config['hic2_reads'] = configIni['HIFIASM'].get('hic2_reads',fallback='')
config['homcov'] = configIni['HIFIASM'].getint('homcov',fallback=0)
config['use_genomescope_homcov'] = configIni['HIFIASM'].getboolean('use_genomescope_homcov',fallback=False)

### QUALITY CONTROL
config['compleasm_download_path'] = configIni['QUALITY_CONTROL'].get('compleasm_download_path',fallback='')
config['compleasm_lineage'] = configIni['QUALITY_CONTROL'].get('compleasm_lineage',fallback='brassicales')
config['compleasm_odb'] = configIni['QUALITY_CONTROL'].get('compleasm_odb',fallback='odb12')

### SCAFFOLDING
config['telomere_motif'] = configIni['SCAFFOLDING'].getboolean('telomere_motif',fallback='TTTAGGG')

### SLURM ARGS
config['mem_mb'] = configIni['SLURM_ARGS'].getint('mem_Mb',fallback=32000)
config['cpus_per_task'] = configIni['SLURM_ARGS'].getint('cpus_per_task',fallback=16)
config['runtime'] = configIni['SLURM_ARGS'].getint('runtime',fallback=120)
config['slurm_partition'] = configIni['SLURM_ARGS'].get('slurm_partition',fallback='pelle')
config['slurm_account'] = configIni['SLURM_ARGS'].get('slurm_account',fallback='slu2026-1-5')
config['tmpdir'] = configIni['SLURM_ARGS'].get('tmpdir',fallback='/scratch')

### set other global variables
config['script_dir'] = os.path.join(workflow.basedir,'scripts')
config['shared_data_dir'] = os.path.join(workflow.basedir,'shared_data')

# include rules

# shared rules
include:"rules/common/common.smk"
include:"rules/common/end_results.smk"

### PREPARING INPUT DATA ###
# bam to fastq
include:"rules/prepare_input_data/bam_to_fastq.smk"

# merge HiC
if HASHIC:
    include:"rules/prepare_input_data/merge_hic.smk"

# input data stats
include:"rules/prepare_input_data/input_data_stats.smk"

### ESTIMATING PLOIDY ###
# genomescope
include:"rules/estimate_ploidy/run_genomescope.smk"

# smudgeplot
include:"rules/estimate_ploidy/run_smudgeplot.smk"

### HIFIASM ###
include:"rules/run_hifiasm/run_hifiasm.smk"
include:"rules/run_hifiasm/gfa_to_fa.smk"
include:"rules/run_hifiasm/run_hifiasm_qc.smk"



# all rule
rule all:
    input:
        set_rule_target()
