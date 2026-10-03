### HIC FILE SORTING ###
def get_HiC1_input_files(wildcards):

    if HASHIC and MERGEHIC[wildcards.sample]:

        return [f for f in os.listdir(os.path.join(config['hic1_reads'],s)) if 'merged' in f][0]

    elif HASHIC and not MERGEHIC[wildcards.sample]:

        return HICFILES[wildcards.sample]["HIC1"][0]

    elif not HASHIC and not MERGEHIC[wildcards.sample]:

        return []

def get_HiC2_input_files(wildcards):

    if HASHIC and MERGEHIC[wildcards.sample]:

        return [f for f in os.listdir(os.path.join(config['hic2_reads'],s)) if 'merged' in f][0]

    elif HASHIC and not MERGEHIC[wildcards.sample]:

        return HICFILES[wildcards.sample]["HIC2"][0]

    elif not HASHIC and not MERGEHIC[wildcards.sample]:

        return []


rule get_input_data_stats:
    input:
        hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        hic1=get_HiC1_input_files,
        hic2=get_HiC2_input_files
    output:
        hifiTsv="{sample}/s2_input_data_qc/Hifi/{sample}.hifi.tsv",
        hicTsv1="{sample}/s2_input_data_qc/HiC1/{sample}.hic1.tsv" if HASHIC else [],
        hicTsv2="{sample}/s2_input_data_qc/HiC2/{sample}.hic2.tsv" if HASHIC else []
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=120,
        mem_mb=10000
    params:
        hic1StatsCommand="seqkit stats {input.hic1}/* --threads {threads} -Ta -o {output.hicTsv1}" if HASHIC else [],
        hic2StatsCommand="seqkit stats {input.hic2}/* --threads {threads} -Ta -o {output.hicTsv2}" if HASHIC else []
    shell:
        '''
        seqkit stats {input.hifiFastq} --threads {threads} -Ta -o {output.hifiTsv}

        # get HiC1 stats
        {params.hic1StatsCommand}

        # get HiC2 stats
        {params.hic2StatsCommand}'''
