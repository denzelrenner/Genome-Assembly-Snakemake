rule get_input_hifi_stats:
    input:
        fastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        
    output:
        tsv="{sample}/s2_input_data_qc/Hifi/{sample}.hifi.tsv",
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=30,
        mem_mb=2000
    shell:
        '''
        seqkit stats {input.fastq} --threads {threads} -Ta -o {output.tsv}
        '''

rule get_input_hic_stats:
    input:
        hic1=lambda w: sampleInfo[w.sample]["all_hic1_files"]
        hic2=lambda w: sampleInfo[w.sample]["all_hic2_files"]
        
    output:
        inputDataStatsComplete="{sample}/s2_input_data_qc/.done.txt"
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=120,
        mem_mb=10000
    shell:
        '''
        seqkit stats {params.hic1} --threads {threads} -Ta -o {wildcards.sample}/s2_input_data_qc/HiC1/{wildcards.sample}.hic1.tsv
        seqkit stats {params.hic2} --threads {threads} -Ta -o {wildcards.sample}/s2_input_data_qc/HiC2/{wildcards.sample}.hic2.tsv
        '''
