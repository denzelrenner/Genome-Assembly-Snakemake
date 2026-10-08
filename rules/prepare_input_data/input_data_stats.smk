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
        hic1=lambda w: sampleInfo[w.sample]["hic1_input_file"],
        hic2=lambda w: sampleInfo[w.sample]["hic2_input_file"],
    output:
        outhic1="{sample}/s2_input_data_qc/HiC1/{sample}.hic1.tsv",
        outhic2="{sample}/s2_input_data_qc/HiC2/{sample}.hic2.tsv",
    params:
        inhic1=lambda w: sampleInfo[w.sample]["all_hic1_files"],
        inhic2=lambda w: sampleInfo[w.sample]["all_hic2_files"],
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=120,
        mem_mb=10000
    shell:
        '''
        seqkit stats {params.inhic1} --threads {threads} -Ta -o {output.outhic1}
        seqkit stats {params.inhic2} --threads {threads} -Ta -o {output.outhic2}
        '''
