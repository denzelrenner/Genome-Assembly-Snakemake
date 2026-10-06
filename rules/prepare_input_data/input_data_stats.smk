rule get_input_data_stats:
    input:
        hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        mergeComplete="{sample}/s1_merge_hic/done.txt",
        
    output:
        hifiTsv="{sample}/s2_input_data_qc/Hifi/{sample}.hifi.tsv",
        inputDataStatsComplete="{sample}/s2_input_data_qc/.done.txt"
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=120,
        mem_mb=10000
    params:
        runHiCStats=lambda w: sampleInfo[w.sample]["use_hic_shell"],
        hic1=lambda w: sampleInfo[w.sample]["all_hic1_files"],
        hic2=lambda w: sampleInfo[w.sample]["all_hic2_files"]
    shell:
        '''
        seqkit stats {input.hifiFastq} --threads {threads} -Ta -o {output.hifiTsv}

        if [ {params.runHiCStats} = true ];then

            seqkit stats {params.hic1} --threads {threads} -Ta -o {wildcards.sample}/s2_input_data_qc/HiC1/{wildcards.sample}.hic1.tsv
            seqkit stats {params.hic2} --threads {threads} -Ta -o {wildcards.sample}/s2_input_data_qc/HiC2/{wildcards.sample}.hic2.tsv

        fi
        
        touch {output.inputDataStatsComplete}'''
