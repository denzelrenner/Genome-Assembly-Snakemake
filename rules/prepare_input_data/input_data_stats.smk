rule get_input_hifi_stats:
    input:
        fastq="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz",
    output:
        tsv="{sample}/s2_input_data_qc/Hifi/{sample}.hifi.tsv",
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=30,
        mem_mb=2000
    priority:
        100
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        seqkit stats {input.fastq} --threads {threads} -Ta -o {output.tsv}
        '''

rule get_input_hic_stats:
    input:
        hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file'],
    output:
        tsv1="{sample}/s2_input_data_qc/HiC1/{sample}.hic1.tsv",
        tsv2="{sample}/s2_input_data_qc/HiC2/{sample}.hic2.tsv",
    params:
        inhic1=lambda w,input: ' '.join([sampleInfo[w.sample]["input_hic1_files"],input.hic1]) if sampleInfo[w.sample]["merge_hic"] else input.hic1,
        inhic2=lambda w,input: ' '.join([sampleInfo[w.sample]["input_hic2_files"],input.hic2]) if sampleInfo[w.sample]["merge_hic"] else input.hic2,
    conda:
        "../../envs/seqkit.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=120,
        mem_mb=10000
    priority:
        100
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        seqkit stats {params.inhic1} --threads {threads} -Ta -o {output.tsv1}
        seqkit stats {params.inhic2} --threads {threads} -Ta -o {output.tsv2}
        '''

rule summarise_input_data_single_sample:
    input:
        hifi="{sample}/s2_input_data_qc/Hifi/{sample}.hifi.tsv",
        hic1=lambda w: expand("{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv",s=[w.sample]) if sampleInfo[w.sample]["useHiC"] else [],
        hic2=lambda w: expand("{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv",s=[w.sample]) if sampleInfo[w.sample]["useHiC"] else [],
    output:
        excel=f'{{sample}}/s2_input_data_qc/qc_results/{{sample}}_input_data_qc_{CURRENTDATE}.xlsx',
    params:
        qcScript=os.path.join(config['script_dir'],'python_scripts','summarise_sequencing_reads.py'),
        outfilePrefix=f'{{sample}}_input_data_qc_{CURRENTDATE}',
        genomesize=lambda w:sampleInfo[w.sample]['haplomeSize'],
    conda:
        "../../envs/pd-tools.yml",
    resources:
        runtime=30,
        mem_mb=1000,
    localrule:
        True
    threads:
        1
    priority:
        100
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        python3 {params.qcScript} -s {wildcards.sample} -g {params.genomesize} -o {params.outfilePrefix} -od {wildcards.sample}/s2_input_data_qc/qc_results
        '''

rule summarise_input_data_all_samples:
    input:
        hifi=expand("{smp}/s2_input_data_qc/Hifi/{smp}.hifi.tsv",smp=SAMPLES),
        hic1=expand("{smp}/s2_input_data_qc/HiC1/{smp}.hic1.tsv",smp=HICSAMPLES),
        hic2=expand("{smp}/s2_input_data_qc/HiC2/{smp}.hic2.tsv",smp=HICSAMPLES),
    output:
        excel=f'summary_results/input_data_qc/input_data_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx',
    params:
        qcScript=os.path.join(config['script_dir'],'python_scripts','summarise_sequencing_reads.py'),
        outfilePrefix=f'input_data_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}',
        inputSamples=SAMPLES,
        genomesizeSamples=GENOMESIZEFORSAMPLES,
    conda:
        "../../envs/pd-tools.yml",
    resources:
        runtime=30,
        mem_mb=1000,
    localrule:
        True
    threads:
        1
    priority:
        100
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH
        
        python3 {params.qcScript} -s {params.inputSamples} -g {params.genomesizeSamples} -o {params.outfilePrefix} -od summary_results/input_data_qc
        '''
