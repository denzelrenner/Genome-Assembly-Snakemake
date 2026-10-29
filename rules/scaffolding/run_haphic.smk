import os
rule align_hic_to_utg:
    input:
        fa='{sample}/s4a_run_hifiasm/{sample}.p_utg.fa',
        hic1="{sample}/s1_get_input_data/HiC1/{sample}.hic_1.fq.gz",
        hic2="{sample}/s1_get_input_data/HiC2/{sample}.hic_2.fq.gz"
    output:
        bam="{sample}/s5_run_haphic/align_hic_to_utg/HiC.bam",
        filteredBam="{sample}/s5_run_haphic/align_hic_to_utg/HiC.filtered.bam",
    params:
        scriptPath=os.path.join(config['script_dir'],"HapHiC","utils"),
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        int(workflow.cores * 0.9)
    resources:
        runtime=1200,
        mem_mb=300000,
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH
        
        # align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
        bwa index {input.fa}
        bwa mem -5SP -t {threads} {input.fa} {input.hic1} {input.hic2} | samblaster | samtools view - -@ {threads} -S -h -b -F 3340 -o {output.bam}

        # filter the alignments with MAPQ 1 (mapping quality M1 and NM 3 edit distance < 3)
        {params.scriptPath}/filter_bam {output.bam} 1 --nm 3 --threads {threads} | samtools view - -b -@ {threads} -o {output.filteredBam}
        '''
        
rule run_haphic:
    input:
        fa='{sample}/s4a_run_hifiasm/{sample}.p_utg.fa',
        bam="{sample}/s5_run_haphic/align_hic_to_utg/HiC.filtered.bam"
    output:
        agp="{sample}/s5_run_haphic/scaffolding/04.build/scaffolds.agp",
        outFasta="{sample}/s5_run_haphic/scaffolding/04.build/scaffolds.fa",
        pdf="{sample}/s5_run_haphic/scaffolding/plot/{sample}_contact_map.pdf"
    params:
        outdir=lambda w: f"{w.sample}/s5_run_haphic/scaffolding",
        nchrs=lambda w: sampleInfo[w.sample]["ploidy"] * sampleInfo[w.sample]["chrom"],
        scriptPath=os.path.join(config['script_dir'],"HapHiC"),
    conda:
        "../../envs/haphic.yml"
    threads:
        int(workflow.cores * 0.10)
    resources:
        runtime=720,
        mem_mb=10000,
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # error handling/check, bug where python from the haphic env is not being used in the shell
        echo $CONDA_PREFIX
        which python
        which python3
        echo $PATH

        # if output dir exists then remove it
        if [ -d {params.outdir} ];then
            rm -r {params.outdir}
        fi

        # make output dir 
        mkdir -p {params.outdir}

        {params.scriptPath}/haphic pipeline {input.fa} {input.bam} {params.nchrs} --threads {threads} --processes {threads} --outdir {params.outdir}

        # make plots dir
        mkdir -p {params.outdir}/plot
        cd {params.outdir}/plot

        # move into plot dir otherwise supplementary files will not be properly placed
        {params.scriptPath}/haphic plot ../../../../{output.agp} ../../../../{input.bam} --prefix {wildcards.sample}_
        '''
        

        

