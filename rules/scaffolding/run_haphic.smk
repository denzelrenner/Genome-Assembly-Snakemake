import os
rule align_hic_to_utg:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.p_utg.fa',
        hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file']
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
        # align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
        bwa index {input.fa}
        bwa mem -5SP -t {threads} {input.fa} {input.hic1} {input.hic2} | samblaster | samtools view - -@ {threads} -S -h -b -F 3340 -o {output.bam}

        # filter the alignments with MAPQ 1 (mapping quality M1 and NM 3 edit distance < 3)
        {params.scriptPath}/filter_bam {output.bam} 1 --nm 3 --threads {threads} | samtools view - -b -@ {threads} -o {output.filteredBam}
        '''
        
rule run_haphic:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.p_utg.fa',
        bam="{sample}/s5_run_haphic/align_hic_to_utg/HiC.filtered.bam"
    output:
        agp="{sample}/s5_run_haphic/04.build/scaffolds.agp",
        fa="{sample}/s5_run_haphic/04.build/scaffolds.fa",
        pdf="{sample}/s5_run_haphic/plot/{sample}_contact_map.pdf"
    params:
        nchrs=lambda w: sampleInfo[w.sample]["ploidy"] * sampleInfo[w.sample]["chrom"],
        scriptPath=os.path.join(config['script_dir'],"HapHiC"),
    conda:
        "../../envs/haphic.yml"
    threads:
        int(workflow.cores * 0.50)
    resources:
        runtime=720,
        mem_mb=100000,
    shell:
        '''
        {params.scriptPath}/haphic pipeline {input.fa} {input.bam} {params.nchrs} --threads {threads} --processes {threads} --outdir {sample}/s5_run_haphic

        {params.scriptPath}/haphic plot {output.agp} {input.bam} --prefix {sample}/s5_run_haphic/plot/{sample}_contact_map
        '''
        
        

