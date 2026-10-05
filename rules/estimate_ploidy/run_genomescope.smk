rule run_genomescope:
    input:
        fastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz"

    output:
        plots=expand("{{sample}}/s3_estimate_ploidy/genomescope/plots/{png}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot'])

    params:
        ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
        outdir="{sample}/s3_estimate_ploidy/genomescope"

    conda:
        "../../envs/genomescope.yml"

    resources:
        mem_mb=200000
    threads:
        int(workflow.cores * 0.5)
    resources:
        runtime=240,
        mem_mb=200000
    shell:
        '''
        mkdir -p {params.outdir}/tmp 
        mkdir -p {params.outdir}/plots
        
        # run kmc
        kmc -k31 -t{threads} -m178 -ci1 -cs2000000 {input.fastq} {params.outdir}/reads {params.outdir}/tmp/
        kmc_tools transform {params.outdir}/reads histogram {params.outdir}/reads.histo -cx2000000

        # run genomescope
        genomescope2 -i {params.outdir}/reads.histo -o {params.outdir} -k 31 -m 2000000 -p {params.ploidy}
        mv {params.outdir}/*.png {params.outdir}/plots
        '''

