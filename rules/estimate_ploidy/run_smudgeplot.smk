rule run_smudgeplot:
    input:
        fastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz"
    output:
        plots=expand("{{sample}}/s3_estimate_ploidy/smudgeplot/plots/{png}.png",png=['output_smudgeplot','output_smudgeplot_log10','output_centralities'])
    params:
        ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
        outdir="{sample}/s3_estimate_ploidy/smudgeplot"
    conda:
        "../../envs/smudgeplot053.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=480,
        mem_mb=200000
    shell:
        '''
        mkdir -p {params.outdir}/plots

        echo 'run fastk'

        FastK -v -t4 -k31 -M36 -T24 {input.fastq} -N{params.outdir}/FastK_Table
            
        echo 'now for smudgeplot'

        # Find all k-mer pairs in the dataset using hetmer module
        smudgeplot hetmers -L 12 -t {threads} -o {params.outdir}/kmerpairs --verbose {params.outdir}/FastK_Table

        # use the .smu file to infer ploidy and create smudgeplot
        smudgeplot all -o {params.outdir}/output {params.outdir}/kmerpairs.smu
        
        mv {params.outdir}/*.png {params.outdir}/plots'''
        
        

