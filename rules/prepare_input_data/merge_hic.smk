rule merge_hic:
    input:
        hic1=lambda wildcards:HICFILES[wildcards.sample]["HIC1"]
        hic2=lambda wildcards:HICFILES[wildcards.sample]["HIC2"]

    output:
        merged1=config['hic1_reads'] + '{sample}.merged_1.fq.gz'
        merged2=config['hic2_reads'] + '{sample}.merged_2.fq.gz'

    params:
        outputDir1=config['hic1_reads']
        outputDir2=config['hic2_reads']

    shell:
        "
        # merge HiC1
        cd {params.outputDir1}

        cat {input.hic1} > {output.merged1}

        # merge HiC2
        cd {params.outputDir2}

        cat {params.input.hic2} > {output.merged2}
        "
