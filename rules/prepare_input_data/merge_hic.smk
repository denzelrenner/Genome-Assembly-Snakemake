import os
rule merge_hic:
    output:
        merged1=os.path.join(config['hic1_reads'],'{sample}','{sample}.merged_1.fq.gz'),
        merged2=os.path.join(config['hic2_reads'],'{sample}','{sample}.merged_2.fq.gz'),
    params:
        hic1=lambda w: sampleInfo[w.sample]["hic1_files"],
        hic2=lambda w: sampleInfo[w.sample]["hic2_files"],
    threads:
        1
    resources:
        runtime=120,
        mem_mb=1000
    shell:
        '''
        # merge HiC1 and HiC2
        cat {params.hic1} > {output.merged1}
        cat {params.hic2} > {output.merged2}
        '''
        
        
