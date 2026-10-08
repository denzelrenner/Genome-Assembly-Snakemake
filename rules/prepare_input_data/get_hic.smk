import os
rule merge_hic:
    output:
        merged1=os.path.join(config['hic1_reads'],'{sample}','{sample}.merged_1.fq.gz'),
        merged2=os.path.join(config['hic2_reads'],'{sample}','{sample}.merged_2.fq.gz'),
    params:
        hic1=lambda w: sampleInfo[w.sample]["input_hic1_files"],
        hic2=lambda w: sampleInfo[w.sample]["input_hic2_files"],
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

rule get_hic:
    output:
        hic1="{sample}/s1_get_input_data/HiC1/{sample}.hic_1.fq.gz".
        hic2="{sample}/s1_get_input_data/HiC2/{sample}.hic_2.fq.gz",
    params:
        in1=lambda w: sampleInfo[w.sample]["input_hic1_files"],
        in2=lambda w: sampleInfo[w.sample]["input_hic2_files"],
        merge=lambda w: sampleInfo[w.sample]["merge_hic"],
    threads:
        1
    resources:
        runtime=120,
        mem_mb=1000
    shell:
        '''

        if [ {params.merge} = true ];then

            # merge HiC1 and HiC2
            cat {params.in1} > {output.hic1}
            cat {params.in2} > {output.hic2}

        elif [ {params.merge} = false ];then

            # copy file to output
            cp {params.in1} {output.hic1}
            cp {params.in2} {output.hic2}

        fi
        '''
        
        
