import os

rule get_hic:
    output:
        hic1="{sample}/s1_get_input_data/HiC1/{sample}.hic_1.fq.gz",
        hic2="{sample}/s1_get_input_data/HiC2/{sample}.hic_2.fq.gz",
    params:
        in1=lambda w: sampleInfo[w.sample]["HiC1"],
        in2=lambda w: sampleInfo[w.sample]["HiC2"],
        merge=lambda w: 'true' if sampleInfo[w.sample]["merge_hic"] else 'false',
    threads:
        1
    localrule:
        True
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
            ln -s {params.in1} {output.hic1}
            ln -s {params.in2} {output.hic2}

        fi
        '''
        
        
