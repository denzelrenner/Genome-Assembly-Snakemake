import os
rule merge_hic:
    input:
        hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz"
    output:
        mergeComplete="{sample}/s1_merge_hic/done.txt"

    params:
        mergeData=lambda w: sampleInfo[w.sample]["merge_hic"],
        merged1=os.path.join(config['hic1_reads'],'{sample}','{sample}.merged_1.fq.gz'),
        merged2=os.path.join(config['hic2_reads'],'{sample}','{sample}.merged_2.fq.gz'),
        hic1=lambda w: sampleInfo[w.sample]["merge_hic1_files"],
        hic2=lambda w: sampleInfo[w.sample]["merge_hic2_files"]
        
    threads:
        1

    resources:
        runtime=120,
        mem_mb=1000

    shell:
        '''
        mergeHiC={params.mergeData}

        if [ $mergeHiC -eq 1 ];then

            # merge HiC1 and HiC2
            cat {params.hic1} > {params.merged1}
            cat {params.hic2} > {params.merged2}
            echo 'HiC was merged' > {output.mergeComplete}

        elif [ $mergeHiC -eq 0 ];then
            echo 'No HiC was merged' > {output.mergeComplete}

        fi
        '''
        
        
