rule run_hifiasm:
    input:
        hifi="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        hic1=lambda w: sampleInfo[w.sample]['hic1_input_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_input_file']
    output:
        complete="{sample}/s4_run_hifiasm/.done.txt"
        utg='{sample}/s4_run_hifiasm/{sample}.p_utg.gfa',
    params:
        useHiCShell=lambda w: sampleInfo[w.sample]["use_hic_shell"],
        ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
        
    conda:
        "../../envs/hifiasm.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=300000
    shell:
        '''
        hifiasm -v # confirm hifiasm version

        if [ {params.useHiCShell} = true ];then

	        hifiasm -o {wildcards.sample}/s4_run_hifiasm/{wildcards.sample} --n-hap {params.ploidy} --h1 {input.hic1} --h2 {input.hic2} -t {threads} {input.hifi}
        
        else
            
            hifiasm -o {wildcards.sample}/s4_run_hifiasm/{wildcards.sample} --n-hap {params.ploidy} -t {threads} {input.hifi}

        fi
    

        # go through all the gfa files
        for gfa in {wildcards.sample}/s4_run_hifiasm/*hap*.p_ctg.gfa; do

            echo $gfa

            if [ -f $gfa ]; then

                # get only file name and not whole path
                file=$(basename "${{gfa}}")

                # find haplome number
                hap=$(awk -F'.' '{{print $(NF-2)}}' <<< $file)

                # change name of gfa file
                mv $gfa {wildcards.sample}/s4_run_hifiasm/{wildcards.sample}."${{hap}}".p_ctg.gfa
                
            fi
        done

        # change name of utg file
        mv {wildcards.sample}/s4_run_hifiasm/*.p_utg.gfa {wildcards.sample}/s4_run_hifiasm/{wildcards.sample}.p_utg.gfa

        echo 'DONE' > {output.complete}
        '''
        
        

