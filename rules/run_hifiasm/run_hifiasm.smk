rule run_hifiasm:
    input:
        hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        smudgeplots=expand("{{sample}}/s3_estimate_ploidy/smudgeplot/plots/{png}.png",png=['output_smudgeplot','output_smudgeplot_log10','output_centralities']),
        genomescopeplots=expand("{{sample}}/s3_estimate_ploidy/genomescope/plots/{png}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        inputdatastats="{sample}/s2_input_data_qc/.done.txt"
    output:
        complete="{sample}/s4_run_hifiasm/.done.txt"
    params:
        useHiCShell=lambda w: sampleInfo[w.sample]["use_hic_shell"],
        ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
        hic1=lambda w: sampleInfo[w.sample]['hic1_input_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_input_file'],
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

        if {params.useHiCShell};then

	        hifiasm -o {wildcards.sample}/s4_run_hifiasm/{wildcards.sample} --n-hap {params.ploidy} --h1 {params.hic1} --h2 {params.hic2} -t {threads} {input.hifiFastq}
        
        else
            
            hifiasm -o {wildcards.sample}/s4_run_hifiasm/{wildcards.sample} --n-hap {params.ploidy} -t {threads} {input.hifiFastq}

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

        echo 'DONE' > {output.complete}
        '''
        
        

