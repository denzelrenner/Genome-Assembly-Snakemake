rule run_hifiasm:
    input:
        hifi="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz",
        hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file']
    output:
        outdir=directory('{sample}/s4a_run_hifiasm'),
        complete='{sample}/s4a_run_hifiasm/.done.txt',
    log:
        stderr='{sample}/s4a_run_hifiasm/hifiasm_run.log',
    params:
        useHiCShell=lambda w: 'true' if sampleInfo[w.sample]["useHiC"] else 'false',
        ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
        hifiasmPloidy=lambda w: sampleInfo[w.sample]["hifiasm_ploidy"],
    conda:
        "../../envs/hifiasm.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=300000
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        hifiasm -v # confirm hifiasm version

        # make outdir
        mkdir -p {output.outdir}

        if [ {params.useHiCShell} = true ];then

	        hifiasm -o {output.outdir}/{wildcards.sample} --n-hap {params.ploidy} --h1 {input.hic1} --h2 {input.hic2} -t {threads} {input.hifi} > {log.stderr} 2>&1
        else
            
            hifiasm -o {output.outdir}/{wildcards.sample} --n-hap {params.ploidy} -t {threads} {input.hifi} > {log.stderr} 2>&1
        fi
    
        # clean up unnecessary files
        rm {output.outdir}/*.bin
        rm {output.outdir}/*.noseq.gfa
        rm {output.outdir}/*.lowQ.bed

        # go through all the gfa files
        for gfa in {output.outdir}/*hap*.p_ctg.gfa; do

            echo $gfa

            if [ -f $gfa ]; then

                # get only file name and not whole path
                file=$(basename "${{gfa}}")

                # find haplome number
                hap=$(awk -F'.' '{{print $(NF-2)}}' <<< $file)

                # change name of gfa file
                mv $gfa {output.outdir}/{wildcards.sample}."${{hap}}".p_ctg.gfa
                
            fi
        done

        # change name of utg file
        mv {output.outdir}/*.p_utg.gfa {output.outdir}/{wildcards.sample}.p_utg.gfa

        # do a sense check for completion by making sure we had the expected number of gfa's
        gfaFiles=({output.outdir}/*hap*.p_ctg.gfa)
        numOfGfaFiles=${{#gfaFiles[@]}}

        # echo complete signal if file number is the same
        if [ "$numOfGfaFiles" -eq "{params.hifiasmPloidy}" ] && [ -f {output.outdir}/{wildcards.sample}.p_utg.gfa ];then

            echo 'DONE' > {output.complete}
        fi

        '''
        
        

