rule run_hifiasm:
    input:
        hifi="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz",
        hic1=lambda w: f"{w.sample}/s1_get_input_data/HiC1/{w.sample}.hic_1.fq.gz" if sampleInfo[w.sample]['has_hic'] else [],
        hic2=lambda w: f"{w.sample}/s1_get_input_data/HiC2/{w.sample}.hic_2.fq.gz" if sampleInfo[w.sample]['has_hic'] else []
    output:
        complete='{sample}/s4a_run_hifiasm/.done.txt',
    log:
        stderr='{sample}/s4a_run_hifiasm/hifiasm_run.log',
    params:
        hasHiCShell=lambda w: 'true' if sampleInfo[w.sample]["has_hic"] else 'false',
        ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
        hifiasmPloidy=lambda w: sampleInfo[w.sample]["hifiasm_ploidy"],
        outdir=lambda w:f'{w.sample}/s4a_run_hifiasm',
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

        # if output dir exists then remove it
        if [ -d {params.outdir} ];then
            rm -r {params.outdir}
        fi

        # make outdir
        mkdir -p {params.outdir}

        if [ {params.hasHiCShell} = true ];then

	        hifiasm -o {params.outdir}/{wildcards.sample} --n-hap {params.ploidy} --h1 {input.hic1} --h2 {input.hic2} -t {threads} {input.hifi} > {log.stderr} 2>&1
        else
            
            hifiasm -o {params.outdir}/{wildcards.sample} --n-hap {params.ploidy} -t {threads} {input.hifi} > {log.stderr} 2>&1
        fi
    
        # clean up unnecessary files
        rm {params.outdir}/*.bin
        rm {params.outdir}/*.noseq.gfa
        rm {params.outdir}/*.lowQ.bed

        # go through all the gfa files
        for gfa in {params.outdir}/*hap*.p_ctg.gfa; do

            echo $gfa

            if [ -f $gfa ]; then

                # get only file name and not whole path
                file=$(basename "${{gfa}}")

                # find haplome number
                hap=$(awk -F'.' '{{print $(NF-2)}}' <<< $file)

                # change name of gfa file
                mv $gfa {params.outdir}/{wildcards.sample}."${{hap}}".p_ctg.gfa
                
            fi
        done

        # change name of utg file
        mv {params.outdir}/*.p_utg.gfa {params.outdir}/{wildcards.sample}.p_utg.gfa

        # do a sense check for completion by making sure we had the expected number of gfa's
        gfaFiles=({params.outdir}/*hap*.p_ctg.gfa)
        numOfGfaFiles=${{#gfaFiles[@]}}

        # echo complete signal if file number is the same
        if [ "$numOfGfaFiles" -eq "{params.hifiasmPloidy}" ] && [ -f {params.outdir}/{wildcards.sample}.p_utg.gfa ];then

            echo 'DONE' > {output.complete}
        fi

        '''
        
        

