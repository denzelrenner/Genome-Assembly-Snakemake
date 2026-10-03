## HIC FILE SORTING ###
def get_HiC1_input_files(wildcards):

    if HASHIC and MERGEHIC[wildcards.sample]:

        return [f for f in os.listdir(os.path.join(config['hic1_reads'],s)) if 'merged' in f][0]

    elif HASHIC and not MERGEHIC[wildcards.sample]:

        return HICFILES[wildcards.sample]["HIC1"][0]

    elif not HASHIC and not MERGEHIC[wildcards.sample]:

        return []

def get_HiC2_input_files(wildcards):

    if HASHIC and MERGEHIC[wildcards.sample]:

        return [f for f in os.listdir(os.path.join(config['hic2_reads'],s)) if 'merged' in f][0]

    elif HASHIC and not MERGEHIC[wildcards.sample]:

        return HICFILES[wildcards.sample]["HIC2"][0]

    elif not HASHIC and not MERGEHIC[wildcards.sample]:

        return []



rule run_hifiasm:
    input:
        hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        hic1=get_HiC1_input_files,
        hic2=get_HiC2_input_files,
        smudgeplots=expand("{{sample}}/s3_estimate_ploidy/smudgeplot/plots/{png}.png",png=['output_smudgeplot','output_smudgeplot_log10','output_centralities']),
        genomescopeplots=expand("{{sample}}/s3_estimate_ploidy/genomescope/plots/{png}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        inputdatastats="{sample}/s2_input_data_qc/Hifi/{sample}.hifi.tsv"
    output:
        expand('{{sample}}/s4_run_hifiasm/{{sample}}.hap{haplomeNumber}.p_ctg.gfa',haplomeNumber=HAPLOMENUMBERS)
    params:
        hicParam="--h1 {input.hic1} --h2 {input.hic2} " if HASHIC else [],
        ploidy=config["ploidy"]
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
	    hifiasm -o {wildcards.sample}/s4_run_hifiasm/{wildcards.sample} --n-hap {params.ploidy} {params.hicParam}-t {threads} {input.hifiFastq}

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
        '''
        
        

