rule run_gphase:
    input:
        fa="{sample}/s4a_run_hifiasm/{sample}.p_utg.fa",
        hifi="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz",
        hic1="{sample}/s1_get_input_data/HiC1/{sample}.hic_1.fq.gz",
        hic2="{sample}/s1_get_input_data/HiC2/{sample}.hic_2.fq.gz",
    output:
        hic_to_utg_sam=temp("{sample}/s5_run_gphase/map.chromap.sam"),
        hic_to_utg_bam=temp("{sample}/s5_run_gphase/map.chromap.bam"),
        scaffoldFasta="{sample}/s5_run_gphase/gphase_output/scaffold_hap/gphase_final.unitig.scaffold.fasta",
    params:
        scriptPath=os.path.join(config['script_dir'],"GPhase"),
        ploidy=lambda w: sampleInfo[w.sample]['ploidy'],
        nchr=lambda w: sampleInfo[w.sample]['chrom'],
        gfa=lambda w: f"{w.sample}/s4a_run_hifiasm/{w.sample}.p_utg.gfa"
    conda:
        "../../envs/gphase.yml"
    threads:
        64
    resources:
        runtime=1200,
        mem_mb=80000,
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH
        export PATH={params.scriptPath}:$PATH

        # remove dir if it exists
        if [ -d {wildcards.sample}/s5_run_gphase ];then
            rm -r {wildcards.sample}/s5_run_gphase 
        fi

        # make out dir and move into it
        mkdir -p {wildcards.sample}/s5_run_gphase 
        cd {wildcards.sample}/s5_run_gphase 

        # map hic to assembly
        chromap -i -r ../../{input.fa} -o index

        chromap --preset hic -x index -r ../../{input.fa} -q 0 -1 ../../{input.hic1} -2 ../../{input.hic2} --remove-pcr-duplicates -t {threads} --SAM -o ../../{output.hic_to_utg_sam}

        samtools view -@ {threads} -bh ../../{output.hic_to_utg_sam} -o ../../{output.hic_to_utg_bam}

        gphase popcnv -f ../../{input.fa} -p output -t {threads} -r ../../{input.hifi}

        gphase pipeline -f ../../{input.fa} -g ../../{params.gfa} -c popcnv/06.genes.round.cn -m ../../{output.hic_to_utg_bam} --n_chr {params.nchr} --n_hap {params.ploidy} --hap_pm 0.6 -p output --rescue --min_len 50

        #chromap -i -r ../../{output.scaffoldFasta} -o reindex_utg

        # Remap Hi-C/Omni-C reads
        #chromap --preset hic -x reindex_utg -r ../../{output.scaffoldFasta} -q 0 -1 {input.hic1} -2 {input.hic2} --remove-pcr-duplicates -t {threads} --SAM -o unitig_mq.sam

        #samtools view -@64 -h unitig_mq.sam | samtools sort -@64 -o unitig_mq.bam

        #samtools view -@64 -h unitig_mq.bam | PretextMap -o unitig.mq.pretext --mapq 0

        '''

rule rename_gphase_fasta:
    input:
        "{sample}/s5_run_gphase/gphase_output/scaffold_hap/gphase_final.unitig.scaffold.fasta"
    output:
        "{sample}/s5_run_gphase/gphase_output/scaffold_hap/scaffolds.fa"
    localrule:
        True
    threads:
        1
    resources:
        runtime=20,
        mem_mb=1000,
    shell:
        '''
        cp {input} {output}
        '''
