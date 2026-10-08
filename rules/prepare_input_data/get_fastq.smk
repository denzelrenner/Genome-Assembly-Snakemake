import os
rule bam_to_fastq:
    input:
        bam=os.path.join(config['hifi_reads'],"{sample}.hifi_reads.bam")
    output:
        fastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz"
    conda:
        "../../envs/minimap2_samtools.yml"
    threads:
        2
    resources:
        runtime=240,
        mem_mb=5000
    shell:
        "samtools fastq --threads {threads} {input.bam} | gzip > {output.fastq}"


rule get_fastq:
    output:
        fastq="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz"
    conda:
        "../../envs/minimap2_samtools.yml"
    threads:
        2
    resources:
        runtime=240,
        mem_mb=5000
    params:
        hifi=lambda w: sampleInfo[w.sample]['hifi']
        operation=lambda w: sampleInfo[w.sample]['hifiOperation']
    shell:
        '''
        if [ {params.operation} = 'convert' ];then

            samtools fastq --threads {threads} {params.hifi} | gzip > {output.fastq}

        elif [ {params.operation} = 'compress' ];then

            gzip {params.hifi} > {output.fastq}

        elif [ {params.operation} = 'copy' ];then

            cp {params.hifi} > {output.fasq}

        fi
        '''


rule bam_to_fastq_with_filter:
    input:
        bam=config['hifi_reads'] + "{sample}.hifi_reads.bam"

    output:
        fastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz"

    params:
        outdir="{sample}/s1_get_fastq",
        DB=config['hifiadapterfiltDB_path']
    
    conda:
        "../../envs/hifiadapterfilt.yml"

    shell:
        '''export PATH=$PATH:{params.DB}

        # sym link to input bam
        ln -s {input.bam} ./

        hifiadapterfilt.sh \
	        -p {sample}.hifi_reads \
        	-l 44 \
        	-m 97 \
        	-t 16 \
        	-o .

        # set variable for filtered reads
        mv {sample}.filt.fastq.gz {sample}.fastq.gz'''
        
