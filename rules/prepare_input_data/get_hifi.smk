import os
rule get_hifi:
    output:
        fastq="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz"
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        8
    resources:
        runtime=240,
        mem_mb=5000,
    params:
        hifi=lambda w: sampleInfo[w.sample]['HiFi'],
        operation=lambda w: sampleInfo[w.sample]['hifi_operation'],
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        if [ {params.operation} = 'convert' ];then

            samtools fastq --threads {threads} {params.hifi} | pigz -p {threads} -6 > {output.fastq}

        elif [ {params.operation} = 'compress' ];then

            piqz -p {threads} -6 -c {params.hifi} > {output.fastq}

        elif [ {params.operation} = 'copy' ];then

            ln -s {params.hifi} {output.fastq}

        fi
        '''

        
