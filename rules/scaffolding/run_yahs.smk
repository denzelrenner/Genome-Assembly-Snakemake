rule align_hic_to_hifiasm_haplome:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa'
    output:
        bam="{sample}/s5_run_yahs/hap{haplomeNumber}/mapped.PT.bam",
        fai='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa.fai',
        ctgSizes='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa.contigsizes'
    params:
        hic1=lambda w: sampleInfo[w.sample]['hic1_input_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_input_file'],
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=300000
    shell:
        '''
        samtools faidx {input.fa} -o {output.fai}

		# create tsv of name and size
		cut -f1,2 {output.fai} > {output.ctgSizes}
	
		# index with bwa
		bwa index {input.fa}

		# map hic reads to fasta from hifiasm
		bwa mem -5SP -T0 -t {threads} {input.fa} {params.hic1} {params.hic2} | \
		pairtools parse --min-mapq 40 --walks-policy 5unique \
		--max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} --chroms-path {output.ctgSizes} | \
		pairtools sort --tmpdir=tmp_pairtools | pairtools dedup --mark-dups --output-stats stats.txt | \
		pairtools split --output-pairs mapped.pairs --output-sam -|samtools view -bS -@ {threads} | \
		samtools sort -@ {threads} -o {output.bam}
        '''


rule run_yahs:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa',
        bam="{sample}/s5_run_yahs/hap{haplomeNumber}/mapped.PT.bam"
    output:
        yahsScaffoldOutput="{sample}/s5_run_yahs/hap{haplomeNumber}/yahs.out_scaffolds_final.fa",
        pretextInput="{sample}/s5_run_yahs/hap{haplomeNumber}/hap{haplomeNumber}_yahs.fa"
    conda:
        "../../envs/yahs.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=300000
    shell:
        '''
        yahs {input.fa} {input.bam}

        # specify hap1 or hap2 for fasta headers
        sed 's/^>/>H{wildcards.haplomeNumber}_/g' {output.yahsScaffoldOutput} > {output.pretextInput}
        '''
        
        
rule combine_yahs_fastas:

    input:
        fa=lambda w:expand("{{sample}}/s5_run_yahs/hap{hap}/hap{hap}_yahs.fa",hap=sampleInfo[w.sample]['ploidyAsArray']),
    output:
        combFasta="{sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa"
    threads:
        1
    resources:
        runtime=20,
        mem_mb=1000
    shell:
        '''
        cat {input.fa} > {output.combFasta}
        '''
