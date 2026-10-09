rule align_hic_to_hifiasm_haplome:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa',
        hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file']
    output:
        bam="{sample}/s5_run_yahs/hap{haplomeNumber}/mapped.PT.bam",
        fai='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa.fai',
        ctgSizes='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa.contigsizes',
        stats="{sample}/s5_run_yahs/hap{haplomeNumber}/stats.txt",
        pairs="{sample}/s5_run_yahs/hap{haplomeNumber}/mapped.pairs",
    params:
        tmpdir="{sample}/s5_run_yahs/hap{haplomeNumber}/tmp_pairtools",
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        int(workflow.cores * 0.9)
    resources:
        runtime=1200,
        mem_mb=300000
    shell:
        '''
        samtools faidx {input.fa} -o {output.fai}

		# create tsv of name and size
		cut -f1,2 {output.fai} > {output.ctgSizes}

                # make pairtools tmp dir
                mkdir -p {params.tmpdir}
	
		# index with bwa
		bwa index {input.fa}

		# map hic reads to fasta from hifiasm
		bwa mem -5SP -T0 -t {threads} {input.fa} {input.hic1} {input.hic2} | \
		pairtools parse --min-mapq 40 --walks-policy 5unique \
		--max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} --chroms-path {output.ctgSizes} | \
		pairtools sort --tmpdir={params.tmpdir} | pairtools dedup --mark-dups --output-stats {output.stats} | \
		pairtools split --output-pairs {output.pairs} --output-sam -|samtools view -bS -@ {threads} | \
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
        int(workflow.cores * 0.1)
    resources:
        runtime=240,
        mem_mb=5000
    shell:
        '''
        yahs -o {wildcards.sample}/s5_run_yahs/hap{wildcards.haplomeNumber}/yahs.out {input.fa} {input.bam}

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
