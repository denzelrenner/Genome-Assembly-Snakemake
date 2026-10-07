rule build_pretextmap:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        
    output:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        hicToAssem="{sample}/s6_build_pretextmap/pretextmap/hic_sorted.bam",
    params:
        hic1=lambda w: sampleInfo[w.sample]['hic1_input_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_input_file'],
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=300000
    shell:
        '''
        echo 'bwa pipe'

        bwa index {input.fa} && \
        bwa mem -t {threads} {input.fa} {params.hic1} {params.hic2} | \
        samtools sort -@{threads} -o {output.hicToAssem} && \
        samtools view -@{threads} -h {output.hicToAssem} | PretextMap -o {input.pretextmap} --mapq 0

        # create index of fasta
        samtools faidx {input.fa} -o {input.fa}.fai
        '''
    
rule create_telomeres_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext"
    output:
        telomereBedgraph="{sample}/s6_build_pretextmap/telomeres/telomeres_telomeric_repeat_windows.bedgraph",
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.25)
    resources:
        runtime=240,
        mem_mb=10000
    shell:
        '''
        echo 'tidk telomeres'
        tidk search --fasta {input.fa} --string TTTAGGG --output telomeres --dir {sample}/s6_build_pretextmap/telomeres --extension bedgraph

        cat {output.telomereBedgraph} | PretextGraph -i {input.pretextmap} -n "telomeres"  
        ''' 

rule create_hic_coverage_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        hicToAssem="{sample}/s6_build_pretextmap/pretextmap/hic_sorted.bam",
        
    output:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",        
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hic.bedgraph",
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=120,
        mem_mb=5000
    shell:
        '''
        echo 'coverage graph hic'

        bedtools genomecov -ibam {input.hicToAssem} -bga > {output.hicCovBedgraph}
        cat {output.hicCovBedgraph} | PretextGraph -i {input.pretextmap} -n "coveragehic"
        ''' 

rule create_hifi_coverage_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
    output:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",        
        hifiCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hifi.bedgraph",
        hifiToAssem="{sample}/s6_build_pretextmap/pretextmap/hifi_sorted.bam",
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=30000
    shell:
        '''
        echo 'minimap mapping hifi reads to scaffold assem'

        minimap2 -ax map-pb -t {threads} {input.fa} {input.hifiFastq} | samtools sort -@{threads} -O BAM -o {output.hifiToAssem} -
        
        echo 'coverage graph hifi'
        
        bedtools genomecov -ibam {output.hifiToAssem} -bga > {output.hifiCovBedgraph}
        cat {output.hifiCovBedgraph} | PretextGraph -i {input.pretextmap} -n "coveragehifi"
        ''' 

rule create_hic_gaps_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hic.bedgraph",
    output:
        covGapsBedgraph="{sample}/s6_build_pretextmap/coverage_gaps/gaps.bedgraph",

    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=30,
        mem_mb=1000
    shell:
        '''
        echo 'gaps graph'

        grep -w 0$ {input.hicCovBedgraph} | sed 's/0$/200/g' > {output.covGapsBedgraph}
        cat {output.covGapsBedgraph} | PretextGraph -i {input.pretextmap} -n "gaps"
        ''' 

rule create_gfastats_gaps_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
    output:
        gfastatsGapsBedgraph="{sample}/s6_build_pretextmap/gfastats_gaps/gfastats_gaps.bedgraph",
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.10)
    resources:
        runtime=30,
        mem_mb=5000
    shell:
        '''
        echo 'gfastats gaps'	

        gfastats ${input.fa} -b gaps | sed "s/$/\t200/" > {output.gfastatsGapsBedgraph}

        cat {output.gfastatsGapsBedgraph} | PretextGraph -i {input.pretextmap} -n "gfastats_gaps"
        ''' 

rule create_barrnap_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
    output:
        barrnapGff="{sample}/s6_build_pretextmap/barrnap/RNA.gff",
        barrnapHits="{sample}/s6_build_pretextmap/barrnap/all_hits.fa",
        barrnapBedgraph="{sample}/s6_build_pretextmap/barrnap/RNA.bedgraph",

    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.10)
    resources:
        runtime=60,
        mem_mb=10000,
    shell:
        '''
        echo 'rDNA'

        barrnap --threads {threads} \
            --kingdom euk \
            --outseq {output.barrnapHits} \
            {input.fa} > {output.barrnapGff}

        grep -v '#' {output.barrnapGff} | cut -f 1,4,5 | sed "s/$/\t200/" > {output.barrnapBedgraph}

        cat {output.barrnapBedgraph} | PretextGraph -i {input.pretextmap} -n "rDNA"
        ''' 

# rule create_trash_bedgraph:
#     input:
#         fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
#        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
#         bam="{sample}/s5_run_yahs/hap{haplomeNumber}/mapped.PT.bam"
#         hifiFastq="{sample}/s1_bam_to_fastq/{sample}.fastq.gz"
#         pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
#         hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hic.bedgraph",
        
#     output:
#         covGapsBedgraph="{sample}/s6_build_pretextmap/coverage_gaps/gaps.bedgraph",
#         gfastatsGapsBedgraph="{sample}/s6_build_pretextmap/gfastats_gaps/gfastats_gaps.bedgraph",
#         barrnapGff="{sample}/s6_build_pretextmap/barrnap/RNA.gff",
#         barrnapHits="{sample}/s6_build_pretextmap/barrnap/all_hits.fa",
#         barrnapBedgraph="{sample}/s6_build_pretextmap/barrnap/RNA.bedgraph"


#     params:
#         ploidy=lambda w: sampleInfo[w.sample]["ploidy"],
#         hic1=lambda w: sampleInfo[w.sample]['hic1_input_file'],
#         hic2=lambda w: sampleInfo[w.sample]['hic2_input_file'],
#     conda:
#         "../../envs/manual-curation.yml"
#     threads:
#         int(workflow.cores * 0.75)
#     resources:
#         runtime=720,
#         mem_mb=300000
#     shell:
#         '''

#         echo 'run trash'

#         TRASH_run.sh {input.fa} --o $yahsdir/s2_pretext/tracks/trash --par 48

#         fasta_base=$(basename "${{input.fa}}")

#         trash_analysis_for_pretext.py -i "Summary.of.repetitive.regions.${fasta_base}.csv"

#         cat trashrepeat_all.bedgraph | PretextGraph -i $premap -n "trash"
#         ''' 