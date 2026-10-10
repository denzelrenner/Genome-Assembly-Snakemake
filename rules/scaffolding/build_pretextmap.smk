rule build_pretextmap:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",    
        hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file'],    
    output:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        bam="{sample}/s6_build_pretextmap/pretextmap/hic_sorted.bam",
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=1200,
        mem_mb=300000
    shell:
        '''
        echo 'bwa pipe'

        bwa index {input.fa} && \
        bwa mem -t {threads} {input.fa} {input.hic1} {input.hic2} | \
        samtools sort -@{threads} -o {output.bam} && \
        samtools view -@{threads} -h {output.bam} | PretextMap -o {output.pretextmap} --mapq 0

        # create index of fasta
        samtools faidx {input.fa} -o {input.fa}.fai
        '''
    
rule create_telomeres_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
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

        ''' 

rule create_hic_coverage_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        hicToAssem="{sample}/s6_build_pretextmap/pretextmap/hic_sorted.bam",
        
    output:
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hic.bedgraph",
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=180,
        mem_mb=10000
    shell:
        '''
        echo 'coverage graph hic'

        bedtools genomecov -ibam {input.hicToAssem} -bga > {output.hicCovBedgraph}
        ''' 

rule create_hifi_coverage_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        hifiFastq="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz",
    output:
        hifiCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hifi.bedgraph",
        bam="{sample}/s6_build_pretextmap/pretextmap/hifi_sorted.bam",
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

        minimap2 -ax map-pb -t {threads} {input.fa} {input.hifiFastq} | samtools sort -@{threads} -O BAM -o {output.bam} -
        
        echo 'coverage graph hifi'
        
        bedtools genomecov -ibam {output.bam} -bga > {output.hifiCovBedgraph}
        ''' 

rule create_hic_gaps_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hic.bedgraph",
    output:
        covGapsBedgraph="{sample}/s6_build_pretextmap/coverage_gaps/gaps.bedgraph",

    conda:
        "../../envs/manual-curation.yml"
    threads:
        2
    resources:
        runtime=30,
        mem_mb=1000
    shell:
        '''
        echo 'gaps graph'

        grep -w 0$ {input.hicCovBedgraph} | sed 's/0$/200/g' > {output.covGapsBedgraph}
        ''' 

rule create_gfastats_gaps_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
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

        gfastats {input.fa} -b gaps | sed "s/$/\t200/" > {output.gfastatsGapsBedgraph}

        ''' 

rule create_barrnap_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
    output:
        barrnapGff="{sample}/s6_build_pretextmap/barrnap/RNA.gff",
        barrnapHits="{sample}/s6_build_pretextmap/barrnap/all_hits.fa",
        barrnapBedgraph="{sample}/s6_build_pretextmap/barrnap/RNA.bedgraph",

    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.10)
    resources:
        runtime=90,
        mem_mb=10000,
    shell:
        '''
        echo 'rDNA'

        barrnap --threads {threads} \
            --kingdom euk \
            --outseq {output.barrnapHits} \
            {input.fa} > {output.barrnapGff}

        grep -v '#' {output.barrnapGff} | cut -f 1,4,5 | sed "s/$/\t200/" > {output.barrnapBedgraph}

        ''' 

rule add_tracks_to_pretext:

    input:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        telomereBedgraph="{sample}/s6_build_pretextmap/telomeres/telomeres_telomeric_repeat_windows.bedgraph",
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hic.bedgraph",
        hifiCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/coverage_output_hifi.bedgraph",
        covGapsBedgraph="{sample}/s6_build_pretextmap/coverage_gaps/gaps.bedgraph",
        gfastatsGapsBedgraph="{sample}/s6_build_pretextmap/gfastats_gaps/gfastats_gaps.bedgraph",
        barrnapBedgraph="{sample}/s6_build_pretextmap/barrnap/RNA.bedgraph",

    output:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.wtracks.map.pretext",
        
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.10)
    resources:
        runtime=120,
        mem_mb=10000,
    shell:
        '''

        cp {input.pretextmap} {output.pretextmap}

        cat {output.telomereBedgraph} | PretextGraph -i {output.pretextmap} -n "telomeres"  
        cat {output.hicCovBedgraph} | PretextGraph -i {output.pretextmap} -n "coveragehic"
        cat {output.hifiCovBedgraph} | PretextGraph -i {output.pretextmap} -n "coveragehifi"
        cat {output.covGapsBedgraph} | PretextGraph -i {output.pretextmap} -n "gaps"
        cat {output.gfastatsGapsBedgraph} | PretextGraph -i {output.pretextmap} -n "gfastats_gaps"
        cat {output.barrnapBedgraph} | PretextGraph -i {output.pretextmap} -n "rDNA"
        ''' 


# rule create_trash_bedgraph:
#     input:
#         fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/yahs_hap1_hap2.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
#        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
#         bam="{sample}/s5_run_yahs/hap{haplomeNumber}/mapped.PT.bam"
#         hifiFastq="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz"
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
#         hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
#         hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file'],
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