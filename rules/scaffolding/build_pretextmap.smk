rule build_pretextmap:
    input:
        bam="{sample}/s6_build_pretextmap/mapping_coverage/hic/hic_sorted.bam",
        telomereBedgraph="{sample}/s6_build_pretextmap/telomeres/telomeres_telomeric_repeat_windows.bedgraph",
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/hic/coverage_output_hic.bedgraph",
        hifiCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/hifi/coverage_output_hifi.bedgraph",
        covGapsBedgraph="{sample}/s6_build_pretextmap/coverage_gaps/gaps.bedgraph",
        gfastatsGapsBedgraph="{sample}/s6_build_pretextmap/gfastats_gaps/gfastats_gaps.bedgraph",
        barrnapBedgraph="{sample}/s6_build_pretextmap/barrnap/RNA.bedgraph",
        trashBedgraph="{sample}/s6_build_pretextmap/trash/trash_filtered_repeat_arrays.bedgraph",

    output:
        pretextmap="{sample}/s6_build_pretextmap/pretextmap/{sample}.map.pretext",
        
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        int(workflow.cores * 0.10)
    resources:
        runtime=1200,
        mem_mb=300000,
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        samtools view -@{threads} -h {input.bam} | PretextMap -o {output.pretextmap} --mapq 0

        cat {input.telomereBedgraph} | PretextGraph -i {output.pretextmap} -n "telomeres"  
        cat {input.hicCovBedgraph} | PretextGraph -i {output.pretextmap} -n "coveragehic"
        cat {input.hifiCovBedgraph} | PretextGraph -i {output.pretextmap} -n "coveragehifi"
        cat {input.covGapsBedgraph} | PretextGraph -i {output.pretextmap} -n "gaps"
        cat {input.gfastatsGapsBedgraph} | PretextGraph -i {output.pretextmap} -n "gfastats_gaps"
        cat {input.barrnapBedgraph} | PretextGraph -i {output.pretextmap} -n "rDNA"
        cat {input.trashBedgraph} | PretextGraph -i {output.pretextmap} -n "trash"
        ''' 
    
rule create_telomeres_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
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
        tidk search --extension bedgraph --string TTTAGGG --output telomeres --dir {wildcards.sample}/s6_build_pretextmap/telomeres {input.fa}

        ''' 

rule create_hic_coverage_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        hic1=lambda w: sampleInfo[w.sample]['hic1_workflow_file'],
        hic2=lambda w: sampleInfo[w.sample]['hic2_workflow_file'], 
    output:
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/hic/coverage_output_hic.bedgraph",
        bam="{sample}/s6_build_pretextmap/mapping_coverage/hic/hic_sorted.bam",
    conda:
        "../../envs/mapping-tools.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=1200,
        mem_mb=300000
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        echo 'map HiC to assembly'

        bwa index {input.fa}

        bwa mem -t {threads} {input.fa} {input.hic1} {input.hic2} | samtools sort -@{threads} -o {output.bam}

        echo 'coverage graph hic'

        bedtools genomecov -ibam {output.bam} -bga > {output.hicCovBedgraph}
        ''' 

rule create_hifi_coverage_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        hifiFastq="{sample}/s1_get_input_data/Hifi/{sample}.fastq.gz",
    output:
        hifiCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/hifi/coverage_output_hifi.bedgraph",
        bam="{sample}/s6_build_pretextmap/mapping_coverage/hifi/hifi_sorted.bam",
    conda:
        "../../envs/manual-curation.yml"
    threads:
        int(workflow.cores * 0.75)
    resources:
        runtime=720,
        mem_mb=100000
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        echo 'minimap mapping hifi reads to scaffold assem'

        minimap2 -ax map-pb -t {threads} {input.fa} {input.hifiFastq} | samtools sort -@{threads} -O BAM -o {output.bam} -
        
        echo 'coverage graph hifi'
        
        bedtools genomecov -ibam {output.bam} -bga > {output.hifiCovBedgraph}
        ''' 

rule create_hic_gaps_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
        hicCovBedgraph="{sample}/s6_build_pretextmap/mapping_coverage/hic/coverage_output_hic.bedgraph",
    output:
        covGapsBedgraph="{sample}/s6_build_pretextmap/coverage_gaps/gaps.bedgraph",
    threads:
        2
    resources:
        runtime=30,
        mem_mb=1000
    localrule:
        True
    shell:
        '''
        echo 'gaps graph'

        grep -w 0$ {input.hicCovBedgraph} | sed 's/0$/200/g' > {output.covGapsBedgraph}
        ''' 

rule create_gfastats_gaps_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
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
        export PATH=$CONDA_PREFIX/bin:$PATH

        echo 'gfastats gaps'	

        gfastats {input.fa} -b gaps | sed "s/$/\t200/" > {output.gfastatsGapsBedgraph}

        ''' 

rule create_barrnap_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
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
        export PATH=$CONDA_PREFIX/bin:$PATH

        echo 'rDNA'

        barrnap --threads {threads} \
            --kingdom euk \
            --outseq {output.barrnapHits} \
            {input.fa} > {output.barrnapGff}

        grep -v '#' {output.barrnapGff} | cut -f 1,4,5 | sed "s/$/\t200/" > {output.barrnapBedgraph}

        ''' 

rule create_trash_bedgraph:
    input:
        fa=lambda w:f"{w.sample}/s5_run_yahs/combined_haps/scaffolds.fa" if sampleInfo[w.sample]['ploidy'] == 2 else f"{w.sample}/s5_run_haphic/04.build/scaffolds.fa",
    output:
        bed="{sample}/s6_build_pretextmap/trash/trash_filtered_repeat_arrays.bed",
        bedgraph="{sample}/s6_build_pretextmap/trash/trash_filtered_repeat_arrays.bedgraph",
    params:
        trash_dir=os.path.join(config['script_dir'],'TRASH_2'),
        trash_src_dir=os.path.join(config['script_dir'],'TRASH_2','src'),
        trash_filter_dir=os.path.join(config['script_dir'],'trash_filter_scripts'),
        outdir=lambda w:f"{w.sample}/s6_build_pretextmap/trash",
    threads:
        16
    resources:
        runtime=240,
        mem_mb=10000,
    conda:
        '../../envs/trash-v2.yml'
    shell:
        '''
        # export PATHs
        export PATH=$CONDA_PREFIX/bin:$PATH
        export PATH={params.trash_dir}:$PATH
        export PATH={params.trash_src_dir}:$PATH
        export PATH={params.trash_filter_dir}:$PATH

        export TRASH_FILTER_SCRIPTS_DIR={params.trash_filter_dir}

        # run TRASH
        TRASH.R -f {input.fa} -o {params.outdir} -p {threads}

        # filter trash output 
        filter_trash.R {params.outdir}/scaffolds.fa_repeats_with_seq.csv {params.outdir}/scaffolds.fa_arrays.csv {params.outdir}/scaffolds.fa_repeats_filtered.csv {params.outdir}/scaffolds.fa_arrays_filtered.csv

        # csv to bedgraph
        awk -F "," -v 'OFS=\\t' '{{ print $3,$1,$2,$11 }}' {params.outdir}/scaffolds.fa_arrays_filtered.csv | sed s/'"'/''/g | sed s/'seqID'/'Chr'/g > {output.bedgraph}

        # create gff
        cut -f 1,2,3 {output.bedgraph} > {output.bed}
        '''
