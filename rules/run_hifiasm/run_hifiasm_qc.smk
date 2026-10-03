rule hifiasm_qc:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa'
    output:
        compleasm='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/compleasm/summary.txt',
        gfastats='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/gfastats/stats.txt',
        qcSummary='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/QC_Summary/output_MasterSpreadsheet.xlsx',
    params:
        odb=config['compleasm_odb'],
        lineage=config['compleasm_lineage'],
        libpath=config['compleasm_download_path'],
        qcScript=config['script_dir'] + '/create_qc_summary_table_v3.py',
    resources:
        runtime=240,
    threads:
        int(workflow.cores * 0.5)
    conda:
        "../../envs/qc-tools.yml"
    shell:
        '''
        # run compleasm
        compleasm run -a {input.fa} -o {wildcards.sample}/s4_run_hifiasm/quality_control/{wildcards.sample}.hap{wildcards.haplomeNumber}/compleasm -l {params.lineage} --odb {params.odb} -L {params.libpath} -t {threads}

        # run gfastats
        gfastats -f {input.fa} -j {threads} > {output.gfastats}

        # run qc summary metrics
        {params.qcScript} --fasta {input.fa} \
                --gfastats {output.gfastats} \
                --compleasm {output.compleasm} \
                -o output --assem_name {wildcards.sample}.hap{wildcards.haplomeNumber} \
                -od {wildcards.sample}/s4_run_hifiasm/quality_control/{wildcards.sample}.hap{wildcards.haplomeNumber}/QC_Summary
        '''
