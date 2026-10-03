rule gather_end_results:
    input:
        '{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/QC_Summary/output_MasterSpreadsheet.xlsx'

    output:
        '{sample}/s0_end_pipe/{sample}.hap{haplomeNumber}/done.txt'
    
    threads:
        1
        
    shell:
        "touch {output}"
