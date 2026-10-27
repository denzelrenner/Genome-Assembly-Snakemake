rule gather_end_results:
    input:
        

    output:
        '{sample}/s0_end_pipe/{sample}.hap{haplomeNumber}/done.txt'
    
    threads:
        1
        
    shell:
        "touch {output}"
