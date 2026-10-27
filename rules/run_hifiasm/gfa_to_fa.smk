rule pctg_gfa_to_fa:
    input:
        '{sample}/s4_run_hifiasm/.done.txt',
    output:
        '{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa'
    threads:
        1
    resources:
        runtime=10,
        mem_mb=1000,
    localrule:
        True
    shell:
        '''
        # convert to fasta file
        awk '/^S/{{print ">"$2;print $3}}' {wildcards.sample}/s4_run_hifiasm/{wildcards.sample}.hap{wildcards.haplomeNumber}.p_ctg.gfa > {output}
        '''
        
rule putg_gfa_to_fa:
    input:
        '{sample}/s4_run_hifiasm/.done.txt',
    output:
        '{sample}/s4_run_hifiasm/{sample}.p_utg.fa',
    threads:
        1
    localrule:
        True
    resources:
        runtime=15,
        mem_mb=1000
    shell:
        '''
        # convert to fasta file
        awk '/^S/{{print ">"$2;print $3}}' {wildcards.sample}/s4_run_hifiasm/{wildcards.sample}.p_utg.gfa' > {output}
        '''