rule pctg_gfa_to_fa:
    input:
        "{sample}/s4_run_hifiasm/.done.txt"
    output:
        '{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa'
    threads:
        1
    resources:
        runtime=10,
        mem_mb=1000
    shell:
        '''
        # go through all the gfa files
        for gfa in {wildcards.sample}/s4_run_hifiasm/*hap*.p_ctg.gfa; do

            echo $gfa

            if [ -f $gfa ]; then

                # get only file name and not whole path
                file=$(basename "${{gfa}}")

                # find haplome number
                hap=$(awk -F'.' '{{print $(NF-2)}}' <<< $file)

                # convert to fasta file
                awk '/^S/{{print ">"$2;print $3}}' $gfa > {wildcards.sample}/s4_run_hifiasm/{wildcards.sample}."${{hap}}".fa

            fi
        done
        '''
        
rule putg_gfa_to_fa:
    input:
        '{sample}/s4_run_hifiasm/{sample}.p_utg.gfa',
    output:
        '{sample}/s4_run_hifiasm/{sample}.p_utg.fa',
    threads:
        1
    resources:
        runtime=15,
        mem_mb=1000
    shell:
        '''
        # convert to fasta file
        awk '/^S/{{print ">"$2;print $3}}' {input} > {output}
        '''