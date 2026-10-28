def return_contig_fasta_file_paths(w):

    # set array of out files
    outputFiles = []

    # get target hifiasm qc file for all individuals
    for smp in SAMPLES:
        
        # get hifiasm qc outfiles
        for i in sampleInfo[smp]['hifiasmploidyAsArray']:
        
            # create files
            contigFasta=f'{smp}/s4a_run_hifiasm/{smp}.hap{i}.fa'

            # set as target files
            outputFiles.append(contigFasta)

    return outputFiles

def return_unitig_fasta_file_paths(w):

    # set array of out files
    outputFiles = []

    # get target hifiasm qc file for all individuals
    for smp in SAMPLES:

        # add unitig fasta
        unitigFasta=f'{smp}/s4a_run_hifiasm/{smp}.p_utg.fa'

        outputFiles.append(unitigFasta)

    return outputFiles

rule pctg_gfa_to_fa:
    input:
        '{sample}/s4a_run_hifiasm/.done.txt',
    output:
        '{sample}/s4a_run_hifiasm/{sample}.hap{haplomeNumber}.fa'
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
        awk '/^S/{{print ">"$2;print $3}}' {wildcards.sample}/s4a_run_hifiasm/{wildcards.sample}.hap{wildcards.haplomeNumber}.p_ctg.gfa > {output}
        '''
        
rule putg_gfa_to_fa:
    input:
        '{sample}/s4a_run_hifiasm/.done.txt',
    output:
        '{sample}/s4a_run_hifiasm/{sample}.p_utg.fa',
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
        awk '/^S/{{print ">"$2;print $3}}' {wildcards.sample}/s4a_run_hifiasm/{wildcards.sample}.p_utg.gfa > {output}
        '''

# this rule is here for now because it is dependent on the top 2 rules actually completing
rule get_hifiasm_fasta_paths:
    input:
        contigs=return_contig_fasta_file_paths,
        unitigs=return_unitig_fasta_file_paths,
    output:
        expand(f'{config['rsync_folder']}/all_samples/hifiasm_file_paths/{ANALYSIS_NAME}_hifiasm_{{fastaType}}_paths_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.txt',fastaType=['unitig','contigs']),
    params:
        workingDir=SNAKEMAKE_OUTPUTDIR,
    threads:
        1
    resources:
        runtime=40,
        mem_mb=1000,
    localrule:
        True
    shell:
        '''
        ls {input.contigs} > {output[0]}
        ls {input.unitigs} > {output[1]}

        # use awk to prepend the working dir to file names
        awk '{{print "{params.workingDir}/"$1}}' {output[0]} > tmp && mv tmp {output[0]}
        awk '{{print "{params.workingDir}/"$1}}' {output[1]} > tmp && mv tmp {output[1]}
        '''