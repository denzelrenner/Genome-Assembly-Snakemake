import os

# get compleasm and gfastats results for all samples, this will become space separated list of input files
def get_hifiasm_qc_files_all_samples(w):

    # set array of out files
    outputFiles = []

    # get target hifiasm qc file for all individuals
    for smp in SAMPLES:
        
        # get hifiasm qc outfiles
        for i in sampleInfo[smp]['hifiasmploidyAsArray']:
        
            # create files
            compleasm=f'{smp}/s4_run_hifiasm/quality_control/{smp}.hap{i}/compleasm/summary.txt'
            gfastats=f'{smp}/s4_run_hifiasm/quality_control/{smp}.hap{i}/gfastats/stats.txt'

            # set as target files
            outputFiles.extend([compleasm,gfastats])

    return outputFiles


rule hifiasm_qc_on_individual_haplome:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa',
    output:
        compleasm='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/compleasm/summary.txt',
        gfastats='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/gfastats/stats.txt',
    params:
        odb=config['compleasm_odb'],
        lineage=config['compleasm_lineage'],
        libpath=config['compleasm_download_path'],
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
        '''

rule summarise_hifiasm_qc_single_sample:
    input:
        compleasm=lambda w:expand('{s}/s4_run_hifiasm/quality_control/{s}.hap{hapNum}/compleasm/summary.txt',s=[w.sample],hapNum=sampleInfo[w.sample]['hifiasmploidyAsArray']),
        gfastats=lambda w:expand('{s}/s4_run_hifiasm/quality_control/{s}.hap{hapNum}/gfastats/stats.txt',s=[w.sample],hapNum=sampleInfo[w.sample]['hifiasmploidyAsArray']),
    output:
        qcSummary=f'{{sample}}/s4_run_hifiasm/quality_control/qc_summary/{{sample}}_hifiasm_qc_{CURRENTDATE}.xlsx',
    params:
        qcScript=os.path.join(config['script_dir'],'python_scripts','create_hifiasm_qc_summary.py'),
        ploidy=lambda w:sampleInfo[w.sample]['hifiasm_ploidy'],
        outfilePrefix=f'{{sample}}_hifiasm_qc_{CURRENTDATE}',
    resources:
        runtime=30,
    threads:
        1
    conda:
        "../../envs/qc-tools.yml"
    shell:
        '''
        # run qc summary excel file
        python3 {params.qcScript} -s {wildcards.sample} -nhap {params.ploidy} -o {params.outfilePrefix} -od {wildcards.sample}/s4_run_hifiasm/quality_control/qc_summary
        '''

rule summarise_hifiasm_qc_all_samples:
    input:
        get_hifiasm_qc_files_all_samples
    output:
        excel=f'summary_results/hifiasm_qc/hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx',
    params:
        qcScript=os.path.join(config['script_dir'],'python_scripts','create_hifiasm_qc_summary.py'),
        outfilePrefix=f'hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}',
        inputSamples=SAMPLES,
        ploidySamples=HIFIASMPLOIDYFORSAMPLES,
    resources:
        runtime=30,
    threads:
        1,
    conda:
        "../../envs/qc-tools.yml"
    shell:
        '''
        python3 {params.qcScript} -s {params.inputSamples} -nhap {params.ploidySamples} -o {params.outfilePrefix} -od summary_results/hifiasm_qc
        '''

