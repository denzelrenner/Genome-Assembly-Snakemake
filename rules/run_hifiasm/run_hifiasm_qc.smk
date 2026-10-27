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
            compleasm=f'{smp}/s4b_run_hifiasm_qc/{smp}.hap{i}/compleasm/summary.txt'
            gfastats=f'{smp}/s4b_run_hifiasm_qc/{smp}.hap{i}/gfastats/stats.txt'

            # set as target files
            outputFiles.extend([compleasm,gfastats])

    return outputFiles

rule download_compleasm_lineages:
    output:
        done=f"{config['compleasm_download_path']}/{{lineage_odb}}.done",
    params:
        libpath=config['compleasm_download_path'],
        lineage=lambda w: w.lineage_odb.split('_')[0],
        odb=lambda w: w.lineage_odb.split('_')[1],
    log:
        stdout=f"{config['compleasm_download_path']}/download_logs/{{lineage_odb}}.out",
        stderr=f"{config['compleasm_download_path']}/download_logs/{{lineage_odb}}.err"
    resources:
        runtime=240,
    threads:
        1
    priority:
        200
    conda:
        "../../envs/qc-tools.yml"
    localrule:
        True
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # make download log file if it does not exist
        mkdir -p {params.libpath}/download_logs

        compleasm download -L {params.libpath} --odb {params.odb} {params.lineage}
        '''



rule hifiasm_qc_on_individual_haplome:
    input:
        fa='{sample}/s4a_run_hifiasm/{sample}.hap{haplomeNumber}.fa',
        compleasm_lineage_downloads=expand(f"{config['compleasm_download_path']}/{{lineage}}.done",lineage=ALLCOMPLEASMLINEAGES)
    output:
        compleasm='{sample}/s4b_run_hifiasm_qc/{sample}.hap{haplomeNumber}/compleasm/summary.txt',
        gfastats='{sample}/s4b_run_hifiasm_qc/{sample}.hap{haplomeNumber}/gfastats/stats.txt',
    params:
        odb=lambda w: sampleInfo[w.sample]['odb'],
        lineage=lambda w: sampleInfo[w.sample]['lineage'],
        libpath=config['compleasm_download_path'],
    resources:
        runtime=240,
    threads:
        int(workflow.cores * 0.5)
    conda:
        "../../envs/qc-tools.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # run compleasm
        compleasm run -a {input.fa} -o {wildcards.sample}/s4b_run_hifiasm_qc/{wildcards.sample}.hap{wildcards.haplomeNumber}/compleasm -l {params.lineage} --odb {params.odb} -L {params.libpath} -t {threads}

        # run gfastats
        gfastats -f {input.fa} -j {threads} > {output.gfastats}
        '''

rule summarise_hifiasm_qc_single_sample:
    input:
        compleasm=lambda w:expand('{s}/s4b_run_hifiasm_qc/{s}.hap{hapNum}/compleasm/summary.txt',s=[w.sample],hapNum=sampleInfo[w.sample]['hifiasmploidyAsArray']),
        gfastats=lambda w:expand('{s}/s4b_run_hifiasm_qc/{s}.hap{hapNum}/gfastats/stats.txt',s=[w.sample],hapNum=sampleInfo[w.sample]['hifiasmploidyAsArray']),
    output:
        qcSummary=f'{{sample}}/s4b_run_hifiasm_qc/qc_summary/{{sample}}_hifiasm_qc_{CURRENTDATE}.xlsx',
    params:
        qcScript=os.path.join(config['script_dir'],'python_scripts','create_hifiasm_qc_summary.py'),
        ploidy=lambda w:sampleInfo[w.sample]['hifiasm_ploidy'],
        outfilePrefix=f'{{sample}}_hifiasm_qc_{CURRENTDATE}',
    resources:
        runtime=30,
        mem_mb=1000,
    localrule:
        True
    threads:
        1
    conda:
        "../../envs/pd-tools.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # run qc summary excel file
        python3 {params.qcScript} -s {wildcards.sample} -nhap {params.ploidy} -o {params.outfilePrefix} -od {wildcards.sample}/s4b_run_hifiasm_qc/qc_summary
        '''

rule summarise_hifiasm_qc_all_samples:
    input:
        get_hifiasm_qc_files_all_samples
    output:
        excel=f'results/all_samples/hifiasm_qc/hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx',
    params:
        qcScript=os.path.join(config['script_dir'],'python_scripts','create_hifiasm_qc_summary.py'),
        outfilePrefix=f'hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}',
        inputSamples=SAMPLES,
        ploidySamples=HIFIASMPLOIDYFORSAMPLES,
    resources:
        runtime=30,
        mem_mb=1000,
    localrule:
        True
    threads:
        1,
    conda:
        "../../envs/pd-tools.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        python3 {params.qcScript} -s {params.inputSamples} -nhap {params.ploidySamples} -o {params.outfilePrefix} -od results/all_samples/hifiasm_qc
        '''

