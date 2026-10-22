import os

rule hifiasm_qc_single_sample:
    input:
        fa='{sample}/s4_run_hifiasm/{sample}.hap{haplomeNumber}.fa'
    output:
        compleasm='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/compleasm/summary.txt',
        gfastats='{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/gfastats/stats.txt',
        qcSummary=f'{{sample}}/s4_run_hifiasm/quality_control/{{sample}}.hap{{haplomeNumber}}/QC_Summary/{{sample}}_hifiasm_qc_{CURRENTDATE}.xlsx',
    params:
        odb=config['compleasm_odb'],
        lineage=config['compleasm_lineage'],
        libpath=config['compleasm_download_path'],
        qcScript=os.path.join(config['script_dir'],'python_scripts','create_qc_summary_table.py'),
        outfilePrefix=f'{{sample}}_hifiasm_qc_{CURRENTDATE}'
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

rule hifiasm_qc_all_samples:
    input:
        compleasm=lambda w:expand('{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/compleasm/summary.txt',zip,sample=SAMPLES,haplomeNumber=HIFIASMPLOIDYFORSAMPLES),
        gfastats=lambda w:expand('{sample}/s4_run_hifiasm/quality_control/{sample}.hap{haplomeNumber}/gfastats/stats.txt',zip,sample=SAMPLES,haplomeNumber=HIFIASMPLOIDYFORSAMPLES)
    output:
        excel=f'shared_results/hifiasm_QC/hifiasm_qc_all_samples_{CURRENTDATE}.xlsx',
    params:
        outfilePrefix=f'hifiasm_qc_all_samples_{CURRENTDATE}'
        inputSamples=SAMPLES
        ploidySamples=HIFIASMPLOIDYFORSAMPLES
    resources:
        runtime=30,
    threads:
        1,
    conda:
        "../../envs/qc-tools.yml"
    shell:
        '''
        # run qc of all the assemblies


        '''

