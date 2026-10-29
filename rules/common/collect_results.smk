expandPatternGenomescope="{{sample}}/genomescope_plots/{{sample}}_{png}"
expandPatternSmudgeplot="{{sample}}/smudgeplots/{{sample}}_{png}"
expandPatternFastaPaths="{{sample}}/hifiasm_file_paths/{{sample}}_hifiasm_{fastaType}_paths"

# no more doublke braces because we are not trying to expand only png and not sample for example
expandPatternGenomescopeZipFolder="{s}/genomescope_plots/{s}_{png}"
expandPatternSmudgeplotZipFolder="{s}/smudgeplots/{s}_{png}"
expandPatternFastaPathsZipFolder="{s}/hifiasm_file_paths/{s}_hifiasm_{fastaType}_paths"


rule build_rsync_folder_single_sample:
    input:
        hifiasm='{sample}/s4b_run_hifiasm_qc/qc_summary/{sample}_hifiasm_qc.xlsx',
        fasta=lambda w:expand(f'{w.sample}/s4b_run_hifiasm_qc/hifiasm_file_paths/{w.sample}_hifiasm_{{fastaType}}_paths.txt',fastaType=['unitig','contigs']),
        genomescope=lambda w:expand(f"{w.sample}/s3_estimate_ploidy/genomescope/{w.sample}_{{png}}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        smudgeplot=lambda w:expand(f"{w.sample}/s3_estimate_ploidy/smudgeplot/{w.sample}_{{png}}.png",png=['smudgeplot','smudgeplot_log10','centralities']),
        seqkit='{sample}/s2_input_data_qc/qc_results/{sample}_input_data_qc.xlsx',
        pretext=lambda w:f"{w.sample}/s6_build_pretextmap/pretextmap/{w.sample}.map.pretext" if sampleInfo[w.sample]['has_hic'] else [],
    output:
        expand(f"{config['rsync_folder']}/individual_samples/{expandPatternSmudgeplot}_{CURRENTDATE}.png",png=['smudgeplot','smudgeplot_log10','centralities']),
        expand(f"{config['rsync_folder']}/individual_samples/{expandPatternGenomescope}_{CURRENTDATE}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        expand(f'{config['rsync_folder']}/individual_samples/{expandPatternFastaPaths}_{CURRENTDATE}.txt',fastaType=['unitig','contigs']),
        hifiasm_out=f'{config['rsync_folder']}/individual_samples/{{sample}}/hifiasm_qc_excel/{{sample}}_hifiasm_qc_{CURRENTDATE}.xlsx',
        seqkit_out=f'{config['rsync_folder']}/individual_samples/{{sample}}/input_reads_qc_excel/{{sample}}_input_data_qc_{CURRENTDATE}.xlsx',
        pretext=f'{config['rsync_folder']}/individual_samples/{{sample}}/pretextmap/{{sample}}_{CURRENTDATE}.map.pretext',
    params:
        inputDate=CURRENTDATE,
        outdir=lambda w:f"{config['rsync_folder']}/individual_samples/{w.sample}",
        hasHiCShell=lambda w: 'true' if sampleInfo[w.sample]["has_hic"] else 'false',
    threads:
        1
    resources:
        runtime=40,
        mem_mb=1000,
    localrule:
        True
    shell:
        '''
        cp {input.hifiasm} {output.hifiasm_out}
        cp {input.seqkit} {output.seqkit_out}

        # get fasta paths
        for f in {input.fasta};do

            file=$(basename $f)
            cp $f {params.outdir}/hifiasm_file_paths/"${{file%.txt}}_{params.inputDate}.txt"

        done
        
        # get genomescope plots
        for f in {input.genomescope};do
            file=$(basename $f)
            cp $f {params.outdir}/genomescope_plots/"${{file%.png}}_{params.inputDate}.png"
        done

        # get smudgeplots
        for f in {input.smudgeplot};do
            file=$(basename $f)
            cp $f {params.outdir}/smudgeplots/"${{file%.png}}_{params.inputDate}.png"
        done

        # get pretextmap if HiC used
        if [ {params.hasHiCShell} = true ];then
            cp {input.pretext} {params.outdir}/pretextmap/{wildcards.sample}_{params.inputDate}.map.pretext
        else
            echo 'No HiC data was given for {wildcards.sample} so I am empty' > {params.outdir}/pretextmap/{wildcards.sample}_{params.inputDate}.map.pretext
        fi

        '''

rule zip_rsync_folder:
    input:
        expand(f"{config['rsync_folder']}/individual_samples/{expandPatternSmudgeplotZipFolder}_{CURRENTDATE}.png",s=SAMPLES,png=['smudgeplot','smudgeplot_log10','centralities']),
        expand(f"{config['rsync_folder']}/individual_samples/{expandPatternGenomescopeZipFolder}_{CURRENTDATE}.png",s=SAMPLES,png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        expand(f'{config['rsync_folder']}/individual_samples/{expandPatternFastaPaths}_{CURRENTDATE}.txt',s=SAMPLES,fastaType=['unitig','contigs']),
        expand(f'{config['rsync_folder']}/individual_samples/{{s}}/hifiasm_qc_excel/{{s}}_hifiasm_qc_{CURRENTDATE}.xlsx',s=SAMPLES),
        expand(f'{config['rsync_folder']}/individual_samples/{{s}}/input_reads_qc_excel/{{s}}_input_data_qc_{CURRENTDATE}.xlsx',s=SAMPLES),
        expand(f'{config['rsync_folder']}/individual_samples/{{s}}/pretextmap/{{s}}_{CURRENTDATE}.map.pretext',s=SAMPLES),
        expand(f'{config['rsync_folder']}/all_samples/hifiasm_file_paths/{ANALYSIS_NAME}_hifiasm_{{fastaType}}_paths_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.txt',fastaType=['unitig','contigs']),
        f'{config['rsync_folder']}/all_samples/input_data_qc_excel/{ANALYSIS_NAME}_input_data_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx',
        f'{config['rsync_folder']}/all_samples/hifiasm_qc_excel/{ANALYSIS_NAME}_hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx',
    output:
        f"{config['rsync_folder']}/{ANALYSIS_NAME}_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.zip",
    params:
        outdir=f"{config['rsync_folder']}",
    localrule:
        True
    threads:
        1
    resources:
        runtime=60,
        mem_mb=1000,
    shell:
        '''
        zip -r {output} {params.outdir}/individual_samples {params.outdir}/all_samples
        '''

