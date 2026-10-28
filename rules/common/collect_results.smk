expandPatternGenomescope="{{sample}}/genomescope_plots/{{sample}}_{png}"
expandPatternSmudgeplot="{{sample}}/smudgeplots/{{sample}}_{png}"

rule build_rsync_folder_single_sample:
    input:
        hifiasm='{sample}/s4b_run_hifiasm_qc/qc_summary/{sample}_hifiasm_qc.xlsx',
        genomescope=lambda w:expand(f"{w.sample}/s3_estimate_ploidy/genomescope/{w.sample}_{{png}}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        smudgeplot=lambda w:expand(f"{w.sample}/s3_estimate_ploidy/smudgeplot/{w.sample}_{{png}}.png",png=['smudgeplot','smudgeplot_log10','centralities']),
        seqkit='{sample}/s2_input_data_qc/qc_results/{sample}_input_data_qc.xlsx',
        pretext=lambda w:f"{w.sample}/s6_build_pretextmap/pretextmap/{w.sample}.map.pretext" if sampleInfo[w.sample]['useHiC'] else [],
    output:
        expand(f"{config['rsync_folder']}/individual_sample/{expandPatternSmudgeplot}_{CURRENTDATE}.png",png=['smudgeplot','smudgeplot_log10','centralities']),
        expand(f"{config['rsync_folder']}/individual_sample/{expandPatternGenomescope}_{CURRENTDATE}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']),
        hifiasm_out=f'{config['rsync_folder']}/individual_sample/{{sample}}/hifiasm_qc_excel/{{sample}}_hifiasm_qc_{CURRENTDATE}.xlsx',
        seqkit_out=f'{config['rsync_folder']}/individual_sample/{{sample}}/input_reads_qc_excel/{{sample}}_input_data_qc_{CURRENTDATE}.xlsx',
        pretext=f'{config['rsync_folder']}/individual_sample/{{sample}}/pretextmap/{{sample}}_{CURRENTDATE}.map.pretext',
    params:
        inputDate=CURRENTDATE,
        outdir=lambda w:f"{config['rsync_folder']}/individual_sample/{w.sample}",
        useHiCShell=lambda w: 'true' if sampleInfo[w.sample]["useHiC"] else 'false',
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
        if [ {params.useHiCShell} = true ];then
            cp {input.pretext} {params.outdir}/pretextmap/{wildcards.sample}_{params.inputDate}.map.pretext
        else
            echo 'No HiC data was given so I am empty' > {params.outdir}/pretextmap/{wildcards.sample}_{params.inputDate}.map.pretext
        fi

        '''

