# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
from datetime import datetime

# wildcard constaints
wildcard_constraints:
    sample='|'.join(SAMPLES), # limit value of sample to our input bams


### GENERAL ###

# get current date which will be used to name combined QC files
CURRENTDATE=datetime.today().strftime('%d_%m_%Y')


def set_rule_target():

    outFiles = []

    # get target files for individuals with HiC
    #outFiles.extend(expand("{smp}/s6_build_pretextmap/{outFile}",smp=HICSAMPLES,outFile=SCAFFOLDINGTARGETFILES))
    outFiles.extend(expand("{smp}/s6_build_pretextmap/pretextmap/{smp}.map.pretext",smp=HICSAMPLES))

    # regardless of HiC or not, add smudgeplot and genomescope to target files, and remove as input for hifiasm
    outFiles.extend(expand("{smp}/s3_estimate_ploidy/smudgeplot/{smp}_{png}.png",smp=SAMPLES,png=['smudgeplot','smudgeplot_log10','centralities']))
    outFiles.extend(expand("{smp}/s3_estimate_ploidy/genomescope/{smp}_{png}.png",smp=SAMPLES,png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']))

    # add hifi stats file for all individuals
    outFiles.extend(expand("{smp}/s2_input_data_qc/Hifi/{smp}.hifi.tsv",smp=SAMPLES))
    
    # add hic stats file for HIC individuals
    outFiles.extend(expand("{smp}/s2_input_data_qc/HiC{pairedEnd}/{smp}.hic{pairedEnd}.tsv",smp=HICSAMPLES,pairedEnd=[1,2]))

    # get target hifiasm qc file for all individuals
    for smp in SAMPLES:
        
        outFiles.extend(expand('{s}/s4b_run_hifiasm_qc/{s}.hap{i}/compleasm/summary.txt',s=[smp],i=sampleInfo[smp]['hifiasmploidyAsArray']))
        outFiles.extend(expand('{s}/s4b_run_hifiasm_qc/{s}.hap{i}/gfastats/stats.txt',s=[smp],i=sampleInfo[smp]['hifiasmploidyAsArray']))
        outFiles.append(f'{smp}/s4b_run_hifiasm_qc/qc_summary/{smp}_hifiasm_qc.xlsx')

    # add excel spreadhseet which summarises hifiasm results for all samples
    outFiles.append(f'{config['rsync_folder']}/all_samples/hifiasm_qc_excel/hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx')

    # get target files for compleasm dbs
    outFiles.extend(expand(f"{config['compleasm_download_path']}/{{lineage}}.done",lineage=ALLCOMPLEASMLINEAGES))

    # get target files for sequencing read qc
    # individual sample files
    outFiles.extend(expand(f'{{s}}/s2_input_data_qc/qc_results/{{s}}_input_data_qc.xlsx',s=SAMPLES))
    
    # all samples file 
    outFiles.append(f'{config['rsync_folder']}/all_samples/input_data_qc_excel/input_data_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx')

    outFiles.extend(expand(f'{config['rsync_folder']}/individual_sample/{{smp}}/hifiasm_qc_excel/{{smp}}_hifiasm_qc_{CURRENTDATE}.xlsx',smp=SAMPLES))
    outFiles.extend(expand(f"{config['rsync_folder']}/individual_sample/{{smp}}/smudgeplots/{{smp}}_{{png}}_{CURRENTDATE}.png",smp=SAMPLES,png=['smudgeplot','smudgeplot_log10','centralities']))
    outFiles.extend(expand(f"{config['rsync_folder']}/individual_sample/{{smp}}/genomescope_plots/{{smp}}_{{png}}_{CURRENTDATE}.png",smp=SAMPLES,png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']))
    outFiles.extend(expand(f'{config['rsync_folder']}/individual_sample/{{smp}}/input_reads_qc_excel/{{smp}}_input_data_qc_{CURRENTDATE}.xlsx',smp=SAMPLES))
    outFiles.extend(expand(f'{config['rsync_folder']}/individual_sample/{{smp}}/pretextmap/{{smp}}_{CURRENTDATE}.map.pretext',smp=HICSAMPLES))

    return outFiles


