# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
from datetime import datetime

# wildcard constaints
wildcard_constraints:
    sample='|'.join(SAMPLES), # limit value of sample to our input bams


### GENERAL ###

# get current date which will be used to name combined QC files
CURRENTDATE=datetime.today().strftime('%d_%m_%Y')

# set list of output files for HiC 
SCAFFOLDINGTARGETFILES=["barrnap/RNA.gff","barrnap/all_hits.fa","barrnap/RNA.bedgraph",
                "gfastats_gaps/gfastats_gaps.bedgraph",
                "coverage_gaps/gaps.bedgraph",
                "pretextmap/hifi_sorted.bam","mapping_coverage/coverage_output_hifi.bedgraph",
                "mapping_coverage/coverage_output_hic.bedgraph",
                "telomeres/telomeres_telomeric_repeat_windows.bedgraph",
                "pretextmap/hic_sorted.bam"]

def set_rule_target():

    outFiles = []

    # get target files for individuals with HiC
    #outFiles.extend(expand("{smp}/s6_build_pretextmap/{outFile}",smp=HICSAMPLES,outFile=SCAFFOLDINGTARGETFILES))
    outFiles.extend(expand("{smp}/s6_build_pretextmap/pretextmap/{smp}.map.pretext",smp=HICSAMPLES))

    # regardless of HiC or not, add smudgeplot and genomescope to target files, and remove as input for hifiasm
    outFiles.extend(expand("{smp}/s3_estimate_ploidy/smudgeplot/plots/{png}.png",smp=SAMPLES,png=['output_smudgeplot','output_smudgeplot_log10','output_centralities']))
    outFiles.extend(expand("{smp}/s3_estimate_ploidy/genomescope/plots/{png}.png",smp=SAMPLES,png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']))

    # add hifi stats file for all individuals
    outFiles.extend(expand("{smp}/s2_input_data_qc/Hifi/{smp}.hifi.tsv",smp=SAMPLES))
    
    # add hic stats file for HIC individuals
    outFiles.extend(expand("{smp}/s2_input_data_qc/HiC{pairedEnd}/{smp}.hic{pairedEnd}.tsv",smp=HICSAMPLES,pairedEnd=[1,2]))

    # get target hifiasm qc file for all individuals
    for smp in SAMPLES:
        
        outFiles.extend(expand('{s}/s4_run_hifiasm/quality_control/{s}.hap{i}/compleasm/summary.txt',s=[smp],i=sampleInfo[smp]['hifiasmploidyAsArray']))
        outFiles.extend(expand('{s}/s4_run_hifiasm/quality_control/{s}.hap{i}/gfastats/stats.txt',s=[smp],i=sampleInfo[smp]['hifiasmploidyAsArray']))
        outFiles.append(f'{smp}/s4_run_hifiasm/quality_control/qc_summary/{smp}_hifiasm_qc_{CURRENTDATE}.xlsx')
        outFiles.append(f'{smp}/s2_input_data_qc/qc_results/{smp}_input_data_qc_{CURRENTDATE}.xlsx')

    # add excel spreadhseet which summarises hifiasm results for all samples
    outFiles.append(f'summary_results/hifiasm_qc/hifiasm_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx')

    # get target files for sequencing read qc
    # individual sample files
    # outFiles.extend(expand(f'{{s}}/s2_input_data_qc/qc_results/{{s}}_input_data_qc_{CURRENTDATE}.xlsx',s=SAMPLES))
    
    # all sampels file
    outFiles.append(f'summary_results/input_data_qc/input_data_qc_{NUMBEROFINPUTSAMPLES}samples_{CURRENTDATE}.xlsx')


    return outFiles


