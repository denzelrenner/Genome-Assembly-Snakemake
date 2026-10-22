# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
from datetime import datetime

# wildcard constaints
wildcard_constraints:
    sample='|'.join(SAMPLES), # limit value of sample to our input bams


### GENERAL ###

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

    # get target qc file for all individuals
    for smp in SAMPLES:
        outFiles.extend(expand(f"{smp}/s4_run_hifiasm/quality_control/{smp}.hap{{hap}}/QC_Summary/output_MasterSpreadsheet.xlsx",hap=[i for i in range(1,sampleInfo[smp]['hifiasm_ploidy']+1)]))

    # regardless of HiC or not, add smudgeplot and genomescope to target files, and remove as input for hifiasm
    outFiles.extend(expand("{smp}/s3_estimate_ploidy/smudgeplot/plots/{png}.png",smp=SAMPLES,png=['output_smudgeplot','output_smudgeplot_log10','output_centralities']))
    outFiles.extend(expand("{smp}/s3_estimate_ploidy/genomescope/plots/{png}.png",smp=SAMPLES,png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot']))

    # add hifi stats file for all individuals
    outFiles.extend(expand("{smp}/s2_input_data_qc/Hifi/{smp}.hifi.tsv",smp=SAMPLES))
    
    # add hic stats file for HIC individuals
    outFiles.extend(expand("{smp}/s2_input_data_qc/HiC{pairedEnd}/{smp}.hic{pairedEnd}.tsv",smp=HICSAMPLES,pairedEnd=[1,2]))

    # get current date which will be used to name combined QC files
    CURRENTDATE=datetime.today().strftime('%d_%m_%Y')
    
    return outFiles


