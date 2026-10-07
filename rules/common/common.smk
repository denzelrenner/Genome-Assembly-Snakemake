# helper function to be used in the actual snakefile, and in other snakemake files as well
import os

# set some global variables
HASHIC = True if (config['hic1_reads'] and config['hic2_reads']) else False

# get samples by reading hifi reads file
#SAMPLES = [f.split('.')[0] for f in os.listdir(config['hifi_reads'])]
# HAPLOMENUMBERS = [str(i) for i in range(1,config['ploidy']+1)] if HASHIC else ['1','2']

# join path of hic reads and sample, then save in dir
HIC1FILES={s:"" for s in SAMPLES}
HIC2FILES={s:"" for s in SAMPLES}
HICFILES={s:{"HIC1":"","HIC2":""} for s in SAMPLES}
MERGEHIC={s:False for s in SAMPLES}
SMUDGEPLOTFILES={s:"" for s in SAMPLES}
GENOMESCOPEFILES={s:"" for s in SAMPLES}

# wildcard constaints
wildcard_constraints:
    sample='|'.join(SAMPLES), # limit value of sample to our input bams


# helper functions

### GENERAL ###

# set list of output files for HiC 
HICTARGETFILES=["barrnap/RNA.gff","barrnap/all_hits.fa","barrnap/RNA.bedgraph",
                "gfastats_gaps/gfastats_gaps.bedgraph",
                "coverage_gaps/gaps.bedgraph",
                "pretextmap/hifi_sorted.bam","mapping_coverage/coverage_output_hifi.bedgraph",
                "mapping_coverage/coverage_output_hic.bedgraph",
                "telomeres/telomeres_telomeric_repeat_windows.bedgraph",
                "pretextmap/hic_sorted.bam"]

def set_rule_target():

    outFiles = []

    # get target files for individuals with HiC
    outFiles.extend(expand("{smp}/s6_build_pretextmap/{outFile}",smp=HICSAMPLES,outFile=HICTARGETFILES))

    # get target file for individuals without HIC
    outFiles.extend(expand("{smp}/s4_run_hifiasm/quality_control/{smp}.hap{hap}/QC_Summary/output_MasterSpreadsheet.xlsx'",smp=NOHICSAMPLES,hap=[1,2]))
    
    return outFiles


