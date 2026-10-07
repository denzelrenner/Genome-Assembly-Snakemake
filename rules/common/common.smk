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

def set_rule_target():

    outFiles = []

    for smpl in sampleInfo:

        for i in range(1,sampleInfo[smpl]['ploidy']+1):
            outFiles.append(f'{smpl}/s0_end_pipe/{smpl}.hap{i}/done.txt')

    return outFiles


