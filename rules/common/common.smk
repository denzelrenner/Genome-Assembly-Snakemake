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



# add hifiasm ploidy column
for s in SAMPLES:

    if not HASHIC:
        break

    
    if config['hic1_reads'] and s in os.listdir(config['hic1_reads']):

        # HIC1FILES[s] = [os.path.join(d,s,f) for f in os.listdir(os.path.join(config['hic1_reads'],s))]
        HICFILES[s]['HIC1'] = [os.path.join(config['hic1_reads'],s,f) for f in os.listdir(os.path.join(config['hic1_reads'],s))]

    if config['hic2_reads'] and s in os.listdir(config['hic2_reads']):

        # HIC2FILES[s] = [os.path.join(d,s,f) for f in os.listdir(os.path.join(config['hic2_reads'],s))]
        HICFILES[s]['HIC2'] = [os.path.join(config['hic2_reads'],s,f) for f in os.listdir(os.path.join(config['hic2_reads'],s))]

        # return ""
    if len(os.listdir(os.path.join(config['hic1_reads'],s))) + len(os.listdir(os.path.join(config['hic2_reads'],s))) > 2:
        MERGEHIC[s] = True

for s in SAMPLES:

    # generate combination of output files for smudgeplot
    smudgeFiles = expand("/s3_estimate_ploidy/smudgeplot/plots/{png}.png",png=['output_smudgeplot','output_smudgeplot_log10','output_centralities'])
    
    SMUDGEPLOTFILES[s] = [s + f for f in smudgeFiles]

    # # generate combination of output files for genomescope
    gscopeFiles = expand("/s3_estimate_ploidy/genomescope/plots/{png}.png",png=['transformed_log_plot','transformed_linear_plot','log_plot','linear_plot'])

    GENOMESCOPEFILES[s] = [s + f for f in gscopeFiles]


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


