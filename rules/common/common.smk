# helper function to be used in the actual snakefile, and in other snakemake files as well
import os

# set some global variables
HASHIC = True if (config['hic1_reads'] and config['hic2_reads']) else False

# get samples by reading hifi reads file
SAMPLES = [f.split('.')[0] for f in os.listdir(config['hifi_reads'])]
HAPLOMENUMBERS = [str(i) for i in range(1,config['ploidy']+1)] if HASHIC else ['1','2']

# join path of hic reads and sample, then save in dir
HIC1FILES={s:"" for s in SAMPLES}
HIC2FILES={s:"" for s in SAMPLES}
HICFILES={s:{"HIC1":"","HIC2":""} for s in SAMPLES}
MERGEHIC={s:False for s in SAMPLES}
SMUDGEPLOTFILES={s:"" for s in SAMPLES}
GENOMESCOPEFILES={s:"" for s in SAMPLES}

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
    haplomeNumber='|'.join(HAPLOMENUMBERS) # limit value of haplome number


# helper functions

### GENERAL ###

def set_rule_target():

    # return [f'{s}/s0_end_pipe/done.txt' for s in SAMPLES]

    return expand('{sample}/s0_end_pipe/{sample}.hap{haplomeNumber}/done.txt',sample=SAMPLES,haplomeNumber=HAPLOMENUMBERS)

def get_directory_array(myDir:str,absPath:bool=True):
    
    if absPath:
        return [os.path.join(myDir,i) for i in os.listdir(myDir)]
    
    elif not absPath:
        return os.listdir(myDir)

### HIFIASM ###
# def get_hifiasm_output_files(wildcards):

#     if HASHIC:

#         # return list of files
#         return expand('{wildcards.sample}/s3_run_hifiasm/{wildcards.sample}.hic.hap{hapNumber}.p_ctg.gfa',hapNumber=HAPLOMENUMBERS)
    
#     else:

#         # return list of files
#         return expand('{wildcards.sample}/s3_run_hifiasm/{wildcards.sample}.hap{hapNumber}.p_ctg.gfa',hapNumber=HAPLOMENUMBERS)

# def get_hifiasm_qc_output_files(wildcards):

#     if HASHIC:

#         # return list of files
#         return expand('{wildcards.sample}/s3_run_hifiasm/quality_control/{wildcards.sample}/{wildcards.sample}.hic.hap{hapNumber}/Compleasm',hapNumber=HAPLOMENUMBERS)
    
#     else:

#         # return list of files
#         return expand('{wildcards.sample}/s3_run_hifiasm/{wildcards.sample}.hap{hapNumber}.fa',hapNumber=HAPLOMENUMBERS)

# def get_hifiasm_input_files(wildcards):

#     if HASHIC:

#         # return list of files
#         return expand('{wildcards.sample}/s3_run_hifiasm/{wildcards.sample}.hic.hap{hapNumber}.p_ctg.gfa',hapNumber=HAPLOMENUMBERS)
    
#     else:

#         # return list of files
#         return expand('{wildcards.sample}/s3_run_hifiasm/{wildcards.sample}.hap{hapNumber}.p_ctg.gfa',hapNumber=HAPLOMENUMBERS)

# # get the HiC files to use for hifiasm
# def get_hifiasm_hic_input(pairedEnd):

#     # return empty string if paths to hic1 and hic2 data not given
#     if not HASHIC:
#         return ''

#     if pairedEnd == 'hic1':

#         # get array of hic files
#         hic1Files = get_directory_array(config['hic1_reads'])

#         # check if the dir has more than one HiC file, in that case there should be a merged HiC file, otherwise if only one HiC file then just return that
#         if len(hic1Files) > 1:
#             return [f for f in hic1Files if 'merged' in f][0]
        
#         elif len(hic1Files) == 1:
#             return hic1Files[0]
    
#     elif pairedEnd == 'hic2':

#         # get array of hic files
#         hic2Files = get_directory_array(config['hic2_reads'])

#         # check if the dir has more than one HiC file, in that case there should be a merged HiC file, otherwise if only one HiC file then just return that
#         if len(hic2Files) > 1:
#             return [f for f in hic2Files if 'merged' in f][0]
        
#         elif len(hic2Files) == 1:
#             return hic2Files[0]

# def set_hifiasm_hic_parameter(wildcards):
#     pass

