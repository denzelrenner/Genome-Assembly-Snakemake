# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
import pandas as pd
import numpy as np

def list_directory(myDir:str,absPath:bool=True):
    
    if absPath:
        return [os.path.join(myDir,i) for i in os.listdir(myDir)]
    
    elif not absPath:
        return os.listdir(myDir)

def return_hic_input_files(hic1_input:str,hic2_input:str):

    hic1Files,hic2Files = [],[]

    # get HiC1 and HiC2 files by checking if a dir or space seprated list of files was given
    if os.path.isdir(hic1_input):
        hic1Files = [f for f in list_directory(hic1_input) if (f.endswith('1.fq.gz') or f.endswith('1.fastq.gz'))]

    # otherwise assume it is a space separated list
    elif not os.path.isdir(hic1_input):
        hic1Files = hic1_input.split()

    if os.path.isdir(hic2_input):
        hic2Files = [f for f in list_directory(hic2_input) if (f.endswith('2.fq.gz') or f.endswith('2.fastq.gz'))]

    # otherwise assume it is a space separated list
    elif not os.path.isdir(hic2_input):
        hic2Files = hic2_input.split()

    return hic1Files,hic2Files


# check whether input sample has HiC
def check_for_hic(row):

    if row['HiC1'] and row['HiC2']:
        hic1,hic2 = return_hic_input_files(row['HiC1'],row['HiC2'])
        if len(hic1) == len(hic2):
            return True
    return False

# reorder HiC files to make sure the concatenating order of the hic data is the same
def fix_hic_file_order(row):

    if row['HiC1'] and row['HiC2']:

        # get hic files
        hic1,hic2 = return_hic_input_files(row['HiC1'],row['HiC2'])

        # set var for output files, HiC1 leads the ordering
        hic1Output = []
        hic2Output = []

        # start looping and find mate file that shares the longest prefix with it
        for f in hic1:

            mateWithLongestSharedPrefix = max(hic2,key=lambda m:len(os.path.commonprefix([f.split('/')[-1],m.split('/')[-1]])))

            hic1Output.append(f)
            hic2Output.append(mateWithLongestSharedPrefix)

        # return space separated string again
        return ' '.join(hic1Output), ' '.join(hic2Output)

    return None,None

def should_hic_be_merged(row):

    # if there is actually any data in that column
    if row['has_hic']:

        # if more than HiC file was given
        hic1,hic2 = return_hic_input_files(row['HiC1'],row['HiC2'])
        if len(hic1) > 1:
            return True        
            
    return False

def set_hifiasm_ploidy(row):

    if row['has_hic']:
        return row['ploidy']
    return 2

def set_scaffolding_tool(row):

    if row['has_hic']:

        if row['ploidy'] > 2:
            return 'haphic'
        elif row['ploidy'] <= 2:
            return 'yahs'

    return 'none'

def separate_lineage_and_odb(row):

    return row['compleasm_lineage_odb'].split('_')

# based on extension we will process the input hifi reads differently
def set_hifi_reads_operation(row):

    if row['HiFi'].endswith('.fastq'):
        return 'compress'

    elif row['HiFi'].endswith('.fastq.gz'):
        return 'copy'

    elif row['HiFi'].endswith('.bam'):
        return 'convert'

     



# process samples

# load in csv
df = pd.read_csv(config['samples_file'],comment='#').set_index("sample", drop=False).fillna(np.nan).replace([np.nan], [None])

# get index values, which are sample names
SAMPLES = df.index.values.tolist()

# add new columns to samples

# hifiasm ploidy, set to 2 automatically if no HiC
df['has_hic'] = df.apply(check_for_hic,axis=1)
df['HiC1'],df['HiC2'] = zip(*df.apply(fix_hic_file_order,axis=1))
df['merge_hic'] = df.apply(should_hic_be_merged,axis=1)
df['hifiasm_ploidy'] = df.apply(set_hifiasm_ploidy,axis=1)
df['ploidyAsArray'] = df.apply(lambda row:[i for i in range(1,row['ploidy']+1)],axis=1)
df['hifiasmploidyAsArray'] = df.apply(lambda row:[i for i in range(1,row['hifiasm_ploidy']+1)],axis=1)
df['scaffolding_tool'] = df.apply(set_scaffolding_tool,axis=1)
df['lineage'],df['odb'] = zip(*df.apply(separate_lineage_and_odb,axis=1))
df['hifi_operation'] = df.apply(set_hifi_reads_operation,axis=1)

# convert df to dictionary, this is why we set sample as index in the pandas read csv, and we also dont drop its comment either
sampleInfo = df.to_dict("index")

# set HiC individuals and non HiC indivuals
HICSAMPLES = [s for s in sampleInfo if sampleInfo[s]['has_hic']]

NOHICSAMPLES = [s for s in sampleInfo if not sampleInfo[s]['has_hic']]

# create list of hifiasm ploidys, order is in the same order as the samples so no misassignemnets
HIFIASMPLOIDYFORSAMPLES = [sampleInfo[s]['hifiasm_ploidy'] for s in SAMPLES]

# create list of genome sizes, order is in the same order as the sample so no missagnments
GENOMESIZEFORSAMPLES = [sampleInfo[s]['haplomeSize'] for s in SAMPLES]

# set number of samples
NUMBEROFINPUTSAMPLES = len(SAMPLES)

# get all lineages to be download
ALLCOMPLEASMLINEAGES = list(set([sampleInfo[s]['compleasm_lineage_odb'] for s in SAMPLES]))

