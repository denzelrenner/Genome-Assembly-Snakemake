# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
import pandas as pd

def list_directory(myDir:str,absPath:bool=True):
    
    if absPath:
        return [os.path.join(myDir,i) for i in os.listdir(myDir)]
    
    elif not absPath:
        return os.listdir(myDir)

# check whether input sample has HiC
def check_for_hic(row):

    if row['HiC1'] and row['HiC2']:
        if len(row['HiC1'].split()) == len(row['HiC2'].split()):
            return True
    return False

# reorder HiC files to make sure the concatenating order of the hic data is the same
def fix_hic_file_order(row):

    if row['HiC1'] and row['HiC2']:

        # get list of each mate file
        hic1Files = row['HiC1'].split()
        hic2Files = row['HiC2'].split()

        # set var for output files, HiC1 leads the ordering
        hic1Output = []
        hic2Output = []

        # start looping and find mate file that shares the longest prefix with it
        for f in hic1Files:

            mateWithLongestSharedPrefix = max(hic2Files,key=lambda m:len(os.path.commonprefix([f,m])))

            hic1Output.append(f)
            hic2Output.append(mateWithLongestSharedPrefix)

        # return space separated string again
        return ' '.join(hic1Output), ' '.join(hic2Output)

    return []

def should_hic_be_merged(row):

    # if there is actually any data in that column
    if row['has_hic']:

        # if more than HiC file was given
        if len(row['HiC1'].split()) > 1:        
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
        sampleInfo[smpl]['hifi_operation'] = 'compress'

    elif row['HiFi'].endswith('.fastq.gz'):
        sampleInfo[smpl]['hifi_operation'] = 'copy'

    elif row['HiFi'].endswith('.bam'):
        sampleInfo[smpl]['hifi_operation'] = 'convert'




# process samples

# load in csv
df = pd.read_csv(config['samples_file'],comment='#').set_index("sample", drop=False)

# get index values, which are sample names
SAMPLES = df.index.values.tolist()

# add new columns to samples

# hifiasm ploidy, set to 2 automatically if no HiC
df['has_hic'] = df.apply(check_for_hic,axis=1)
df['HiC'],df['HiC2'] = zip(*df.apply(fix_hic_file_order,axis=1))
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

print(sampleInfo)
print(df)