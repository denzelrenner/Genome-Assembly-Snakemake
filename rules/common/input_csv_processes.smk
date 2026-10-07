# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
import pandas as pd

def get_directory_array(myDir:str,absPath:bool=True):
    
    if absPath:
        return [os.path.join(myDir,i) for i in os.listdir(myDir)]
    
    elif not absPath:
        return os.listdir(myDir)

# functions
def set_hifiasm_ploidy(row):

    if row['useHiC']:
        return row['ploidy']
    return 2

# the total hic files in the dir
def set_hic1_files_to_merge(row):

    if row['useHiC']:

        # get list of hic files
        sampleHic1Dir = os.path.join(config['hic1_reads'],row['sample'])

        hic1Files = get_directory_array(sampleHic1Dir)

        # get input file name for downstream tools, this is dependent on how many files are in the HiC dir because thta determine if the file is merged or not
        hic1OutputFile = []

        if len(hic1Files) == 1:
            hic1OutputFile = hic1Files[0]

        elif len(hic1Files) > 1:
            hic1OutputFile = os.path.join(sampleHic1Dir,f"{row['sample']}.merged_1.fq.gz")

        else:
            hic1OutputFile = []

        return [' '.join(hic1Files),hic1OutputFile]
    
    return [[],[]]

# the total hic files in the dir
def set_hic2_files_to_merge(row):

    if row['useHiC']:

        # get list of hic files
        sampleHic2Dir = os.path.join(config['hic2_reads'],row['sample'])

        hic2Files = get_directory_array(sampleHic2Dir)

        # get input file name for downstream tools, this is dependent on how many files are in the HiC dir because thta determine if the file is merged or not
        hic2OutputFile = []

        if len(hic2Files) == 1:
            hic1OutputFile = hic1Files[0]

        elif len(hic2Files) > 1:
            hic2OutputFile = os.path.join(sampleHic2Dir,f"{row['sample']}.merged_2.fq.gz")

        else:
            hic2OutputFile = []

        return [' '.join(hic2Files),hic2OutputFile]

    return [[],[]]

def set_scaffolding_tool(row):

    if row['ploidy'] > 2:
        return 'haphic'

    elif row['ploidy'] <= 2:
        return 'yahs'




# process samples

# load in csv
df = pd.read_csv(config['samples_file']).set_index("sample", drop=False)

# get index values, which are sample names
SAMPLES = df.index.values.tolist()

# make HiC directory for each sample if it does not exist already
for sample in SAMPLES:

    os.makedirs(os.path.join(config['hic1_reads'],sample),exist_ok=True)
    os.makedirs(os.path.join(config['hic2_reads'],sample),exist_ok=True)


# add new columns to samples

# hifiasm ploidy, set to 2 automatically if no HiC
df['hifiasm_ploidy'] = df.apply(set_hifiasm_ploidy,axis=1)
df['hic1_files'],df['hic1_input_file'] = zip(*df.apply(set_hic1_files_to_merge,axis=1))
df['hic2_files'],df['hic2_input_file'] = zip(*df.apply(set_hic2_files_to_merge,axis=1))
df['scaffolding_tool'] = df.apply(set_scaffolding_tool,axis=1)

# convert df to diction
sampleInfo = df.to_dict("index")

# add all HiC files we expect to have in input dir
for smpl in sampleInfo:

    # if we are using HiC which implies there is a HiC file that exists
    if sampleInfo[smpl]['useHiC']:

        if (sampleInfo[smpl]['hic1_input_file'] in sampleInfo[smpl]['hic1_files']):
            sampleInfo[smpl]['all_hic1_files'] = sampleInfo[smpl]['hic1_files']
            sampleInfo[smpl]['all_hic2_files'] = sampleInfo[smpl]['hic2_files']

        elif (sampleInfo[smpl]['hic1_input_file'] not in sampleInfo[smpl]['hic1_files']):
            sampleInfo[smpl]['all_hic1_files'] = ' '.join([sampleInfo[smpl]['hic1_input_file'], sampleInfo[smpl]['hic1_files']])
            sampleInfo[smpl]['all_hic2_files'] = ' '.join([sampleInfo[smpl]['hic2_input_file'], sampleInfo[smpl]['hic2_files']])
    
    # if we are not using HiC then we dont even include it 
    else:
        sampleInfo[smpl]['all_hic1_files'] = []
        sampleInfo[smpl]['all_hic2_files'] = []

# add merge HiC column
for smpl in sampleInfo:

    sampleInfo[smpl]['merge_hic'] = 'false'

    if (sampleInfo[smpl]['hic1_input_file'] not in sampleInfo[smpl]['hic1_files']) and sampleInfo[smpl]['useHiC']:
        sampleInfo[smpl]['merge_hic'] = 'true'

# add has bash true false or use HiC
for smpl in sampleInfo:
    
    sampleInfo[smpl]['use_hic_shell'] = 'false'

    if sampleInfo[smpl]['useHiC']:
        sampleInfo[smpl]['use_hic_shell'] = 'true'
