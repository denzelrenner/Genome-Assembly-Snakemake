#!/usr/bin/env python3

import os
import argparse
import logging
import pandas as pd
from PowerhouseOfTheCell import *

# Input:
# name of a sample or space separated list of samples
# an expected haplome size
# an outdir and out prefix

# Output:
# if only a single sample is given an excel spreadsheet of information in sample/s2_input_data_qc/excel_summary/file.xlsx
# if a space separated list of samples is given an excel spreadsheet of information in collective_results/s2_input_data_qc/excel_summary/file.xlsx

vars = argparse.ArgumentParser(description='Create summary statistics from seqkit analysis for an individual sample or many samples')

vars.add_argument('-s','--samples',dest='sample',nargs='+', type=str, required=True, help='sample(s) to process stats for')

vars.add_argument('-o','--output_prefix',dest='outputPrefix', type=str, required=False,default='output', help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=True,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('-g','--genome_size',dest='genomeSize', nargs='+',type=int, required=True, help='expected haplome size (in bp), order should follow that of samples')

args = vars.parse_args()

# function to create data frame from seqkit tsv input
def seqkit_tsv_to_df(inputFilePath:str,haplomeSize:int,sampleAsCol:bool=False,combineRows:bool=False):

    # open tsv file
    df = pd.read_csv(inputFilePath,sep='\t',header=0)

    # add column for file size, divide number of bytes by how much is in one gigabit (GiB), and also gigabyte (Gb)
    df['file_size (in Gib)'] = round(df['file'].apply(os.path.getsize)/(1024**3),2)
    df['file_size (in Gb)'] = round(df['file'].apply(os.path.getsize)/(10**9),2)

    # add column for total sequenced dna in gigabases
    df['sum_len (in Gigabases)'] = df['sum_len']//(10**9)

    # add column for sequencing depth
    df['sequencing_depth'] = round(df['sum_len']/haplomeSize,2)

    # add , to everything that is not a string, i.e the first 3 columns of output
    numcols = list(set(df.columns.values) - set(['file','format','type']))
    df[numcols] = df[numcols].map(lambda x:f'{x:,}')

    # change index to be sample name
    if sampleAsCol:
        df['sample'] = [inputFilePath.split('/')[0]] * len(df)

    # combine multiple rows so they are one big result
    if combineRows:
        df = df.groupby('sample').agg(lambda x:'\n\n\n'.join(x)).reset_index()

    # add last 5 columns to start of column list
    cols = list(df.columns)
    cols = cols[-5:][::-1] + cols[:-5]

    # insert file column in first position
    cols.remove('file')
    cols.insert(1,'file')

    # subset the df using the new column arrangeemnt
    df = df[cols]

    return df

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# create dictionary of sample and haplome size
sample_haplomeSize = dict(zip(args.sample,args.genomeSize))

# run different analysis depending on if we have a single sample or many
if len(sample_haplomeSize) == 1:

    for s,hapSize in sample_haplomeSize.items():

        # first check if we have hifi, hic1, and hic2 files
        hifi = (seqkit_tsv_to_df(f'{s}/s2_input_data_qc/Hifi/{s}.hifi.tsv',haplomeSize=hapSize,sampleAsCol=True),'hifi') if os.path.exists(f'{s}/s2_input_data_qc/Hifi/{s}.hifi.tsv') else None
        hic1 = (seqkit_tsv_to_df(f'{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv',haplomeSize=hapSize,sampleAsCol=True),'hic1') if os.path.exists(f'{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv') else None
        hic2 = (seqkit_tsv_to_df(f'{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv',haplomeSize=hapSize,sampleAsCol=True),'hic2') if os.path.exists(f'{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv') else None

        # array to store tuples of df and sheet name
        dataframes_sheetnames = [i for i in [hifi,hic1,hic2] if i]

        many_df_to_excel(dataframes_sheetnames,outputDir=outdir,outputFilePrefix=args.outputPrefix)

# if many samples given we want summary stats for all of them
elif len(args.sample) > 1:

    # array to store tuples of df and sheet name
    dataframes_sheetnames = []

    # create a list to store the dfs from each sample and then merge them and add to dataframes_sheetname
    toMergeHifi = []
    toMergeHic1 = []
    toMergeHic2 = []

    # create sample index, this is important for merging cells in output, counter to help build start and stop, and stats from 1 because of 1 based indexing
    sample_index = {}
    counter = 1
    
    # go through each of the samples
    for s,hapSize in sample_haplomeSize.items():

        # create dfs for each input data type
        hifi = seqkit_tsv_to_df(f'{s}/s2_input_data_qc/Hifi/{s}.hifi.tsv',haplomeSize=hapSize,sampleAsCol=True) if os.path.exists(f'{s}/s2_input_data_qc/Hifi/{s}.hifi.tsv') else None
        hic1 = seqkit_tsv_to_df(f'{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv',haplomeSize=hapSize,sampleAsCol=True) if os.path.exists(f'{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv') else None
        hic2 = seqkit_tsv_to_df(f'{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv',haplomeSize=hapSize,sampleAsCol=True) if os.path.exists(f'{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv') else None

        # add to list of dfs to be merged if file was found
        if isinstance(hifi,pd.DataFrame):
            toMergeHifi.append(hifi)

        if isinstance(hic1,pd.DataFrame):
            toMergeHic1.append(hic1)

        if isinstance(hic2,pd.DataFrame):
            toMergeHic2.append(hic2)

        # hic1 and hic2 should be the same so we can use one to build the dict of sample index
        if isinstance(hic1,pd.DataFrame) or isinstance(hic2,pd.DataFrame):
                
            # get from samples column from concat df
            rows = hic1['sample'].to_list()

            # add start index
            sample_index[s] = {'start':counter,'stop':0}

            # add to the counter
            counter += len(rows)

            # add stop index
            sample_index[s]['stop'] = counter - 1


    # merge and add to dataframe sheetname, explode to separate cells that are multiline, remove sample from array of columns be different size to the others
    if toMergeHifi:
        dataframes_sheetnames.append((concat_multiple_dataframes(dfs=toMergeHifi),'hifi'))

    if toMergeHic1:
        dataframes_sheetnames.append((concat_multiple_dataframes(dfs=toMergeHic2),'hic1'))

    if toMergeHic2:
        dataframes_sheetnames.append((concat_multiple_dataframes(dfs=toMergeHic2),'hic2'))


    # create excel out 
    with pd.ExcelWriter(f'{outdir}/{args.outputPrefix}.xlsx',engine='xlsxwriter') as writer:

        for df,sheetName in dataframes_sheetnames:

            # write to excel out
            df.to_excel(writer,sheet_name=sheetName,index=False)

            # get workbook and current sheet
            workbook = writer.book
            worksheet = writer.sheets[sheetName]

            # set format
            format = workbook.add_format({'text_wrap': True})

            # from col 1 to the end setting the column width to 200 characters, and adding text wrap format
            worksheet.set_column(0,len(df.columns.tolist()), 50, format)

            # merge columns for hic1 and hic2, we dont do this for hifi because there will only ever be one hifi file
            # we will always have more than hic1/hic2 file if its available because we have the file that is copied to the working dir, and the original file
            if sheetName in ['hic1','hic2']:

                # merge rows for each sample
                for s in sample_index:

                    worksheet.merge_range(sample_index[s]['start'], 0, sample_index[s]['stop'], 0, s, format)
                    