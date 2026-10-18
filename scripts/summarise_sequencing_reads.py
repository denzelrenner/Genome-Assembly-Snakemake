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

vars.add_argument('-s','--sample',dest='sample',nargs='+', type=str, required=True, help='sample(s) to process stats for')

vars.add_argument('-o','--output_prefix',dest='outputPrefix', type=str, required=False,default='output', help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=True,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('-g','--genome_size',dest='genomeSize', type=int, required=True,default=210000000, help='expected haplome size (in bp)')

args = vars.parse_args()

# function to merge multiple dataframes, input is list of dataframe 
def conact_multiple_dataframes(dfs:list):
    return pd.concat(dfs) if len(dfs) > 1 else dfs[0]

# write df to excel and tsv, only single sheet in excel file
def write_to_tsv_and_excel(df:pd.DataFrame,outputDir:str,outputFilePrefix:str,excel=True,csv=True,seperator='\t',sheetName='Sheet1'):

    # make dirs if it does not exist
    os.makedirs(outputDir,exist_ok=True)

    # write to tsv and excel file
    if csv:
        df.to_csv(f'{outputDir}/{outputFilePrefix}.tsv',sep=seperator,index=False)

    if excel:
        df.to_excel(f'{outputDir}/{outputFilePrefix}.xlsx',sheet_name=sheetName,index=False)

# write multuple dfs to an excel spreadsheet
# df_sheetname is an array containing a pandas dataframe and its sheet name
def many_df_to_excel(df_sheetName:list,outputDir:str,outputFilePrefix:str):

    with pd.ExcelWriter(f'{outputDir}/{outputFilePrefix}.xlsx',engine='xlsxwriter') as writer:

        for df,sheetName in df_sheetName:

            df.to_excel(writer,sheet_name=sheetName,index=False)

            # get workbook and current sheet
            workbook = writer.book
            worksheet = writer.sheets[sheetName]

            # set format
            format = workbook.add_format({'text_wrap': True})

            # from col 1 to the end setting the column width to 200 characters, and adding format
            worksheet.set_column(1,len(list(df.columns.values)), 50, format)


# function to create data frame from seqkit tsv input
def create_seqkit_df(inputFilePath:str,sampleAsCol:bool=False,combineRows:bool=False):

    # open tsv file
    df = pd.read_csv(inputFilePath,sep='\t',header=0)

    # add column for sequencing depth
    df['sequencing_depth'] = round(df['sum_len']/args.genomeSize,2)

    # add column for total sequenced dna in gigabases
    df['sum_len (in Gigabases)'] = df['sum_len']//(10**9)

    # add column for file size, divide number of bytes by how much is in one gigabit (GiB), and also gigabyte (Gb)
    df['file_size (in Gib)'] = round(df['file'].apply(os.path.getsize)/(1024**3),2)
    df['file_size (in Gb)'] = round(df['file'].apply(os.path.getsize)/(10**9),2)

    # add , to everything that is not a string, i.e the first 3 columns of output
    numcols = list(set(df.columns.values) - set(['file','format','type']))
    df[numcols] = df[numcols].map(lambda x:f'{x:,}')

    # change index to be sample name
    if sampleAsCol:
        df['sample'] = [inputFilePath.split('/')[0]] * len(df)

    # combine multiple rows so they are one big result
    if combineRows:
        df = df.groupby('sample').agg(lambda x:'\n\n\n'.join(x)).reset_index()

    return df

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# set col width of excel spreadsheet cells
pd.set_option('display.max_colwidth', 0)

# run different analysis depending on if we have a single sample or many
if len(args.sample) == 1:

    # first check if we have hifi, hic1, and hic2 files
    hifi = (create_seqkit_df(f'{args.sample[0]}/s2_input_data_qc/Hifi/{args.sample[0]}.hifi.tsv'),'hifi') if os.path.exists(f'{args.sample[0]}/s2_input_data_qc/Hifi/{args.sample[0]}.hifi.tsv') else None
    hic1 = (create_seqkit_df(f'{args.sample[0]}/s2_input_data_qc/HiC1/{args.sample[0]}.hic1.tsv'),'hic1') if os.path.exists(f'{args.sample[0]}/s2_input_data_qc/HiC1/{args.sample[0]}.hic1.tsv') else None
    hic2 = (create_seqkit_df(f'{args.sample[0]}/s2_input_data_qc/HiC2/{args.sample[0]}.hic2.tsv'),'hic2') if os.path.exists(f'{args.sample[0]}/s2_input_data_qc/HiC2/{args.sample[0]}.hic2.tsv') else None

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

    # go through each of the samples
    for s in args.sample:

        # create dfs for each input data type
        hifi = create_seqkit_df(f'{s}/s2_input_data_qc/Hifi/{s}.hifi.tsv',sampleAsCol=True,combineRows=True) if os.path.exists(f'{s}/s2_input_data_qc/Hifi/{s}.hifi.tsv') else None
        hic1 = create_seqkit_df(f'{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv',sampleAsCol=True,combineRows=True) if os.path.exists(f'{s}/s2_input_data_qc/HiC1/{s}.hic1.tsv') else None
        hic2 = create_seqkit_df(f'{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv',sampleAsCol=True,combineRows=True) if os.path.exists(f'{s}/s2_input_data_qc/HiC2/{s}.hic2.tsv') else None

        # add to list of dfs to be merged if file was found
        if isinstance(hifi,pd.DataFrame):
            toMergeHifi.append(hifi)

        if isinstance(hic1,pd.DataFrame):
            toMergeHic1.append(hic1)

        if isinstance(hic2,pd.DataFrame):
            toMergeHic2.append(hic2)

    # merge and add to dataframe sheetname, explode to separate cells that are multiline, remove sample from array of columns be different size to the others
    if toMergeHifi:
        dataframes_sheetnames.append((conact_multiple_dataframes(dfs=toMergeHifi),'hifi'))

    if toMergeHic1:

        dataframes_sheetnames.append((conact_multiple_dataframes(dfs=toMergeHic1),'hic1'))

    if toMergeHic2:

        dataframes_sheetnames.append((conact_multiple_dataframes(dfs=toMergeHic2),'hic2'))

    # create excel out 
    many_df_to_excel(dataframes_sheetnames,outputDir=outdir,outputFilePrefix=args.outputPrefix)